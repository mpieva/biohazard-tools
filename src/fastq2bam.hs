import Bio.Bam
import Bio.Bam.Evan ( removeWarts )
import Bio.Iteratee.ZLib
import Bio.Prelude
import Paths_biohazard_tools ( version )
import System.Console.GetOpt
import System.IO

import qualified Data.ByteString       as B
import qualified Data.ByteString.Char8 as S
import qualified Data.Vector.Generic   as V

data Opts = Opts { output  :: BamMeta -> Iteratee [BamRec] IO ()
                 , inputs  :: [Input]
                 , verbose :: Bool
                 , merge   :: Maybe Int
                 , lowqual :: Int }

defaultOpts :: Opts
defaultOpts = Opts { output  = protectTerm . pipeBamOutput
                   , inputs  = []
                   , verbose = False
                   , merge   = Nothing
                   , lowqual = 20 }

data Input = Input { _read1  :: FilePath         -- ^ file with first read (or other stuff)
                   ,  read2  :: Maybe FilePath   -- ^ optional file with second read
                   , index1  :: Maybe FilePath   -- ^ optional file with first index
                   , index2  :: Maybe FilePath   -- ^ optional file with second index
                   , lindex1 :: Int }           -- ^ length of first index contained in first read
  deriving Show

plainInput :: FilePath -> Input
plainInput fn = Input fn Nothing Nothing Nothing 0

getopts :: [String] -> ([Opts -> IO Opts], [String], [String])
getopts = getOpt (ReturnInOrder add_read1) options
  where
    options =
        [ Option "o" ["output"]         (ReqArg set_output "FILE") "Write output to FILE"
        , Option "1" ["read-one"]        (ReqArg add_read1 "FILE") "Parse FILE for anything"
        , Option "2" ["read-two"]        (ReqArg add_read2 "FILE") "Parse FILE for second mates"
        , Option "I" ["index-one"]        (ReqArg add_idx1 "FILE") "Parse FILE for first index"
        , Option "J" ["index-two"]        (ReqArg add_idx2 "FILE") "Parse FILE for second index"
        , Option "l" ["length-index-one"] (ReqArg set_lidx1 "NUM") "Read 1 ends on NUM index bases"
        , Option "m" ["merge-overlap"]   (OptArg set_merge "QUAL") "Attempt to merge or trim reads"
        , Option "q" ["merge-qual"]       (ReqArg set_qual "QUAL") "Minimum quality for merge is QUAL"
        , Option "v" ["verbose"]               (NoArg set_verbose) "Print progress information"
        , Option "h?" ["help","usage"]               (NoArg usage) "Print this helpful message" ]

    set_output "-" c = return $ c { output  = pipeBamOutput }
    set_output  fn c = return $ c { output  = writeBamFile fn }
    set_verbose    c = return $ c { verbose = True }

    set_merge Nothing  c =                    return $ c { merge   = Just 200 }
    set_merge (Just a) c = readIO a >>= \m -> return $ c { merge   = Just  m  }
    set_qual        a  c = readIO a >>= \q -> return $ c { lowqual =       q  }

    add_read1 fn c = return $ c { inputs = plainInput fn : inputs c }
    add_read2 fn c = return $ c { inputs = at_head (\i -> i { read2  = Just fn }) (inputs c) }
    add_idx1  fn c = return $ c { inputs = at_head (\i -> i { index1 = Just fn }) (inputs c) }
    add_idx2  fn c = return $ c { inputs = at_head (\i -> i { index2 = Just fn }) (inputs c) }

    set_lidx1  a c = readIO a >>= \n -> return $  c { inputs = at_head (\i -> i { lindex1 = n}) (inputs c) }

    at_head f [    ] = [ f $ plainInput "-" ]
    at_head f (i:is) = f i : is

    usage _ = do pn <- getProgName
                 let t = "Usage: " ++ pn ++ " [OPTION...]\n" ++
                         "Reads multiple FastA or FastQ files and converts them to BAM.  See manpage for details."
                 hPutStrLn stderr $ usageInfo t options
                 exitSuccess


main :: IO ()
main = do (opts, [], errors) <- getopts `fmap` getArgs
          unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
          conf <- foldl (>>=) (return defaultOpts) opts
          pgm <- addPG (Just version)

          let eff_inputs = if null (inputs conf) then [ plainInput "-" ] else inputs conf
          when (verbose conf) $ mapM_ (hPrint stderr) eff_inputs

          foldr ((>=>) . enum_input) run (reverse eff_inputs) $
                joinI $ progress (verbose conf) $
                joinI $ concatMapStream (maybe concatDuals (mergeDuals $ lowqual conf) $ merge conf) $
                output conf (pgm mempty)


type UpToTwo a = (a, Maybe a)

one :: a -> UpToTwo a
one a = (a, Nothing)

two :: a -> a -> UpToTwo a
two a b = (a, Just b)

mapU2 :: (a -> b) -> UpToTwo a -> UpToTwo b
mapU2 f (a,b) = (f a, fmap f b)

concatDuals :: UpToTwo a -> [a]
concatDuals (a,Just  b) = [a,b]
concatDuals (a,Nothing) = [ a ]

mergeDuals :: Int -> Int -> UpToTwo BamRec -> [BamRec]
mergeDuals lowq highq (r1, Just r2) = mergeBam lowq highq default_fwd_adapters default_rev_adapters r1 r2
mergeDuals lowq highq (r1, Nothing) = trimBam lowq highq default_fwd_adapters r1

-- Enumerates a file.  Sequence and quality end up in b_seq and b_qual.
fromFastq :: (MonadIO m, MonadMask m) => FilePath -> Enumerator [BamRec] m a
fromFastq fp = enumAny fp $= enumInflateAny $= parseFastqCassava $= mapStream removeWarts
  where
    enumAny "-" = enumHandle defaultBufSize stdin
    enumAny  f  = enumFile   defaultBufSize   f

enum_input :: (MonadIO m, MonadMask m) => Input -> Enumerator [UpToTwo BamRec] m a
enum_input (Input r1 mr2 mi1 mi2 il1) = enum $= mapStream (addIdx il1)
  where
    enum = withIndex mi1 "XI" "YI" $ withIndex mi2 "XJ" "YJ" $
           maybe (fromFastq r1 $= mapStream one) (enumDual r1) mr2

addIdx :: Int -> UpToTwo BamRec -> UpToTwo BamRec
addIdx 0 brs = brs
addIdx l (br1, mbr2) = ( doext br1', fmap doext mbr2 )
  where
    l' = V.length (b_seq br1) - l

    br1'     = br1 { b_seq  = V.take l' (b_seq br1), b_qual = V.take l' (b_qual br1) }
    doext br = br  { b_exts = updateE "XI" (Text xi) $ updateE "YI" (Text yi) $ b_exts br }

    xi = S.pack . map showNucleotides . V.toList . V.drop l' $ b_seq  br1
    yi = B.pack . map ((+) 33 . unQ)  . V.toList . V.drop l' $ b_qual br1

-- Given an enumerator and maybe a filename, read index sequences from
-- the file and merge them with the numerator.
withIndex :: (MonadIO m, MonadMask m)
          => Maybe FilePath -> BamKey -> BamKey
          -> Enumerator [UpToTwo BamRec] m a -> Enumerator [UpToTwo BamRec] m a
withIndex Nothing      _    _ enum = enum
withIndex (Just fp) tagi tagq enum = mergeEnums enum (fromFastq fp) (convStream combine)
  where
    combine = do seqrecs <- lift headStream
                 idxrec  <- headStream
                 when (b_qname (fst seqrecs) /= b_qname idxrec) . error $
                        "read names do not match: " ++ shows (b_qname (fst seqrecs)) " & " ++ show (b_qname idxrec)

                 let idxseq  = S.pack $ map showNucleotides $ V.toList $ b_seq idxrec
                     idxqual = B.pack $ map   ((+33) . unQ) $ V.toList $ b_qual idxrec
                 return [ flip mapU2 seqrecs $
                        \r -> r { b_exts = (if B.null idxqual then id else insertE tagq (Text idxqual))
                                         $ insertE tagi (Text idxseq) $ b_exts r } ]

-- Enumerate dual files.  We read two FastQ files and match them up.  We
-- must make sure the names match, and we will flag everything as
-- 1st/2nd mate, no matter if the syntactic warts were present in the
-- files themselves.
enumDual :: (MonadIO m, MonadMask m)
         => FilePath -> FilePath -> Enumerator [UpToTwo BamRec] m a
enumDual f1 f2 = mergeEnums (fromFastq f1 $= mapStream one) (fromFastq f2) (convStream combine)
  where
    combine = do (firstMate, Nothing) <- lift headStream
                 secondMate           <- headStream

                 when (b_qname firstMate /= b_qname secondMate) . error $
                        "read names do not match: " ++ shows (b_qname firstMate) " & " ++ show (b_qname secondMate)

                 let qc = (b_flag firstMate .|. b_flag secondMate) .&. flagFailsQC
                     addx k = maybe id (updateE k) $ maybe (lookup k (b_exts secondMate)) Just $ lookup k (b_exts firstMate)
                     add_indexes = addx "XI" . addx "XJ" . addx "YI" . addx "YJ"

                 return [ two (firstMate  { b_flag = qc .|.  flagFirstMate .|. flagPaired .|. b_flag firstMate .&. complement flagSecondMate
                                          , b_exts = add_indexes $ b_exts firstMate })
                              (secondMate { b_flag = qc .|. flagSecondMate .|. flagPaired .|. b_flag secondMate .&. complement flagFirstMate
                                          , b_exts = add_indexes $ b_exts secondMate }) ]


progress :: MonadIO m => Bool -> Enumeratee [UpToTwo BamRec] [UpToTwo BamRec] m b
progress False = mapChunks id
progress True  = eneeCheckIfDonePass (icont . go 0 0)
  where
    go !_ !_ k (EOF         mx) = idone (liftI k) (EOF mx)
    go !l !n k (Chunk    [   ]) = liftI $ go l n k
    go !l !n k (Chunk as@(a:_)) = do
        let !n' = n + length as
            !nm = b_qname (fst a)
            !l' = l `max` S.length nm
        when (n .&. complement 0x1fff /= n' .&. complement 0x1fff) $ liftIO $ do
            hPutStr stderr $ "\27[K" ++
                replicate (l' - S.length nm) ' '
                ++ S.unpack nm ++ ", "
                ++ shows n' " records processed\n"
            hFlush stderr
        eneeCheckIfDonePass (icont . go l' n') . k $ Chunk as


