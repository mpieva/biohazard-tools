{-# LANGUAGE CPP, Rank2Types #-}
{-
This is a validator/fixup for paired end BAM files, that is more
efficient than 'samtools sort -n' followed by 'samtools fixmate'.

We want both: to quickly join separate mates together again from the
information about the mate's mapping coordinate, but at the same time
deal with broken files where that doesn't actually work.  Whenever we
join mates, we also check if the flags are consistent and fix them if
they aren't.

In the end, the code will work...

 - splendidly, if mates are already adjacent, in which case everything
   is streamed.
 - well, if the input is sorted properly, in which case most reads
   stream, but improper pairs need to queue until the mate is reached.
 - reasonably, if there are occasional widows, which will be queued
   to the very end and sorted by hashed-qname before they are recognized
   and repaired.
 - awkwardly, if sorting is violated, flags are wrong or widows are
   the rule, because then it degenerates to a full sort by qname.
-}

import Bio.Bam                           hiding ( mergeInputs, combineCoordinates )
import Bio.Prelude                       hiding ( yield )
import Bio.Util.Numeric                         ( showNum )
import Control.Concurrent.Async
import Control.Concurrent.STM.TQueue
import Control.Concurrent.STM.TVar
import Paths_biohazard_tools                    ( version )
import PriorityQueue
import System.Console.GetOpt
import System.Process
#if MIN_VERSION_process(1,2,1)
                                         hiding ( createPipe )
#endif

import qualified Data.ByteString as S
import qualified Data.ByteString.Builder as B
import qualified Data.Vector.Generic as V

data Verbosity = Silent | Errors | Warnings | Notices deriving (Eq, Ord)
data KillMode  = KillNone | KillUu | KillAll deriving (Eq, Ord)

data Config = CF { report_mrnm  :: !Bool
                 , report_mpos  :: !Bool
                 , report_isize :: !Bool
                 , report_flags :: !Bool
                 , report_fflag :: !Bool
                 , report_ixs   :: !Bool
                 , verbosity    :: Verbosity
                 , killmode     :: KillMode
                 , infilter     :: BamPair -> Bool
                 , output       :: BamMeta -> Iteratee [BamRec] IO ExitCode
                 , pqconf       :: PQ_Conf
                 , fixsven      :: Maybe Int }

config0 :: IO Config
config0 = return $ CF True True False True False True Errors KillNone
                      (const True) (fmap (const ExitSuccess) . protectTerm . pipeBamOutput)
                      (PQ_Conf 1000 200 "/var/tmp/" undefined) Nothing

options :: [OptDescr (Config -> IO Config)]
options = [
    Option "o" ["output"]                                  (ReqArg set_output "FILE") "Write output to FILE",
    Option "X" ["exec"]                                    (NoArg             return) "Send FastQ output to a program",
    Option "n" ["dry-run","validate"]                      (NoArg       set_validate) "No output, validate only",
    Option "k" ["kill-widows"]     (NoArg (\c -> return $ c { killmode =  KillAll })) "Delete all widows",
    Option "u" ["kill-unmapped"]   (NoArg (\c -> return $ c { killmode =   KillUu })) "Delete unmapped widows",
    Option [ ] ["kill-none"]       (NoArg (\c -> return $ c { killmode = KillNone })) "Never delete widows (default)",

    Option "v" ["verbose"]        (NoArg (\c -> return $ c { verbosity = Notices  })) "Print informational messages",
    Option "w" ["warnings"]       (NoArg (\c -> return $ c { verbosity = Warnings })) "Print warnings and errors",
    Option [ ] ["errors"]         (NoArg (\c -> return $ c { verbosity = Errors   })) "Print only errors (default)",
    Option "q" ["quiet"]          (NoArg (\c -> return $ c { verbosity = Silent   })) "Print only fatal errors",

    Option "" ["report-mrnm"]     (NoArg (\c -> return $ c { report_mrnm  =  True })) "Report wrong mate reference name (default yes)",
    Option "" ["report-mpos"]     (NoArg (\c -> return $ c { report_mpos  =  True })) "Report wrong mate position (default yes)",
    Option "" ["report-isize"]    (NoArg (\c -> return $ c { report_isize =  True })) "Report wrong insert size (default no)",
    Option "" ["report-flags"]    (NoArg (\c -> return $ c { report_flags =  True })) "Report wrong flags (default yes)",
    Option "" ["report-fflag"]    (NoArg (\c -> return $ c { report_fflag =  True })) "Report commonly inconsistent flags (default no)",
    Option "" ["report-ixs"]      (NoArg (\c -> return $ c { report_ixs   = False })) "Report mismatched index fields (default yes)",

    Option "" ["no-report-mrnm"]  (NoArg (\c -> return $ c { report_mrnm  = False })) "Do not report wrong mate reference name",
    Option "" ["no-report-mpos"]  (NoArg (\c -> return $ c { report_mpos  = False })) "Do not report wrong mate position",
    Option "" ["no-report-isize"] (NoArg (\c -> return $ c { report_isize = False })) "Do not report wrong insert size",
    Option "" ["no-report-flags"] (NoArg (\c -> return $ c { report_flags = False })) "Do not report wrong flags",
    Option "" ["no-report-fflag"] (NoArg (\c -> return $ c { report_fflag = False })) "Do not report commonly inconsistent flags",
    Option "" ["no-report-ixs"]   (NoArg (\c -> return $ c { report_ixs   = False })) "Do not report mismatched index fields",

    Option "" ["only-mapped"]   (NoArg (\c -> return $ c { infilter = mapped_only })) "Ignore totally unmapped input",
    Option "" ["fix-sven"]                                (ReqArg set_fixsven "QUAL") "Trim 3' ends of avg qual lower than QUAL",

    Option "M" ["max-memory"]                           (ReqArg set_max_mem     "MB") "Use at most MB megabytes per queue",
    Option "L" ["max-merge"]                            (ReqArg set_max_merge  "NUM") "Merge at most NUM files at a time",
    Option "T" ["temp-path"]                            (ReqArg set_temp_path "PATH") "Store temporary files in PATH",

    Option "h?" ["help","usage"] (NoArg usage) "Print this helpful message and exit",
    Option "V"  ["version"]      (NoArg  vrsn) "Print version number and exit" ]
  where
    usage _ = do pn <- getProgName
                 let blah = "Usage: " ++ pn ++ " [OPTION...] [FILE...]\n" ++
                            "Merge BAM files, rearrange them to move mate pairs together, \n" ++
                            "output a file with consistent mate pair information.  Alternatively, \n" ++
                            "pipe multiple FastQ streams in lock step to a program."
                 hPutStrLn stderr $ usageInfo blah options
                 exitSuccess

    vrsn _ = do pn <- getProgName
                hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
                exitSuccess

    set_output "-" c = return $ c { output = fmap (const ExitSuccess) . pipeBamOutput }
    set_output  f  c = return $ c { output = fmap (const ExitSuccess) . writeBamFile f }
    set_validate   c = return $ c { output = \_ -> ExitSuccess <$ skipToEof }
    set_fixsven  a c = readIO a >>= \q -> return $ c { fixsven = Just q }

    set_max_mem   a c = readIO a >>= \x -> return $ c { pqconf = (pqconf c) { max_mb    = x } }
    set_max_merge a c = readIO a >>= \x -> return $ c { pqconf = (pqconf c) { max_merge = x } }
    set_temp_path a c =                    return $ c { pqconf = (pqconf c) { temp_path = a } }


mapped_only :: BamPair -> Bool
mapped_only p = case p of
        Singleton a -> okay a
        LoneMate  a -> okay a
        Pair    a b -> okay a || okay b
  where
    okay = (\r -> not (isUnmapped r) || (isPaired r && not (isMateUnmapped r))) . unpackBam

pipe_to :: FilePath -> [String] -> ([Config -> IO Config], t1, t2) -> ([Config -> IO Config], t1, t2)
pipe_to cmd args (opts, errs, fs) = (mkout : opts, errs, fs)
  where
    mk1out key test (as, flush, qs, vs, ps, rfds)
        | all (/= key) as = return (as, flush, qs, vs, ps, rfds)
        | otherwise = do
            (pout, pin) <- createPipe
            setFdOption pin CloseOnExec True
            queue <- newTQueueIO
            vnum <- newTVarIO (0::Int)
            pid <- async $ flush_fastq queue vnum pin
            link pid

            return ( map (\a -> if a == key then "/dev/fd/" ++ show pout else a) as
                   , \br -> when (test br) (modifyTVar' vnum succ >> writeTQueue queue (Just br)) >> flush br
                   , queue : qs
                   , vnum : vs
                   , pid : ps
                   , pout : rfds )

    mkout cfg = do
        (args', flush_bam, queues, vars, pids, rfds) <- mk1out "CLOWNS" isFirstMate =<<
                                                        mk1out "JOKERS" isSecondMate =<<
                                                        mk1out "MIDDLE" (not . isPaired)
                                                          (args, const (return ()), [], [], [], [])
        pid_cmd <- spawnProcess cmd args'
        mapM_ closeFd rfds

        return $ cfg { output = \_ -> do
            mapStreamM_ (\br -> atomically $ do ns <- mapM readTVar vars
                                                when (minimum ns > 64) retry
                                                flush_bam br)
            liftIO $ do atomically $ mapM_ (flip writeTQueue Nothing) queues
                        mapM_ wait pids
                        waitForProcess pid_cmd }

    flush_fastq qq nn fd = do
            mbr <- atomically $ readTQueue qq <* modifyTVar' nn pred
            case mbr of
                Just br -> do fdPutLazy fd . B.toLazyByteString $
                                    B.char8 '@' <> B.byteString (b_qname br) <>
                                    (if isFirstMate  br then B.char8 '/' <> B.char8 '1' else mempty) <>
                                    (if isSecondMate br then B.char8 '/' <> B.char8 '2' else mempty) <>
                                    B.char8 '\n' <> V.foldr ((<>) . B.char8 . showNucleotides) mempty (b_seq br) <>
                                    B.char8 '\n' <> B.char8 '+' <> B.char8 '\n' <>
                                    V.foldr ((<>) . B.word8 . (+) 33 . unQ) (B.char8 '\n') (b_qual br)
                              flush_fastq qq nn fd
                Nothing -> return ()


main :: IO ExitCode
main = do (args,cmd) <- break (`elem` ["-X","--exec"]) `fmap` getArgs
          let (opts, files, errors) = (case cmd of _:cmd':args' -> pipe_to cmd' args' ; _ -> id)
                                      $ getOpt Permute options args

          unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
          config <- foldl (>>=) config0 opts
          add_pg <- addPG $ Just version
          mergeInputs files >=> run                          $ \hdr ->
            filterStream (infilter config)                  =$
            re_pair config (meta_refs hdr)                  =$
            mapChunks (maybe id do_trim (fixsven config))   =$
            (output config) (add_pg hdr)


-- | Fix a pair of reads.  Right now fixes their order and checks that
-- one is 1st mate, the other 2nd mate.  More fixes to come.

fixmate :: MonadIO m => BamRaw -> BamRaw -> Mating r m [BamRec]
fixmate r s | isFirstMate (unpackBam r) && isSecondMate (unpackBam s) = sequence [go r s, go s r]
            | isSecondMate (unpackBam r) && isFirstMate (unpackBam s) = sequence [go s r, go r s]
            | otherwise = liftIO $ do hPutStrLn stderr $ "Names match, but 1st mate / 2nd mate flags do not: "
                                                       ++ unpack (b_qname (unpackBam r))
                                      hPutStrLn stderr $ "There is no clear way to fix this file.  Giving up."
                                      exitFailure
  where
    -- position of 5' end
    pos5 a = if isReversed a then b_pos a + alignedLength (b_cigar a) else b_pos a

    -- transfer info from b to a
    go p q | null problems = return a
           | otherwise = do infos <- filter (not . null) `fmap` sequence [ m | (_,_,m) <- problems ]
                            unless (null infos) $ liftIO $ hPutStrLn stderr $ message infos
                            return $ foldr (\(_,m,_) -> m) a problems
      where
        a = unpackBam p
        b = unpackBam q

        problems = filter (\(x,_,_) -> not x) checks
        checks = [ (b_mrnm a  == b_rname b,     \x -> x { b_mrnm  = b_rname b },     count_mrnm)
                 , (b_mpos a  == b_pos b,       \x -> x { b_mpos  = b_pos b },       count_mpos)
                 , (b_isize a == computedIsize, \x -> x { b_isize = computedIsize }, count_isize)
                 , (b_flag a === computedFlag,  \x -> x { b_flag  = computedFlag },  count_flags)
                 , (b_flag a =!= computedFlag,  \x -> x { b_flag  = computedFlag },  count_fflag)
                 , (b_indices a == common_indices, setIndices common_indices, count_ixs) ]

        message infos = "fixing " ++ shows (b_qname a `S.append` if isFirstMate a then "/1" else "/2")
                        ": \t" ++ intercalate ", " infos
        !computedFlag' = (if b_rname a == invalidRefseq then (.|. flagUnmapped) else id) .
                         (if b_rname b == invalidRefseq then (.|. flagMateUnmapped) else id) .
                         (if isReversed b then (.|. flagMateReversed) else (.&. complement flagMateReversed)) .
                         (if isUnmapped b then (.|. flagMateUnmapped) else (.&. complement flagMateUnmapped)) .
                         (if isFailsQC  b then (.|. flagFailsQC) else id) $
                         b_flag a

        !properly_paired = computedFlag' .&. (flagUnmapped .|. flagMateUnmapped) == 0 && b_rname a == b_rname b
        !computedFlag    = if properly_paired then computedFlag' else computedFlag' .&. complement flagProperlyPaired
        !computedIsize   = if properly_paired then pos5 b - pos5 a else 0

        reduce f | f .&. flagMateUnmapped == 0 = f .&. complement flagFailsQC
                 | otherwise = f .&. complement (flagFailsQC .|. flagMateReversed)

        f1 === f2 = reduce f1 == reduce f2
        f1 =!= f2 = f1 /= f2 && f1 === f2

        onlyIf f m = (\z -> if z then m else "") `fmap` tells f

        count_mrnm  = do modify $ \c -> c { num_mrnm = 1 + num_mrnm c }
                         let ra = unRefseq (b_mrnm a); rb = unRefseq (b_rname b)
                         onlyIf report_mrnm $ printf "MRNM %d is wrong (%d)" ra rb

        count_mpos  = do modify $ \c -> c { num_mpos = 1 + num_mpos c }
                         onlyIf report_mpos $ printf "MPOS %d is wrong (%d)" (b_mpos a) (b_pos b)

        count_isize = do modify $ \c -> c { num_isize = 1 + num_isize c }
                         onlyIf report_isize $ printf "ISIZE %d is wrong (%d)" (b_isize a) computedIsize

        count_flags = do modify $ \c -> c { num_flags = 1 + num_flags c }
                         onlyIf report_flags $ printf "FLAG %03X is wrong (+%03X,-%03X)" (b_flag a) fp fm

        count_fflag = do modify $ \c -> c { num_fflag = 1 + num_fflag c }
                         onlyIf report_fflag $ printf "FLAG %03X is technically wrong (+%03X,-%03X)" (b_flag a) fp fm

        count_ixs = do modify $ \c -> c { num_ixs = 1 + num_ixs c }
                       onlyIf report_ixs $ printf "Index fields %s are wrong (%s)" (show $ b_indices a) (show common_indices)

        fp = computedFlag .&. complement (b_flag a)
        fm = complement computedFlag .&. b_flag a

        index_fields = [ "XI", "XJ", "YI", "YJ", "RG", "BC" ]
        b_indices x = [ extAsString key x | key <- index_fields ]
        common_indices = zipWith max (b_indices a) (b_indices b)

        setIndices is x = x { b_exts = add_new . remove_old $ b_exts x }
          where
            add_new y = foldr (:) y $ zip index_fields $ map Text is
            remove_old y = foldr deleteE y index_fields

-- | Turns a widow into a single.  Basically removes the pairing
-- related flags and clear the information concerning the mate.
divorce :: BamRec -> BamRec
divorce b = b { b_flag = b_flag b .&. complement pair_flags
              , b_mrnm = invalidRefseq
              , b_mpos = invalidPos
              , b_isize = 0 }
  where
    pair_flags = flagPaired .|. flagProperlyPaired .|.
                 flagFirstMate .|. flagSecondMate .|.
                 flagMateUnmapped .|. flagMateReversed

-- I think this can work with priority queues alone:
--
-- - One contains incomplete pairs ordered by mate position.  When we
--   reach a given position and find the 2nd mate, the minimum in this
--   queue must be the 1st mate (or another 1st mate matching another
--   read we'll find here).
--
-- - One contains incomplete pairs ordered by (hash of) qname.  This one
--   is only used if we missed a mate for some reason.  After we read
--   the whole input, all remaining pairs can be pulled off this queue
--   in order of increasing (hash of) qname.
--
-- - At any given position, we will have a number of 1st mates that have
--   been waiting in the queue and a number of 2nd mates that are coming
--   in from the input.  We dump both sets into a queue by qname, then
--   pull them out in pairs.  Stuff that comes off as anything else than
--   a pair gets queued up again.

data MatingStats = MS { total_in   :: !Int
                      , total_out  :: !Int
                      , singletons :: !Int
                      , widows     :: !Int
                      , num_mrnm :: !Int
                      , num_mpos :: !Int
                      , num_isize :: !Int
                      , num_flags :: !Int
                      , num_fflag :: !Int
                      , num_ixs :: !Int }

report_stats :: MatingStats -> String
report_stats ms = unlines [
    "number of records read:          " ++ showNum (total_in ms),
    "number of records written:       " ++ showNum (total_out ms),
    "number of true singletons:       " ++ showNum (singletons ms),
    "number of widows:                " ++ showNum (widows ms),
    "number of repaired MRNM values:  " ++ showNum (num_mrnm ms),
    "number of repaired MPOS values:  " ++ showNum (num_mpos ms),
    "number of repaired ISIZE values: " ++ showNum (num_isize ms),
    "number of repaired FLAGS values: " ++ showNum (num_flags ms),
    "number of common FLAGS problems: " ++ showNum (num_fflag ms),
    "number of index field problems:  " ++ showNum (num_ixs ms) ]

data Queues = QS { right_here :: !(PQ ByQName)
                 , in_order   :: !(PQ ByMatePos)
                 , messed_up  :: !(PQ ByQName) }

type Lens s a = forall f. Functor f => (a -> f a) -> s -> f s

right_here' :: Lens Queues (PQ ByQName)
right_here' f q = (\p -> q { right_here = p }) <$> f (right_here q)

in_order' :: Lens Queues (PQ ByMatePos)
in_order'   f q = (\p -> q { in_order   = p }) <$> f (in_order   q)

messed_up' :: Lens Queues (PQ ByQName)
messed_up'  f q = (\p -> q { messed_up  = p }) <$> f (messed_up  q)

ms0 :: MatingStats
ms0 = MS 0 0 0 0 0 0 0 0 0 0

getSize :: (Ord a, Sized a) => (Queues -> PQ a) -> Mating r m Int
getSize sel = getq $ sizePQ . sel

enqueue :: (MonadIO m, Ord a, Sized a) => a -> Lens Queues (PQ a) -> Mating r m ()
enqueue a sel = Mating $ \k s o q c r -> liftIO (sel (enqueuePQ (pqconf c) a) q) >>= \q' -> k () s o q' c r

peekMin :: (Ord a, Sized a) => (Queues -> PQ a) -> Mating r m (Maybe a)
peekMin sel = getq $ fmap fst . viewMinPQ . sel

fetchMin :: (Ord a, Sized a) => Lens Queues (PQ a) -> Mating r m (Maybe a)
fetchMin sel = modq . sel $ \q -> case viewMinPQ q of Nothing     -> (Nothing, q)
                                                      Just (a,q') -> (Just  a, q')

discardMin :: (Ord a, Sized a) => Lens Queues (PQ a) -> Mating r m ()
discardMin sel = () <$ fetchMin sel


note, warn, err :: MonadIO m => String -> Mating r m ()
note msg = do v <- tells verbosity ; unless (v < Notices)  $ liftIO $ hPutStrLn stderr $ "[fixpair] info:    " ++ msg
warn msg = do v <- tells verbosity ; unless (v < Warnings) $ liftIO $ hPutStrLn stderr $ "[fixpair] warning: " ++ msg
err  msg = do v <- tells verbosity ; unless (v < Errors)   $ liftIO $ hPutStrLn stderr $ "[fixpair] error:   " ++ msg

report' :: MonadIO m => Mating r m ()
report' = do o <- gets total_out
             when (o `mod` 0x40000 == 0) $ do
                     ms <- getSize messed_up
                     note $ printf "out: %d, mess: %d" o ms

report :: MonadIO m => BamRaw -> Mating r m ()
report br = do i <- gets total_in
               o <- gets total_out
               when (i `mod` 0x100000 == 0) $ do
                     hs <- getSize right_here
                     os <- getSize in_order
                     ms <- getSize messed_up
                     rr <- getRefseqs
                     let BamRec{..} = unpackBam br
                         rn = unpack . sq_name $ getRef rr b_rname
                         at = if b_rname == invalidRefseq || b_pos == invalidPos
                              then "" else printf "@%s/%d,\t" rn b_pos
                         off = fromIntegral $ b_virtual_offset `shiftR` 36
                     note $ printf "%+5dMB, %sin: %9d, out: %9d, here: %6d, wait: %6d, mess: %6d"
                                   (off::Int) (at::String) i o hs os ms

no_mate_here :: MonadIO m => String -> BamRaw -> Mating r m ()
no_mate_here l br = do note $ let b = unpackBam br
                              in "[" ++ l ++ "] record "
                                 ++ shows (b_qname b) (if isFirstMate b then "/1" else "/2")
                                 ++ " did not have a mate at the right location."
                       let !br' = br_copy br
                       enqueue (byQName br') messed_up'

no_mate_ever :: MonadIO m => BamRaw -> Mating r m ()
no_mate_ever b = do let b' = unpackBam b
                    err $ "record " ++ shows (b_qname b') " (" ++
                          shows (extAsInt 1 "XI" b') ") did not have a mate at all."
                    modify $ \c -> c { widows = 1 + widows c }
                    kill <- tells killmode
                    case kill of
                        KillAll  -> return ()
                        KillUu   -> unless (isUnmapped b') $ yield [divorce b']
                        KillNone -> yield [divorce b']

-- Basically the CPS version of the State Monad.  CPS is necessary to be
-- able to call 'eneeCheckIfDone' in the middle, and that fixes the
-- underlying monad to an 'Iteratee' and the ultimate return type to an
-- 'Iteratee', too.  Pretty to work with, not pretty to look at.
type Sink r m = Stream [BamRec] -> Iteratee [BamRec] m r
newtype Mating r m a = Mating { runMating ::
    (a -> MatingStats -> Sink r m -> Queues -> Config -> Refs -> Iteratee [BamPair] m (Iteratee [BamRec] m r))
       -> MatingStats -> Sink r m -> Queues -> Config -> Refs -> Iteratee [BamPair] m (Iteratee [BamRec] m r) }

instance Functor (Mating r m) where
    fmap f m = Mating $ \k -> runMating m (k . f)

instance Applicative (Mating r m) where
    pure a = Mating $ \k -> k a
    u <*> v = Mating $ \k -> runMating u (\a -> runMating v (k . a))

instance Monad (Mating r m) where
    return a = Mating $ \k -> k a
    m >>=  k = Mating $ \k2 -> runMating m (\a -> runMating (k a) k2)

instance MonadIO m => MonadIO (Mating r m) where
    liftIO f = Mating $ \k s o q c r -> liftIO f >>= \a -> k a s o q c r


lift'it :: Monad m => Iteratee [BamPair] m a -> Mating r m a
lift'it m = Mating $ \k s o q c r -> m >>= \a -> k a s o q c r

tells :: (Config -> a) -> Mating r m a
tells f = Mating $ \k s o q c -> k (f c) s o q c

gets :: (MatingStats -> a) -> Mating r m a
gets f = Mating $ \k s -> k (f s) s

getq :: (Queues -> a) -> Mating r m a
getq f = Mating $ \k s o q -> k (f q) s o q

modq :: (Queues -> (b, Queues)) -> Mating r m b
modq f = Mating $ \k s o q -> case f q of (b, q') -> k b s o q'

modify :: (MatingStats -> MatingStats) -> Mating r m ()
modify f = Mating $ \k s -> (k () $! f s)

getRefseqs :: Mating r m Refs
getRefseqs = Mating $ \k s o q c r -> k r s o q c r

fetchNext :: MonadIO m => Mating r m (Maybe BamPair)
fetchNext = do r <- lift'it tryHead
               case r of Nothing -> return ()
                         Just (Singleton x) -> do modify $ \s -> s { total_in = 1 + total_in s } ; report x
                         Just (Pair    _ x) -> do modify $ \s -> s { total_in = 2 + total_in s } ; report x
                         Just (LoneMate  x) -> do modify $ \s -> s { total_in = 1 + total_in s } ; report x
               return r

yield :: MonadIO m => [BamRec] -> Mating r m ()
yield rs = Mating $ \k s o q c r -> let !s' = s { total_out = length rs + total_out s }
                                    in eneeCheckIfDone (\o' -> k () s' o' q c r) . o $ Chunk rs

-- To ensure proper cleanup, we require the priority queues to be created
-- outside.  Since one is continually reused, it is important that a PQ
-- that is emptied no longer holds on to files on disk.
re_pair :: MonadIO m => Config -> Refs -> Enumeratee [BamPair] [BamRec] m a
re_pair cf rs = eneeCheckIfDone $ \out -> runMating go finish ms0 out (QS makePQ makePQ makePQ) cf' rs
   where
    cf' = cf { pqconf = (pqconf cf) { croak = unless (verbosity cf < Notices) . hPutStrLn stderr . (++) "[fixpair] info:    " } }

    go = fetchNext >>= go'

    -- At EOF, flush everything.
    go' Nothing = peekMin right_here >>= \mm -> case mm of
            Just (ByQName _ _ qq) -> do complete_here (br_self_pos qq)
                                        flush_here Nothing  -- flush_here loops back here
            Nothing               -> flush_in_order  -- this ends the whole operation

    -- Single read?  Pass through and go on.
    -- Paired read?  Does it belong 'here'?
    go' (Just (Singleton x)) = modify (\c -> c { singletons = 1 + singletons c }) >> yield [unpackBam x] >> go
    go' (Just (Pair    x y)) = fixmate x y >>= yield >> go
    go' (Just (LoneMate  r)) = peekMin right_here >>= \mm -> case mm of

            -- there's nothing else here, so here becomes redefined
            Nothing             -> enqueueThis r >> go

            Just (ByQName _ _ qq) -> case compare (br_self_pos r) (br_self_pos qq) of
                -- nope, r is out of order and goes to 'messed_up'
                LT -> do warn $ "record " ++ show (br_qname r) ++ " is out of order."
                         let !r' = br_copy r
                         enqueue (byQName r') messed_up'
                         go

                -- nope, r comes later.  we need to finish our business here
                GT -> do complete_here (br_self_pos qq)
                         flush_here (Just (LoneMate r))

                -- it belongs here or there is nothing else here
                EQ -> enqueueThis r >> go


    -- lonely guy, belongs either here or needs to wait for the mate
    enqueueThis r | br_self_pos r >= br_mate_pos r = enqueue (byQName r) right_here'
                  | otherwise             = r' `seq` enqueue (ByMatePos r') in_order'
        where r' = br_copy r

    -- Flush the in_order queue to messed_up, since those didn't find
    -- their mate the ordinary way.  Afterwards, flush the messed_up
    -- queue.
    flush_in_order = fetchMin in_order' >>= \zz -> case zz of
        Just (ByMatePos b) -> no_mate_here "flush_in_order" b >> flush_in_order
        Nothing            -> flush_messed_up

    -- Flush the messed up queue.  Everything should come off in pairs,
    -- unless something is broken.
    flush_messed_up = fetchMin messed_up' >>= flush_mess1

    flush_mess1 Nothing                 = return ()
    flush_mess1 (Just (ByQName _ ai a)) = fetchMin messed_up' >>= flush_mess2 ai a

    flush_mess2  _ a Nothing = no_mate_ever a

    flush_mess2 ai a b'@(Just (ByQName _ bi b))
        | ai /= bi || br_qname a /= br_qname b = no_mate_ever a >> report' >> flush_mess1 b'
        | otherwise                            = fixmate a b    >>= yield >> report' >> flush_messed_up


    -- Flush the right_here queue.  Everything should come off in pairs,
    -- if not, it goes to messed_up.  When done, loop back to 'go'
    flush_here  r = fetchMin right_here' >>= flush_here1 r

    flush_here1 r Nothing = go' r
    flush_here1 r (Just a) = fetchMin right_here' >>= flush_here2 r a

    flush_here2 r (ByQName _ _ a) Nothing = do no_mate_here "flush_here2/Nothing" a
                                               flush_here r

    flush_here2 r (ByQName _ ai a) b'@(Just (ByQName _ bi b))
        | ai /= bi || br_qname a /= br_qname b = no_mate_here "flush_here2/Just" a >> flush_here1 r b'
        | otherwise                            = fixmate a b >>= yield >> flush_here r


    -- add stuff coming from 'in_order' to 'right_here'
    complete_here pivot = do
            zz <- peekMin in_order
            case zz of
                Nothing -> return ()
                Just (ByMatePos b)
                       | pivot  > br_mate_pos b -> do discardMin in_order'
                                                      no_mate_here "complete_here" b
                                                      complete_here pivot

                       | pivot == br_mate_pos b -> do discardMin in_order'
                                                      enqueue (byQName b) right_here'
                                                      complete_here pivot

                       | otherwise -> return ()

    finish () st o (QS x y z) _cf _rs = liftIO $ do
        closePQ x ; closePQ y ; closePQ z
        hPutStrLn stderr $ report_stats st
        return (liftI o)


-- | To catch pairs whose mates are adjacent (either because the file
-- has never been sorted or because it has been group-sorted), we apply
-- preprocessing.  The idea is that if we can catch these pairs early,
-- the priority queues never fill up and we save a ton of processing.
-- Now to make the re-pair algorithm work well, we need to merge-sort
-- inputs.  But after that, the pairs have been separated.  So we apply
-- the preprocessing to each input file, then merge then, then run
-- re-pair.

data BamPair = Singleton BamRaw | Pair BamRaw BamRaw | LoneMate BamRaw


mergeInputs :: (MonadIO m, MonadMask m) => [FilePath] -> Enumerator' BamMeta [BamPair] m a
mergeInputs = go0
  where
    go0 [        ] = enumG $ enumFd defaultBufSize stdInput
    go0 (fp0:fps0) = go fp0 fps0

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

    enum1 "-" = enumG $ enumFd   defaultBufSize stdInput
    enum1  fp = enumG $ enumFile defaultBufSize    fp

    enumG ee k = ee >=> run $ joinI $ decodeAnyBam $ \h -> quick_pair (k h)


quick_pair :: Monad m => Enumeratee [BamRaw] [BamPair] m a
quick_pair = eneeCheckIfDone go0
  where
    go0 k = tryHead >>= maybe (return $ liftI k) (\x -> go1 x k)

    go1 x k | not (isPaired (unpackBam x)) = eneeCheckIfDone go0 . k $ Chunk [Singleton x]
            | otherwise                    = tryHead >>= maybe (return . k $ Chunk [LoneMate x]) (\y -> go2 x y k)

    go2 x y k | b_qname (unpackBam x) == b_qname (unpackBam y) = eneeCheckIfDone go0 . k $ Chunk [Pair x y]
              | otherwise                                      = eneeCheckIfDone (go1 y) . k $ Chunk [LoneMate x]


combineCoordinates :: Monad m => BamMeta -> Enumeratee [BamPair] [BamPair] (Iteratee [BamPair] m) a
combineCoordinates _ = mergeSortStreams (?)
  where u ? v = if (bp_rname u, bp_pos u) < (bp_rname v, bp_pos v) then Less else NotLess

bp_rname :: BamPair -> Refseq
bp_rname (Singleton u) = b_rname $ unpackBam u
bp_rname (Pair    u _) = b_rname $ unpackBam u
bp_rname (LoneMate  u) = b_rname $ unpackBam u

bp_pos :: BamPair -> Int
bp_pos (Singleton u) = b_pos $ unpackBam u
bp_pos (Pair    u _) = b_pos $ unpackBam u
bp_pos (LoneMate  u) = b_pos $ unpackBam u


do_trim :: Int -> [BamRec] -> [BamRec]
do_trim q = scan_empties . map trim1
  where
    trim1 b = case [ l | l <- [0 .. V.length (b_qual b) -1], avquallow (V.drop l qs) ] of
                [ ] -> b
                l:_ -> trim_3 l b
      where
        qs | isReversed b = V.reverse (b_qual b)
           | otherwise    =            b_qual b

    scan_empties (x:y:z)
        | b_qname x == b_qname y
            = if V.null (b_qual x) || V.null (b_qual y)
                then scan_empties z
                else x : y : scan_empties z

    scan_empties (x:z)
        = if V.null (b_qual x)
           then scan_empties z
           else x : scan_empties z

    scan_empties [] = []

    avquallow vec = V.sum (V.map (fromIntegral . unQ) vec) <= q * V.length vec

