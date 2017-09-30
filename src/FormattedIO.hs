{-# LANGUAGE LambdaCase #-}
module FormattedIO where

import Bio.Bam                   hiding ( Region(..) )
import Bio.Prelude               hiding ( bracket )
import Control.Monad.Catch              ( bracket )
import Data.Functor.Identity
import System.FilePath
import System.IO

import qualified Data.ByteString                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.HashMap.Strict                as M

import Diet ( addI )
import Symtab

data AnnoSet = Hits [Bytes] | NearMiss Bytes Int deriving Show
type Annotation = ( L.ByteString, AnnoSet )

print_header :: Handle -> IO ()
print_header hdl = hPutStr hdl "#RegionName\tID\n"

make_table :: Handle -> Iteratee [Annotation] IO ()
make_table hdl = mapStreamM_ format
  where
    format (_, Hits []) = return ()
    format (n, Hits hs) = do L.hPut hdl n ; hPutChar hdl '\t'
                             S.hPut hdl (head hs)
                             forM_ (tail hs) $ \h -> hPutChar hdl ':' >> S.hPut hdl h
                             hPutChar hdl '\n'
    format (n, NearMiss a d) = do L.hPut hdl n ; hPutChar hdl '\t'
                                  S.hPut hdl a ; hPutChar hdl '@'
                                  hPutStrLn hdl (show d)

make_summary :: a -> Iteratee [(t, AnnoSet)] IO (a, HashMap Bytes Int)
make_summary f = (,) f <$> foldStream count_annos M.empty
  where
    count_annos m (_, Hits hs) = foldl' count1 m hs
    count_annos m _ = m
    count1 m h = (M.insert h $! 1 + M.lookupDefault (0::Int) h m) m

print_summary :: [(FilePath, HashMap Bytes Int)] -> IO ()
print_summary mm = do
    let keymap = M.unions $ map (M.map (const ()) . snd) mm

    putStr "ID"
    forM_ mm $ \(f,_) -> do putChar '\t' ; putStr f
    putChar '\n'

    sequence_ [ do S.putStr anno
                   forM_ mm $ \(_,counts) -> do let v = M.lookupDefault 0 anno counts
                                                putChar '\t' ; putStr (show v)
                   putChar '\n'
              | ( anno, () ) <- M.toList keymap ]

write_file_pattern :: FilePath -> FilePath -> Iteratee [Annotation] IO ()
write_file_pattern pat fname =
    bracket (liftIO $ openFile (make_oname pat) WriteMode) (liftIO . hClose)
            (\hdl -> liftIO (print_header hdl) >> make_table hdl)
  where
    (qname, ext) = splitExtension fname
    (path, name) = splitFileName qname

    make_oname ('%':'%':cs) = '%' : make_oname cs
    make_oname ('%':'s':cs) = qname ++ make_oname cs
    make_oname ('%':'p':cs) = path ++ make_oname cs
    make_oname ('%':'n':cs) = name ++ make_oname cs
    make_oname ('%':'e':cs) = ext ++ make_oname cs
    make_oname (      c:cs) = c : make_oname cs
    make_oname [          ] = [ ]


readGeneTable :: Iteratee [Region] CPS RevSymtab
readGeneTable = do
    mapStreamM_ $ \(gi, crm, strands, s, e) -> do
        val <- findSymbol (L.toStrict gi)
        sequence_ $ withSenses strands $ \str ->
            findDiet (str crm) >>=
                liftIO . addI (min s e) (max s e) (fromIntegral val)
    invertTab <$> lift get_syms

readChromTable :: L.ByteString -> ChromTable
readChromTable = foldl' add M.empty . map L.words . L.lines
  where
    add m (x:_) | L.head x == '#' = m
    add m (ens:ucsc:ofs:_) =
        let k = L.fromChunks [ S.copy (L.toStrict  ens) ]
            x = L.fromChunks [ S.copy (L.toStrict ucsc) ]
            y = ereadInt (L.toStrict ofs)
        in x `seq` y `seq` M.insert k (x,y) m
    add _ (ens:_) = error $ "too few fields in translation of chromosome " ++ show (L.unpack ens)
    add m [     ] = m

readMyBam :: Monad m => BamMeta -> Enumeratee [BamRec] [Region] m b
readMyBam hdr = filterStream (not . isUnmapped) ><> mapStream xlate
  where
    getref r = sq_name $ getRef (meta_refs hdr) r
    str r = if b_flag r .&. 0x10 == 0 then Forward else Reverse

    xlate r = ( L.fromStrict $ b_qname r
              , L.fromStrict $ getref (b_rname r)
              , str r
              , b_pos r
              , b_pos r + alignedLength (b_cigar r) )

-- 'MonadIO' is needed to support gzip.
readInput :: MonadIO m => Enumeratee Bytes [Region] m b
readInput it = do
    isBam >>= \case
        Just  d -> joinI $ d $ \h -> mapStream unpackBam =$ readMyBam h it      -- sort-of Bam
        Nothing -> do
            dropWhileStreamBS (\c -> c == 32 || c == 9 || c == 10 || c == 13 || c == 160)
            peekStreamBS >>= \case
                Nothing               -> return it                              -- empty?!
                Just c | c == c2w '>' -> readMyFasta it                         -- Fasta?
                       | c == c2w '@' -> joinI $ decodeSam $ flip readMyBam it  -- probably Sam
                       | otherwise    -> do
                         (l1,c1) <- iLookAhead $ (,) <$> takeWhileBS isPrinting <*> peekStreamBS
                         case trySam l1 of
                            _ | c1 /= Just 10 && c1 /= Just 13 && c1 /= Nothing ->
                                    error "this doesn't look like either Bed or Sam"
                            Right _   -> joinI $ decodeSam $ flip readMyBam it  -- headerless Sam
                            Left  e   -> (e::SomeException) `seq` readMyBed it  -- whatever, probably Bed
  where
    isPrinting c = (c >= 32 && c < 127) || c == 9
    trySam l1 = runIdentity $ enumPure1Chunk l1 >=> tryRun $ decodeSam' noRefs =$ stream2stream

takeWhileBS :: (Word8 -> Bool) -> Iteratee Bytes m Bytes
takeWhileBS p = icont (step mempty) Nothing
  where
    step bfr (Chunk str)
      | B.null str     = icont (step bfr) Nothing
      | otherwise      = case B.span p str of
        (hd, tl)
          | B.null tl -> icont (step (mappend bfr str)) Nothing
          | otherwise -> idone (mappend bfr hd) (Chunk tl)
    step bfr stream    = idone bfr stream


readMyFasta :: Monad m => Enumeratee Bytes [Region] m b
readMyFasta = enumLinesBS ><> filterStream (S.isPrefixOf ">") ><> concatMapStream proc_line
  where
    proc_line = map parse_label . S.split ';' . S.drop 1
    parse_label label = (L.fromStrict label, L.fromStrict crm, read_strand [str], lft, rht)
        where
            (crm, l1) = S.break (':'==) label
            (str, l2) = S.break (':'==) $ S.drop 1 l1
            (lft, l3) = fromMaybe (error "parse error in FASTA header") $ S.readInt $ S.drop 1 l2
            (rht,  _) = fromMaybe (error "parse error in FASTA header") $ S.readInt $ S.drop 1 l3

readMyBed :: Monad m => Enumeratee Bytes [Region] m b
readMyBed = enumLinesBS ><> concatMapStream (go . S.split '\t')
  where go [] = []
        go (w:_) | w == "track" || S.head w == '#' = []
        go (crm:s:e:gi:strand) = [(L.fromStrict gi, L.fromStrict crm, read_strand (drop 1 strand), ereadInt s, ereadInt e)]
        go ws = error $ "not enough fields in BED file" ++ case ws of
                    [w] | length (S.words w) >= 4 -> " (whitespace instead of tabs?)"
                    _                             -> ""

read5col :: Monad m => Maybe ChromTable -> Enumeratee Bytes [Region] m b
read5col mct = enumLinesBS ><> concatMapStream (maybe id (map . xlate_chrom) mct . go . S.words)
  where go [   ]                   = []
        go (w:_) | S.head w == '#' = []
        go (gi:crm:s:e:strand)     = [(L.fromStrict gi, L.fromStrict crm, read_strand strand, ereadInt s, ereadInt e)]
        go _ = error $ "too few fields in legacy annotation file"

read_strand :: [Bytes] -> Senses
read_strand [   ]                   = Both
read_strand (s:_) | S.null s        = Both
                  | S.head s == '-' = Reverse
                  | s == "0"        = Both
                  | otherwise       = Forward

xlate_chrom :: ChromTable -> Region -> Region
xlate_chrom ct (n, crm, ss, s, e) = case M.lookup crm ct of
        Just (c,o) -> (n, c, ss, s+o-1, e+o)                        -- hit, so translate
        Nothing -> error $ "chromosome unknown: " ++ L.unpack crm   -- translate failed

ereadInt :: Bytes -> Int
ereadInt s = case S.readInt s of Just (x,_) -> x
                                 Nothing -> error $ "couldn't parse Int at " ++ S.unpack s

myReadFilesWith :: ( FilePath -> L.ByteString -> IO a ) -> [ FilePath ] -> IO [a]
myReadFilesWith k [] = fmap (:[]) $ L.getContents >>= k "stdin"
myReadFilesWith k fs = forM fs $ \f -> if f == "-" then L.getContents >>= k "stdin" else L.readFile f >>= k f

