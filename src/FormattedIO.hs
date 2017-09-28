{-# LANGUAGE LambdaCase #-}
module FormattedIO where

import Bio.Bam                      hiding ( Region(..) )
import Bio.Prelude
import System.FilePath
import System.IO

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

{- make_summary :: [( FilePath, Iteratee [(t, AnnoSet)] m (HashMap Bytes Int)
                      -> IO (Iteratee s IO (HashMap Bytes Integer)))] -> IO ()
make_summary fs = do
    mm <- forM fs $ \(_,enum) -> enum >=> run $ foldStream count_annos M.empty

    -- let mm = [ foldl' count_annos M.empty annos
             -- | ( _, annos ) <- annotated_files ]

    -- XXX I should convert to set, shouldn't I?
    let keymap = M.unions $ map (M.map (const ())) mm

    putStr "ID"
    forM_ fs $ \(f,_) -> do putChar '\t' ; putStr f
    putChar '\n'

    sequence_ [ do S.putStr anno
                   forM_ mm $ \counts -> do let v = M.lookupDefault 0 anno counts
                                            putChar '\t' ; putStr (show v)
                   putChar '\n'
              | ( anno, () ) <- M.toList keymap ] -}

    {- make_summary :: [( FilePath, [Annotation] )] -> IO ()
    make_summary annotated_files = do
        let mm = [ foldl' count_annos M.empty annos
                 | ( _, annos ) <- annotated_files ]

        let keymap = M.unions $ map (M.map (const ())) mm

        putStr "ID"
        sequence_ [ do putChar '\t' ; putStr f | (f,_) <- annotated_files ]
        putChar '\n'

        sequence_ [ do S.putStr anno
                       sequence_ [ do let v = M.lookupDefault 0 anno counts
                                      putChar '\t' ; putStr (show v)
                                 | counts <- mm ]
                       putChar '\n'
                  | ( anno, () ) <- M.toList keymap ] -}

{- write_file_pattern :: FilePath -> [( FilePath, [(L.ByteString,AnnoSet)] )] -> IO ()
write_file_pattern pat = mapM_ go
  where
    go (fname, recs) = withFile (make_oname pat) WriteMode $ \h -> make_table h recs
      where
        (qname, ext) = splitExtension fname
        (path, name) = splitFileName qname

        make_oname ('%':'%':cs) = '%' : make_oname cs
        make_oname ('%':'s':cs) = qname ++ make_oname cs
        make_oname ('%':'p':cs) = path ++ make_oname cs
        make_oname ('%':'n':cs) = name ++ make_oname cs
        make_oname ('%':'e':cs) = ext ++ make_oname cs
        make_oname (      c:cs) = c : make_oname cs
        make_oname [] = [] -}

readGeneTable :: (String -> IO ()) -> [Region] -> CPS r RevSymtab
readGeneTable report inp = do
    add_all inp (1::Int)
    symtab <- get_syms
    return $! invertTab symtab
  where
    report' n c = io . report $ "reading annotations: " ++ shows n "\27[K" ++ [c]

    add_all [] n = report' n '\n'
    add_all ((gi, crm, strands, s, e):xs) n = do
        when (n `mod` 1024 == 0) (report' n '\r')
        val <- findSymbol gi
        sequence_ $ withSenses strands $ \str ->
                io . addI (min s e) (max s e) (fromIntegral val)
                    =<< findDiet (str crm)
        add_all xs $! n+1


readChromTable :: L.ByteString -> ChromTable
readChromTable = foldl' add M.empty . map L.words . L.lines
  where
    add m (x:_) | L.head x == '#' = m
    add m (ens:ucsc:ofs:_) =
        let k = L.fromChunks [ S.copy (L.toStrict  ens) ]
            x = L.fromChunks [ S.copy (L.toStrict ucsc) ]
            y = ereadInt ofs
        in x `seq` y `seq` M.insert k (x,y) m
    add _ (ens:_) = error $ "too few fields in translation of chromosome " ++ show (L.unpack ens)
    add m [     ] = m

readMySam :: [L.ByteString] -> [Region]
readMySam = xlate . map (L.split '\t')
  where
    xlate ((qname : flag0 : rname : pos0 : _mapq : cigar : _) : rs)
        | flag .&. 4 == 0 = ( qname, rname, str, pos-1, pos+len-1 ) : xlate rs
      where
        flag = ereadInt flag0
        pos = ereadInt pos0
        len = fromMaybe (error "parse error in CIGAR field") (textual_cigar_to_len cigar)
        str = if flag .&. 0x10 == 0 then Forward else Reverse

    xlate (_:rs) = xlate rs
    xlate [] = []

textual_cigar_to_len :: L.ByteString -> Maybe Int
textual_cigar_to_len s0 = go s0 0
  where
    go s acc | L.null s = return acc
    go s acc = do (n,s') <- L.readInt s
                  (op,t) <- L.uncons s'
                  go t $! acc + w n op

    w n 'M' = n ; w n 'D' = n ; w n 'N' = n ; w _ _ = 0

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
                       | c == c2w '@' -> joinI $ decodeSam $ flip readMyBam it  -- probably sam
                -- XXX  take printable (!) chars until line feed
                -- stopped at other than line feed? -> junk
                -- first line has >11 fields and #5 llooks like a cigar? -> -- Sam?
                       | otherwise    -> readMyBed it                           -- Bed?

{- case dropWhile (L.all isSpace) (L.lines s) of
    -- _        | isBam s                       -> readMyBam s
    ls@(l:_) | L.head l == '>'               -> readMyFasta ls
             | may_well_be_sam l             -> readMySam ls
    ls                                       -> readMyBed ls
  where
    may_well_be_sam l = let ws = L.split '\t' l
                        in L.head l == '@' || (length ws >= 11 && isJust (textual_cigar_to_len (ws !! 5))) -}


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
        go (crm:s:e:gi:strand) = [(L.fromStrict gi, L.fromStrict crm, read_strand (drop 1 strand), ereadIntS s, ereadIntS e)]
        go ws = error $ "not enough fields in BED file" ++ case ws of
                    [w] | length (S.words w) >= 4 -> " (whitespace instead of tabs?)"
                    _                             -> ""

read5cFile :: [L.ByteString] -> [Region]
read5cFile = foldr go [] . map L.words
  where go [] = id
        go (w:_) | L.head w == '#' = id
        go (gi:crm:s:e:strand) = (:) (gi, crm, read_strand (map L.toStrict strand), ereadInt s, ereadInt e)
        go _ = error $ "too few fields in legacy annotation file"

read_strand :: [Bytes] -> Senses
read_strand [   ]                   = Both
read_strand (s:_) | S.null s        = Both
                  | S.head s == '-' = Reverse
                  | s == "0"        = Both
                  | otherwise       = Forward

xlate_chrom :: Maybe ChromTable -> Region -> Region
xlate_chrom Nothing reg = reg                                       -- no table, no translation
xlate_chrom (Just ct) (n, crm, ss, s, e) = case M.lookup crm ct of
        Just (c,o) -> (n, c, ss, s+o-1, e+o)                        -- hit, so translate
        Nothing -> error $ "chromosome unknown: " ++ L.unpack crm   -- translate failed


ereadInt :: L.ByteString -> Int
ereadInt s = case L.readInt s of Just (x,_) -> x
                                 Nothing -> error $ "couldn't parse Int at " ++ L.unpack s

ereadIntS :: Bytes -> Int
ereadIntS s = case S.readInt s of Just (x,_) -> x
                                  Nothing -> error $ "couldn't parse Int at " ++ S.unpack s


myReadFilesWith :: ( FilePath -> L.ByteString -> IO a ) -> [ FilePath ] -> IO [a]
myReadFilesWith k [] = fmap (:[]) $ L.getContents >>= k "stdin"
myReadFilesWith k fs = forM fs $ \f -> if f == "-" then L.getContents >>= k "stdin" else L.readFile f >>= k f

