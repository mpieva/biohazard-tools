module FormattedIO where

import Bio.Prelude
import System.FilePath
import System.IO

import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.Map                           as M

import Diet
import Symtab

data AnnoSet = Hits [S.ByteString] | NearMiss S.ByteString Int deriving Show
type Annotation = ( L.ByteString, AnnoSet )

make_table :: Handle -> [Annotation] -> IO ()
make_table hdl anns = do hPutStr hdl "#RegionName\tID\n" ; mapM_ format anns
  where
    format (_, Hits []) = return ()
    format (n, Hits hs) = do L.hPut hdl n ; hPutChar hdl '\t'
                             S.hPut hdl (head hs)
                             forM_ (tail hs) $ \h -> hPutChar hdl ':' >> S.hPut hdl h
                             hPutChar hdl '\n'
    format (n, NearMiss a d) = do L.hPut hdl n ; hPutChar hdl '\t'
                                  S.hPut hdl a ; hPutChar hdl '@'
                                  hPutStrLn hdl (show d)

make_one_table :: Handle -> [( FilePath, [Annotation] )] -> IO ()
make_one_table hdl = make_table hdl . concat . map snd

make_summary :: [( FilePath, [Annotation] )] -> IO ()
make_summary annotated_files = do
    let mm = [ foldl' count_annos M.empty annos
             | ( _, annos ) <- annotated_files ]

    let keymap = M.unions $ map (M.map (const ())) mm

    putStr "ID"
    sequence_ [ do putChar '\t' ; putStr f | (f,_) <- annotated_files ]
    putChar '\n'

    sequence_ [ do S.putStr anno
                   sequence_ [ do let v = M.findWithDefault 0 anno counts
                                  putChar '\t' ; putStr (show v)
                             | counts <- mm ]
                   putChar '\n'
              | ( anno, () ) <- M.toList keymap ]
  where
    count_annos m (_, Hits hs) = foldl' count1 m hs
    count_annos m _ = m
    count1 m h = (M.insert h $! 1 + M.findWithDefault (0::Int) h m) m

write_file_pattern :: FilePath -> [( FilePath, [(L.ByteString,AnnoSet)] )] -> IO ()
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
        make_oname [] = []

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
        forM_ strands $ \str ->
                io . addI (min s e) (max s e) (fromIntegral val) =<< findDiet (crm,str)
        add_all xs $! n+1


readChromTable :: L.ByteString -> ChromTable
readChromTable = foldl' add M.empty . map L.words . L.lines
  where
    add m (x:_) | L.head x == '#' = m
    add m (ens:ucsc:ofs:_) =
        let x = ucsc -- XXX L.fromChunks [shelve ucsc]
            y = ereadInt ofs
        in x `seq` y `seq` M.insert (ens {-L.fromChunks [shelve ens]-}) (x,y) m
    add _ (ens:_) = error $ "too few fields in translation of chromosome " ++ show (L.unpack ens)
    add m [     ] = m

readMySam :: [L.ByteString] -> [Region]
readMySam = xlate . map (L.split '\t')
  where
    xlate ((qname : flag0 : rname : pos0 : _mapq : cigar : _) : rs)
        | flag .&. 4 == 0 = ( qname, rname, [str], pos-1, pos+len-1 ) : xlate rs
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

{-
readMyBam :: L.ByteString -> [Region]
readMyBam s = map xlate $ filter (\r -> b_flag r .&. 4 == 0) recs
  where
    Bam _ refs recs = readBam s
    getref r = L.fromChunks [fst $ refs ! r]
    xlate r = (b_qname r, getref (b_rname r), [str], b_pos r, b_pos r + cigarToAlnLen (b_cigar r))
       where str = if b_flag r .&. 0x10 == 0 then Forward else Reverse
-}

readInput :: L.ByteString -> [Region]
readInput s = case dropWhile (L.all isSpace) (L.lines s) of
    -- _        | isBam s                       -> readMyBam s
    ls@(l:_) | L.head l == '>'               -> readMyFasta ls
             | may_well_be_sam l             -> readMySam ls
    ls                                       -> readMyBed ls
  where
    may_well_be_sam l = let ws = L.split '\t' l
                        in L.head l == '@' || (length ws >= 11 && isJust (textual_cigar_to_len (ws !! 5)))

readMyFasta :: [L.ByteString] -> [Region]
readMyFasta = concatMap proc_line . filter (\s -> not (L.null s) && L.head s == '>')
  where
    proc_line = map parse_label . L.split ';' . L.tail
    parse_label label = (label, crm, read_strand [str], lft, rht)
        where
            (crm, l1) = L.break (':'==) label
            (str, l2) = L.break (':'==) $ L.drop 1 l1
            (lft, l3) = fromMaybe (error "parse error in FASTA header") $ L.readInt $ L.drop 1 l2
            (rht, _) = fromMaybe (error "parse error in FASTA header") $ L.readInt $ L.drop 1 l3


readMyBed :: [L.ByteString] -> [Region]
readMyBed = foldr go [] . map (L.split '\t')
  where go [] = id
        go (w:_) | w == L.pack "track" || L.head w == '#' = id
        go (crm:s:e:gi:strand) = (:) (gi, crm, read_strand (drop 1 strand), ereadInt s, ereadInt e)
        go ws = error $ "not enough fields in BED file" ++ case ws of
                    [w] | length (L.words w) >= 4 -> " (whitespace instead of tabs?)"
                    _                             -> ""

read5cFile :: [L.ByteString] -> [Region]
read5cFile = foldr go [] . map L.words
  where go [] = id
        go (w:_) | L.head w == '#' = id
        go (gi:crm:s:e:strand) = (:) (gi, crm, read_strand strand, ereadInt s, ereadInt e)
        go _ = error $ "too few fields in legacy annotation file"

read_strand :: [L.ByteString] -> [Sense]
read_strand [   ]                      = [Forward, Reverse]
read_strand (s:_) | L.null s           = [Forward, Reverse]
                  | L.head s == '-'    = [Reverse]
                  | s == L.pack "0"    = [Reverse, Forward]
                  | otherwise          = [Forward]

xlate_chrom :: Maybe ChromTable -> Region -> Region
xlate_chrom Nothing reg = reg                                       -- no table, no translation
xlate_chrom (Just ct) (n, crm, ss, s, e) = case M.lookup crm ct of
        Just (c,o) -> (n, c, ss, s+o-1, e+o)                        -- hit, so translate
        Nothing -> error $ "chromosome unknown: " ++ L.unpack crm   -- translate failed


ereadInt :: L.ByteString -> Int
ereadInt s = case L.readInt s of Just (x,_) -> x
                                 Nothing -> error $ "couldn't parse Int at " ++ L.unpack s


myReadFilesWith :: ( FilePath -> L.ByteString -> IO a ) -> [ FilePath ] -> IO [a]
myReadFilesWith k [] = fmap (:[]) $ L.getContents >>= k "stdin"
myReadFilesWith k fs = forM fs $ \f -> if f == "-" then L.getContents >>= k "stdin" else L.readFile f >>= k f


