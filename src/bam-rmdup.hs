import Bio.Bam
import Bio.Bam.Rmdup
import Bio.Prelude
import Bio.Util.Numeric                         ( showNum, showOOM, estimateComplexity )
import Data.ByteString.Builder
import Paths_biohazard_tools                    ( version )
import PriorityQueue
import System.Console.GetOpt
import System.IO                                ( withFile, IOMode(..) )

import qualified Data.ByteString.Char8          as S
import qualified Data.HashMap.Strict            as M
import qualified Data.IntMap                    as I
import qualified Data.Sequence                  as Z
import qualified Data.Vector.Generic            as V
import qualified Data.Vector.Unboxed            as U hiding ( replicate )
import qualified Data.Vector.Unboxed.Mutable    as U ( replicate, read, write )

data Conf = Conf {
    output :: Maybe ((BamRec -> Seqid) -> BamMeta -> Iteratee [BamRec] IO ()),
    strand_preserved :: Bool,
    collapse :: Bool -> Collapse,
    clean_multimap :: BamRec -> IO (Maybe BamRec),
    keep_all :: Bool,
    keep_unaligned :: Bool,
    keep_improper :: Bool,
    transform :: BamRec -> Maybe BamRec,
    min_len :: Int,
    min_qual :: Qual,
    get_label :: HashMap Seqid Seqid -> BamRec -> Seqid,
    putResult :: String -> IO (),
    putHistogram :: U.Vector Int -> IO (),
    debug :: String -> IO (),
    which :: Which,
    circulars :: Refs -> IO (IntMap (Seqid,Int), Refs),
    pqconf :: PQ_Conf }

-- | Which reference sequences to scan
data Which = Allrefs | Some Refseq Refseq | Unaln deriving Show

defaults :: Conf
defaults = Conf { output = Nothing
                , strand_preserved = True
                , collapse = cons_collapse' (Q 60)
                , clean_multimap = check_flags
                , keep_all = False
                , keep_unaligned = False
                , keep_improper = False
                , transform = Just
                , min_len = 0
                , min_qual = Q 0
                , get_label = get_library
                , putResult = putStr
                , putHistogram = \_ -> return ()
                , debug = \_ -> return ()
                , which = Allrefs
                , circulars = \rs -> return (I.empty, rs)
                , pqconf = PQ_Conf 1000 200 "/var/tmp/" undefined }

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option  "o" ["output"]         (ReqArg set_output "FILE") "Write to FILE (default: no output, count only)",
    Option  "O" ["output-lib"]     (ReqArg set_lib_out "PAT") "Write each lib to file named following PAT",
    Option  "H" ["histogram"]      (ReqArg set_hist   "FILE") "Write histogram to FILE",
    Option  [ ] ["debug"]          (NoArg  set_debug_out)     "Write textual debugging output",
    Option  "z" ["circular"]  (ReqArg add_circular "CHR:LEN") "Refseq CHR is circular with length LEN",
    Option  "R" ["refseq"]         (ReqArg set_range "RANGE") "Read only range of reference sequences",
    Option  "p" ["improper-pairs"] (NoArg  set_improper)      "Include improper pairs",
    Option  "u" ["unaligned"]      (NoArg  set_unaligned)     "Include unaligned reads and pairs",
    Option  "1" ["single-read"]    (NoArg  set_single)        "Pretend there is no second mate",
    Option  "m" ["multimappers"]   (NoArg  set_multi)         "Process multi-mappers (by dropping secondary alignments)",
    Option  "c" ["cheap"]          (NoArg  set_cheap)         "Cheap computation: skip the consensus calling",
    Option  "k" ["keep","mark-only"] (NoArg set_keep)         "Mark duplicates, but include them in output",
    Option  "Q" ["max-qual"]       (ReqArg set_qual "QUAL")   "Set maximum quality after consensus call to QUAL",
    Option  "l" ["min-length"]     (ReqArg set_len "LEN")     "Discard reads shorter than LEN",
    Option  "q" ["min-mapq"]       (ReqArg set_mapq "QUAL")   "Discard reads with map quality lower than QUAL",
    Option  "s" ["no-strand"]      (NoArg  set_no_strand)     "Strand of alignments is uninformative",
    Option  "r" ["ignore-rg"]      (NoArg  set_no_rg)         "Ignore read groups when looking for duplicates",
    Option  "M" ["max-memory"]     (ReqArg set_max_mem "MB")  "Use at most MB megabytes per queue",
    Option  "L" ["max-merge"]      (ReqArg set_max_mrg "NUM") "Merge at most NUM files at a time",
    Option  "T" ["temp-path"]    (ReqArg set_temp_pth "PATH") "Store temporary files in PATH",
    Option  "v" ["verbose"]        (NoArg  set_verbose)       "Print more diagnostics",
    Option "h?" ["help","usage"]   (NoArg  (const usage))     "Display this message",
    Option "V"  ["version"]        (NoArg  (const vrsn))      "Display version number and exit" ]

  where
    set_output "-" c =                    return $ c { output = Just $ \_ -> pipeBamOutput, putResult = hPutStr stderr }
    set_output   f c =                    return $ c { output = Just $ \_ -> writeBamFile f }
    set_lib_out  f c =                    return $ c { output = Just $       writeLibBamFiles f }
    set_debug_out  c =                    return $ c { output = Just $ \_ -> pipeSamOutput, putResult = hPutStr stderr }
    set_qual     n c = readIO n >>= \q -> return $ c { collapse = cons_collapse' (Q q) }
    set_no_strand  c =                    return $ c { strand_preserved = False }
    set_verbose    c =                    return $ c { debug = hPutStr stderr }
    set_improper   c =                    return $ c { keep_improper = True }
    set_single     c =                    return $ c { transform = make_single }
    set_cheap      c =                    return $ c { collapse = cheap_collapse' }
    set_keep       c =                    return $ c { keep_all = True }
    set_unaligned  c =                    return $ c { keep_unaligned = True }
    set_len      n c = readIO n >>= \l -> return $ c { min_len = l }
    set_mapq     n c = readIO n >>= \q -> return $ c { min_qual = Q q }
    set_no_rg      c =                    return $ c { get_label = get_no_library }
    set_multi      c =                    return $ c { clean_multimap = clean_multi_flags }
    set_hist    fp c =                    return $ c { putHistogram = writeHistogramFile fp }
    set_max_mem  a c = readIO a >>= \x -> return $ c { pqconf = (pqconf c) { max_mb    = x } }
    set_max_mrg  a c = readIO a >>= \x -> return $ c { pqconf = (pqconf c) { max_merge = x } }
    set_temp_pth a c =                    return $ c { pqconf = (pqconf c) { temp_path = a } }

    set_range    a c
        | a == "A" || a == "a" = return $ c { which = Allrefs }
        | a == "U" || a == "u" = return $ c { which = Unaln }
        | otherwise = case reads a of
                [ (x,"")    ] -> return $ c { which = Some (Refseq $ x-1) (Refseq $ x-1) }
                [ (x,'-':b) ] -> readIO b >>= \y ->
                                 return $ c { which = Some (Refseq $ x-1) (Refseq $ y-1) }
                _ -> fail $ "parse error in " ++ show a

    add_circular a c = case break ((==) ':') a of
        (nm,':':r) -> case reads r of
            [(l,[])] | l > 0 -> return $ c { circulars = add_circular' (fromString nm) l (circulars c) }
            _ -> fail $ "couldn't parse length " ++ show r ++ " for " ++ show nm
        _ -> fail $ "couldn't parse \"circular\" argument " ++ show a

    add_circular' nm l io refs = do
        (m1, refs') <- io refs
        case filter (S.isPrefixOf nm . sq_name . snd) $ zip [0..] $ toList refs' of
            [(k,a)] | sq_length a >= l -> let m2     = I.insert k (sq_name a,l) m1
                                              refs'' = Z.update k (a { sq_length = l }) refs'
                                          in return (m2, refs'')
                    | otherwise -> fail $ "cannot wrap " ++ show nm ++ " to " ++ show l
                                       ++ ", which is more than the original " ++ show (sq_length a)
            [] -> fail $ "no match for target sequence " ++ show nm
            _ -> fail $ "target sequence " ++ show nm ++ " is ambiguous"

vrsn :: IO a
vrsn = do pn <- getProgName
          hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
          exitSuccess

usage :: IO a
usage = do p <- getProgName
           hPutStrLn stderr $ "Usage: " ++ usageInfo (p ++ info) options
           exitSuccess
  where
    info = " [option...] [bam-file...]\n\
           \Removes PCR duplicates from BAM files and calls a consensus for each duplicate set.  \n\
           \Input files must be sorted by coordinate and are merged on the fly.  Options are:"

cons_collapse' :: Qual -> Bool -> Collapse
cons_collapse' m False = cons_collapse m
cons_collapse' m True  = cons_collapse_keep m

cheap_collapse' :: Bool -> Collapse
cheap_collapse'  False = cheap_collapse
cheap_collapse'  True  = cheap_collapse_keep

-- | Get library from BAM record.
-- This gets the read group from a bam record, then the library for read
-- group.  This will work correctly if and only if the RG-LB field is
-- the name of the "Ur-Library", the common one before the first
-- amplification.
--
-- If no RG-LB field is present, RG-SM is used instead.  This will work
-- if and only if no libraries were aliquotted and then pooled again.
--
-- Else the RG-ID field is used.  This will work if and only if read
-- groups correspond directly to libraries.
--
-- If no RG is present, the empty string is returned.  This serves as
-- fall-back.

get_library, get_no_library :: M.HashMap Seqid Seqid -> BamRec -> Seqid
get_library  tbl br = M.lookupDefault rg rg tbl where rg = extAsString "RG" br
get_no_library _  _ = S.empty

mk_rg_tbl :: BamMeta -> M.HashMap Seqid Seqid
mk_rg_tbl hdr = M.fromList
    [ (rg_id, rg_lb)
    | ("RG",fields) <- meta_other_shit hdr
    , rg_id <- take 1   [ i | ("ID",i) <- fields ]
    , rg_lb <- take 1 $ [ l | ("LB",l) <- fields ]
                     ++ [ s | ("SM",s) <- fields ]
                     ++ [ rg_id ] ]

data Counts = Counts { tin          :: !Int
                     , tout         :: !Int
                     , good_singles :: !Int
                     , good_total   :: !Int }

main :: IO ()
main = do
    args <- getArgs
    when (null args) usage
    let (opts, files, errors) = getOpt Permute options args
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Conf{..} <- foldr (>=>) return opts defaults

    let pqconf1 = pqconf { croak = debug, max_mb = 1 +      max_mb pqconf `div` 16 }
        pqconf2 = pqconf { croak = debug, max_mb = 1 + 15 * max_mb pqconf `div` 16 }

    add_pg <- addPG $ Just version
    mergeInputRanges which files >=> run $ \hdr -> do
       (circtable_, refs') <- liftIO $ circulars (meta_refs hdr)
       let tbl = mk_rg_tbl hdr
           circtable = I.map (\(nm,l) -> (nm, l, maybe l sq_length $ find ((==) nm . sq_name) (meta_refs hdr))) circtable_

       unless (M.null tbl) $ liftIO $ do
                debug "mapping of read groups to libraries:\n"
                mapM_ debug [ unpack k ++ " --> " ++ unpack v ++ "\n" | (k,v) <- M.toList tbl ]

       let cleanup = cleanup2 . transform . unpackBam

           cleanup2  Nothing  = return []
           cleanup2 (Just  b) = cleanup3 <$> clean_multimap b

           cleanup3 (Just  b)
                | keep_unaligned || is_aligned b
                , keep_improper || is_proper b
                , eff_len b >= min_len              = [b]
           cleanup3 _                               = [ ]

       (ct,hist,ou) <- progressBam "Rmdup at" refs' 0x8000 debug                                          =$
                       takeWhileE is_halfway_aligned                                                      =$
                       concatMapStreamM cleanup                                                           =$
                       ( let norm (nm,lnat,lpad) r = [ normalizeFwd nm lnat lpad r ]
                         in mapSortAtGroups pqconf1 Preserved circtable norm )                             =$
                       filterStream (\b -> b_mapq b >= min_qual)                                          =$
                       case output of
                            Nothing -> rmdup (get_label tbl) strand_preserved (cheap_collapse' keep_all)  =$
                                       zipStreams3 (mapStream snd =$ count_all (get_label tbl))
                                                   (mapStream fst =$ collect_histogram)
                                                   (skipToEof <$ skipToEof)

                            Just  o -> rmdup (get_label tbl) strand_preserved (collapse keep_all)         =$
                                       zipStreams3 (mapStream snd =$ count_all (get_label tbl))
                                                   (mapStream fst =$ collect_histogram)
                                                   (mapStream snd =$
                                                    (let wrap (_nm,ln,_lp) = map (either id id) . wrapTo ln
                                                     in mapSortAtGroups pqconf2 Destroyed circtable wrap $
                                                     o (get_label tbl) (add_pg hdr { meta_refs = refs' })))

       let do_copy = progressBam "Copying junk at" refs' 0x8000 debug ><>
                     concatMapStreamM cleanup

       case which of Unaln              -> lift (run $ do_copy =$ ou)
                     _ | keep_unaligned -> lift (run $ do_copy =$ ou)
                     _                  -> lift (run ou)

       liftIO $ do putHistogram hist
                   putResult . unlines $
                       "\27[K#RG\tin\tout\tin@MQ20\tsingle@MQ20\tunseen\ttotal\t%unique\t%exhausted"
                       : map (uncurry do_report) (M.toList ct)


do_report :: Seqid -> Counts -> String
do_report lbl Counts{..} = intercalate "\t" fs
  where
    fs = label : showNum tin : showNum tout : showNum good_total : showNum good_singles :
         report_estimate (estimateComplexity good_total good_singles)

    label = if S.null lbl then "--" else unpack lbl

    report_estimate  Nothing                = [ "N/A" ]
    report_estimate (Just good_grand_total) =
            [ showOOM (grand_total - fromIntegral tout)
            , showOOM grand_total
            , showFFloat (Just 1) rate []
            , showFFloat (Just 1) exhaustion [] ]
      where
        grand_total = good_grand_total * fromIntegral tout / fromIntegral good_total
        exhaustion  = 100 * fromIntegral good_total / good_grand_total
        rate        = 100 * fromIntegral tout / fromIntegral tin :: Double


collect_histogram :: MonadIO m => Iteratee [Int] m (U.Vector Int)
collect_histogram = do
    vec <- liftIO $ U.replicate 1024 0
    mapStreamM_ $ \i -> when (i<1024) . liftIO $ U.read vec i >>= U.write vec i . succ
    liftIO $ U.unsafeFreeze vec

writeHistogramFile :: FilePath -> U.Vector Int -> IO ()
writeHistogramFile fp v = withFile fp WriteMode $ \hdl ->
    hPutBuilder hdl .
    U.ifoldr (\i f b -> intDec i <> char7 '\t' <> intDec f <> char7 '\n' <> b) mempty .
    U.reverse . U.dropWhile (==0) . U.reverse $ v

-- | Counting reads:  we count total read in (ti), total reads out (to),
-- good (confidently mapped) singletons out (gs), total good
-- (confidently mapped) reads out (gt).  Paired reads count 1, unpaired
-- reads count 2, and at the end we divide by 2.  This ensures that we
-- don't double count mate pairs, while still working mostly sensibly in
-- the presence of broken BAM files.

count_all :: Monad m => (BamRec -> Seqid) -> Iteratee [BamRec] m (M.HashMap Seqid Counts)
count_all lbl = M.map fixup `liftM` foldStream plus M.empty
  where
    plus m b = M.insert (lbl b) cs m
      where
        !cs = plus1 (M.lookupDefault (Counts 0 0 0 0) (lbl b) m) b

    plus1 (Counts ti to gs gt) b = Counts ti' to' gs' gt'
      where
        !w   = if isPaired b then 1 else 2
        !ti' = ti + w * extAsInt 1 "XP" b
        !to' = to + w
        !gs' = if b_mapq b >= Q 20 && extAsInt 1 "XP" b == 1 then gs + w else gs
        !gt' = if b_mapq b >= Q 20                           then gt + w else gt

    fixup (Counts ti to gs gt) = Counts (div ti 2) (div to 2) (div gs 2) (div gt 2)

eff_len :: BamRec -> Int
eff_len b | isProperlyPaired b = abs $ b_isize b
          | otherwise          = V.length $ b_seq b

is_halfway_aligned :: BamRaw -> Bool
is_halfway_aligned = isValidRefseq . b_rname . unpackBam

is_aligned :: BamRec -> Bool
is_aligned b = not (isUnmapped b && isMateUnmapped b) && isValidRefseq (b_rname b)

is_proper :: BamRec -> Bool
is_proper b = not (isPaired b) || (isMateUnmapped b == isUnmapped b && isProperlyPaired b)

make_single :: BamRec -> Maybe BamRec
make_single b | isPaired b && isSecondMate b = Nothing
              | isUnmapped b                 = Nothing
              | not (isPaired b)             = Just b
              | otherwise = Just $ b { b_flag = b_flag b .&. complement pair_flags
                                     , b_mrnm = invalidRefseq
                                     , b_mpos = invalidPos
                                     , b_isize = 0 }
  where
    pair_flags = flagPaired .|. flagProperlyPaired .|.
                 flagFirstMate .|. flagSecondMate .|.
                 flagMateUnmapped


mergeInputRanges :: (MonadIO m, MonadMask m)
                 => Which -> [FilePath] -> Enumerator' BamMeta [BamRaw] m a
mergeInputRanges Allrefs  fps   = mergeInputs combineCoordinates fps
mergeInputRanges  _  [        ] = \k -> return $ k mempty
mergeInputRanges rng (fp0:fps0) = go fp0 fps0
  where
    enum1  fp k1 = case rng of Allrefs  -> decodeAnyBamFile                 fp k1
                               Some x y -> decodeBamFileRange           x y fp k1
                               Unaln    -> decodeWithIndex eneeBamUnaligned fp k1

    go fp [       ] = enum1 fp
    go fp (fp1:fps) = mergeEnums' (go fp1 fps) (enum1 fp) combineCoordinates

    decodeBamFileRange x y = decodeWithIndex $
            \idx -> foldr ((>=>) . eneeBamRefseq idx) return [x..y]


decodeWithIndex :: (MonadIO m, MonadMask m)
                => (BamIndex () -> Enumeratee [BamRaw] [BamRaw] m a)
                -> FilePath -> (BamMeta -> Iteratee [BamRaw] m a)
                -> m (Iteratee [BamRaw] m a)
decodeWithIndex enum fp k0 = do
    idx <- liftIO $ readBamIndex fp
    decodeAnyBamFile fp >=> run $ enum idx . k0


writeLibBamFiles :: MonadIO m
                 => FilePath -> (BamRec -> Seqid) -> BamMeta -> Iteratee [BamRec] m ()
writeLibBamFiles fp lbl hdr = tryHead >>= go M.empty
  where
    go m  Nothing  = liftIO . mapM_ run $ M.elems m
    go m (Just br) = do
        let !l = lbl br
        let !it = M.lookupDefault (writeBamFile (fp `subst` l) hdr) l m
        it' <- liftIO $ enumPure1Chunk [br] it
        let !m' = M.insert l it' m
        tryHead >>= go m'

    subst [            ] _ = []
    subst ('%':'s':rest) l = unpack l ++ subst rest l
    subst ('%':'%':rest) l = '%' : subst rest l
    subst ( c :    rest) l =  c  : subst rest l


check_flags :: Monad m => BamRec -> m (Maybe BamRec)
check_flags b | extAsInt 1 "HI" b /= 1 = fail "cannot deal with HI /= 1"
              | extAsInt 1 "IH" b /= 1 = fail "cannot deal with IH /= 1"
              | extAsInt 1 "NH" b /= 1 = fail "cannot deal with NH /= 1"
              | otherwise              = return $ Just b

clean_multi_flags :: Monad m => BamRec -> m (Maybe BamRec)
clean_multi_flags b = return $ if extAsInt 1 "HI" b /= 1 then Nothing else Just b'
  where
    b' = b { b_exts = deleteE "HI" $ deleteE "IH" $ deleteE "NH" $ b_exts b }

-- | Groups sorted records and selectively maps a function over some of
-- these groups, while preserving sorting.
--
-- This function is applied twice, once for the normalization phase and
-- once for the wrapping around phase.  Normalization only ever
-- increases coordinates, so ordering is preserved and we can stream
-- most of the reads;  this is done if the 'order' parameter is
-- 'Preserved'.
--
-- Wrapping sometimes decreases coordinates, so we need to buffer
-- everything to restore sorting.  We simply run everything through a
-- 'PriorityQueue' to get full sorting and controlled external
-- buffering.  While some microoptimizations are possible, this is easy,
-- robust, and tunable.

data Order = Destroyed | Preserved

mapSortAtGroups :: MonadIO m => PQ_Conf -> Order -> IntMap a -> (a -> BamRec -> [BamRec]) -> Enumeratee [BamRec] [BamRec] m b
mapSortAtGroups cf order m f = eneeCheckIfDonePass no_group
  where
    no_group k (Just e) = idone (liftI k) $ EOF (Just e)
    no_group k Nothing  = tryHead >>= maybe (idone (liftI k) $ EOF Nothing) (\a -> no_group_1 a k Nothing)

    no_group_1 _ k (Just e) = idone (liftI k) $ EOF (Just e)
    no_group_1 a k Nothing  =
        case I.lookup (b_rname_int a) m of
            Nothing  -> eneeCheckIfDonePass no_group . k $ Chunk [a]
            Just arg -> do pq <- liftIO . pushTo makePQ $ f arg a
                           cont_group (b_rname a) arg pq k Nothing

    cont_group _rn _arg _pq k (Just  e) = idone (liftI k) $ EOF (Just e)
    cont_group  rn  arg  pq k  Nothing  = tryHead >>= maybe flush_eof check1
      where
        flush_eof  = enumPQ (const $ return                            ) (Refseq maxBound, maxBound) pq k
        flush_go a = enumPQ (const $ eneeCheckIfDonePass (no_group_1 a)) (Refseq maxBound, maxBound) pq k
        check1 a | b_rname a == rn = do pq' <- liftIO . pushTo pq $ f arg a
                                        case order of
                                          Preserved -> enumPQ (eneeCheckIfDonePass . cont_group rn arg) (Refseq 0, 0) pq' k
                                          Destroyed -> cont_group rn arg pq' k Nothing
                 | otherwise       = flush_go a

    b_rname_int = fromIntegral . unRefseq . b_rname

    pushTo :: PQ BySelfPos -> [BamRec] -> IO (PQ BySelfPos)
    pushTo = foldM (\q a -> packBam a >>= \b -> enqueuePQ cf (BySelfPos b) q)

    -- Flushes from the PQ anything with position less than 'pos'.
    enumPQ :: Monad m
           => (PQ BySelfPos -> Iteratee [BamRec] m b
                    -> Iteratee [BamRec] m (Iteratee [BamRec] m b))
           -> (Refseq, Int) -> PQ BySelfPos
           -> (Stream  [BamRec] -> Iteratee [BamRec] m b)
           -> Iteratee [BamRec] m (Iteratee [BamRec] m b)
    enumPQ ke pos pq k = case viewMinPQ pq of
        Just (BySelfPos r, pq')
            | br_self_pos r < pos
                -> eneeCheckIfDone (enumPQ ke pos pq') . k $ Chunk [unpackBam r]
        _       -> ke pq (liftI k)


-- | Normalize a read's alignment to fall into the "forward canonical
-- region" of [e..e+l].  Takes the name of the reference sequence, its
-- length and its padded length.  This normalization is ensured to work
-- (doesn't need information from the mate) and always increases
-- coordinates.
normalizeFwd :: Seqid -> Int -> Int -> BamRec -> BamRec
normalizeFwd nm lnat lpad b = b { b_pos  = norm $ b_pos  b
                                , b_mpos = norm $ b_mpos b
                                , b_mapq = if dups_are_fine then Q 37 else b_mapq b
                                , b_exts = if dups_are_fine then deleteE "XA" (b_exts b) else b_exts b }
  where
    dups_are_fine  = all_match_XA (extAsString "XA" b)
    all_match_XA s = case S.split ';' s of [xa1, xa2] | S.null xa2 -> one_match_XA xa1
                                           [xa1]                   -> one_match_XA xa1
                                           _                       -> False
    one_match_XA s = case S.split ',' s of (sq:pos:_) | sq == nm   -> pos_match_XA pos ; _ -> False
    pos_match_XA s = case S.readInt s   of Just (p,z) | S.null z   -> int_match_XA p ;   _ -> False
    int_match_XA p | p >= 0    = norm  (p-1) == norm (b_pos b) && not (isReversed b)
                   | otherwise = norm (-p-1) == norm (b_pos b) && isReversed b

    norm p = o + (p-o) `mod` lnat where o = lpad - lnat

