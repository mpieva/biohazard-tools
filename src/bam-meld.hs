-- Reads multiple BAM files, melds them by keeping the best hit for
-- every entry.  All input files must be parallel (same reads, same
-- order, no omissions).  The best hit and the new mapq are calculated
-- by combining appropriate optional fields.  Presets exist for common
-- types of aligners, other schemes can be configured flexibly.
--
-- Paired end support is easy:  since all input files are either
-- unsorted (and strictly parallel) or sorted by read name, pairs are
-- also sorted together.  So all we have to do is (maybe) exchange first
-- and seocnd mate.

import Bio.Bam
import Bio.Prelude
import Paths_biohazard_tools                    ( version )
import System.Console.GetOpt

import qualified Data.ByteString.Char8 as S
import qualified Data.Sequence         as Z

data Conf = Conf {
    c_score  :: Maybe (BamPair -> Int),
    c_output :: BamMeta -> Iteratee [BamRec] IO (),
    c_merge  :: Enumeratee [BamPair] [[BamPair]] (Iteratee [[BamPair]] IO) () }

defaultConf :: Conf
defaultConf = Conf Nothing (protectTerm . pipeBamOutput) iter_transpose

defaultScore :: BamPair -> Int
defaultScore r = 30 * getExt "XM" r + 45 * getExt "XO" r + 15 * getExt "XG" r

getExt :: BamKey -> BamPair -> Int
getExt k (Single a) = extAsInt 0 k a
getExt k (Pair a b) = extAsInt 0 k a + extAsInt 0 k b


-- | Enumerates a list of BAM files.  Meta records are merged sensibly,
-- records are merged using the supplied "merging Enumeratee".  Results
-- in something close to an Enumerator (not quite, because the merged
-- headers need to be passed along).
enum_bam_files :: (MonadIO m, MonadMask m)
               => Enumeratee [BamPair] [[BamPair]] (Iteratee [[BamPair]] m) a
               -> [ FilePath ]
               -> Enumerator' BamMeta [[BamPair]] m a
enum_bam_files _etee [    ] = return . ($ mempty)
enum_bam_files  etee (f1:fs1) = go (decodeAnyBamOrSamFile f1 $== find_pairs $== mapStream (:[])) fs1
  where
    go e1 [    ] k = e1 k
    go e1 (f:fs) k = go e1 fs $
                          \h1 -> (decodeAnyBamOrSamFile f $== adjust h1 $== find_pairs)
                         (\h2 -> joinI . etee $ ilift lift (k $ h1 `mappend` h2)) >>= run

      -- How to merge?  We keep the stream from e1 as is, the refids in
      -- e2 are shifted down by the number of refseqs in h1.  Headers
      -- are merged by concatenating the reference lists and appending
      -- the headers using mappend.  The Monoid instance does
      -- everything.

    adjust h = let o    = Z.length . meta_refs $ h
                   f br = br { b_rname = b_rname br `plus` o
                             , b_mrnm  = b_mrnm  br `plus` o }
               in mapStream f

    r        `plus` _ | r == invalidRefseq = r
    Refseq r `plus` o                      = Refseq (r + fromIntegral o)

data BamPair = Single BamRec | Pair BamRec BamRec

find_pairs :: Enumeratee [BamRec] [BamPair] m a
find_pairs = mapStream Single

unpair :: Enumeratee [BamPair] [BamRec] m a
unpair = mapChunks (concatMap unpair1)
  where
    unpair1 (Single a) = [a]
    unpair1 (Pair a b) = [a,b]

p_qname :: BamPair -> Seqid
p_qname (Single a) = b_qname a
p_qname (Pair a _) = b_qname a

p_mapq :: BamPair -> Qual
p_mapq (Single a) = b_mapq a
p_mapq (Pair a _) = b_mapq a

p_is_unmapped :: BamPair -> Bool
p_is_unmapped (Single a) = isUnmapped a
p_is_unmapped (Pair a b) = isUnmapped a && isUnmapped b

set_mapq :: BamPair -> Qual -> BamPair
set_mapq (Single a) q = Single (a { b_mapq = q })
set_mapq (Pair a b) q = Pair (a { b_mapq = q }) (b { b_mapq = q })

meld :: BamMeta -> (BamPair -> Int) -> [BamPair] -> BamPair
meld hdr score rs | all p_is_unmapped rs = head rs
                  | all_equal (map p_qname rs) = set_mapq best' mapq
                  | otherwise = error $ "BAMs are not in the same order or sequences are missing: "
                                     ++ show (map p_qname rs)
  where
    all_equal [] = error "no input (not supposed to happen)"
    all_equal (x:xs) = all (x ==) xs

    ( best : rs' ) = sortBy (\a b -> score a `compare` score b) $ filter (not . p_is_unmapped) rs
    mapq = case rs' of [    ] -> p_mapq best
                       (r2:_) -> Q . fromIntegral $ fromIntegral (unQ (p_mapq best))
                                              `min` (score r2 - score best)


    split_xa br = let s = extAsString "XA" br in if S.null s then id else (++) (S.split ';' s)

    get_xas (Single a) (one,two) = (split_xa a one, two)
    get_xas (Pair a b) (one,two) = (split_xa a one, split_xa b two)

    (xas1, xas2) = foldr enc_xas (foldr get_xas ([],[]) (best:rs')) rs'

    add_xas xas b = b { b_exts = updateE "XA" (Text (S.intercalate (S.singleton ';') xas)) (b_exts b) }

    best' = case best of Single a -> Single (add_xas xas1 a)
                         Pair a b -> Pair (add_xas xas1 a) (add_xas xas2 b)

    enc_xas (Single a) (one,two) = (encode a one,two)
    enc_xas (Pair a b) (one,two) = (encode a one,encode b two)

    encode b xas | isUnmapped b = xas
                 | otherwise = S.intercalate (S.singleton ',') [ rnm, pos, cig, nm ] : xas
      where
        nm =  fromString $ show $ extAsInt 0 "NM" b
        cig = fromString $ show $ b_cigar b
        pos = fromString $ (if isReversed b then '-' else '+') : show (b_pos b)
        rnm = sq_name $ getRef (meta_refs hdr) (b_rname b)


options :: [OptDescr (Conf -> IO Conf)]
options =
    [ Option "o" ["output"]   (ReqArg set_output "FILE") "Send output to FILE"
    , Option "u" ["unsorted"] (NoArg  set_unsorted)      "Input is unsorted"
    , Option "s" ["sorted"]   (NoArg  set_sorted)        "Input is sorted by name"
    , Option "w" ["weight"]   (ReqArg set_weight "XX:Y") "Set the badness of field XX to Y"
    , Option [ ] ["bwa"]      (NoArg  set_bwa)           "Preset for alignments from 'bwa' (uses XM, XO, XG)"
    , Option [ ] ["anfo"]     (NoArg  set_anfo)          "Preset for alignments from 'anfo' (uses UQ, PQ)"
    , Option [ ] ["blast"]    (NoArg  set_blast)         "Preset for alignments from 'blast' (uses AS)"
    , Option [ ] ["blat"]     (NoArg  set_blat)          "Preset for alignments from 'blat' (uses NM)"
    , Option "h?" ["help","usage"] (NoArg usage)         "Display this information and exit"
    , Option "V"  ["version"]      (NoArg  vrsn)         "Display version number and exit" ]

vrsn :: Conf -> IO Conf
vrsn _ = do pn <- getProgName
            hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
            exitSuccess

usage :: Conf -> IO Conf
usage _ = putStrLn (usageInfo blurb options) >> exitSuccess
  where
    blurb = "Merges multiple bam files containing the same sequences, keeping only\n\
            \the best hit for each.  Attempts to be configurable to bam files from\n\
            \various sources and attempts to calculate a sensible map quality.\n"

set_output :: String -> Conf -> IO Conf
set_output "-" c = return $ c { c_output = pipeBamOutput }
set_output  fn c = return $ c { c_output = writeBamFile fn }

set_unsorted :: Conf -> IO Conf
set_unsorted c = return $ c { c_merge = iter_transpose }

set_sorted :: Conf -> IO Conf
set_sorted c = return $ c { c_merge = merge_by_name }

set_weight :: String -> Conf -> IO Conf
set_weight (a:b:':':rest) c = do
    w <- readIO rest
    let f r = getExt (fromString [a,b]) r * w + maybe 0 ($ r) (c_score c)
    return $ c { c_score = Just f }
set_weight s _ = error $ "illegal weight specification " ++ show s

set_bwa, set_anfo, set_blat, set_blast :: Conf -> IO Conf
set_bwa c = return $ c { c_score = Just defaultScore }
set_anfo c = return $ c { c_score = Just $ \r -> getExt "UQ" r }
set_blat c = return $ c { c_score = Just $ \r -> getExt "NM" r * 30 }
set_blast c = return $ c { c_score = Just $ \r -> getExt "AS" r * (-3) }

main :: IO ()
main = do
    ( opts, files, errors ) <- getOpt Permute options `fmap` getArgs
    conf <- foldM (flip id) defaultConf opts

    let errors' | null files = "no input files" : errors
                | otherwise  = errors

    unless (null errors') $ do
        mapM_ (hPutStrLn stderr) errors'
        exitFailure

    add_pg <- addPG (Just version)
    enum_bam_files (c_merge conf) files >=> run                              $ \hdr ->
        joinI $ mapStream (meld hdr $ fromMaybe defaultScore $ c_score conf) $
        joinI $ unpair $ c_output conf (add_pg hdr)


iter_transpose :: Monad m => Enumeratee [BamPair] [[BamPair]] (Iteratee [[BamPair]] m) a
iter_transpose = eneeCheckIfDone step
  where
    step k = do mx <- tryHead ; my <- lift tryHead ; step' k mx my

    step' k Nothing Nothing = idone (liftI k) $ EOF Nothing
    step' k (Just x) (Just ys) | p_qname x == p_qname (head ys) = iter_transpose . k $ Chunk [x:ys]
    step' _ _ _ = error "files do not contain the same query records"

merge_by_name :: Monad m => Enumeratee [BamPair] [[BamPair]] (Iteratee [[BamPair]] m) a
merge_by_name = ensure_sorting ><> merge'
  where
    merge'     = eneeCheckIfDone (\k -> tryHead >>= \mx -> lift tryHead >>= \my -> merge''' k mx my)
    merge'x my = eneeCheckIfDone (\k -> tryHead >>= \mx ->                         merge''' k mx my)
    merge'y mx = eneeCheckIfDone (\k ->                    lift tryHead >>= \my -> merge''' k mx my)

    merge''' k  Nothing   Nothing  = idone (liftI k) $ EOF Nothing
    merge''' k  Nothing  (Just ys) = merge'y Nothing . k $ Chunk [ys]
    merge''' k (Just  x)  Nothing  = merge'x Nothing . k $ Chunk [[x]]
    merge''' k (Just  x) (Just ys) = case p_qname x `compareNames` p_qname (head ys) of
            LT -> merge'x (Just ys) . k $ Chunk [[  x ]]
            EQ -> merge'            . k $ Chunk [ x:ys ]
            GT -> merge'y (Just  x) . k $ Chunk [   ys ]

ensure_sorting :: Monad m => Enumeratee [BamPair] [BamPair] m a
ensure_sorting = eneeCheckIfDonePass (icont . step)
  where
    step k (EOF       mx) = idone (liftI k) $ EOF mx
    step k (Chunk [    ]) = liftI $ step k
    step k (Chunk (x:xs)) = step' x k $ Chunk xs

    step' x1 k (EOF   mx) = idone (k $ Chunk [ x1 ]) $ EOF mx
    step' x1 k (Chunk []) = liftI $ step' x1 k
    step' x1 k (Chunk (x2:xs)) = case p_qname x1 `compareNames` p_qname x2 of
            GT -> error "input is not sorted by qname"
            _  -> eneeCheckIfDone (\k' -> step' x2 k' (Chunk xs)) . k $ Chunk [ x1 ]

