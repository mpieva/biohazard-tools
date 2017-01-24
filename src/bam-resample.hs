-- Resample m out of n `virtual' BAM records.
--
-- Strategy for fair down sampling:  we first count the number of
-- records, then scan again to sample.  Input must be grouped by QNAME
-- (sorted by QNAME is fine).
--
-- Usage: resample [NUM] [FILE...]

import Bio.Bam
import Bio.Prelude
import Paths_biohazard_tools ( version )
import System.Random

main :: IO ()
main = do
    args <- getArgs
    case args of
        []             -> complain
        [_num]         -> complain
        num_ : files   -> case reads num_ of
            [(num,"")] -> main' num files
            _          -> complain

complain :: IO ()
complain = do pn <- getProgName
              hPutStr stderr $ pn ++ ", version " ++ showVersion version
                            ++ "\nUsage: " ++ pn ++ " <num> [file...]\n"
              exitFailure

main' :: Int -> [String] -> IO ()
main' num files = do
    hPutStr stderr "counting... "
    total <- enumInputs files >=> run $
             joinI $ decodeAnyBam $ \_hdr ->
             joinI $ groupOn (b_qname . unpackBam) $
             foldStream (\a _ -> 1+a) 0
    hPutStr stderr $ shows total " records.\n"

    add_pg <- addPG (Just version)
    enumInputs files >=> run $
             joinI $ decodeAnyBam $
             joinI . groupOn (b_qname . unpackBam) .
             joinI . resample num total .
             protectTerm . pipeBamOutput . add_pg


resample :: MonadIO m => Int -> Int -> Enumeratee [[BamRaw]] [BamRaw] m a
resample m0 n0 | m0 > n0 = error "upsampling requested"
resample m0 n0 = eneeCheckIfDone (go m0 n0)
  where
    go  !m !n k = tryHead >>= maybe (return (liftI k)) (go' m n k)
    go' !m !n k a = do r <- liftIO $ randomRIO (0,n-1)
                       if r < m
                         then eneeCheckIfDone (go (m-1) (n-1)) . k $ Chunk a
                         else go m (n-1) k

groupOn :: (Monad m, Eq b) => (a -> b) -> Enumeratee [a] [[a]] m c
groupOn f = eneeCheckIfDone (\k -> tryHead >>= maybe (return $ liftI k) (\a -> go k [a] (f a)))
  where
    go  k acc fa = tryHead >>= maybe (return . k $ Chunk [reverse acc]) (go' k acc fa)
    go' k acc fa b | fa == f b = go k (b:acc) fa
                   | otherwise = eneeCheckIfDone (\k' -> go k' [b] (f b)) . k $ Chunk [reverse acc]
