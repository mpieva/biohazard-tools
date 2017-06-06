-- This sorts a number of Bam files by coordinate.  It's not exactly
-- fast and usually provides little benefit compared to samtools.
-- However, it serves as a stress test for the externally backed
-- PriorityQueue module; it should work under all circumstances (lots of
-- input, tight memory, spill to thousands of files), and never run out
-- of memory or file descriptors or slow down to impractical levels.
--
-- Since this is not generally useful, we only compile this as a test
-- suite and don't bother installing it.
--
-- XXX  One way this could be useful is if it provided a sort that was
--      efficient when applied to the "almost sorted" output of
--      bam-fixpair (or a similar program, its postprocessed output).
--      The idea would be to copy everything to a temporary file that
--      follows the sort order.  Those records that don't fit the normal
--      order go to a PriorityQueue; at the end, we unfold the
--      PriorityQueue and merge with the temporary stream.  If this
--      avoids most of the merging work (and keeps the number of streams
--      to be read from low), it should be faster than a general sort.
--      It still revert to the general sort if the input is unsorted.

import Bio.Prelude
import Bio.Bam
import PriorityQueue
import System.IO

pqconf :: PQ_Conf
pqconf = PQ_Conf 300 30 "/var/tmp/" $ hPutStrLn stderr

main :: IO ()
main = getArgs >>= main'

main' :: [String] -> IO ()
main' [           ] = exitSuccess
main' (out : files) =
    withFile out WriteMode $ \hdl ->
        concatInputs files >=> run $ \hdr -> do
            pq <- progressNum "read" 1000000 (hPutStrLn stderr) =$
                  foldStreamM (\pq br -> enqueuePQ pqconf (BySelfPos br) pq) makePQ
            liftIO $ hPutStrLn stderr "Input phase done, writing output."
            enumPQ pq >=> run $
                progressNum "wrote" 1000000 (hPutStrLn stderr) =$
                writeBamHandle hdl hdr

enumPQ :: Monad m => PQ BySelfPos -> Enumerator [BamRaw] m b
enumPQ pq it = case viewMinPQ pq of
    Nothing -> return it
    Just (BySelfPos r, pq') -> enumPure1Chunk [r] it >>= enumPQ pq'
