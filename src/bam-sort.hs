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
--
--      In a quick test, the opportunistic sort is almost twice as fast
--      as the general version on the output of bam-fixpair (with few
--      actual paired-end reads).  Is that worth it?

import Bio.Prelude
import Bio.Bam
import PriorityQueue
import System.IO

pqconf :: PQ_Conf
pqconf = PQ_Conf 300 30 "/var/tmp/" $ hPutStrLn stderr

main :: IO ()
main = getArgs >>= main'

main' :: [String] -> IO ()
main' ("1" : out : files) =
    withFile out WriteMode $ \hdl ->
        concatInputs files >=> run $ \hdr -> do
            pq <- progressNum "read" 1000000 (hPutStrLn stderr) =$
                  foldStreamM (\pq br -> enqueuePQ pqconf (BySelfPos br) pq) makePQ
            liftIO $ hPutStrLn stderr "Input phase done, writing output."
            enumPQ pq >=> run $
                progressNum "wrote" 1000000 (hPutStrLn stderr) =$
                writeBamHandle hdl hdr

-- This is ugly.  We should be able to rewind the temporary handle
-- instead of reopening the file.  Somehow, that yielded an empty stream
-- :-(  Not sure if any of this is worth the trouble, though.
main' ("2" : out : files) = do
    (hdr, pq) <- withFile "/var/tmp/tmp.bam" WriteMode $ \htmp ->
                     mergeInputs combineCoordinates files >=> run $ \hdr ->
                         (,) hdr <$> ( progressNum "read" 1000000 (hPutStrLn stderr)
                                    =$ isolate_sorted_part (writeBamHandle htmp hdr) )

    hPutStrLn stderr "Input phase done, writing output."
            -- hSeek htmp AbsoluteSeek 0
    withFile out WriteMode $ \hdl ->
        decodeAnyBamFile "/var/tmp/tmp.bam" >=> run $ \_ ->
        -- enumHandle defaultBufSize htmp >=> run $
            -- joinI $ decodeAnyBam $ \_ ->
            joinI $ eneePQ pq $
            joinI $ progressNum "wrote" 1000000 (hPutStrLn stderr) $
            writeBamHandle hdl hdr

main' _ = exitSuccess


enumPQ :: Monad m => PQ BySelfPos -> Enumerator [BamRaw] m b
enumPQ pq it = case viewMinPQ pq of
    Nothing -> return it
    Just (BySelfPos r, pq') -> enumPure1Chunk [r] it >>= enumPQ pq'

eneePQ :: Monad m => PQ BySelfPos -> Enumeratee [BamRaw] [BamRaw] m b
eneePQ pq0 it0 = tryHead >>= go pq0 it0
  where
    go pq it Nothing = lift $ enumPure1Chunk [r|BySelfPos r <- unfoldr viewMinPQ pq] it
    go pq it (Just a) = case viewMinPQ pq of
            Nothing -> lift (enumPure1Chunk [a] it) >>= mapChunks id
            Just (BySelfPos b,pq')
                | br_self_pos a <= br_self_pos b
                    -> lift (enumPure1Chunk [a] it) >>= eneePQ pq
                | otherwise
                    -> lift (enumPure1Chunk [b] it) >>= \it' -> go pq' it' (Just a)


-- How to make it opportunistic?  We need a PriorityQueue and a stream.
-- It's probably easiest if the stream is a temporary Bam file.  What
-- comes in goes into the stream unless its coordinate is smaller than
-- that of the previous record.  Else it goes to the PQ.
--
-- We take an 'Iteratee' and send the sorted subset to it.  Everything
-- else is accumulated in the PQ.  At EOF, we 'run' the 'Iteratee', but
-- return the PQ.

isolate_sorted_part :: MonadIO m => Iteratee [BamRaw] m () -> Iteratee [BamRaw] m (PQ BySelfPos)
isolate_sorted_part it0 = tryHead >>= ini
  where
    ini  Nothing  = makePQ <$ lift (run it0)
    ini (Just br) = do it' <- lift $ enumPure1Chunk [br] it0
                       (itZ, pqZ, _) <- foldStreamM go (it', makePQ, br_self_pos br)
                       lift $ run itZ
                       return pqZ

    go :: MonadIO m => (Iteratee [BamRaw] m (), PQ BySelfPos, (Refseq, Int))
       -> BamRaw ->  m (Iteratee [BamRaw] m (), PQ BySelfPos, (Refseq, Int))
    go (it, pq, pos) br
        | br_self_pos br >= pos = do
            it' <- enumPure1Chunk [br] it
            return (it', pq, br_self_pos br)
        | otherwise = do
            pq' <- liftIO $ enqueuePQ pqconf (BySelfPos br) pq
            return (it, pq', pos)



