{-# LANGUAGE BangPatterns, RecordWildCards, NoImplicitPrelude #-}
module PriorityQueue (
        Sized(..),
        PQ_Conf(..),
        MonadPQ,

        PQ,
        makePQ,
        enqueuePQ,
        viewMinPQ,
        sizePQ,
        listToPQ,
        pqToList
) where

import Bio.Prelude
import Codec.Compression.Snappy         ( compress, decompress )
import Control.Monad.Trans.Class        ( lift )
import Control.Monad.Trans.State hiding ( get, put )
import Data.Binary                      ( Binary, Put, get, put )
import Data.Binary.Get                  ( runGetOrFail )
import Data.Binary.Put                  ( runPut )
import Foreign.C.Error                  ( throwErrnoIfMinus1 )
import Foreign.C.String                 ( withCAString )
import Foreign.C.Types                  ( CChar(..), CInt(..), CSize(..) )
import System.IO

import qualified Data.ByteString                as S
import qualified Data.ByteString.Internal       as S
import qualified Data.ByteString.Lazy           as L
import qualified Data.ByteString.Lazy.Internal  as L

-- | A Priority Queue that can fall back to external storage.
--
-- Note that such a Priority Queue automatically gives rise to an
-- external sorting algorithm:  enqueue everything, dequeue until empty.
--
-- Whatever is to be stored in this queue needs to be in 'Binary',
-- because it may need to be moved to external storage on demand.  We
-- also need a way to estimate the memory consumption of an enqueued
-- object.  When constructing the queue, the maximum amount of RAM to
-- consume is set.  Note that open input streams use memory for
-- buffering, too.
--
-- Enqueued objects are kept in an in memory heap until the memory
-- consumption becomes too high.  At that point, the whole heap is
-- sorted and dumped to external storage.  If necessary, the file to do
-- so is created and kept open.  The newly created stream is added to a
-- heap so that dequeueing objects amounts to performing a merge sort on
-- multiple external streams.  To conserve on file descriptors, we
-- concatenate multiple streams into a single file, then use pread(2) on
-- that as appropriate.  If too many streams are open (how do we set
-- that limit?), we do exactly that:  merge-sort all streams and the
-- in-memory heap into a single new stream.  One file is created for
-- each generation of streams, so that mergind handles streams of
-- roughly equal length.
--
--
-- XXX  Need a way to decide when too many streams are open.  That point
--      is reached when seeking takes about as much time as reading
--      (which depends on buffer size and system characteristics), so
--      that an additional merge pass becomes economical.
--
-- XXX  Need to fClose everything at some point.  No GC for 'Fd's!

data PQ_Conf = PQ_Conf {
        max_mb :: Int,          -- ^ memory limit
        temp_path :: FilePath } -- ^ path to temporary files
        -- functions to report progress go here
        -- limit on streams per gen goes here

type MonadPQ = StateT PQ_Conf IO

data PQ a = PQ { heap       :: !(SkewHeap a)        -- generation 0 (in RAM)
               , heap_size  :: !Int                 -- approx. size of gen 0
               , sizePQ     :: !Int                 -- number of items in total
               , spills     :: [( Fd, SkewHeap (Stream a) )] }   -- gen 1 (XXX)
  deriving Show

-- | Things that have a (conservative) size estimate.  Doesn't need to
-- be accurate, it's used to decide when to spill to external memory.
class Binary a => Sized a where usedBytes :: a -> Int

-- For testing purposes.
instance Sized Int where usedBytes _ = 8

-- | Creates a priority queue.  Note that the priority queue creates
-- files, which will only be cleaned up if deletePQ is called.
makePQ :: (Ord a, Sized a) => PQ a
makePQ = PQ Empty 0 0 []
{- makePQ = gets temp_path >>= \path -> lift $ do
            pn <- getProgName
            PQ Empty 0 0 [] <$>
                withCAString (path ++ "/" ++ pn ++ "-XXXXXX") (\p ->
                    throwErrnoIfMinus1 "mkstemp" (mkstemp p) <* unlink p) -}

foreign import ccall unsafe "stdlib.h mkstemp" mkstemp :: Ptr CChar -> IO Fd
foreign import ccall unsafe "unistd.h unlink"  unlink  :: Ptr CChar -> IO CInt

-- | Enqueues an element.
-- This operation may result in the creation of a file or in an enormous
-- merge of already created files.
enqueuePQ :: (Ord a, Sized a) => a -> PQ a -> MonadPQ (PQ a)
enqueuePQ a (PQ p s c ss) = do
    msz <- (*) 1000000 <$> gets max_mb
    let !s' = s + usedBytes a
        !pq' = PQ (insertH a p) s' (succ c) ss
    if s' >= msz
      then flushPQ pq'
      else return pq'

listToPQ :: (Ord a, Sized a) => [a] -> MonadPQ (PQ a)
listToPQ = foldM (flip enqueuePQ) makePQ


-- | Takes 99% of the in-memory portion of the PQ and flushes it to
-- generation one in external storage.

flushPQ :: (Ord a, Sized a) => PQ a -> MonadPQ (PQ a)
flushPQ PQ{..} = do
    msz <- (*) 10000 <$> gets max_mb
    let (sheap, ssize, lheap) = splitH msz Empty 0 heap
    spills' <- spillTo (heapToList lheap) spills
    -- str <- lift . externalize extmem . runPut . foldMap put . heapToList $ lheap
    return $! PQ sheap ssize sizePQ spills' -- (insertH (decodeStream str) streams) extmem
  where
    -- We take up to 1% of the maximum size off generation 0 and make that the
    -- new generation 0.  The remainder is flushed.  The idea is that the next
    -- items may be extracted soon, so it wouldn't pay to flush them.
    splitH msz sh ssz h = case viewMin h of
        Just (b,h') | usedBytes b + ssz < msz ->
            splitH msz (insertH b sh) (usedBytes b + ssz) h'
        _ -> (sh, ssz, h)

-- Spill a list of thingies to external storage.  A list of generations
-- of external spills is supplied.  We always spill to the first
-- generation by appending to the 'Fd' and creating a new stream that is
-- added to the 'SkewHeap'.  If the list is empty, we have to create the
-- first generation spill.
--
-- XXX  After spilling, we check if the first generation overflowed.  If
-- so, we merge by spilling it to the tail of the list and recreating a
-- new, empty first generation.  (Have to remember to close the 'Fd's.)

spillTo :: (Ord a, Sized a) => [a] -> [( Fd, SkewHeap (Stream a) )] -> MonadPQ [( Fd, SkewHeap (Stream a) )]
spillTo as [                         ] = mkEmptySpill >>= spillTo as . (:[])
spillTo as (( fd, streams ) : spills ) = do
    str <- lift . externalize fd . runPut . foldMap put $ as
    return (( fd, insertH (decodeStream str) streams ) : spills )

mkEmptySpill :: MonadPQ ( Fd, SkewHeap (Stream a) )
mkEmptySpill = do path <- gets temp_path
                  fd <- lift $ do pn <- getProgName
                                  withCAString (path ++ "/" ++ pn ++ "-XXXXXX")
                                       (\p -> throwErrnoIfMinus1 "mkstemp" (mkstemp p) <* unlink p)
                  return ( fd, Empty )



-- | Creates a disk backed stream from a heap.  The heap is unravelled
-- in ascending order and written to a new segment in a temporary file.
-- That segment is then converted back to a lazy bytestring, which is
-- returned.  (Yes, evil lazy IO.)
externalize :: Fd -> L.ByteString -> IO L.ByteString
externalize fd s = do
        pos0 <- fdSeek fd RelativeSeek 0
        fdPutLazy fd $ snappy s
        pos1 <- fdSeek fd RelativeSeek 0
        unsnappy <$> fdGetLazy fd pos0 pos1

fdGetLazy :: Fd -> FileOffset -> FileOffset -> IO L.ByteString
fdGetLazy fd p0 p1
    | p0 == p1 = return L.empty
    | otherwise = do
        let len = fromIntegral (p1-p0) `min` 32768
        buf <- mallocForeignPtrBytes len
        len' <- withForeignPtr buf $ \p ->
            throwErrnoIfMinus1 "pread" $
                pread fd p (fromIntegral len) p0
        L.Chunk (S.PS buf 0 $ fromIntegral len') <$>
            unsafeInterleaveIO (fdGetLazy fd (p0 + fromIntegral len') p1)

foreign import ccall unsafe "unistd.h pread"
    pread :: Fd -> Ptr a -> CSize -> FileOffset -> IO CSsize


testList :: Ord a => [a] -> [a]
testList (a:b:cs) | a <= b    = a : testList (b:cs)
                  | otherwise = error "sorting violated?!"
testList [a] = [a]
testList [ ] = [ ]

decodeStream :: Binary a => L.ByteString -> Stream a
decodeStream = go
  where
    go s | L.null s                     = Nil
    go s = case runGetOrFail get s of
        Left ( _rest, _consumed, err ) -> error err
        Right ( rest, _consumed,  !a ) -> Cons a (decodeStream rest)

-- XXX Borked!!
viewMinPQ :: (Ord a, Sized a, Show a) => PQ a -> Maybe (a, PQ a)
viewMinPQ PQ{..} =
    case viewMinL spills of

        Just (a, ss') -> case viewMin heap of

            Just (b, h') | b <= a -> do
                Just (b, PQ h' (heap_size - usedBytes b) (pred sizePQ) spills)

            _ -> Just (a, PQ heap heap_size (pred sizePQ) ss')
                    -- (case s' of Nil      -> ss'
                    --            Cons _ _ -> insertH s' ss'))

        -- Just (Nil, _) -> error "WTF?! (2)"

        Nothing -> case viewMin heap of

            Just (b ,h') -> do
                Just (b, PQ h' (heap_size - usedBytes b) (pred sizePQ) spills)

            Nothing -> Nothing

viewMinL :: Ord a => [( Fd, SkewHeap (Stream a) )] -> Maybe (a, [( Fd, SkewHeap (Stream a) )])
viewMinL = minimumM . each
  where
    each [              ] = [ ]
    each ( (fd,hp) : xs ) =
        case viewMin hp of
            Nothing                -> k
            Just (Cons a Nil, hp') -> (a, (fd,          hp') : xs) : k
            Just (Cons a   s, hp') -> (a, (fd,insertH s hp') : xs) : k
      where
        k = [ (a, (fd,hp):xs') | (a,xs') <- each xs ]

    minimumM [] = Nothing
    minimumM xs = Just $ minimumBy (comparing fst) xs


pqToList :: (Ord a, Sized a, Show a) => PQ a -> [a]
pqToList = testList . unfoldr viewMinPQ

-- We need an in-memory priority queue.  Here's a skew heap.
data SkewHeap a = Empty | Node a (SkewHeap a) (SkewHeap a) deriving Show

singleton :: a -> SkewHeap a
singleton x = Node x Empty Empty

unionH :: Ord a => SkewHeap a -> SkewHeap a -> SkewHeap a
Empty              `unionH` t2                 = t2
t1                 `unionH` Empty              = t1
t1@(Node x1 l1 r1) `unionH` t2@(Node x2 l2 r2)
   | x1 <= x2                                 = Node x1 (t2 `unionH` r1) l1
   | otherwise                                = Node x2 (t1 `unionH` r2) l2

insertH :: Ord a => a -> SkewHeap a -> SkewHeap a
insertH x heap = singleton x `unionH` heap

viewMin :: Ord a => SkewHeap a -> Maybe (a, SkewHeap a)
viewMin Empty        = Nothing
viewMin (Node x l r) = Just (x, unionH l r)

sizeH :: SkewHeap a -> Int
sizeH  Empty       = 0
sizeH (Node _ l r) = sizeH l + sizeH r

heapToList :: (Binary a, Ord a) => SkewHeap a -> [a]
heapToList = unfoldr viewMin

data Stream a = Nil | Cons !a (Stream a) deriving Show

-- Streams are ordered by looking at just the first item.
instance Eq a => Eq (Stream a) where
    Cons a _ == Cons b _ = a == b
    Nil      == Nil      = True
    _        == _        = False

instance Ord a => Ord (Stream a) where
    Nil      `compare` Nil      = EQ
    Nil      `compare` Cons _ _ = LT
    Cons _ _ `compare` Nil      = GT
    Cons a _ `compare` Cons b _ = compare a b


-- This could be more efficient, but that requires access to the
-- internals of snappy.  One Chunk per compressed segment would be
-- ideal.
snappy :: L.ByteString -> L.ByteString
snappy = L.concat . map snap . L.toChunks
  where
    snap s | S.null s = L.empty
    snap s = let t = compress $ S.take 0xffff s
             in fromIntegral (S.length t .&. 0xff) `L.cons`
                fromIntegral (S.length t `shiftR` 8) `L.cons`
                L.fromStrict t `L.append` snap (S.drop 0xffff s)

unsnappy :: L.ByteString -> L.ByteString
unsnappy s | L.null s = L.empty
unsnappy s = let l = fromIntegral (L.index s 0) .|.
                     fromIntegral (L.index s 1) `shiftL` 8
             in L.Chunk (decompress (L.toStrict (L.take l (L.drop 2 s))))
                        (unsnappy (L.drop (l+2) s))



