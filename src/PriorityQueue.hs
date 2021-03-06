module PriorityQueue (
        Sized(..),
        PQ_Conf(..),

        PQ,
        makePQ,
        closePQ,
        enqueuePQ,
        viewMinPQ,
        sizePQ,
        listToPQ,
        pqToList,

        ByQName(..),
        byQName,
        ByMatePos(..),
        BySelfPos(..),
        br_qname,
        br_mate_pos,
        br_self_pos,
        br_copy
) where

import Bio.Prelude
import Bio.Bam                   hiding ( Stream )
import Data.Binary                      ( Binary, get, put )
import Data.Binary.Get                  ( runGetOrFail )
import Data.Binary.Put                  ( execPut )
import Foreign.C.Error                  ( throwErrnoIfMinus1 )
import Foreign.C.String                 ( withCAString )
import Foreign.C.Types                  ( CChar(..), CInt(..), CSize(..) )
import System.IO                        ( SeekMode(RelativeSeek) )

import qualified Data.ByteString                  as S
import qualified Data.ByteString.Builder.Internal as B
import qualified Data.ByteString.Internal         as S ( createAndTrim )
import qualified Data.ByteString.Unsafe           as S ( unsafeUseAsCStringLen )
import qualified Data.ByteString.Lazy             as L
import qualified Data.ByteString.Lazy.Internal    as L ( ByteString(..) )

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
-- buffering, too, this should be taken into account, at least
-- approximately.
--
-- Enqueued objects are kept in an in-memory heap until the memory
-- consumption becomes too high.  At that point, the whole heap is
-- sorted and dumped to external storage.  If necessary, the file to do
-- so is created and kept open.  The newly created stream is added to a
-- heap so that dequeueing objects amounts to performing a merge sort on
-- multiple external streams.  If too many streams are open, we do
-- exactly that:  merge-sort all streams in one generation into a single
-- stream in the next generation.  The discrete generations ensure that
-- merging handles streams of roughly equal length.  To conserve on file
-- descriptors, we concatenate multiple streams into a single file, then
-- use pread(2) on that as appropriate.  This way, we need only one file
-- descriptor per generation.
--
-- A generic way to decide when too many streams are open or how much
-- memory we can afford to use would be nice.  That first condition is
-- reached when seeking takes about as much time as reading (which
-- depends on buffer size and system characteristics), so that an
-- additional merge pass becomes economical.  The second condition not
-- only depends on how much memory we have, but also on how much we will
-- use for the merging logic.  Right now, both are just a parameter.

data PQ_Conf = PQ_Conf {
        max_mb :: Int,          -- ^ memory limit
        max_merge :: Int,       -- ^ fan-in limit
        temp_path :: FilePath,  -- ^ path to temporary files
        croak :: String -> IO () }

data PQ a = PQ { heap       :: !(SkewHeap a)        -- generation 0 (in RAM)
               , heap_size  :: !Int                 -- approx. size of gen 0
               , sizePQ     :: !Int                 -- number of items in total
               , spills     :: [( Fd, SkewHeap (Stream a) )] }
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

-- | Closes open file descriptors.  There is only one way file
-- descriptors become detectably unused, which is when we migrate one
-- generation of temporary streams to the next.  In this case, the 'Fd'
-- is closed.  Therefore, when we are done with a PQ, we have to close
-- these file handles manually.  (Naturally, the 'PQ' can no longer be
-- used after calling this.)
closePQ :: PQ a -> IO ()
closePQ = mapM_ (closeFd . fst) . spills

-- | Enqueues an element.
-- This operation may result in the creation of a file or in an enormous
-- merge of already created files.
enqueuePQ :: (Ord a, Sized a) => PQ_Conf -> a -> PQ a -> IO (PQ a)
enqueuePQ cfg a (PQ p s c ss) = do
    let !msz = 1000000 * max_mb cfg
        !s' = s + usedBytes a
        !pq' = PQ (insertH a p) s' (succ c) ss
    if s' >= msz
      then flushPQ cfg pq'
      else return pq'

listToPQ :: (Ord a, Sized a) => PQ_Conf -> [a] -> IO (PQ a)
listToPQ cfg = foldM (flip (enqueuePQ cfg)) makePQ


-- | Takes 99% of the in-memory portion of the PQ and flushes it to
-- generation one in external storage.

flushPQ :: (Ord a, Sized a) => PQ_Conf -> PQ a -> IO (PQ a)
flushPQ cfg PQ{..} = do
    croak cfg "Spilling memory to gen 0."
    let msz = 10000 * max_mb cfg
        (sheap, ssize, lheap) = splitH msz Empty 0 heap
    spills' <- spillTo cfg 0 (heapToList lheap) spills
    return $! PQ sheap ssize sizePQ spills'
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
-- After spilling, we check if the first generation overflowed.  If
-- so, we merge by spilling it to the tail of the list and recreating a
-- new, empty first generation.

spillTo :: (Ord a, Sized a) => PQ_Conf -> Int -> [a] -> [( Fd, SkewHeap (Stream a) )] -> IO [( Fd, SkewHeap (Stream a) )]
spillTo cfg g as [                         ] = mkEmptySpill cfg >>= spillTo cfg g as . (:[])
spillTo cfg g as (( fd, streams ) : spills ) = do
    str <- externalize fd . foldMap (execPut . put) $ as
    let streams' = insertH (decodeStream str) streams
    pos <- fdSeek fd RelativeSeek 0
    croak cfg $ "Gen " ++ shows g " has " ++ shows (sizeH streams') " streams (fd "
             ++ shows fd ", " ++ shows (pos `shiftR` 20) "MB)."
    if sizeH streams' == max_merge cfg
      then do croak cfg $ "Spilling gen " ++ shows g " to gen " ++ shows (succ g) "."
              (:) <$> mkEmptySpill cfg
                  <*> spillTo cfg (succ g) (unfoldr viewMinS streams') spills
                  <*  closeFd fd
      else return (( fd, streams' ) : spills )

mkEmptySpill :: PQ_Conf -> IO ( Fd, SkewHeap (Stream a) )
mkEmptySpill cfg = do pn <- getProgName
                      fd <- withCAString (temp_path cfg ++ "/" ++ pn ++ "-XXXXXX") $
                                \p -> throwErrnoIfMinus1 "mkstemp" (mkstemp p) <* unlink p
                      return ( fd, Empty )


-- | Creates a disk backed stream from a heap.  The heap is unravelled
-- in ascending order and written to a new segment in a temporary file.
-- That segment is then converted back to a lazy bytestring, which is
-- returned.  (Yes, evil lazy IO.)
externalize :: Fd -> B.Builder -> IO L.ByteString
externalize fd s = do
        pos0 <- fdSeek fd RelativeSeek 0
        fdPutLazy fd $ lz4 $ B.toLazyByteStringWith
            (B.untrimmedStrategy lz4_bsize lz4_bsize) L.empty s
        pos1 <- fdSeek fd RelativeSeek 0
        unLz4 <$> fdGetLazy fd pos0 pos1

fdGetLazy :: Fd -> FileOffset -> FileOffset -> IO L.ByteString
fdGetLazy fd p0 p1
    | p0 == p1 = return L.empty
    | otherwise = do
        let len = (p1-p0) `min` (1024*1024)
        str <- S.createAndTrim (fromIntegral len) $ \p ->
                    fmap fromIntegral $
                        throwErrnoIfMinus1 "pread" $
                            pread fd p (fromIntegral len) p0
        L.Chunk str <$>
            unsafeInterleaveIO (fdGetLazy fd (p0 + fromIntegral (S.length str)) p1)

foreign import ccall unsafe "stdlib.h mkstemp" mkstemp :: Ptr CChar -> IO Fd
foreign import ccall unsafe "unistd.h unlink"  unlink  :: Ptr CChar -> IO CInt
foreign import ccall unsafe "unistd.h pread"   pread   :: Fd -> Ptr a -> CSize -> FileOffset -> IO CSsize


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

viewMinPQ :: (Ord a, Sized a) => PQ a -> Maybe (a, PQ a)
viewMinPQ PQ{..} =
    case viewMinL spills of

        Just (a, ss') -> case viewMin heap of

            Just (b, h') | b <= a ->
                Just (b, PQ h' (heap_size - usedBytes b) (pred sizePQ) spills)

            _ -> Just (a, PQ heap heap_size (pred sizePQ) ss')

        Nothing -> case viewMin heap of

            Just (b ,h') ->
                Just (b, PQ h' (heap_size - usedBytes b) (pred sizePQ) spills)

            Nothing -> Nothing

viewMinL :: Ord a => [( Fd, SkewHeap (Stream a) )] -> Maybe (a, [( Fd, SkewHeap (Stream a) )])
viewMinL = minimumM . each
  where
    each [              ] = [ ]
    each ( (fd,hp) : xs ) =
        case viewMinS hp of
            Nothing       -> k
            Just (a, hp') -> (a, (fd, hp') : xs) : k
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

viewMinS :: Ord a => SkewHeap (Stream a) -> Maybe (a, SkewHeap (Stream a))
viewMinS h = case viewMin h of
    Nothing               -> Nothing
    Just (Nil       ,  _) -> error "WTF?!"
    Just (Cons a Nil, h') -> Just (a, h')
    Just (Cons a   s, h') -> Just (a, insertH s h')

sizeH :: SkewHeap a -> Int
sizeH  Empty       = 0
sizeH (Node _ l r) = sizeH l + sizeH r + 1

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


lz4_bsize, lz4_csize :: Int
lz4_bsize = 16*1024*1024 - 128
lz4_csize = 16*1000*1000


-- We just compress and frame each 'Chunk' individually.
-- 'toLazyByteStringWith' should give us 'Chunk's of reasonable size;
-- but we enforce a maximum size to make decompression easier.
lz4 :: L.ByteString -> L.ByteString
lz4 = L.foldrChunks go L.empty
  where
    go s | S.null s               = id
         | S.length s > lz4_csize = go (S.take lz4_csize s) . go (S.drop lz4_csize s)
         | otherwise              = L.Chunk . unsafePerformIO $
                S.createAndTrim lz4_bsize $ \pdest ->
                    S.unsafeUseAsCStringLen s $ \(psrc,srcLen) -> do
                        dLen <- c_lz4 psrc (plusPtr pdest 4)
                                      (fromIntegral srcLen) (fromIntegral lz4_bsize)
                        -- place compressed size at beginning (slow and portable)
                        pokeByteOff pdest 0 (fromIntegral (dLen `shiftR`  0) :: Word8)
                        pokeByteOff pdest 1 (fromIntegral (dLen `shiftR`  8) :: Word8)
                        pokeByteOff pdest 2 (fromIntegral (dLen `shiftR` 16) :: Word8)
                        pokeByteOff pdest 3 (fromIntegral (dLen `shiftR` 24) :: Word8)
                        return $ fromIntegral dLen + 4

foreign import ccall unsafe "lz4.h LZ4_compress_default"
    c_lz4 :: Ptr CChar  -- source
          -> Ptr CChar  -- dest
          -> CInt       -- sourceSize
          -> CInt       -- maxDestSize
          -> IO CInt    -- number of bytes written


unLz4 :: L.ByteString -> L.ByteString
unLz4 s | L.null  s = L.empty
        | otherwise = L.Chunk ck' . unLz4 $ L.drop (l+4) s
  where
    -- size (slow and portable)
    l   = fromIntegral (L.index s 0) `shiftL`  0 .|.
          fromIntegral (L.index s 1) `shiftL`  8 .|.
          fromIntegral (L.index s 2) `shiftL` 16 .|.
          fromIntegral (L.index s 3) `shiftL` 24

    ck  = L.toStrict $ L.take (l+4) s
    ck' = unsafePerformIO $
             S.createAndTrim lz4_bsize $ \pdst ->
                S.unsafeUseAsCStringLen ck $ \(psrc,srcLen) ->
                    fromIntegral <$> c_unlz4 (plusPtr psrc 4) pdst
                            (fromIntegral srcLen - 4) (fromIntegral lz4_bsize)

foreign import ccall unsafe "lz4.h LZ4_decompress_safe"
    c_unlz4 :: Ptr CChar  -- source
            -> Ptr Word8  -- dest
            -> CInt       -- sourceSize
            -> CInt       -- maxDestSize
            -> IO CInt    -- number of bytes written


data ByQName = ByQName { _bq_hash :: !Int
                       , _bq_alnid :: !Int
                       , _bq_rec :: !BamRaw }

-- Note on XI:  this is *not* the index read (MPI EVAN convention), but
-- the number of the alignment (Segemehl convention, I believe).
-- Fortunately, the index is never numeric, so this is reasonably safe.
byQName :: BamRaw -> ByQName
byQName b = ByQName (hash $ br_qname b) (extAsInt 0 "XI" $ unpackBam b) b

instance Eq ByQName where
    ByQName ah ai a == ByQName bh bi b =
        (ah, ai, br_qname a) == (bh, bi, br_qname b)

instance Ord ByQName where
    ByQName ah ai a `compare` ByQName bh bi b =
        (ah, ai, b_qname (unpackBam a)) `compare` (bh, bi, b_qname (unpackBam b))

newtype ByMatePos = ByMatePos BamRaw

instance Eq ByMatePos where
    ByMatePos a == ByMatePos b =
        br_mate_pos a == br_mate_pos b

instance Ord ByMatePos where
    ByMatePos a `compare` ByMatePos b =
        br_mate_pos a `compare` br_mate_pos b

newtype BySelfPos = BySelfPos BamRaw

instance Eq BySelfPos where
    BySelfPos a == BySelfPos b =
        br_self_pos a == br_self_pos b

instance Ord BySelfPos where
    BySelfPos a `compare` BySelfPos b =
        br_self_pos a `compare` br_self_pos b

instance Binary ByQName   where put (ByQName _ _ r) = put (raw_data r) ; get = byQName   . bamRaw 0 <$> get
instance Binary ByMatePos where put (ByMatePos   r) = put (raw_data r) ; get = ByMatePos . bamRaw 0 <$> get
instance Binary BySelfPos where put (BySelfPos   r) = put (raw_data r) ; get = BySelfPos . bamRaw 0 <$> get

instance Sized ByQName   where usedBytes (ByQName _ _ r) = S.length (raw_data r) + 80
instance Sized ByMatePos where usedBytes (ByMatePos   r) = S.length (raw_data r) + 64
instance Sized BySelfPos where usedBytes (BySelfPos   r) = S.length (raw_data r) + 64

br_qname :: BamRaw -> Seqid
br_qname = b_qname . unpackBam

br_mate_pos :: BamRaw -> (Refseq, Int)
br_mate_pos = (b_mrnm &&& b_mpos) . unpackBam

br_self_pos :: BamRaw -> (Refseq, Int)
br_self_pos = (b_rname &&& b_pos) . unpackBam

br_copy :: BamRaw -> BamRaw
br_copy br = bamRaw (virt_offset br) $! S.copy (raw_data br)

