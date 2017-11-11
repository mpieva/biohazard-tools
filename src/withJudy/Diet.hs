{-# LANGUAGE ForeignFunctionInterface, EmptyDataDecls #-}
module Diet (
    MDiet, IDiet, Word, newDiet, addI, lookupI,
    lookupIs, lookupICovered, lookupIPartial,
    lookupLeft, lookupRight, unsafeFreezeDiet,
    unions
) where

import Bio.Prelude       hiding ( pi, Word, loop )
import Foreign.C.Types          ( CULong(..), CInt(..) )
import Foreign.ForeignPtr       ( newForeignPtr, withForeignPtr, ForeignPtr )
import Foreign.Marshal.Alloc    ( malloc )
import Foreign.Marshal.Utils    ( with )
import Foreign.Ptr              ( Ptr, FunPtr, nullPtr )
import Foreign.Storable         ( peek, poke )

-- | "Discrete Interval Encoding Tree:"
-- Stores a mapping from Int to a set of Ints.  Consecutive intervals
-- with identical values are stored compactly as follows:  We use a
-- JudyL map, and for each contiguous interval, store the mapped value
-- for the key just outside its right edge and store the old value for
-- the key at its left edge.  Lookup at some position is done by looking
-- up the smallest key bigger(!) than the interesting position.  Lookup
-- for an interval is done the same way, values are concatenated.  If
-- lookup fails, the result is an empty list.  To store lists, we
-- actually store a Judy1 array, which is converted to a list on demand.
--
-- Note:  There are data structures specialized to the storage of
-- regions, e.g. priority search queues.  However, these are only
-- advantageous if regions tend to overlap a lot, and we don't expect
-- that to happen.  Instead, we choose to benefit from the substantial
-- tuning that went into Judy.

-- This is the immutable version of the data structure, see 'MDiet' and
-- 'unsafeFreezeDiet'.

newtype IDiet = IDiet FPJudyL deriving Show

-- | Mutable version of 'IDiet', see 'unsafeFreezeDiet'.
newtype MDiet = MDiet FPJudyL deriving Show

data JudyLTag
type JudyL = Ptr JudyLTag
type PJudyL = Ptr JudyL
type FPJudyL = ForeignPtr JudyL

data Judy1Tag
type Judy1 = Ptr Judy1Tag
type PJudy1 = Ptr Judy1

type PError = Ptr ()

type Word = CULong
type PWord = Ptr CULong

foreign import ccall unsafe "&finalize_diet" finalizeDiet :: FunPtr (PJudyL -> IO())
foreign import ccall unsafe "Judy.h JudyLFirst" judyLFirst ::  JudyL -> PWord -> PError -> IO PJudy1
foreign import ccall unsafe "Judy.h JudyLNext"  judyLNext  ::  JudyL -> PWord -> PError -> IO PJudy1
foreign import ccall unsafe "Judy.h JudyLPrev"  judyLPrev  ::  JudyL -> PWord -> PError -> IO PJudy1
foreign import ccall unsafe "Judy.h JudyLIns"   judyLIns   :: PJudyL ->  Word -> PError -> IO PJudy1

foreign import ccall unsafe "Judy.h Judy1First" judy1First ::  Judy1 -> PWord -> PError -> IO CInt
foreign import ccall unsafe "Judy.h Judy1Next"  judy1Next  ::  Judy1 -> PWord -> PError -> IO CInt
foreign import ccall unsafe "Judy.h Judy1Set"   judy1Set   :: PJudy1 ->  Word -> PError -> IO CInt

-- | Creates a new, empty DIET.
newDiet :: IO MDiet
newDiet = do p <- malloc
             poke p nullPtr
             MDiet `fmap` newForeignPtr finalizeDiet p

-- | Turns an 'MDiet' into an 'IDiet' without making a copy.  The
-- original 'MDiet' must not be changed afterwards.
unsafeFreezeDiet :: MDiet -> IO IDiet
unsafeFreezeDiet (MDiet d) = return $! IDiet d


-- | Reads a list of integers from a pointer to a 'Judy1' array.  If the
-- pointer is null, an empty list is returned.  This function is
-- dangerous if the array is modified while it is traversed lazily.  No
-- crashed should result, but completely unpredictable results.
pJudy1ToLazyAscList :: ForeignPtr a -> PJudy1 -> [Word]
pJudy1ToLazyAscList  _ pj1 | pj1 == nullPtr = []
pJudy1ToLazyAscList fp pj1 = judy1ToLazyAscList fp $ unsafeDupablePerformIO (peek pj1)

judy1ToLazyAscList :: ForeignPtr a -> Judy1 -> [Word]
judy1ToLazyAscList fp j1 = loop 0 judy1First
  where
    loop i act = unsafeDupablePerformIO $
        withForeignPtr fp $ \_ ->
        with i $ \pi -> do
        rc <- act j1 pi nullPtr
        if rc == 0
            then return []
            else do i' <- peek pi
                    return $ i' : loop i' judy1Next


-- | Look an interval up.  This is equivalent to looking up each
-- contained position, but much more efficient.  A list of lists of
-- integers is returned, each corresponds to the stored integers for
-- some position in the query interval that differs from the previous
-- one.  Empty intervals (where 'begin' is greater or equal than 'end')
-- return the empty list.
--
-- Note use of 'judyLNext' (not 'judyLFirst'): this gets whatever
-- follows(!) the start position, which is correct, since we store
-- values for intervals at the position past their end.

lookupIs :: Int -> Int -> IDiet -> [[Word]]
lookupIs begin end _ | begin >= end = []
lookupIs begin end (IDiet fpjl) = loop (fromIntegral begin)
  where
    loop i = unsafeDupablePerformIO $
        withForeignPtr fpjl $ \pjl ->
        with i $ \pi -> do
        jl <- peek pjl
        pj1 <- judyLNext jl pi nullPtr
        if pj1 == nullPtr
            then return []
            -- scan until we're past the end, but the last result is still returned!
            else do let xs = pJudy1ToLazyAscList fpjl pj1
                    i' <- peek pi
                    if i' < fromIntegral end then return $ xs : loop i'
                                             else return [xs]

-- | Lookup of union of annotations.  Values mapped from any covered
-- position are returned.
lookupI :: Int -> Int -> IDiet -> [Word]
lookupI a b = unions . lookupIs a b

-- | Lookup intersection of annotations.  Only value mapped from every
-- covered position are returned.
lookupICovered :: Int -> Int -> IDiet -> [Word]
lookupICovered a b = intersections . lookupIs a b

-- | Lookup partial annotations.  Only values mapped from some, but not
-- all of the covered positions are returned.
lookupIPartial :: Int -> Int -> IDiet -> [Word]
lookupIPartial a b =  edges_only . lookupIs a b
  where edges_only xs = unions xs `setminus` intersections xs

-- | Add an annotation for an interval to an existing map.
addI :: Int -> Int -> Word -> MDiet -> IO ()
addI begin end val (MDiet fpjl) | begin == end = return ()
                                | begin > end  = addI end begin val (MDiet fpjl)
                                | otherwise    = withForeignPtr fpjl $ \pjl -> do
    -- make sure 'begin' and 'end' have an entry of their own.  If the
    -- entry is empty, copy the old annotation; in the case of end, the
    -- new annotation gets added later.
    fillInAnno begin
    fillInAnno end

    -- interior handling: add annotation everywhere where some
    -- annotation already exists; judyLNext for lookup ensures we don't
    -- overwrite the beginning
    jl <- peek pjl
    with (fromIntegral begin) $ \pi ->
         let loop = do pj1 <- judyLNext jl pi nullPtr
                       when (pj1 /= nullPtr) $ do
                            i <- peek pi
                            when (i <= fromIntegral end) $
                                judy1Set pj1 val nullPtr >> loop
         in loop

  where
    -- Make sure a given position has an explicit annotation.  First
    -- lookup from the position (inclusive), if nothing is found or
    -- something is found, but not at the exact position, we insert an
    -- annotation (initially empty).  If something was found, it's then
    -- copied.
    fillInAnno k = withForeignPtr fpjl $ \pjl ->
                   with (fromIntegral k) $ \pi -> do
        pp <- peek pjl >>= \jl -> judyLFirst jl pi nullPtr
        if pp == nullPtr
          -- found nothing at all: insert an empty set
          then do pj1 <- judyLIns pjl (fromIntegral k) nullPtr
                  poke pj1 nullPtr
          else do
            i <- peek pi
            -- found something, but not at the exact location
            -- insert an empty set, then lookup the other one again
            -- (necessary, because the pointer probably changed!) and
            -- copy it
            when (i /= fromIntegral k) $ do
                j1 <- peek pp
                pnew <- judyLIns pjl (fromIntegral k) nullPtr
                poke pnew nullPtr
                judy1Copy j1 pnew

judy1Copy :: Judy1 -> PJudy1 -> IO ()
judy1Copy j1 pj2 =
    with 0 $ \pi -> do
    let loop act = do rc <- act j1 pi nullPtr
                      if rc == 0
                          then return ()
                          else do i <- peek pi
                                  _ <- judy1Set pj2 i nullPtr
                                  loop judy1Next
    loop judy1First



-- | Union of many sets with sets being represented by sorted lists.
unions :: [[Word]] -> [Word]
unions = foldr (~~) []
  where
    [] ~~ xs = xs
    xs ~~ [] = xs
    (x:xs) ~~ (y:ys) = case compare x y of
        LT -> x : (   xs  ~~ (y:ys))
        EQ -> x : (   xs  ~~    ys )
        GT -> y : ((x:xs) ~~    ys )

-- | Intersection of many sets with sets being represented by sorted
-- lists.
intersections :: [[Word]] -> [Word]
intersections [    ] = []
intersections (a:as) = foldr (~~) a as
  where
    [] ~~ _ = []
    _ ~~ [] = []
    (x:xs) ~~ (y:ys) = case compare x y of
        LT ->         xs  ~~ (y:ys)
        EQ -> x : (   xs  ~~    ys )
        GT ->      (x:xs) ~~    ys

-- | Difference between sets represented as sorted lists.  @setminus as
-- bs@ returns every element that is in @as@ but not in @bs@.
setminus :: [Word] -> [Word] -> [Word]
[] `setminus`  _ = []
xs `setminus` [] = xs
(x:xs) `setminus` (y:ys) = case compare x y of
    LT -> x :    xs  `setminus` (y:ys)
    EQ ->        xs  `setminus`    ys
    GT ->     (x:xs) `setminus`    ys


-- | Get the closest set of annotations to the left of an interval.  We
-- look at the start and to its left.  Returns the distance and the list
-- of annotations.
lookupLeft :: Int -> IDiet -> ( Int, [Word] )
lookupLeft s (IDiet fpjl) = unsafeDupablePerformIO $
    withForeignPtr fpjl $ \pjl ->
    with (fromIntegral s) $ \pi -> do
    jl <- peek pjl

    pj_r <- judyLNext jl pi nullPtr
    pj_l <- judyLPrev jl pi nullPtr
    i <- fromIntegral `fmap` peek pi

    let open  = pJudy1ToLazyAscList fpjl pj_r
        close = pJudy1ToLazyAscList fpjl pj_l
    return ( i-s, close `setminus` open )


-- | Get the closest set of annotations to the right of an interval.  To
-- this end, we look up the first and second stored positions truly
-- right of the interval.  The first contains intervals already open in
-- the query, the second contains the interesting stuff.
lookupRight :: Int -> IDiet -> ( Int, [Word] )
lookupRight e (IDiet fpjl) = unsafeDupablePerformIO $
    withForeignPtr fpjl $ \pjl ->
    with (fromIntegral e) $ \pi -> do
    jl <- peek pjl
    pj1 <- judyLFirst jl pi nullPtr
    pj2 <- judyLNext  jl pi nullPtr
    i <- fromIntegral `fmap` peek pi

    let open  = pJudy1ToLazyAscList fpjl pj1
        close = pJudy1ToLazyAscList fpjl pj2
    return ( i-e+1, close `setminus` open )

