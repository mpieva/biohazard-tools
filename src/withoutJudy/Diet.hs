module Diet (
    MDiet, IDiet, Word, newDiet, addI, lookupI,
    lookupIs, lookupICovered, lookupIPartial,
    lookupLeft, lookupRight, unsafeFreezeDiet,
    unions
) where

-- | "Discrete Interval Encoding Tree", but without the use of Judy.
-- The API looks as if the underlying data structure is mutable, but
-- only for compatibility reasons.  Compare the version using Judy to
-- make sense of this.

import Bio.Prelude       hiding ( loop, Word )

import qualified Data.IntMap as M
import qualified Data.IntSet as S

type Word = M.Key

-- | A Discrete Interval Encoding Tree.
newtype IDiet = IDiet (M.IntMap S.IntSet) deriving Show

-- | Mutable version of 'IDiet', see 'unsafeFreezeDiet'.
newtype MDiet = MDiet (IORef IDiet)

-- | Creates a new, empty DIET.
newDiet :: IO MDiet
newDiet = MDiet <$> newIORef (IDiet M.empty)

-- | Turns an 'MDiet' into an 'IDiet' without making a copy.  The
-- original 'MDiet' must not be changed afterwards.
unsafeFreezeDiet :: MDiet -> IO IDiet
unsafeFreezeDiet (MDiet v) = readIORef v

-- | Look an interval up.  This is equivalent to looking up each
-- contained position, but much more efficient.  A list of lists of
-- integers is returned, each corresponds to the stored integers for
-- some position in the query interval that differs from the previous
-- one.  Empty intervals (where 'begin' is greater or equal than 'end')
-- return the empty list.
--
-- Note use of 'lookupGT' (not 'lookupGE'): this gets whatever
-- follows(!) the start position, which is correct, since we store
-- values for intervals at the position past their end.

lookupIs :: Int -> Int -> IDiet -> [[Word]]
lookupIs begin end _ | begin >= end = []
lookupIs begin end (IDiet diet) = loop begin
  where
    loop !i =
        case M.lookupGT i diet of
            Nothing       -> []
            -- scan until we're past the end, but the last result is still returned!
            Just (i', xs) ->
                S.toAscList xs : if i' < fromIntegral end then loop i' else []

-- | Lookup of union of annotations.  Values mapped from any covered
-- position are returned.
lookupI :: Int -> Int -> IDiet -> [Word]
lookupI a b = unions . lookupIs a b

-- | Lookup intersection of annotations.  Only values mapped from every
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
addI begin end val (MDiet  ref)
    | begin == end = return ()
    | begin > end  = addI end begin val (MDiet ref)
    | otherwise    =
        modifyIORef ref $ \(IDiet d) ->
            IDiet .
            -- add interior annotations where necessary
            interior begin .
            -- make sure 'begin' and 'end' have an entry of their own.  If the
            -- entry is empty, copy the old annotation
            fillInAnno end .
            fillInAnno begin $ d
  where
    -- interior handling: add annotation with 'val' everywhere where
    -- some annotation already exists; 'lookupGT' for lookup ensures we
    -- don't overwrite the beginning
    interior :: M.Key -> M.IntMap S.IntSet -> M.IntMap S.IntSet
    interior !i !diet =
        case M.lookupGT i diet of
            Nothing         -> diet
            Just (i',set)
                | i > end   -> diet
                | otherwise ->
                    let !set' = S.insert val set
                    in interior i' (M.insert i' set' diet)

    -- Make sure a given position has an explicit annotation.  First
    -- lookup from the position (inclusive), if nothing is found or
    -- something is found, but not at the exact position, we insert an
    -- annotation (initially empty).  If something was found, it's then
    -- copied.
    fillInAnno :: M.Key -> M.IntMap S.IntSet -> M.IntMap S.IntSet
    fillInAnno k diet =
        case M.lookupGE k diet of
            Nothing         -> M.insert k S.empty diet
            Just (i,set)
                -- found something at the right location: all set
                | i == k    -> diet
                -- found something somewhere to the right:
                -- insert a copy at the exact location
                | otherwise -> M.insert k set diet


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
lookupLeft s (IDiet diet) =
    let open = case M.lookupGT s diet of
                    Nothing -> []
                    Just (_,set) -> S.toAscList set

    in case M.lookupLE s diet of
        Nothing      -> (  -s, [] )
        Just (i,set) -> ( i-s, S.toAscList set `setminus` open )


-- | Get the closest set of annotations to the right of an interval.  To
-- this end, we look up the first and second stored positions truly
-- right of the interval.  The first contains intervals already open in
-- the query, the second contains the interesting stuff.
lookupRight :: Int -> IDiet -> ( Int, [Word] )
lookupRight e (IDiet diet) =
    case M.lookupGE e diet of
        Nothing -> ( maxBound-e+1, [] )
        Just (k, open) ->
            case M.lookupGT k diet of
                Nothing -> ( maxBound-e+1, [] )
                Just (i, close) ->
                    ( i-e+1, S.toAscList close `setminus` S.toAscList open )
