{-# LANGUAGE Rank2Types #-}
module Symtab where

import qualified Data.HashMap.Strict                as M
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.ByteString.Char8              as S
import qualified Data.Vector                        as V

import Bio.Prelude
import Control.Monad.IO.Class
import Diet

type Gene       = S.ByteString
type Chrom      = L.ByteString
type ChromTable = M.HashMap L.ByteString (Chrom, Int)
type Start      = Int
type End        = Int
type Region     = ( L.ByteString, Chrom, Senses, Start, End )

data Senses = None | Forward | Reverse | Both deriving ( Show, Enum )

withSenses :: Senses -> ((a -> Either a a) -> b) -> [b]
withSenses None    _ = []
withSenses Forward k = [k Right]
withSenses Reverse k = [k Left]
withSenses Both    k = [k Left, k Right]

-- We use `Right` for forward, and `Left` for reverse,
-- just like when driving on roads.
type MAnnotab   = M.HashMap (Either Chrom Chrom) MDiet
type IAnnotab   = M.HashMap (Either Chrom Chrom) IDiet
type Symtab     = M.HashMap Gene Int
type RevSymtab  = V.Vector Gene

findSymbol :: Bytes -> CPS Int
findSymbol s = do
    let !k = S.copy s
    t <- get_syms
    case M.lookup k t of
        Just x  -> return x
        Nothing -> do let !l = M.size t
                      modify_syms $ M.insert k l
                      return l

findDiet :: Either Chrom Chrom -> CPS MDiet
findDiet k = do
    m <- get_anno
    case M.lookup k m of
        Just d  -> return d
        Nothing -> do d <- liftIO newDiet
                      modify_anno $ M.insert k d
                      return d


invertTab :: Symtab -> RevSymtab
invertTab t = V.replicate (M.size t) "" V.// [ (y,x) | (x,y ) <- M.toList t ]

eraseStrand :: Region -> Region
eraseStrand ( n, c, _, s, e ) = ( n, c, Both, s, e )

-- CPS ~~ StateT (Symtab, Annotab) IO
newtype CPS a = CPS { runCPS :: forall r . (a -> Symtab -> MAnnotab -> IO r) -> Symtab -> MAnnotab -> IO r }

instance Functor CPS where
    fmap f m = CPS $ \k -> runCPS m (k . f)

instance Applicative CPS where
    pure a = CPS $ \k -> k a
    a <*> b = CPS $ \k -> runCPS a (\f -> runCPS b (k . f))

instance Monad CPS where
    return a = CPS $ \k -> k a
    m >>= k = CPS $ \c -> runCPS m (\a -> runCPS (k a) c)

instance MonadIO CPS where
    liftIO m = CPS $ \k s t -> m >>= \a -> k a s t

get_syms :: CPS Symtab
get_syms = CPS $ \k s -> k s s

modify_syms :: (Symtab -> Symtab) -> CPS ()
modify_syms f = CPS $ \k s -> k () $! f s

get_anno :: CPS MAnnotab
get_anno = CPS $ \k s t -> k t s t

modify_anno :: (MAnnotab -> MAnnotab) -> CPS ()
modify_anno f = CPS $ \k s t -> k () s $! f t

execCPS :: CPS a -> IO a
execCPS m = runCPS m (\a _ _ -> return a) M.empty M.empty
