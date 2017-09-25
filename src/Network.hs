module Network where

import Bio.Prelude hiding ( loop )
import Network.Socket

import qualified Data.Binary.Get                    as B
import qualified Data.Binary.Put                    as B
import qualified Data.ByteString.Char8              as S
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Network.Socket.ByteString          as N ( recv )
import qualified Network.Socket.ByteString.Lazy     as N ( sendAll )

import FormattedIO
import Symtab

data Which   = Covered | Fragment | Inclusive | Nearest deriving (Enum, Show)

data Request = StartAnno S.ByteString
             | AddAnno Region
             | EndAnno
             | Anno Which Region deriving Show

data Response = UnknownAnno | KnownAnno | Result L.ByteString AnnoSet deriving Show

putWord :: ( Bits a, Integral a ) => a -> B.Put
putWord i | i < 0x80  = do B.putWord8 (fromIntegral i)
          | otherwise = do B.putWord8 $ (fromIntegral i .&. 0x7f) .|. 0x80
                           putWord $ i `shiftR` 7

getWord :: ( Bits a, Integral a ) => B.Get a
getWord = do
    w <- B.getWord8
    if w < 0x80 then return (fromIntegral w)
                else do w' <- getWord
                        return $ (w' `shiftL` 7) .|. (fromIntegral w .&. 0x7f)

getRequest :: B.Get Request
getRequest = do key <- getWord
                case key of
                    0 -> StartAnno <$> getString
                    1 -> AddAnno <$> ( (,,,,) <$> getLazyString
                                              <*> getLazyString
                                              <*> (toSense <$> getWord)
                                              <*> getWord
                                              <*> getWord )
                    2 -> return EndAnno
                    _ -> pack (toEnum (key - 3)) <$> getLazyString
                                                 <*> getLazyString
                                                 <*> (toSense <$> getWord)
                                                 <*> getWord
                                                 <*> getWord
  where
    pack w n m s u v = Anno w (n,m,s,u,v)

    toSense :: Int -> [Sense]
    toSense 0 = []
    toSense 1 = [Forward]
    toSense 2 = [Reverse]
    toSense _ = [Forward, Reverse]


putRequest :: Request -> B.Put
putRequest (StartAnno name)      = do putWord (0::Int)
                                      putString name
putRequest (AddAnno (n,m,s,u,v)) = do putWord (1::Int)
                                      putLazyString n
                                      putLazyString m
                                      putWord $ fromSense s
                                      putWord $ u
                                      putWord $ v
putRequest (EndAnno)             = do putWord (2::Int)
putRequest (Anno w (n,m,s,u,v))  = do putWord $ fromEnum w + 3
                                      putLazyString n
                                      putLazyString m
                                      putWord $ fromSense s
                                      putWord $ u
                                      putWord $ v

fromSense :: [Sense] -> Int
fromSense []        = 0
fromSense [Forward] = 1
fromSense [Reverse] = 2
fromSense _         = 3


getResponse :: B.Get Response
getResponse = do key <- getWord
                 case key :: Int of
                    0 -> return UnknownAnno
                    1 -> return KnownAnno
                    _ -> do n <- getLazyString
                            l <- getWord
                            case l of 0 -> do a <- getString
                                              d <- getWord
                                              return (Result n (NearMiss a d))
                                      1 -> do a <- getString
                                              d <- getWord
                                              return (Result n (NearMiss a (-d)))
                                      _ -> do hs <- replicateM (l-2) getString
                                              return (Result n (Hits hs))

putResponse :: Response -> B.Put
putResponse (UnknownAnno) = putWord (0::Int)
putResponse (KnownAnno)   = putWord (1::Int)
putResponse (Result n as) = do putWord (2::Int)
                               putLazyString n
                               case as of NearMiss a d -> do putWord (if d >= 0 then 0 else 1 ::Int) ; putString a ; putWord (abs d)
                                          Hits hs      -> do putWord (length hs + 2) ; mapM_ putString hs


getString :: B.Get S.ByteString
getString = getWord >>= B.getByteString

getLazyString :: B.Get L.ByteString
getLazyString = L.fromChunks . (:[]) <$> getString

putString :: S.ByteString -> B.Put
putString s = do putWord $ S.length s ; B.putByteString s

putLazyString :: L.ByteString -> B.Put
putLazyString s = do putWord $ L.length s ; B.putLazyByteString s



service :: B.Get a -> (b -> B.Put) -> (s -> a -> IO (s, Maybe b)) -> Socket -> s -> IO s
service g p f sk = loop (B.runGetIncremental g)
  where
    loop (B.Done rest _ a) s = do (s',b) <- f s a
                                  maybe (return ()) (N.sendAll sk . B.runPut . p) b
                                  loop (B.runGetIncremental g `B.pushChunk` rest) s'
    loop (B.Fail    _ _ m) _ = do close sk
                                  fail m
    loop (B.Partial     k) s = do inp <- N.recv sk 4000
                                  if S.null inp
                                      then s <$ close sk
                                      else loop (k $ Just inp) s

