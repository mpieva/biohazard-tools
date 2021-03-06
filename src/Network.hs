module Network where

import Bio.Prelude                    hiding ( loop )
import Control.Monad.IO.Class
import Control.Monad.Trans.RWS.Strict
import Network.Socket

import qualified Data.Binary.Get                    as B
import qualified Data.Binary.Put                    as B
import qualified Data.ByteString                    as S
import qualified Network.Socket.ByteString          as N ( recv )
import qualified Network.Socket.ByteString.Lazy     as N ( sendAll )

import FormattedIO
import Symtab

data Which   = Covered | Fragment | Inclusive | Nearest deriving (Enum, Show)

--  The protocol, initially:
--  -> StartAnno name
--          Open the annotation set named "name".
--  <- KnownAnno
--          The server already has this set and brought it into scope or
--  <- UnknownAnno
--          The server doesn't know about this set and created it.
--
--  With an open annotation set:
--  -> AddAnno region
--          Add to the currently open annotation set or
--  -> EndAnno
--          The set is complete now and to be brought into scope.
--
--  At any time:
--  -> Anno which region
--          Request annotation of "region" using lookup strategy
--          "which" from all sets in scope.
--  <- Result name annos
--          Resulting annotation for the region named "name".

data Request = StartAnno Bytes
             | AddAnno Region
             | EndAnno
             | StartFile Bytes
             | Anno Which Region
             | EndStream
    deriving Show

data Response = UnknownAnno
              | KnownAnno
              | StartedFile Bytes
              | Result Bytes AnnoSet
              | EndedStream
    deriving Show

putWord :: ( Bits a, Integral a ) => a -> B.Put
putWord i | i < 0x80  =    B.putWord8 (fromIntegral i)
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
                    1 -> AddAnno <$> getRegion
                    2 -> return EndAnno
                    3 -> StartFile <$> getString
                    4 -> return EndStream
                    _ -> Anno (toEnum (key - 5)) <$> getRegion

getRegion :: B.Get Region
getRegion = (,,,,) <$> getString
                   <*> getString
                   <*> fmap toEnum getWord
                   <*> getWord
                   <*> getWord

putRegion :: Region -> B.Put
putRegion (n,m,s,u,v) = do putString n
                           putString m
                           putWord $ fromEnum s
                           putWord u
                           putWord v

putRequest :: Request -> B.Put
putRequest (StartAnno name) = do putWord (0::Int)
                                 putString name
putRequest (AddAnno region) = do putWord (1::Int)
                                 putRegion region
putRequest  EndAnno         =    putWord (2::Int)
putRequest (StartFile name) = do putWord (3::Int)
                                 putString name
putRequest  EndStream       =    putWord (4::Int)
putRequest (Anno w region)  = do putWord $ fromEnum w + 5
                                 putRegion region

getResponse :: B.Get Response
getResponse = do key <- getWord
                 case key :: Int of
                    0 -> return UnknownAnno
                    1 -> return KnownAnno
                    3 -> StartedFile <$> getString
                    4 -> return EndedStream
                    _ -> do n <- getString
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
putResponse  UnknownAnno  = putWord (0::Int)
putResponse  KnownAnno    = putWord (1::Int)
putResponse (Result n as) = do putWord (2::Int)
                               putString n
                               case as of NearMiss a d -> do putWord (if d >= 0 then 0 else 1 ::Int)
                                                             putString a
                                                             putWord (abs d)
                                          Hits hs      -> do putWord (length hs + 2)
                                                             mapM_ putString hs
putResponse (StartedFile nm) = putWord (3::Int) >> putString nm
putResponse  EndedStream     = putWord (4::Int)


getString :: B.Get Bytes
getString = getWord >>= B.getByteString

putString :: Bytes -> B.Put
putString s = do putWord $ S.length s ; B.putByteString s


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


type Client = RWST Socket () Bytes IO

receive :: B.Get a -> Client a
receive getter = do
    sock <- ask
    incoming <- get
    loop sock (B.runGetIncremental getter `B.pushChunk` incoming)
  where
    loop _sock (B.Fail    _ _ m) = fail m
    loop _sock (B.Done rest _ a) = put rest >> return a
    loop  sock (B.Partial     k) = do inp <- liftIO $ N.recv sock 4000
                                      if S.null inp then liftIO (close sock) >> fail "connection closed unexpectedly"
                                                    else loop sock $ k $ Just inp

send :: B.Put -> Client ()
send putter = do sock <- ask ; liftIO $ N.sendAll sock $ B.runPut putter
