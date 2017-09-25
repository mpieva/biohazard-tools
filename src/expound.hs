{- This is now supposed to read both Ensembl- and UCSC-style input, and
 - it had a very weird heuristic to convert them into each other.  We
 - drop this shit now, the rules are:  following the specification, all
 - coordinates are 0-based half open intervals for BED, BAM,
 - weird-FastA, 5col; 1-based closed intervals for SAM.  Annotations in
 - 5col format (and only in 5col format!) can be rewritten using a
 - translation table.  If the coordinate conventions change, the offsets
 - in the translation table need to reflect that.
 -
 - Strand information is optional now, too.  We handle it this way:
 -
 - * If the strand field starts with '-', only the reverse strand is meant.
 - * If the strand field is "0" or absent, both strands are meant.
 - * Else only the forward strand is meant.
 - * Every annotation is reported at most once.
 -
 - Internally, we *always* use zero-based, half-open intervals.  We also
 - use UCSC-names whereever possible.  The only exception is if no
 - translation table is given, then we operate with Ensembl-names and
 - the program blows up if it hits a UCSC name.  In that case, the
 - occasional off-by-one error is easily possible.
 -}

import Bio.Prelude           hiding ( Word, loop )
import Data.Vector                  ( (!) )
import Network.Socket        hiding ( send )
import Paths_biohazard_tools        ( getDataFileName )
import System.Console.GetOpt
import System.Directory             ( doesFileExist )
import System.IO                    ( withFile, IOMode(..) )
import System.IO.Unsafe             ( unsafeInterleaveIO )

import qualified Data.Binary.Get                    as B
import qualified Data.Binary.Put                    as B
import qualified Data.ByteString.Lazy.Char8         as L
import qualified Data.ByteString.Char8              as S
import qualified Data.Map                           as M
import qualified Network.Socket.ByteString          as N ( recv )
import qualified Network.Socket.ByteString.Lazy     as N ( sendAll )

import Diet
import FormattedIO
import Network
import Symtab

data Options = Options
     { optOutput      :: [( FilePath, [Annotation] )] -> IO ()
     , optWhich       :: Which
     , optAnnotations :: [ (Maybe ChromTable -> [L.ByteString] -> [Region], FilePath) ]
     , optChromTable  :: Maybe ChromTable
     , optXForm       :: Region -> Region
     , optProgress    :: String -> IO ()
     , optPort        :: Maybe String
     , optHost        :: Maybe String
     }

defaultOptions :: Options
defaultOptions = Options
     { optOutput      = make_one_table stdout
     , optWhich       = Inclusive
     , optAnnotations = []
     , optChromTable  = Nothing
     , optXForm       = id
     , optProgress    = hPutStr stderr
     , optPort        = Nothing
     , optHost        = Nothing
     }

options, common_options :: [OptDescr (Options -> IO Options)]
options =
     [ Option ['o']     ["output"]
         (ReqArg (\f opts -> return $ opts { optOutput = if f == "-" then make_one_table stdout else \anns -> withFile f WriteMode (\h -> make_one_table h anns)})
                 "FILE")
         "write output to FILE"
     , Option ['O']     ["output-pattern"]
         (ReqArg (\f opts -> return $ opts { optOutput = write_file_pattern f })
                 "PAT")
         "write output to files PAT, replace %s with input file basename"
     , Option []        ["summarize"]
         (NoArg (\   opts -> return $ opts { optOutput = make_summary }))
         "summarize annotations to stdout, don't generate a table"
     , Option ['a']     ["annotation"]
         (ReqArg (\f opts -> return $ opts { optAnnotations = (const readMyBed,f) : optAnnotations opts })
                 "FILE")
         "read annotations from FILE in BED format"
     , Option ['A']     ["legacy-annotation"]
         (ReqArg (\f opts -> do hPutStr stderr legacy_warning
                                return $ opts { optAnnotations = (\ct -> map (xlate_chrom ct) . read5cFile,f) : optAnnotations opts })
                 "FILE")
         "read annotations from FILE in legacy five column format"
     , Option ['c']     ["chroms"]
         (ReqArg ((\ f opts -> readDataFile f >>= \s -> return $
                               opts { optChromTable = Just (readChromTable s) }) )
                 "FILE")
         "read chr translation table from FILE"
     , Option ['s']     ["nostrand"]
         (NoArg (\ opts -> return $ opts { optXForm = eraseStrand }))
         "ignore strand information in region files"
     , Option []        ["covered"]
         (NoArg (\ opts -> return $ opts { optWhich = Covered }))
         "output only annotations that cover the query interval"
     , Option []        ["partial"]
         (NoArg (\ opts -> return $ opts { optWhich = Fragment }))
         "output only annotations that cover only part of the query"
     , Option []        ["inclusive"]
         (NoArg (\ opts -> return $ opts { optWhich = Inclusive }))
         "output annotations that overlap at least part of the query"
     , Option []        ["nearest"]
         (NoArg (\ opts -> return $ opts { optWhich = Nearest }))
         "output the closest annotation if none overlaps"
     ]

common_options =
     [ Option ['p']     ["port"]
         (ReqArg (\ p opts -> return $ opts { optPort = Just p } ) "PORT")
         "listen on PORT or connect to PORT"
     , Option ['H']     ["host"]
         (ReqArg (\ h opts -> return $ opts { optHost = Just h }) "HOST")
         "connect to HOST"
     , Option ['q']     ["quiet"]
         (NoArg (\ opts -> return $ opts { optProgress = \_ -> return () }))
         "do not output progress reports"
     , Option ['h','?'] ["help","usage"]
         (NoArg (\ _ -> do putStrLn (usageInfo header $ options ++ common_options) ; exitSuccess ))
         "display this information"
     ]


legacy_warning :: String
legacy_warning = unlines [ "Warning: The five column format is not a standard format.  Consider"
                         , "         converting your annotations to the standard BED format." ]
header :: String
header = unlines [ "Usage: coord2anno [OPTION...] files...", ""
                 , "Annotates genomic regions given by coordinates.  Input files"
                 , "and annotation files are BED or BAM files.  Can operate in"
                 , "client-server-mode.  Options are:" ]

zipWithLM :: (a -> b -> IO c) -> [a] -> [b] -> IO [c]
zipWithLM f (x:xs) (y:ys) = do z <- f x y
                               zs <- unsafeInterleaveIO $ zipWithLM f xs ys
                               return (z:zs)
zipWithLM _ _ _ = return []


type LookupFn = RevSymtab -> Int -> Int -> [IDiet] -> AnnoSet

standardLookup :: (Int -> Int -> IDiet -> [Word]) -> LookupFn
standardLookup which syms s e diets =
    Hits . map (syms!) . map fromIntegral . unions $ map (which s e) diets

lookupNearest :: LookupFn
lookupNearest syms s e diets =
    let overlaps = unions $ map (lookupI s e) diets
        closeby  = sort $ filter (not . null . snd) $ map (lookupLeft s) diets ++ map (lookupRight e) diets
    in case (overlaps, closeby) of
        ([], [    ]) -> Hits []
        ([], (as:_)) -> NearMiss (syms ! fromIntegral a) d where (d,(a:_)) = as
        (is,      _) -> Hits . map (syms!) $ map fromIntegral is


interpretWhich :: Which -> RevSymtab -> Int -> Int -> [IDiet] -> AnnoSet
interpretWhich Covered   = standardLookup lookupICovered
interpretWhich Fragment  = standardLookup lookupIPartial
interpretWhich Inclusive = standardLookup lookupI
interpretWhich Nearest   = lookupNearest

readDataFile :: FilePath -> IO L.ByteString
readDataFile fp | any ('/' ==) fp = L.readFile fp
                | otherwise       = do e <- doesFileExist fp
                                       if e then L.readFile fp
                                            else getDataFileName fp >>= L.readFile

main :: IO ()
main = do
    args <- getArgs
    when (null args) $ putStrLn (usageInfo header $ options ++ common_options) >> exitSuccess

    let (os, files, errors) = getOpt Permute (options ++ common_options) args
    opts <- foldr (>=>) return os defaultOptions
    unless (null errors) $ mapM_ (hPutStr stderr) errors >> exitFailure

    case (optHost opts, optPort opts) of
        (Nothing, Just  p) -> do let (os', files', errors') = getOpt Permute common_options args
                                 opts' <- foldr (>=>) return os' defaultOptions
                                 unless (null errors') $ mapM_ (hPutStr stderr . (++) "in server mode: ") errors' >> exitFailure
                                 unless (null files') $ hPutStrLn stderr "no file arguments allowed in server mode" >> exitFailure
                                 run_server p opts'

        (Nothing, Nothing) -> run_standalone opts files
        (Just  h, Just  p) -> run_client h p opts files
        (Just  _, Nothing) -> hPutStrLn stderr "--host given, but no --port" >> exitFailure


run_standalone :: Options -> [FilePath] -> IO ()
run_standalone opts files = execCPS $ do
    emptyD <- io $ unsafeFreezeDiet =<< newDiet
    syms <- readGeneTable (optProgress opts) . concat =<< forM (optAnnotations opts)
            (\( reader, filepath ) ->
                reader (optChromTable opts). L.lines
                <$> io (if filepath == "-" then L.getContents else L.readFile filepath))
    anno <- io . sequenceA . M.map unsafeFreezeDiet =<< get_anno

    let annotate1 filename which num (name,c,strs,s,e) = do
                when ((1+num) `mod` 16384 == 0) . optProgress opts .
                    shows filename . (++) " (" . shows num . (++) "): " .
                    (++) (L.unpack name) $ "\27[K\r"
                return (name, which syms s e [ M.findWithDefault emptyD (c,str) anno | str <- strs ])

        annotate which f s = (,) f <$> zipWithLM (annotate1 f which) [0::Int ..]
                                                 (map (optXForm opts) $ readInput s)

    io $ optOutput opts =<< myReadFilesWith (annotate $ interpretWhich $ optWhich opts) files
    io $ optProgress opts "\n"


run_client :: HostName -> ServiceName -> Options -> [FilePath] -> IO ()
run_client h p opts files = do
    let hints = defaultHints { addrFlags = [AI_ADDRCONFIG, AI_CANONNAME] }
    addrs <- getAddrInfo (Just hints) (Just h) (Just p)
    let addr = head addrs
    sock <- socket (addrFamily addr) (addrSocketType addr) (addrProtocol addr)
    connect sock (addrAddress addr)

    (_, incoming) <- runClient (forM (optAnnotations opts) $ \( reader, filepath ) ->
                                    sendGeneTable filepath . reader (optChromTable opts) . L.lines
                                    =<< lift (if filepath == "-" then L.getContents else L.readFile filepath))
                               sock S.empty

    optOutput opts =<< myReadFilesWith (annotate sock incoming) files
    close sock

  where
    sendGeneTable fp regions = do
        send (putRequest (StartAnno (S.pack fp)))
        resp <- receive getResponse
        case resp of KnownAnno -> return ()
                     UnknownAnno -> do let loop _ [] = return ()
                                           loop i (r:rs) = do
                                                when ((i+1) `mod` 65556 == 0) $ lift $ optProgress opts $
                                                    "uploading " ++ shows fp ": " ++ shows i "\27[K\r"
                                                send (putRequest (AddAnno r))
                                                (loop $! i+1) rs
                                       loop (0::Int) regions
                                       send (putRequest EndAnno)
                                       return ()
                     _ -> fail "WTF?!"

    annotate sock incoming fp s = do
        done <- newMVar False
        in_flight <- newMVar (0::Int)

        _ <- forkIO $ do let inp = map (optXForm opts) $ readInput s
                         forM_ inp $ \rgn -> do takeMVar in_flight >>= \f -> putMVar in_flight $! succ f
                                                N.sendAll sock $ B.runPut $ putRequest (Anno (optWhich opts) rgn)
                         _ <- takeMVar done
                         putMVar done True

        let loop i (B.Done   rest _ (Result name annoset)) = do
                    ff <- takeMVar in_flight
                    putMVar in_flight $! ff-1
                    f <- if ff == 1 && S.null rest
                           then readMVar done
                           else return False
                    rs <- if f then return []
                               else unsafeInterleaveIO $ loop (1+i)
                                        (B.runGetIncremental getResponse `B.pushChunk` rest)
                    when ((1+i) `mod` 16384 == 0) . optProgress opts .
                        shows fp . (++) " (" . shows i . (++) "): " .
                        (++) (L.unpack name) $ "\27[K\r"
                    return ((name,annoset):rs)

            loop _ (B.Done _ _ _) = do close sock
                                       fail "unexpected answer from server"
            loop _ (B.Fail _ _ m) = do close sock
                                       fail m
            loop i (B.Partial  k) = do inp <- N.recv sock 4000
                                       if S.null inp then close sock >>
                                                          fail "connection closed unexpectedly"
                                                     else loop i (k $ Just inp)

        (,) fp <$> loop (0::Int) (B.runGetIncremental getResponse `B.pushChunk` incoming)


type HalfTable = ( S.ByteString, MAnnotab, Symtab )
type FullTable = ( IAnnotab, Symtab, RevSymtab )
type FullTables = M.Map S.ByteString FullTable

run_server :: ServiceName -> Options -> IO ()
run_server svname opts = do
    full_tables <- newMVar M.empty
    listener <- socket AF_INET Stream defaultProtocol
    port <- fromIntegral <$> (readIO svname :: IO Int)
    bind listener $ SockAddrInet port iNADDR_ANY
    listen listener 1
    optProgress opts $ "waiting for connections on port " ++ shows port "\27[K\n"
    let loop = do (sk, addr) <- accept listener
                  optProgress opts $ "connection from " ++ shows addr "\27[K\n"
                  let pr m = optProgress opts $ shows addr ": " ++ m
                  _ <- forkIO $ service getRequest putResponse (serverfn pr full_tables) sk (Nothing,[]) >> return ()
                  loop
    loop

serverfn :: (String -> IO ())
         -> MVar FullTables
         -> ( Maybe HalfTable, [FullTable] )
         -> Request
         -> IO ( ( Maybe HalfTable, [FullTable] ), Maybe Response)
serverfn _ _fulltables (Just _, _) (StartAnno _) = fail "client tried to open a second annoset"
serverfn pr fulltables (Nothing, fts) (StartAnno name) = do
    tablemap <- readMVar fulltables
    case M.lookup name tablemap of
        Just m  -> do pr $ "brought known annotation " ++ shows name " in scope\n"
                      return ( (Nothing, m:fts), Just KnownAnno)
        Nothing -> do pr $ "starting new annotation set " ++ shows name "\n"
                      return ( (Just (name, M.empty, M.empty), fts), Just UnknownAnno )

serverfn _ _ (Nothing,_) (AddAnno _) = fail "client tried growing an annotation set without opening it"
serverfn _ _ (Just (name, annotab, symtab), st) (AddAnno (gi, chrom, strs, s, e)) =
    runCPS (findSymbol gi >>= \val ->
            forM_ strs $ \str -> findDiet (chrom,str) >>=
                                 io . addI (min s e) (max s e) (fromIntegral val))
           (\_ symtab' annotab' -> return ((Just (name, annotab', symtab'), st), Nothing))
           symtab annotab

serverfn _ _ (Nothing,_) EndAnno = fail "client tried ending an annotation set without opening it"
serverfn pr fulltables (Just (name, annotab, symtab), s) EndAnno = do
    annotab' <- sequenceA $ M.map unsafeFreezeDiet annotab
    let anno' = (annotab', symtab, invertTab symtab)
    tablemap <- takeMVar fulltables
    putMVar fulltables $! M.insert name anno' tablemap
    pr $ "memoized new annotation set " ++ shows name "\n"
    return ((Nothing, anno':s), Nothing)

serverfn _ _ st@(_,annosets) (Anno which (name,c,strs,s,e)) = do
    emptyD <- unsafeFreezeDiet =<< newDiet
    let result = Result name . foldr combine_anno_sets (Hits []) $
           [ interpretWhich which syms s e [ M.findWithDefault emptyD (c,str) anno | str <- strs ]
           | (anno, _, syms) <- annosets ]
    return (st, Just result)

combine_anno_sets :: AnnoSet -> AnnoSet -> AnnoSet
combine_anno_sets (Hits []) x = x
combine_anno_sets x (Hits []) = x
combine_anno_sets (Hits xs) (NearMiss _ _) = Hits xs
combine_anno_sets (NearMiss _ _) (Hits ys) = Hits ys
combine_anno_sets (Hits xs) (Hits ys)      = Hits (xs++ys)
combine_anno_sets (NearMiss a d) (NearMiss b d')
    | abs d <= abs d' = NearMiss a d
    | otherwise       = NearMiss b d'


receive :: B.Get a -> Client a
receive getter = Client $ \sock incoming -> loop sock (B.runGetIncremental getter `B.pushChunk` incoming)
  where
    loop _sock (B.Fail    _ _ m) = fail m
    loop _sock (B.Done rest _ a) = return (a,rest)
    loop  sock (B.Partial     k) = do inp <- N.recv sock 4000
                                      if S.null inp then close sock >> fail "connection closed unexpectedly"
                                                    else loop sock $ k $ Just inp

send :: B.Put -> Client ()
send putter = Client $ \sock incoming -> do N.sendAll sock $ B.runPut putter ; return ((), incoming)


newtype Client a = Client { runClient :: Socket -> S.ByteString -> IO (a, S.ByteString) }

instance Functor Client where
    fmap f c = Client $ \s t -> runClient c s t >>= \(a, t') -> return (f a, t')

instance Applicative Client where
    pure a = Client $ \_ incoming -> return (a, incoming)
    a <*> b = Client $ \s t -> runClient a s t >>= \(f, t') ->
                               runClient b s t' >>= \(x, t'') ->
                               return (f x, t'')
instance Monad Client where
    return a = Client $ \_ incoming -> return (a, incoming)
    m >>= k = Client $ \sock incoming -> runClient m sock incoming >>= \(a, incoming') ->
                                         runClient (k a) sock incoming'

lift :: IO a -> Client a
lift k = Client $ \_ incoming -> k >>= \a -> return (a,incoming)


