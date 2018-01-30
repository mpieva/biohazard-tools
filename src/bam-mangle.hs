{-# LANGUAGE RecordWildCards #-}
-- Generic bam mangler.  Mostly filters, some mutating operations, some
-- extraction of data.  Simple expression language, comparable to Unix
-- 'find', but not to 'awk'.

import Bio.Adna
import Bio.Bam                           hiding ( ParseError )
import Bio.Prelude                       hiding ( try, many )
import Control.Monad.Trans.RWS.Strict
import Data.HashMap.Strict                      ( lookupDefault, fromList )
import GHC.Float                                ( float2Double )
import Paths_biohazard_tools                    ( version )
import System.Console.GetOpt
import Text.Parsec                       hiding ( (<|>) )
import Text.Parsec.Language                     ( emptyDef )
import Text.Regex.Posix                         ( (=~) )

import qualified Data.ByteString                  as B
import qualified Data.Vector.Generic              as V
import qualified Text.Parsec.Token                as P

data Conf = Conf
    { conf_output :: BamMeta -> Iteratee [BamRaw] IO ()
    , conf_expr   :: Maybe String
    , conf_limit  :: Maybe Int }


defaultConf :: Conf
defaultConf = Conf (protectTerm . pipeBamOutput) Nothing Nothing

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o"  ["output"]         (ReqArg set_output "FILE") "Send output to FILE",
    Option "S"  ["sam-pipe"]       (NoArg     set_sam_output) "Send sam to stdout",
    Option "e"  ["expr","filter"]  (ReqArg set_expr   "EXPR") "Use EXPR as filter",
    Option "n"  ["numout","limit"] (ReqArg set_limit   "NUM") "Output at most NUM records",
    Option "p"  ["print"]          (ReqArg set_print  "EXPR") "Print result of EXPR for each record",
    Option "s"  ["sum"]            (ReqArg set_sum    "EXPR") "Sum results of EXPR",
    Option "h?" ["help","usage"]   (NoArg         disp_usage) "Display this message",
    Option "V"  ["version"]        (NoArg       disp_version) "Display version and exit" ]
  where
    set_output "-" c = return $ c { conf_output = pipeBamOutput }
    set_output  fp c = return $ c { conf_output = writeBamFile fp }
    set_sam_output c = return $ c { conf_output = joinI . mapStream unpackBam . pipeSamOutput }
    set_expr     s c = return $ c { conf_expr   = Just s }
    set_limit    a c =   (\l -> c { conf_limit  = Just l }) <$> readIO a

    disp_version _ = do pn <- getProgName
                        hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
                        exitSuccess

    disp_usage _ = do p <- getProgName
                      hPutStrLn stderr $ "Usage: " ++ usageInfo (p ++ info) options
                      exitSuccess

    info = " [option...] [bam-file...]\n" ++
           "Filters a set of BAM files by applying an expression to each record.  Options are:"

    set_sum e c = case parse_str p_num_atom e of
        Left   err  -> hPrint stderr err >> exitFailure
        Right f_num -> return $ c { conf_output = \m -> foldStream (p m) 0 >>= liftIO . print  }
            where p m acc = (+) acc . fst . evalExpr f_num m

    set_print e c = case parse_str (Left <$> p_string_atom <|> Right <$> p_num_atom) e of
        Left           err  -> hPrint stderr err >> exitFailure
        Right (Left  f_str) -> return $ c { conf_output = mapStreamM_ . p }
            where p m = B.hPut stdout . (`B.snoc` 10) . fst . evalExpr f_str m
        Right (Right f_num) -> return $ c { conf_output = mapStreamM_ . p }
            where p m = print . fst . evalExpr f_num m

type Parser a = Parsec String () a

-- The plan:  We take a filter expression and a number of input files.  Inputs
-- are merged, if that doesn't make sense, use only one.  The expression
-- evaluates to a 'Bool', we filter on the result.  We can take a second
-- expression to extract some data, which we may want to sum up.  Where
-- do we put the few modifying functions?

parse_str :: Parser a -> String -> Either ParseError a
parse_str p = parse (spaces *> p <* eof) ""

main :: IO ()
main = do
    (opts, files0, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    conf <- foldr (>=>) return opts defaultConf
    add_pg <- addPG $ Just version

    let (expr, files) = maybe (case files0 of e:fs -> (e,fs)
                                              [  ] -> error "Expression expected."
                              ) (\e -> (e,files0)) (conf_expr conf)

    case parse_str p_bool_expr expr of
        Left err -> hPrint stderr err >> exitFailure
        Right f  -> concatInputs files >=> run $ \meta ->
                    joinI $ mapMaybeStream (runExpr f meta) $
                    joinI $ maybe (mapChunks id) takeStream (conf_limit conf) $
                    conf_output conf (add_pg meta)


-- | We compile an expression of type \"a\" to a Haskell function that
-- takes a BamHeader, a BamRaw, and returns an \"a\".
type Expr = RWS (RGData,Refs) () BamRaw
type RGData = HashMap Bytes (Bytes,Bytes)

runExpr :: Expr Bool -> BamMeta -> BamRaw -> Maybe BamRaw
runExpr e m = unp . evalExpr e m
  where
    unp (True,  s) = Just s
    unp (False, _) = Nothing


evalExpr :: Expr a -> BamMeta -> BamRaw -> (a, BamRaw)
evalExpr e m = unp . runRWS e (rgs, meta_refs m)
  where
    unp (a, s, ()) = (a, s)

    !rgs = fromList [ (rg,(sm,lb))
                    | ("RG", stuff) <- meta_other_shit m
                    , ("ID", rg) <- stuff
                    , ("LB", lb) <- stuff
                    , ("SM", sm) <- stuff ]

p_bool_expr :: Parser (Expr Bool)
p_bool_expr = chainl1 p_bool_conj (orE <$  operator "||")
  where
    orE x y = x >>= \e -> if e then pure True else y

p_bool_conj :: Parser (Expr Bool)
p_bool_conj = chainl1 p_bool_atom (andE <$ operator "&&")
  where
    andE x y = x >>= \e -> if e then y else pure False

-- | Boolean atoms:
-- - string comparisons
-- - numeric comparisons
-- - definedness/isnum/isstring of a tag
-- - negation
-- - true and false
-- - the flags

p_bool_atom :: Parser (Expr Bool)
p_bool_atom = choice
    [ char '(' *> spaces *> p_bool_expr <* char ')' <* spaces
    , fmap not     <$> (operator "!"          *> p_bool_atom)
    , is_defined   <$> (reserved "defined"    *> p_field_name)
    , is_num       <$> (reserved "isnum"      *> p_field_name)
    , is_string    <$> (reserved "isstring"   *> p_field_name)
    , fmap not     <$> (reserved "not"        *> p_bool_atom)
    , isDeaminated <$> (reserved "deaminated" *> natural)

    , try $ (\x o y -> liftA2 o x y) <$> p_string_atom <*> str_op <*> p_string_atom
    , try $ (\x o y -> liftA2 o x y) <$> p_num_atom    <*> num_op <*> p_num_atom
    , named_predicate =<< identifier ]
  where
    getF f = return $ gets $ f . unpackBam

    named_predicate "true"  = return $ pure True
    named_predicate "false" = return $ pure False

    named_predicate "paired"       = getF isPaired
    named_predicate "properly"     = getF isProperlyPaired
    named_predicate "unmapped"     = getF isUnmapped
    named_predicate "mapped"       = getF $ not . isUnmapped
    named_predicate "mate-unmapped"= getF isMateUnmapped
    named_predicate "mate-mapped"  = getF $ not . isMateUnmapped
    named_predicate "reversed"     = getF isReversed
    named_predicate "mate-reversed"= getF isMateReversed
    named_predicate "first-mate"   = getF isFirstMate
    named_predicate "second-mate"  = getF isSecondMate
    named_predicate "auxillary"    = getF isAuxillary
    named_predicate "failed"       = getF isFailsQC
    named_predicate "good"         = getF $ not . isFailsQC
    named_predicate "duplicate"    = getF isDuplicate
    named_predicate "trimmed"      = getF isTrimmed
    named_predicate "merged"       = getF isMerged
    named_predicate "alternative"  = getF isAlternative
    named_predicate "exact-index"  = getF isExactIndex

    named_predicate "clear-failed" = return $ do_clearF
    named_predicate "set-failed"   = return $ do_setF
    named_predicate "set-trimmed"  = return $ setFF 1
    named_predicate "set-merged"   = return $ setFF 2

    named_predicate key = unexpected $ shows key ": not a known predicate"

    num_op :: Parser (Double -> Double -> Bool)
    num_op = choice [ (==) <$ operator "==", (<=) <$ operator "<=", (>=) <$ operator ">="
                    , (/=) <$ operator "!=", (<)  <$ operator  "<", (>)  <$ operator  ">" ]

    str_op :: Parser (Bytes -> Bytes -> Bool)
    str_op = choice [ (=~) <$ operator "~", (.) not . (=~) <$ operator "!~" ]

    is_num key = gets $ \br -> case lookup key (b_exts (unpackBam br)) of
            Just (Int   _) -> True
            Just (Float _) -> True
            _              -> False

    is_string key = gets $ \br -> case lookup key (b_exts (unpackBam br)) of
            Just (Text _) -> True
            Just (Bin  _) -> True
            Just (Char _) -> True
            _             -> False

    is_defined key = gets $ isJust . lookup key . b_exts . unpackBam

    do_clearF = True <$ modify (\br ->
                    let f_hi = B.index (raw_data br) 15
                        rd'  = B.concat [ B.take 15 (raw_data br)
                                        , B.singleton (clearBit f_hi 1)
                                        , B.drop 16 (raw_data br) ]
                    in if testBit f_hi 1 then br { raw_data = rd' } else br)

    do_setF   = True <$ modify (\br ->
                    let f_hi = B.index (raw_data br) 15
                        rd'  = B.concat [ B.take 15 (raw_data br)
                                        , B.singleton (setBit f_hi 1)
                                        , B.drop 16 (raw_data br) ]
                    in if testBit f_hi 1 then br else br { raw_data = rd' })

    setFF f   = True <$ modify (\br -> let b  = unpackBam br
                                           ff = extAsInt 0 "FF" b
                                           b' = b { b_exts = updateE "FF" (Int (ff .|. f)) $ b_exts b }
                                       in if ff .|. f == ff then br else unsafePerformIO (packBam b'))

isDeaminated :: Int -> Expr Bool
isDeaminated n = gets $ \br ->
    case unpackBam br of { b ->
    case alnFromMd (b_seq b) (b_cigar b) <$> getMd b of
        Nothing  -> False
        Just aln -> V.any (== cg) (V.take    n  aln) ||
                    V.any (== cg) (V.drop (l-n) aln) ||
                    V.any (== ga) (V.drop (l-n) aln)
          where
            l  = V.length aln
            cg = npair nucsC nucsT
            ga = npair nucsG nucsA  }


-- | Numeric-valued atoms (float or int, everything is converted to
-- double).  This includes:
-- - tagged fields (default to 0 if missing or wrong type)
-- - predefined fields: RNAME, POS, MRNM, MPOS, ISIZE, LENGTH, MAPQ
-- - mnemonic constants: wrongness, unknownness, rgquality
-- - literals

p_num_atom :: Parser (Expr Double)
p_num_atom = ( pure . either fromIntegral id <$> naturalOrFloat )
         <|> ( numeric_field                 =<< identifier )
  where
    getF f = return $ gets $ fromIntegral . f . unpackBam

    numeric_field "RNAME"       = getF $ unRefseq . b_rname
    numeric_field "POS"         = getF $              b_pos
    numeric_field "MRNM"        = getF $ unRefseq .  b_mrnm
    numeric_field "RNEXT"       = getF $ unRefseq .  b_mrnm
    numeric_field "MPOS"        = getF $             b_mpos
    numeric_field "PNEXT"       = getF $             b_mpos
    numeric_field "ISIZE"       = getF $            b_isize
    numeric_field "TLEN"        = getF $            b_isize
    numeric_field "MAPQ"        = getF $ unQ .       b_mapq
    numeric_field "LENGTH"      = getF $ V.length .   b_seq
    numeric_field "LEN"         = getF $ V.length .   b_seq

    numeric_field "unknownness" = getF $ extAsInt    0 "Z0"
    numeric_field "rgquality"   = getF $ extAsInt 9999 "Z1"
    numeric_field "wrongness"   = getF $ extAsInt    0 "Z2"

    numeric_field key
        | length key /= 2       = unexpected $ shows key ": not a known numeric field"
        | otherwise             = return $ gets $ \br ->
            case lookup (fromString key) (b_exts (unpackBam br)) of
                Just (Int   x) -> fromIntegral x
                Just (Float x) -> float2Double x
                _              -> 0



-- | Parses name of a tagged field.  This is an alphabetic character
-- followed by an alphanumeric character.
p_field_name :: Parser BamKey
p_field_name = do nm <- identifier
                  guard (length nm == 2)
                  return $ fromString nm
    <?> "2-letter tag"

-- | String-valued atoms.  This includes:
-- - tagged fields (default to "" if missing or wrong type)
-- - predefined fields: RNAME, MRNM
-- - convenience functions:  library, sample
-- - literals
p_string_atom :: Parser (Expr Bytes)
p_string_atom = choice
    [ lookupRef b_rname <$  reserved "RNAME"
    , lookupRef  b_mrnm <$  reserved "MRNM"
    , lookupRef  b_mrnm <$  reserved "RNEXT"
    , get_library       <$  reserved "library"
    , get_sample        <$  reserved "sample"
    , from_string_field <$> p_field_name
    , pure . fromString <$> stringLiteral ]
  where
    from_string_field key = gets $ \br -> case lookup key (b_exts (unpackBam br)) of
            Just (Text x) -> x
            Just (Bin  x) -> x
            Just (Char x) -> B.singleton x
            _             -> B.empty

    lookupRef f = asks snd >>= \m -> gets $ sq_name . getRef m . f . unpackBam

    lookup_rg b = asks $ lookupDefault (B.empty,B.empty) (extAsString "RG" b) . fst

    get_library :: Expr Bytes
    get_library = do b <- gets unpackBam
                     let lb = extAsString "LB" b
                     if B.null lb then snd <$> lookup_rg b else return lb

    get_sample :: Expr Bytes
    get_sample = gets unpackBam >>= fmap fst . lookup_rg


operator :: String -> Parser ()
operator nm = try (do cs <- many1 (oneOf ":!#$%&*+./<=>?@\\^|-~")
                      guard $ cs == nm
                      spaces)
              <?> show nm

reserved :: String -> Parser ()
reserved nm = try (do c <- letter
                      cs <- many (alphaNum <|> char '-')
                      guard $ c:cs == nm
                      spaces)
              <?> show nm

identifier :: Parser String
identifier = (:) <$> letter <*> many (alphaNum <|> char '-') <* spaces

natural :: Parser Int
natural = foldl' (\a c -> 10*a + digitToInt c) 0 <$> many1 digit <* spaces <?> "natural number"

naturalOrFloat :: Parser (Either Integer Double)
naturalOrFloat = P.naturalOrFloat (P.makeTokenParser emptyDef) <?> "number"

stringLiteral :: Parser String
stringLiteral  = P.stringLiteral  (P.makeTokenParser emptyDef) <?> "string"
