{-# LANGUAGE RecordWildCards #-}
-- Generic bam mangler.  Mostly filters, some mutating operations, some
-- extraction of data.  Simple expression language, comparable to Unix
-- 'find', but not to 'awk'.

import Bio.Adna                                 ( alnFromMd )
import Bio.Bam                           hiding ( ParseError )
import Bio.Prelude                       hiding ( try )
import Control.Monad.Trans.RWS.Strict
import Data.Functor.Identity                    ( Identity )
import Data.HashMap.Strict                      ( lookupDefault, fromList )
import GHC.Float                                ( float2Double )
import Paths_biohazard_tools                    ( version )
import System.Console.GetOpt
import Text.Parsec                       hiding ( (<|>) )
import Text.Parsec.Token
import Text.Regex.Posix                         ( (=~) )

import qualified Data.ByteString                  as B
import qualified Data.Vector.Generic              as V

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
    set_limit    a c = readIO a >>= \l -> return $ c { conf_limit = Just l }

    disp_version _ = do pn <- getProgName
                        hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
                        exitSuccess

    disp_usage _ = do p <- getProgName
                      hPutStrLn stderr $ "Usage: " ++ usageInfo (p ++ info) options
                      exitSuccess

    info = " [option...] [bam-file...]\n\
           \Filters a set of BAM files by applying an expression to each record.  Options are:"

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
-- evaluates to a bool, we filter on the result.  We can take a second
-- expression to extract some data, which we may want to sum up.  Where
-- do we put the few modifying functions?

parse_str :: Parser a -> String -> Either ParseError a
parse_str p = parse (whiteSpace tp *> p <* eof) ""

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


-- | We compile an expression of type a to a Haskell function that takes
-- a BamHeader, a BamRaw, and returns an a.
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
p_bool_expr = (>>=)
    <$> p_bool_conj
    <*> option return
          ((\y x -> if x then pure x else y) <$> (reservedOp tp "||" *> p_bool_expr))

p_bool_conj :: Parser (Expr Bool)
p_bool_conj = (>>=)
    <$> p_bool_atom
    <*> option return
          ((\y x -> if x then y else pure x) <$> (reservedOp tp "&&" *> p_bool_conj))

-- | Boolean atoms:
-- - string comparisons
-- - numeric comparisons
-- - definedness/isnum/isstring of a tag
-- - negation
-- - true and false
-- - the flags

p_bool_atom :: Parser (Expr Bool)
p_bool_atom = choice
    [ parens tp p_bool_expr
    , fmap not     <$> (reservedOp tp "!"        *> p_bool_atom)
    , is_defined   <$> (reserved tp "defined"    *> p_field_name)
    , is_num       <$> (reserved tp "isnum"      *> p_field_name)
    , is_string    <$> (reserved tp "isstring"   *> p_field_name)
    , fmap not     <$> (reserved tp "not"        *> p_bool_atom)
    , isDeaminated <$> (reserved tp "deaminated" *> natural tp)

    , try $ (\x o y -> liftA2 o x y) <$> p_string_atom <*> str_op <*> p_string_atom
    , try $ (\x o y -> liftA2 o x y) <$> p_num_atom    <*> num_op <*> p_num_atom
    , named_predicate =<< identifier tp ]
  where
    named_predicate "true"  = return $ pure True
    named_predicate "false" = return $ pure False

    named_predicate "paired"       = return $ gets (isPaired         . unpackBam)
    named_predicate "properly"     = return $ gets (isProperlyPaired . unpackBam)
    named_predicate "unmapped"     = return $ gets (isUnmapped       . unpackBam)
    named_predicate "mate-unmapped"= return $ gets (isMateUnmapped   . unpackBam)
    named_predicate "reversed"     = return $ gets (isReversed       . unpackBam)
    named_predicate "mate-reversed"= return $ gets (isMateReversed   . unpackBam)
    named_predicate "first-mate"   = return $ gets (isFirstMate      . unpackBam)
    named_predicate "second-mate"  = return $ gets (isSecondMate     . unpackBam)
    named_predicate "auxillary"    = return $ gets (isAuxillary      . unpackBam)
    named_predicate "failed"       = return $ gets (isFailsQC        . unpackBam)
    named_predicate "duplicate"    = return $ gets (isDuplicate      . unpackBam)
    named_predicate "trimmed"      = return $ gets (isTrimmed        . unpackBam)
    named_predicate "merged"       = return $ gets (isMerged         . unpackBam)
    named_predicate "alternative"  = return $ gets (isAlternative    . unpackBam)
    named_predicate "exact-index"  = return $ gets (isExactIndex     . unpackBam)

    named_predicate "clear-failed" = return $ do_clearF
    named_predicate "set-failed"   = return $ do_setF
    named_predicate "set-trimmed"  = return $ setFF 1
    named_predicate "set-merged"   = return $ setFF 2

    named_predicate key = unexpected $ shows key ": not a known predicate"

    num_op :: Parser (Double -> Double -> Bool)
    num_op = choice [ (==) <$ reservedOp tp "==", (<=) <$ reservedOp tp "<=", (>=) <$ reservedOp tp ">="
                    , (/=) <$ reservedOp tp "!=", (<)  <$ reservedOp tp  "<", (>)  <$ reservedOp tp  ">" ]

    str_op :: Parser (Bytes -> Bytes -> Bool)
    str_op = choice [ (=~) <$ reservedOp tp "~", (.) not . (=~) <$ reservedOp tp "!~" ]

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

isDeaminated :: Integer -> Expr Bool
isDeaminated nn = gets $ \br ->
    case unpackBam br of { b ->
    case alnFromMd (b_seq b) (b_cigar b) <$> getMd b of
        Nothing  -> False
        Just aln -> V.any (== cg) (V.take n aln) ||
                    V.any (== cg) (V.drop (l-n) aln) ||
                    V.any (== ga) (V.drop (l-n) aln)
          where
            l  = V.length aln
            n  = fromIntegral nn
            cg = ( nucsC, nucsT )
            ga = ( nucsG, nucsA ) }


-- | Numeric-valued atoms (float or int, everything is converted to
-- double).  This includes:
-- - tagged fields (default to 0 if missing or wrong type)
-- - predefined fields: RNAME, POS, MRNM, MPOS, ISIZE, LENGTH, MAPQ
-- - mnemonic constants: wrongness, unknownness, rgquality
-- - literals

p_num_atom :: Parser (Expr Double)
p_num_atom = ( pure . either fromIntegral id <$> naturalOrFloat tp )
         <|> ( numeric_field                 =<< identifier tp )
  where
    numeric_field "RNAME"       = return $ gets $ fromIntegral . unRefseq . b_rname . unpackBam
    numeric_field "POS"         = return $ gets $ fromIntegral .              b_pos . unpackBam
    numeric_field "MRNM"        = return $ gets $ fromIntegral . unRefseq .  b_mrnm . unpackBam
    numeric_field "RNEXT"       = return $ gets $ fromIntegral . unRefseq .  b_mrnm . unpackBam
    numeric_field "MPOS"        = return $ gets $ fromIntegral .             b_mpos . unpackBam
    numeric_field "PNEXT"       = return $ gets $ fromIntegral .             b_mpos . unpackBam
    numeric_field "ISIZE"       = return $ gets $ fromIntegral .            b_isize . unpackBam
    numeric_field "TLEN"        = return $ gets $ fromIntegral .            b_isize . unpackBam
    numeric_field "MAPQ"        = return $ gets $ fromIntegral . unQ .       b_mapq . unpackBam
    numeric_field "LENGTH"      = return $ gets $ fromIntegral . V.length .   b_seq . unpackBam
    numeric_field "LEN"         = return $ gets $ fromIntegral . V.length .   b_seq . unpackBam

    numeric_field "unknownness" = return $ gets $ fromIntegral . extAsInt    0 "Z0" . unpackBam
    numeric_field "rgquality"   = return $ gets $ fromIntegral . extAsInt 9999 "Z1" . unpackBam
    numeric_field "wrongness"   = return $ gets $ fromIntegral . extAsInt    0 "Z2" . unpackBam

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
p_field_name = do nm <- identifier tp -- a <- letter
                  guard (length nm == 2)
                  return $ fromString nm
    <?> "2-letter tag"

-- | String-valued atoms.  This includes:
-- - tagged fields (default to "" if missing or wrong type)
-- - predefined fields: RNAME, MRNM
-- - convenience functions:  library?
-- - literals
p_string_atom :: Parser (Expr Bytes)
p_string_atom = choice
    [ from_string_field <$> p_field_name
    , lookupRef b_rname <$  reserved tp "RNAME"
    , lookupRef  b_mrnm <$  reserved tp "MRNM"
    , lookupRef  b_mrnm <$  reserved tp "RNEXT"
    , get_library       <$  reserved tp "library"
    , get_sample        <$  reserved tp "sample"
    , pure . fromString <$> stringLiteral tp ]
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


tp :: GenTokenParser String u Identity
tp = makeTokenParser $ LanguageDef
    { commentStart    = ""
    , commentEnd      = ""
    , commentLine     = ""
    , nestedComments  = False
    , identStart      = letter
    , identLetter     = alphaNum <|> char '-'
    , opStart         = oneOf ":!#$%&*+./<=>?@\\^|-~"
    , opLetter        = oneOf ":!#$%&*+./<=>?@\\^|-~"
    , reservedNames   = words "not defined isnum isstring deaminated"
    , reservedOpNames = words "|| && > < >= <= == != ! ~ !~"
    , caseSensitive   = True }

