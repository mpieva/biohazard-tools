-- Generic bam mangler.  Mostly filters, some mutating operations, some
-- extraction of data.  Simple expression language, comparable to Unix
-- 'find', but not to 'awk'.

import Bio.Adna                                 ( alnFromMd, NPair )
import Bio.Bam
import Bio.Prelude
import Control.Monad.Trans.RWS.Strict
import Data.Attoparsec.ByteString.Char8         ( (<?>) )
import GHC.Float                                ( float2Double )
import Paths_biohazard_tools                    ( version )
import System.Console.GetOpt

import qualified Data.Attoparsec.ByteString.Char8 as A
import qualified Data.ByteString                  as B
import qualified Data.HashMap.Strict              as H
import qualified Data.Vector.Generic              as V

data Conf = Conf
    { conf_output :: BamMeta -> Iteratee [BamRaw] IO ()
    , conf_expr   :: Maybe String
    , conf_limit  :: Maybe Int }


defaultConf :: Conf
defaultConf = Conf (protectTerm . pipeBamOutput) Nothing Nothing

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "o" ["output"] (ReqArg set_output "FILE") "Send output to FILE",
    Option "S" ["sam-pipe"] (NoArg set_sam_output) "Send sam to stdout",
    Option "e" ["expr","filter"] (ReqArg set_expr "EXPR") "Use EXPR as filter",
    Option "n"  ["numout","limit"] (ReqArg set_limit "NUM") "Output at most NUM records",
    Option "h?" ["help","usage"] (NoArg disp_usage) "Display this message",
    Option "V" ["version"] (NoArg disp_version) "Display version and exit" ]
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



-- The plan:  We take a filter expression and a number of input files.  Inputs
-- are merged, if that doesn't make sense, use only one.  The expression
-- evaluates to a bool, we filter on the result.  We can take a second
-- expression to extract some data, which we may want to sum up.  Where
-- do we put the few modifying functions?


main :: IO ()
main = do
    (opts, files0, errors) <- getOpt Permute options <$> getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    conf <- foldr (>=>) return opts defaultConf
    add_pg <- addPG $ Just version

    let (expr, files) = maybe (case files0 of e:fs -> (e,fs)
                                              [  ] -> error "Expression expected."
                              ) (\e -> (e,files0)) (conf_expr conf)

    case A.parseOnly (A.skipSpace *> p_bool_expr <* A.endOfInput) (fromString expr) of
        Left err -> hPutStrLn stderr err >> exitFailure
        Right f  -> mergeInputs combineCoordinates files >=> run $ \meta ->
                    joinI $ mapMaybeStream (runExpr f meta) $
                    joinI $ maybe (mapChunks id) takeStream (conf_limit conf) $
                    conf_output conf (add_pg meta)


runExpr :: Expr Bool -> BamMeta -> BamRaw -> Maybe BamRaw
runExpr e m = unp . runRWS e (rgs, meta_refs m)
  where
    unp (True,  s, ()) = Just s
    unp (False, _, ()) = Nothing

    !rgs = H.fromList [ (rg,(sm,lb))
                      | ("RG", stuff) <- meta_other_shit m
                      , ("ID", rg) <- stuff
                      , ("LB", lb) <- stuff
                      , ("SM", sm) <- stuff ]


-- | We compile an expression of type a to a Haskell function that takes
-- a BamHeader, a BamRaw, and returns an a.
type Expr = RWS (RGData,Refs) () BamRaw
type RGData = H.HashMap Bytes (Bytes,Bytes)


p_bool_expr :: A.Parser (Expr Bool)
p_bool_expr = (>>=)
    <$> p_bool_conj
    <*> A.option return
          ((\y x -> if x then pure x else y) <$> (literal "||" *> p_bool_expr))


p_bool_conj :: A.Parser (Expr Bool)
p_bool_conj = (>>=)
    <$> p_bool_atom
    <*> A.option return
          ((\y x -> if x then y else pure x) <$> (literal "&&" *> p_bool_conj))


-- | Boolean atoms:
-- - string comparisons
-- - numeric comparisons
-- - definedness/isnum/isstring of a tag
-- - negation
-- - true and false
-- Missing:  tah flags

p_bool_atom :: A.Parser (Expr Bool)
p_bool_atom = A.choice
    [ A.char '(' *> A.skipSpace *> p_bool_expr <* A.char ')' <* A.skipSpace
    , pure True <$ literal "true"
    , pure False <$ literal "false"
    , is_defined <$> (literal "defined"  *> p_field_name)
    , is_num     <$> (literal "isnum"    *> p_field_name)
    , is_string  <$> (literal "isstring" *> p_field_name)
    , fmap not <$> (literal "not" *> p_bool_atom)
    , fmap not <$> (literal "!" *> p_bool_atom)

    , gets (isPaired         . unpackBam) <$ literal "paired"
    , gets (isProperlyPaired . unpackBam) <$ literal "properly"
    , gets (isUnmapped       . unpackBam) <$ literal "unmapped"
    , gets (isMateUnmapped   . unpackBam) <$ literal "mate-unmapped"
    , gets (isReversed       . unpackBam) <$ literal "reversed"
    , gets (isMateReversed   . unpackBam) <$ literal "mate-reversed"
    , gets (isFirstMate      . unpackBam) <$ literal "first-mate"
    , gets (isSecondMate     . unpackBam) <$ literal "second-mate"
    , gets (isAuxillary      . unpackBam) <$ literal "auxillary"
    , gets (isFailsQC        . unpackBam) <$ literal "failed"
    , gets (isDuplicate      . unpackBam) <$ literal "duplicate"
    , gets (isTrimmed        . unpackBam) <$ literal "trimmed"
    , gets (isMerged         . unpackBam) <$ literal "merged"
    , gets (isVestigial      . unpackBam) <$ literal "vestigial"
    , gets (isDeaminated     . unpackBam) <$ literal "deaminated"

    , do_clearF <$ literal "clear-failed"
    , do_setF   <$ literal "set-failed"
    , setFF 1   <$ literal "set-trimmed"
    , setFF 2   <$ literal "set-merged"

    , A.try $ (\x o y -> liftA2 o x y) <$> p_string_atom <*> str_op <*> p_string_atom
    , A.try $ (\x o y -> liftA2 o x y) <$> p_num_atom    <*> num_op <*> p_num_atom ]
  where
    num_op :: A.Parser (Double -> Double -> Bool)
    num_op = A.choice [ (==) <$ literal "==", (<=) <$ literal "<=", (>=) <$ literal ">="
                      , (/=) <$ literal "!=", (<)  <$ literal  "<", (>)  <$ literal  ">" ]

    str_op :: A.Parser (Bytes -> Bytes -> Bool)
    str_op = A.choice [ (==) <$ literal "~", (/=) <$ literal "!~" ]

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

isDeaminated :: BamRec -> Bool
isDeaminated b | V.length (b_seq b) < 4 = False      -- stupid corner case, don't want to crash
isDeaminated b =
    case alnFromMd (b_seq b) (b_cigar b) <$> getMd b of
        Nothing  -> False
        Just aln -> V.head aln == cg || V.head (V.tail aln) == cg ||
                    V.last aln == cg || V.last (V.init aln) == cg ||
                    V.last aln == ga || V.last (V.init aln) == ga
  where
    cg, ga :: NPair
    cg = ( nucsC, nucsT )
    ga = ( nucsG, nucsA )


-- | Numeric-valued atoms (float or int, everything is converted to
-- double).  This includes:
-- - tagged fields (default to 0 if missing or wrong type)
-- - predefined fields: RNAME, POS, MRNM, MPOS, ISIZE, LENGTH, MAPQ
-- - mnemonic constants: wrongness, unexpectedness, rgquality
-- - literals

p_num_atom :: A.Parser (Expr Double)
p_num_atom = A.choice
    [ from_num_field                                       <$> p_field_name
    , gets (fromIntegral . unRefseq . b_rname . unpackBam) <$  literal "RNAME"
    , gets (fromIntegral .              b_pos . unpackBam) <$  literal "POS"
    , gets (fromIntegral . unRefseq .  b_mrnm . unpackBam) <$  literal "MRNM"
    , gets (fromIntegral . unRefseq .  b_mrnm . unpackBam) <$  literal "RNEXT"
    , gets (fromIntegral .             b_mpos . unpackBam) <$  literal "MPOS"
    , gets (fromIntegral .             b_mpos . unpackBam) <$  literal "PNEXT"
    , gets (fromIntegral .            b_isize . unpackBam) <$  literal "ISIZE"
    , gets (fromIntegral .            b_isize . unpackBam) <$  literal "TLEN"
    , gets (fromIntegral . unQ .       b_mapq . unpackBam) <$  literal "MAPQ"
    , gets (fromIntegral . V.length .   b_seq . unpackBam) <$  literal "LENGTH"

    , gets (fromIntegral . extAsInt    0 "Z0" . unpackBam) <$  literal "unknownness"
    , gets (fromIntegral . extAsInt 9999 "Z1" . unpackBam) <$  literal "rgquality"
    , gets (fromIntegral . extAsInt    0 "Z2" . unpackBam) <$  literal "wrongness"
    , pure                                                 <$> A.double <* A.skipSpace ]
  where
    from_num_field key = gets $ \br -> case lookup key (b_exts (unpackBam br)) of
            Just (Int   x) -> fromIntegral x
            Just (Float x) -> float2Double x
            _              -> 0

literal :: Bytes -> A.Parser ()
literal key = A.string key *> A.skipSpace <?> unpack key

-- | Parses name of a tagged field.  This is an alphabetic character
-- followed by an alphanumeric character.
p_field_name :: A.Parser BamKey
p_field_name = A.try $ do a <- A.letter_ascii
                          b <- A.letter_ascii <|> A.digit
                          c <- A.peekChar
                          guard $ maybe True (\z -> not (A.isDigit z || A.isAlpha_ascii z)) c
                          A.skipSpace
                          return $ fromString [a,b]

-- | String-valued atoms.  This includes:
-- - tagged fields (default to "" if missing or wrong type)
-- - predefined fields: RNAME, MRNM
-- - convenience functions:  library?
-- - literals
p_string_atom :: A.Parser (Expr Bytes)
p_string_atom = A.choice
    [ from_string_field <$> p_field_name
    , lookupRef b_rname <$  literal "RNAME"
    , lookupRef  b_mrnm <$  literal "MRNM"
    , lookupRef  b_mrnm <$  literal "RNEXT"
    , get_library       <$  literal "library"
    , get_sample        <$  literal "sample"
    , pure              <$> p_string_lit ]
  where
    from_string_field key = gets $ \br -> case lookup key (b_exts (unpackBam br)) of
            Just (Text x) -> x
            Just (Bin  x) -> x
            Just (Char x) -> B.singleton x
            _             -> B.empty

    lookupRef f = asks snd >>= \m -> gets $ sq_name . getRef m . f . unpackBam

    lookup_rg b = asks $ H.lookupDefault (B.empty,B.empty) (extAsString "RG" b) . fst

    get_library :: Expr Bytes
    get_library = do b <- gets unpackBam
                     let lb = extAsString "LB" b
                     if B.null lb then snd <$> lookup_rg b else return lb

    get_sample :: Expr Bytes
    get_sample = gets unpackBam >>= fmap fst . lookup_rg


-- | String literal.  We'll need escapes, but for now this should do.
p_string_lit :: A.Parser Bytes
p_string_lit =     A.char '"' *> A.takeWhile (/='"') <* A.char '"' <* A.skipSpace
               <|> A.char '\'' *> A.takeWhile (/='\'') <* A.char '\'' <* A.skipSpace

