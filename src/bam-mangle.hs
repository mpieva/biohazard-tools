-- Generic bam mangler.  Mostly filters, some mutating operations, some
-- extraction of data.  Simple expression language, comparable to Unix
-- 'find', but not to 'awk'.

import Bio.Bam
import Bio.Prelude
import GHC.Float ( float2Double )
import qualified Data.Attoparsec.ByteString.Char8 as A
import Paths_biohazard_tools ( version )
import System.Console.GetOpt
import qualified Data.ByteString as B
import qualified Data.Vector.Generic as V
import qualified Data.HashMap.Strict as H

-- data Conf = Conf
    -- { conf_output :: Iteratee Bytes IO ()
            -- ^ Output.  Depending on the expression supplied, we might
            -- send either BAM or some sort of text.

-- The plan:  We take a filter expression and a number of input files.  Inputs
-- are merged, if that doesn't make sense, use only one.  The expression
-- evaluates to a bool, we filter on the result.  We can take a second
-- expression to extract some data, which we may want to sum up.  Where
-- do we put the few modifying functions?


main :: IO ()
main = do
    expr : files <- getArgs
    case A.parseOnly (A.skipSpace *> p_bool_expr <* A.endOfInput) (fromString expr) of
        Left err -> hPutStrLn stderr err >> exitFailure
        Right f  -> go f files
  where
    go f files = mergeInputs combineCoordinates files >=> run $ \meta ->
                 joinI $ filterStream (f meta) $
                 joinI $ mapStream unpackBam $
                 pipeBamOutput meta


-- | We compile an expression of type a to a Haskell function that takes
-- a BamHeader, a BamRaw, and returns an a.
type Expr a = BamMeta -> BamRaw -> a

p_bool_expr :: A.Parser (Expr Bool)
p_bool_expr = (\a f h b -> f h b $ a h b)
    <$> p_bool_conj
    <*> A.option (\_ _ a -> a)
          ((\x h b a -> a || x h b) <$> (literal "||" *> p_bool_expr))


p_bool_conj :: A.Parser (Expr Bool)
p_bool_conj = (\a f h b -> f h b $ a h b)
    <$> p_bool_atom
    <*> A.option (\_ _ a -> a)
          ((\x h b a -> a && x h b) <$> (literal "&&" *> p_bool_conj))


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
    , const . const True <$ literal "true"
    , const . const False <$ literal "false"
    , is_defined <$> (literal "defined"  *> p_field_name)
    , is_num     <$> (literal "isnum"    *> p_field_name)
    , is_string  <$> (literal "isstring" *> p_field_name)
    , (\e h b -> not (e h b)) <$> (literal "not" *> p_bool_atom)
    , (\e h b -> not (e h b)) <$> (literal "!" *> p_bool_atom)

    , const (isPaired         . unpackBam) <$ literal "paired"
    , const (isProperlyPaired . unpackBam) <$ literal "properly"
    , const (isUnmapped       . unpackBam) <$ literal "unmapped"
    , const (isMateUnmapped   . unpackBam) <$ literal "mate-unmapped"
    , const (isReversed       . unpackBam) <$ literal "reversed"
    , const (isMateReversed   . unpackBam) <$ literal "mate-reversed"
    , const (isFirstMate      . unpackBam) <$ literal "first-mate"
    , const (isSecondMate     . unpackBam) <$ literal "second-mate"
    , const (isAuxillary      . unpackBam) <$ literal "auxillary"
    , const (isFailsQC        . unpackBam) <$ literal "failed"
    , const (isDuplicate      . unpackBam) <$ literal "duplicate"
    , const (isTrimmed        . unpackBam) <$ literal "trimmed"
    , const (isMerged         . unpackBam) <$ literal "merged"
    , const (isVestigial      . unpackBam) <$ literal "vestigial"

    , A.try $ (\x o y h b -> x h b `o` y h b) <$>
        p_string_atom <*> ((==) <$ literal "eq" <|> (/=) <$ literal "neq") <*> p_string_atom
    , A.try $ (\x o y h b -> x h b `o` y h b) <$> p_num_atom <*> num_op <*> p_num_atom ]
  where
    num_op :: A.Parser (Double -> Double -> Bool)
    num_op = A.choice
        [ (==) <$ literal "=="
        , (/=) <$ literal "/="
        , (<=) <$ literal "<="
        , (>=) <$ literal ">="
        , (<) <$ literal "<"
        , (>) <$ literal ">" ]

    is_num key _ br = case lookup key (b_exts (unpackBam br)) of
            Just (Int   _) -> True
            Just (Float _) -> True
            _              -> False

    is_string key _ br = case lookup key (b_exts (unpackBam br)) of
            Just (Text _) -> True
            Just (Bin  _) -> True
            Just (Char _) -> True
            _             -> False

    is_defined key _ = isJust . lookup key . b_exts . unpackBam


-- | Numeric-valued atoms (float or int, everything is converted to
-- double).  This includes:
-- - tagged fields (default to 0 if missing or wrong type)
-- - predefined fields: RNAME, POS, MRNM, MPOS, ISIZE, LENGTH, MAPQ
-- - mnemonic constants: wrongness, unexpectedness, rgquality
-- - literals

p_num_atom :: A.Parser (Expr Double)
p_num_atom = A.choice
    [ from_num_field                                        <$> p_field_name
    , const (fromIntegral . unRefseq . b_rname . unpackBam) <$  literal "RNAME"
    , const (fromIntegral .              b_pos . unpackBam) <$  literal "POS"
    , const (fromIntegral . unRefseq .  b_mrnm . unpackBam) <$  literal "MRNM"
    , const (fromIntegral . unRefseq .  b_mrnm . unpackBam) <$  literal "RNEXT"
    , const (fromIntegral .             b_mpos . unpackBam) <$  literal "MPOS"
    , const (fromIntegral .             b_mpos . unpackBam) <$  literal "PNEXT"
    , const (fromIntegral .            b_isize . unpackBam) <$  literal "ISIZE"
    , const (fromIntegral .            b_isize . unpackBam) <$  literal "TLEN"
    , const (fromIntegral . unQ .       b_mapq . unpackBam) <$  literal "MAPQ"
    , const (fromIntegral . V.length .   b_seq . unpackBam) <$  literal "LENGTH"

    , const (fromIntegral . extAsInt    0 "Z0" . unpackBam) <$  literal "unknownness"
    , const (fromIntegral . extAsInt 9999 "Z1" . unpackBam) <$  literal "rgquality"
    , const (fromIntegral . extAsInt    0 "Z2" . unpackBam) <$  literal "wrongness"
    , const . const                                         <$> A.double <* A.skipSpace ]
  where
    from_num_field key _ br = case lookup key (b_exts (unpackBam br)) of
            Just (Int   x) -> fromIntegral x
            Just (Float x) -> float2Double x
            _              -> 0

literal :: Bytes -> A.Parser ()
literal key = A.string key *> A.skipSpace

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
    , const . const     <$> p_string_lit ]
  where
    from_string_field key _ br = case lookup key (b_exts (unpackBam br)) of
            Just (Text x) -> x
            Just (Bin  x) -> x
            Just (Char x) -> B.singleton x
            _             -> B.empty

    lookupRef f m b = sq_name . getRef (meta_refs m) . f . unpackBam $ b

    get_library m = \b -> let lb = extAsString "LB" $ unpackBam b
                          in if B.null lb then H.lookupDefault B.empty (extAsString "RG" $ unpackBam b) lbs else lb
      where
        !lbs = H.fromList [ (rg,lb) | ("RG", stuff) <- meta_other_shit m
                                    , ("LB", lb) <- stuff
                                    , ("ID", rg) <- stuff ]

    get_sample m = \b -> H.lookupDefault B.empty (extAsString "RG" $ unpackBam b) sms
      where
        !sms = H.fromList [ (rg,sm) | ("RG", stuff) <- meta_other_shit m
                                    , ("SM", sm) <- stuff
                                    , ("ID", rg) <- stuff ]


-- | String literal.  We'll need escapes, but for now this should do.
p_string_lit :: A.Parser Bytes
p_string_lit =     A.char '"' *> A.takeWhile (/='"') <* A.char '"' <* A.skipSpace
               <|> A.char '\'' *> A.takeWhile (/='\'') <* A.char '\'' <* A.skipSpace

