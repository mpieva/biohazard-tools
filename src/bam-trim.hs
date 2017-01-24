import Bio.Bam
import Bio.Prelude
import Paths_biohazard_tools                ( version )
import System.Console.GetOpt

data Conf = Conf { c_trim_pred :: [Nucleotides] -> [Qual] -> Bool
                 , c_pass_pred :: BamRec -> Bool }

options :: [OptDescr (Conf -> IO Conf)]
options = [ Option "q"  ["minq"]    (ReqArg set_minq "Q") "Trim where quality is below Q"
          , Option "m"  ["mapped"]  (NoArg set_monly)     "Trim only mapped sequences"
          , Option "h?" ["help"]    (NoArg usage)         "Display this text and exit"
          , Option "V"  ["version"] (NoArg vrsn)          "Display version number and exit" ]

set_minq :: String -> Conf -> IO Conf
set_minq s c = readIO s >>= \q -> return $ c { c_trim_pred = trim_low_quality (Q q) }

set_monly :: Conf -> IO Conf
set_monly c = return $ c { c_pass_pred = \r -> isMerged r || isUnmapped r }

vrsn :: Conf -> IO Conf
vrsn _ = do pn <- getProgName
            hPutStrLn stderr $ pn ++ ", version " ++ showVersion version
            exitSuccess

usage :: Conf -> IO Conf
usage _ = do hPutStrLn stderr $ usageInfo info options ; exitSuccess
  where info = "Simple trimming of sequences in Bam files.  Reads a Bam file from stdin,\n\
               \trims sequences of low quality, writes Bam to stdout.  Does not trim\n\
               \merged reads."


main :: IO ()
main = do
    (opts, files, errors) <- getOpt Permute options `fmap` getArgs

    unless (null errors) $ mapM_ (hPutStrLn stderr) errors
    c <- foldM (flip id) (Conf (trim_low_quality (Q 20)) isMerged) opts
    unless (null errors && null files) exitFailure

    let do_trim r | c_pass_pred c r' = Left r
                  | otherwise        = Right $ trim_3' (c_trim_pred c) r'
            where r' = unpackBam r

    add_pg <- addPG (Just version)
    concatDefaultInputs >=> run $
        joinI . mapStream do_trim .
        protectTerm . pipeBamOutput . add_pg

