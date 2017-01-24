import Control.Exception                    ( try, IOException )
import Distribution.PackageDescription      ( PackageDescription(..) )
import Distribution.Simple
import Distribution.Simple.InstallDirs      ( mandir, CopyDest (NoCopyDest) )
import Distribution.Simple.LocalBuildInfo   ( LocalBuildInfo(..), absoluteInstallDirs )
import Distribution.Simple.Setup            ( copyDest, copyVerbosity, fromFlag, installVerbosity )
import Distribution.Simple.Utils
import Distribution.Verbosity               ( Verbosity )
import System.Exit                          ( exitSuccess )
import System.FilePath                      ( splitDirectories, joinPath, (</>) )

main :: IO ()
main = do
  defaultMainWithHooks $ simpleUserHooks
    { postCopy = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ copyVerbosity flags) (fromFlag $ copyDest flags)

    , postInst = \ _ flags pkg lbi ->
         installManpages pkg lbi (fromFlag $ installVerbosity flags) NoCopyDest
    }
  exitSuccess

installManpages :: PackageDescription -> LocalBuildInfo -> Verbosity -> CopyDest -> IO ()
installManpages pkg lbi verbosity copy = do
    installOrdinaryFiles verbosity (mandir (absoluteInstallDirs pkg lbi copy))
        [ ("man", joinPath mp) | ("man":mp) <- map splitDirectories $ extraSrcFiles pkg ]

installOrdinaryFiles' :: Verbosity -> FilePath -> [(FilePath, FilePath)] -> IO ()
installOrdinaryFiles' verb dest = mapM_ go
  where
    go :: (FilePath, FilePath) -> IO (Either IOException ())
    go (base,src) = try $ installOrdinaryFile verb (base </> src) (dest </> src)

