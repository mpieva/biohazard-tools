set -e

cabal sandbox init

echo -e "\e[91mTesting with GHC 7.8\e[0m"
cabal install -w ghc-7.8.4 --only-dep --enable-tests --force-reinstalls
cabal configure -w ghc-7.8.4 --enable-library-profiling --disable-executable-profiling -O0
cabal build
cabal test
cabal clean
echo -e "\e[32mDone with GHC 7.8\e[0m"

echo -e "\e[91mTesting with GHC 7.10\e[0m"
cabal install -w ghc-7.10.1 --only-dep --enable-tests --force-reinstalls
cabal configure -w ghc-7.10.1 --disable-library-profiling --disable-executable-profiling -O0
cabal build
cabal test
cabal clean
echo -e "\e[32mDone with GHC 7.10\e[0m"

echo -e "\e[91mTesting with GHC 8.0\e[0m"
cabal install -w ghc-8.0.1 --only-dep --enable-tests --force-reinstalls
cabal configure -w ghc-8.0.1 --disable-library-profiling --disable-executable-profiling -O0
cabal build
cabal test
cabal haddock
cabal clean
echo -e "\e[32mDone with GHC 8.0\e[0m"
