biohazard-tools
===============

This is a collection of command line utilities that do useful stuff
involving BAM files for Next Generation Sequencing data:

* `bam-fixpair`: brings mates from paired end runs together and fixes
  their flags

* `bam-meld`: melds multiple BAM files together, retaining only the best
  alignment for each read

* `bam-rewrap`: wraps alignments around the origin of a circular reference.
  (Use this to get sensible alignments to mitochondria.  `bam-rmdup` 
  includes similar functionality.)

* `bam-rmdup`: removes PCR duplicates from BAM files and computes a
  consensus sequence for each cluster of dups

* `fastq2bam`: converts FastQ to BAM.  It handles many weird
  combinations of one to four input files and optionally trims reads and
  merges read pairs.

Installation
------------

`biohazard-tools` uses Cabal, the standard installation mechanism for
Haskell.  It depends on the `biohazard` library and additional stuff
from Hackage.  To install, follow these steps:

* install a useable Haskell environment, either
 * install the Haskell Platform or have someone install it (see
   http://haskell.org/platform/ for details), or
 * install GHC (see http://haskell.org/ghc) and bootstrap Cabal (see
   http://haskell.org/cabal), and
 * run `cabal update` (takes a while to download the current package list),
* run `cabal install
  https://bitbucket.org/ustenzel/biohazard/get/0.6.10.tar.gz
  https://bitbucket.org/ustenzel/biohazard-tools/get/0.1.tar.gz`
  (takes even longer).

When done, on an unmodified Cabal setup, you will find the binaries in 
`${HOME}/cabal/bin`.  Cabal can install them in a different place, please 
refer to the Cabal documentation at http://www.haskell.org/cabal/ if 
you need that.  Sometimes, repeated installations and re-installations can result 
in a thoroughly unusable state of the Cabal package collection.  If you get error 
messages that just don't make sense anymore, please refer to 
http://www.vex.net/~trebla/haskell/sicp.xhtml; among other useful things, it 
tells you how to wipe a package database without causing more destruction.
