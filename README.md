biohazard-tools
===============

This is a collection of command line utilities that do useful stuff
involving BAM files for Next Generation Sequencing data:

* `bam-fixpair`: brings mates from paired end runs together and fixes
  their flags.  Can also pipe into programs that expect _two_ or 
  _three_ FastQ files as input.

* `bam-mangle`: filters BAM files through a flexible boolean expression 
  language.  A little bit like `find` (but not quite like `awk`) for BAM.

* `bam-meld`: melds multiple BAM files together, retaining only the best
  alignment for each read.

* `bam-rewrap`: wraps alignments around the origin of a circular reference.
  (Use this to get sensible alignments to mitochondria.  `bam-rmdup` 
  includes similar functionality.)

* `bam-rmdup`: removes PCR duplicates from BAM files and computes a
  consensus sequence for each cluster of replicates.

* `expound`: annotates regions, similar to what `bedtools` does.
  `expound` has a neat client/server mode so it doesn't have to read the
  same annotation file over and over again.

* `fastq2bam`: converts FastQ to BAM.  It handles many weird
  combinations of one to four input files and optionally trims reads and
  merges read pairs.

Installation
------------

`biohazard-tools` uses Cabal, the standard installation mechanism for
Haskell.  It depends on the `biohazard` library and additional stuff
from Hackage.  To install, follow these steps:

* install GHC (see http://haskell.org/ghc) and Cabal (see
   http://haskell.org/cabal),
* run `cabal update` (takes a while to download the current package list),
* install Judy, either from http://judy.sourceforge.net, or through your
  package manager (it's named `libjudy-dev` or similar),
* run `cabal install
  https://bitbucket.org/ustenzel/biohazard/get/master.tar.gz
  https://bitbucket.org/ustenzel/biohazard-tools/get/master.tar.gz`
  (takes even longer).
* If you don't have or don't like Judy, you can install without it by
  running `cabal install -f -judy
  https://bitbucket.org/ustenzel/biohazard/get/master.tar.gz
  https://bitbucket.org/ustenzel/biohazard-tools/get/master.tar.gz`.
  The performance of `expound` may suffer a bit.


When done, on an unmodified Cabal setup, you will find the binaries in 
`${HOME}/cabal/bin`.  Cabal can install them in a different place, please 
refer to the Cabal documentation at http://www.haskell.org/cabal/ if 
you need that.  Sometimes, repeated installations and re-installations can result 
in a thoroughly unusable state of the Cabal package collection.  If you get error 
messages that just don't make sense anymore, please refer to 
http://www.vex.net/~trebla/haskell/sicp.xhtml; among other useful things, it 
tells you how to wipe a package database without causing further destruction.


Expound
=======

This is a reimagination of old code.  We modernize the code, use a
broadcast protocol to find the server, and thereby make C/S mode
completely transparent, implicit and optional.  Once it runs, it's
integrated into bam-mangle to make it more flexible.  But for now, the
code is barely kept alive.

To do:

* support Bam again
 * the 'Region' type might work better with strict ByteStrings now
 * maybe lazy ByteStrings become unnecessary
* see if C/S still works at all
* add the broadcast discovery
* integrate with `bam-mangle`
* document --output-patterns properly
