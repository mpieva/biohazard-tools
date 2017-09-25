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
* run `cabal install
  https://bitbucket.org/ustenzel/biohazard/get/master.tar.gz
  https://bitbucket.org/ustenzel/biohazard-tools/get/master.tar.gz`
  (takes even longer).

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

This is old code that annotates stuff in weird legacy formats (formerly
called `coord2anno`).  It suffered bitrot and is arguably best replaced
with `Bedtools`.

The only reason this could ever be useful is because it has a
client-server mode, which is handy if many small files need to be
processed and reading the annotation file becomes the bottleneck.
(Annotation files are frequently big and rather expensive to parse.)

A reimagination would modernize the code, use a broadcast protocol to
find the server, and make C/S mode completely transparent, implicit and
optional.  The lookup could even be integrated into bam-mangle to make
it more flexible.  But for now, the code is barely kept alive.

* make `expound` compile again
 * short byte strings as keys in maps
 * hashmaps?
 * support Bam again
* edit `expound.1` to reflect renaming, redesign, whatever
* do something about the silly translation tables
* make Judy optional, when in doubt by not compiling `expound`
