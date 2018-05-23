biohazard-tools
===============

This is a collection of command line utilities that do useful stuff
involving BAM files for Next Generation Sequencing data:

* `bam-fixpair`: brings mates from paired end runs together and fixes
  their flags

* `bam-mangle`: filters BAM files through a flexible boolean expression 
  language.  A little bit like `awk` for BAM.

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

Building the "known good" version.

Install the following prerequisites:

- [Haskell Stack](https://haskellstack.org).
- The [snappy](https://github.com/google/snappy) library from Google.
  (This may well be in the operating system's packages.)

Clone the `biohazard` and `biohazard-tools` repositories into the same
directory:

    git clone https://vcs.eva.mpg.de/git_bioinformatics/biohazard.git
    git clone https://vcs.eva.mpg.de/git_bioinformatics/biohazard-tools.git

Change into the `biohazard` directory and check out the `0.6.15` tag:

    git checkout 0.6.15

Change into the `biohazard-tools` directory and use Stack to build:

    stack build
