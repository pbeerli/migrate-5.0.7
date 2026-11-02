# Migrate 6.1.79

Released November 2026

Migrate estimates effective population sizes, migration rates, divergence
times and related parameters of n populations from genetic data. It uses
coalescent theory and integrates over genealogies and mutation histories
with Markov chain Monte Carlo (Bayesian inference). Results are written to a
text outfile and to a PDF report with tables, posterior histograms and
skyline plots; the PDF is tagged for accessibility (PDF/UA-1).

Content
-------
- Quick installation (source code)
- Data types
- Models and analyses
- Running the examples
- Parallel runs (MPI)
- Converting data
- History, documentation, guidance
- Contact
- License


Quick installation (source code)
--------------------------------
Requirements: a C compiler, autoconf, automake and pkg-config; optional
OpenMPI (for migrate-n-mpi). zlib and libharu (PDF output) are bundled and
used when no system version is found.

    git clone https://github.com/pbeerli/migrate-5.0.7.git
    cd migrate-5.0.7
    autoreconf -i
    ./configure            # or ./configure --enable-debug
    make

This creates src/migrate-n and, when mpicc is found, src/migrate-n-mpi.

Installing the build tools:

    sudo apt-get install autoconf automake pkg-config   # Ubuntu/Debian
    brew install autoconf automake pkg-config           # macOS (Homebrew)

To put the binaries on your search path:

    sudo make install          # into /usr/local/bin

or, without administrator rights:

    ./configure --prefix=$HOME
    make
    make install               # into $HOME/bin; make sure it is on your PATH

If the build fails, please contact me (see below).


Data types
----------
- DNA sequences: finite-sites models from JC69 to Tamura-Nei (JC69, K2P,
  F81, F84, HKY, TN93), with gamma-distributed rate variation among sites
  and optional estimation of sequencing error.
- SNP data (linked SNPs derived from sequences; for resequencing data the
  full sequences are usually the better choice).
- Microsatellites: Brownian motion approximation or stepwise mutation
  model (the stepwise model is slow); repeat numbers or fragment lengths.
- Electrophoretic/allelic data (infinite allele model).
- Unphased diploid data: haplotype phases can be integrated over
  (haplotyping; use heating).


Models and analyses
-------------------
- Structured coalescent: mutation-scaled population sizes (Theta) and
  immigration rates (M) for a full or custom migration matrix (fixed,
  symmetric or tied parameters), n-island models and geographic distance.
- Population divergence (split times) with or without migration.
- Population growth.
- Fractional (Mittag-Leffler) coalescent: one alpha per population for
  environmentally heterogeneous reproduction; time scale chosen with
  mittag-leffler-scale=SITE|LOCUS|INVARIANT.
- Skyline plots of parameters through time; skyline=PARAM estimates
  piecewise-constant multipliers of Theta and M.
- Dated samples (tip dates) with a mutation rate per generation and
  generations per year.
- Rate variation among loci, inheritance scalars (e.g. mtDNA, X, Y).
- Assignment of individuals to populations.
- Multi-locus combination: 'All' (product of per-locus posteriors) and
  'Joint' (joint combination with Monte Carlo error).
- Marginal likelihoods for Bayes factor model comparison (thermodynamic
  integration with Bezier correction, stepping stone, harmonic mean); use
  static heating (heating=YES:1:LINEAR:N or POWER:N) for model comparison.
- Convergence diagnostics: effective sample sizes and Gelman-Rubin
  statistic across replicates.
- Relabeling/merging of populations and random subsets of individuals.


Running the examples
--------------------
The example directory contains data sets with matching parmfiles
(twoswisstowns, assign, gap, hapmap, msat, snp). For example:

    cd example
    ../src/migrate-n parmfile.twoswisstowns -nomenu

Without -nomenu migrate starts with an interactive menu. Start with a small
data set and a simple model before running large analyses (see GUIDANCE).


Parallel runs (MPI)
-------------------
migrate-n-mpi distributes loci (and replicates) over MPI processes:

    mpirun -np 11 migrate-n-mpi parmfile -nomenu

One process is the master; more workers than loci times replicates are not
used. On clusters, compile with the cluster's mpicc
(./configure MPICC=/path/to/mpicc).


Converting data
---------------
The contribution directory has conversion scripts, among them
vcf2mig.py (VCF to migrate format; see python vcf2mig.py --help),
structure2mig.py, mig2snp.py, fasta2genotype.py and stacks2migrate.py.


History, documentation, guidance
--------------------------------
- HISTORY lists bug fixes and new features by version.
- GUIDANCE gives advice when migrate does not deliver the answers you
  expect.
- The manual (migratedoc.pdf) and tutorials are available from
  https://peterbeerli.com/migrate-html5/index.html


Contact
-------
Peter Beerli
beerli@fsu.edu
or migrate-support@googlegroups.com


License
-------
Copyright 1997-2026 Peter Beerli

MIT open source license:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies
or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Last update: November 2, 2026
