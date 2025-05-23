Migrate 5.0.7
--------------
Released Summer 2025

Quick installation guide 
------------------------
- Unpack the compressed distribution file
- in the top directory do:
  ./configure
  or
  ./configure --enable-debug
  and then do
  make
  this should create the binary and
  if you have openmpi installed: it should create also the parallel version migrate-n-mpi
  sudo make install
  will move the binaries to /usr/local/bin,
  if you do not have access to that use
  ./configure --prefix=/pathtoyourhomedirectory
  make
  
  if I have time I will distribute binaries for mac and windows, for unix (and mac) the source install should
  work best.
  

Content
-------
- Overview
- Analyses summary
- Computer systems
- History of MIGRATE
- Distribution
- Installation
- Documentation
- Example folder
- Troubleshoothing
- Disclaimer


Migrate estimates population parameters, effective population sizes
migration rates and divergence times of n populations, using genetic data.  
It uses a coalescent theory approach taking into account history of 
mutations and uncertainty of the genealogy. 

The estimates of the parameter values are achieved by Bayesian inference (BI).
The output is presented in an TEXT file _and_ in a PDF file. The PDF file
contains all currently main tables, posterior histograms (BI) and 
skylineplots are supported in the PDF.

If you fail to see all plots in the outfile PDF, consider to use a 
non-adobe PDF reader. Recent Adobe Reader fails on some output files 
reporting empty histograms; other PDF-readers such as Preview.app on 
mac or Nitro PDF reader are fine.


 
Currently the following data types are supported:
-DNA sequence data
 - finite sites model: Tamura-Nei and all models that can be 
   expressed as a subset of Tamura-Nei
 - finite sites model + rate variation among sites: TN + Gamma

-SNP data (single nucleotide polymorphism)
 - SNP are derived from sampled sequences and are completely linked
   except that we know that the sites are variable, with resequencing project data
   I suggest NOT to use SNPS but the full sequences, the Bayes factor methods allows
   to break down into independent loci easily. 

-Microsatellite data 
 - Brownian motion model: a continuous approximation to the stepwise mutation model.
 - stepwise mutation model (if you want to finish your work in this century do not use this))

-Electrophoretic marker data (infinite allele model).


Analyses
--------
[+ working and described
 - experimental and needs further description]

+ Estimation of population sizes and migration rates of a migration matrix
  model, or arbitrarily subsets of a migration matrix model, or an n-island
  model. Allowing for a geographic distance matrix so that geographic effects 
  can be removed out of the analysis, Estimation of mutation-scaled divergence
  time among populations.
+ Marginal likelihood calculation to assist calculation of Bayes Factor (BI)
- Allows a variable mutation rate AMONG loci estimated from the data.
- For sequences: allows a variable substitution rate among sites.
+ For microsatellites: allows the definition of repeatnumber and use of fragment length as input
+ Facilitates analyses of multimodal search space distributions with heating
  scheme and/or multi-run analyses.
- Histogram of events over time 
- Plot of expected parameters through time (skyline plots) for all parameters. 
- Dated samples
+ Relabeling and merging of populations
+ Random subset of individuals per population
- Assignment of indiviudals to populations (population number needs to be known)
- Haplotype assignment (for small numbers of different haplotypes) and integration over all
  potential haplotype


Computer systems
----------------
You can fetch Migrate from the website http://popgen.sc.fsu.edu
or http://peterbeerli.com/migrate-html5
as source code or binary executables. Currently I supply binaries for

- Macintosh:    migrate-n
- Windows: migrate-n.exe 

The source code should compile on all platforms (windows may be tricky)
- single CPU version (with or without using AVX or threading (use configure and make)
- parallel version using either openpmpi or mpich2 (compiles to migrate-n-mpi) 


The file is compressed as tar.gz or as zip file.
 	
The documentation contains information about
how to compile and use a parallelized version of migrate 
so that it can run concurrently on computer clusters
(using MPI [preferrably OpenMPI]).


History about bug fixes and new features
----------------------------------------
read the HISTORY file.

Distribution
------------
Migrate can be fetched from the www-site 
http://popgen.sc.fsu.edu/ 
or https://peterbeerli.com/migrate-html5/index.html

Installation
------------
(a) Binaries
Unpack the compressed archive, open the directory migrate-4.4.x
- Mac: from the binary distribution copy migrate-n to /usr/local/bin
  and then use the Terminal.app 
- Windows: the preferred way to use migrate-n.exe is through the commandline 
  environment.
- UNIX: I suggest to use the source distribution to make the best use of your system, then place
  the binary in /usr/local/bin or some other directory that is searched and then call migrate-n from
  the commandline.






Compiling the Source
---------------------
Download the program and read the documentation and try the program
on a small [!] data set.
 
For UNIX systems
the binary can go to standard directories (e.g. /usr/local/bin),
the rudimentary man page can go to the /usr/local/man/man1.

(b) Source (UNIX and MACs)
1. unpack the distribution file
2. cd migrate-VERSION/src
3. type "./configure"
   This will create the Makefile 
   [./configure --help for more options]
 	
4. type "make" (please report warnings and especially errors).
   If you have a multiprocessor machine (non-Macintosh computer!) you perhaps want to try 
   "make thread" (this allows parallel execution of chains when 
   using the heating scheme). 
   If you compile on a Macintosh with INTEL CPU and MACOS 10.10+ try
   "make". On Macs configure will use automatically use the built-in 
   threading framework (GrandCentral) that will be faster than the standard threading.

   The result of the compilation should be an executable 
   "migrate-n" in the current directory [it is called "migrate-n" because
   on some computer system there is a system program called "migrate"]

5. sudo make install 
   This will install the programs and man-page into usr/local/bin, 
   /usr/local/man/man1
   [you need to be root or administrator to do this; this step is NOT necessary, 
   to use the program, but it would be convenient for all users
   of your system]
   if this does not work: 
   move migrate-n to /usr/local/bin or $HOME/bin or some 
   other convenient place.

6. change directory to "example"
   run "migrate-n parmfile.testbayes", 
   on my mac this takes a few minutes.
   If this test fails, please let me know!


Documentation
-------------
You need to download it separately from
http://popgen.sc.fsu.edu/Downloads.html

It is a PDF file and called  migratedoc.pdf.
The pdf file can be viewed and printed using Acrobat 
or any other PDF viewer. New versions of Acrobat have issues reading the 
PDF output file from migrate, if you see unfinished graphs in your output PDF
use another viewer (for example NITRO PDF). On macs use Preview.app (it is anyway better than Acrobat Reader).


Examples
--------
In the directory "example" you can find some example data sets.
You might wan to try the two parmfiles.testml and parmfile.testbayes
Use the Terminal.app (on macosx), or xterm (Linux), or cmd (Windows),
change directory to the example directory and then execute
for mac and unix: ../migrate-n parmfile.testbayes
for windows: ..\migrate-n parmfile.testbayes
or the ml version.

Guidance
--------
Look at the file GUIDANCE for some help when migrate is not delivering the answers you think it should deliver.

Fan-mail, complaints, questions and error-reports
-------------------------------------------------
Peter Beerli
beerli@fsu.edu
or
migrate-support@googlegroups.com


Disclaimers
-----------
Copyright 1997-2021 Peter Beerli
 
** MIT OPENSOURCE LICENSE*****************************************************************
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
*  of the Software, and to permit persons to whom the Software is furnished to do
*  so, subject to the following conditions:
* 
*  The above copyright notice and this permission notice shall be included in all copies
*  or substantial portions of the Software.
* 
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
*  INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
*  PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
*  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
*  CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
*  OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
* 

Last update:
March 13  2021









