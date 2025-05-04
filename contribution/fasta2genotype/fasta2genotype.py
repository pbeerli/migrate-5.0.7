#!/usr/bin/env python
#
# (c) Paul Maier
# May 2, 2016
# fasta2genotype.py
# V 1.9
# Written for Python 2.7.5
#
# This program takes a fasta file listing all alleles of all individuals at all loci
# as well as a list of individuals/populations and list of variable loci
# then outputs sequence data in one of five formats:
# (1) migrate-n, (2) Arlequin, (3) DIYabc, (4) LFMM, (5) Phylip, (6) G-Phocs, or (8) Treemix
# (8) Additionally, the data can be coded as unique sequence integers (haplotypes)
#     in Structure/Genepop/SamBada/Bayescan/Arlequin/GenAlEx format
#     or summarized as allele frequencies by population
#
# Execute program in the following way:
# python fasta2genotype.py [fasta file] [whitelist file] [population file] [VCF file] [output name]
#
# The fasta headers should be formatted exactly the way done by populations.pl in Stacks v. 1.12
# e.g. >CLocus_#_Sample_#_Locus_#_Allele_#
# Where the CLocus is the locus number in the catalog produced by cstacks
# Sample is the individual number from denovo_map.pl, in the same order specified with -s flags
# Locus refers to the locus number called for that individual; this information is not used
# Allele must be 0 or 1; this script assumes diploid organisms; homozygote is assumed if only 0
# Restriction enzyme sites can be removed if requested, for either single- or double-digest setups.
# Simply select this option and provide the 5' sequence(s). It might help to examine the fasta file.
#
# The list of "white" loci is a column vector of catalog loci that are variable and should thus be retained
# This can be extracted from many of the stacks output files, including the vcf and structure files
# If all loci are desired, include numbers for all loci and select option to retain all.
#
# The populations file should have three columns:
# [1] Sample ID from stacks, [2] Individual ID, [3] Population ID
# These can be named anything, e.g.
# SampleID	IndID	        PopID
# 1	        Y120037	        Wawona
# 2	        Y120048	        Westfall
# etc...
# Use tabs, don't use any spaces.
# Individual IDs must be 9 characters or shorter for [1] Migrate option.
# Population IDs must start with an alphabetical character (not a number) for [1] Migrate option.
#
# Loci can be pruned using a coverage threshold, or missing data threshold.
# A .vcf file from the SAME run of populations.pl must be provided in this case.
# If no pruning is desired, simply use 'NA' as the input argument.
# e.g. python fasta2genotype.py ./fasta ./whitelist ./populations NA MyOutput
# If pruning is desired, IndividualID spelling in VCF file must match that in populations file.
#       This also means there should be the same list of IndividualIDs in both files.
# This program assumes the same VCF format found in Stacks v. 1.12
#       Data should consist of three elements separated by colons: genotype, coverage, likelihood
# Note: Coverage info is not provided for consensus (non-variable) sequences by STACKS
#       Therefore, these loci are removed if coverage pruning is chosen.
#       This will be true even if "All Loci" is chosen.
# Additionally, the LFMM format relies on SNP genotypes found in the VCF file.
#       A VCF file must be included for this option.
#
# In addition to coverage filtering, several other quality control measures can be selected.
#       Alleles under a locus-wide frequency can be removed
#       Alleles represented under a given frequency of populations can be removed (e.g. 4 pops of 16, frequency of 0.25)
#       Monomorphic loci can be removed.
#       Loci under a given locus-wide missing data threshold can be removed.
#       Loci can be removed from populations where those populations fall under a missing data threshold for those loci
#       Individuals missing a threshold proportion of loci can be removed.
#
# TO DO:
# Speed up 'remove by coverage' part of Seqs ... way too slow
# Allow LFMM format to operate free of VCF file, (this has been accomplished already for Phylip format)
# Allow Phylip format to summarize at the population, instead of individual, level
# Allow order of individuals and pops to be specified from input files
#
#
#
#
#
# What follows is a brief demo of a single sample dataset in all available output formats.
#
# Migrate file format:
#
# <Number of populations> <number of loci>
# <number of sites for locus1> <number of sites for locus 2> ...
# <Number of gene copies> <title for population>
# Ind1a <locus 1 gene copy 1 sequence>
# Ind1b <locus 1 gene copy 2 sequence>
# Ind2a <locus 1 gene copy 1 sequence>
# Ind2b <locus 1 gene copy 2 sequence>
# Ind1a <locus 2 gene copy 1 sequence>
# Ind1b <locus 2 gene copy 2 sequence>
# Ind2a <locus 2 gene copy 1 sequence>
# Ind2b <locus 2 gene copy 2 sequence>
# <Number of gene copies> <title for population>
# Ind1a <locus 1 gene copy 1 sequence>
# Ind1b <locus 1 gene copy 2 sequence>
# Ind2a <locus 1 gene copy 1 sequence>
# Ind2b <locus 1 gene copy 2 sequence>
# Ind1a <locus 2 gene copy 1 sequence>
# Ind1b <locus 2 gene copy 2 sequence>
# Ind2a <locus 2 gene copy 1 ????????>
# Ind2b <locus 2 gene copy 2 ????????>
# ...   ...
#
#
# Arlequin file format:
#
# [Profile]
#       "ProjectName"
#              Project
#              NbSamples=#
#              GenotypicData=1
#              GameticPhase=0
#              DataType=DNA
#              LocusSeparator=TAB
#              MissingData="?"
# [Data]
#       [[Samples]]
#               SampleName="PopID"
#               SampleSize=#
#               SampleData={
# Ind1  1       <Locus 1 copy 1>        <Locus 2 copy 1> ...
#               <Locus 1 copy 2>        <Locus 2 copy 2> ...
# Ind1  2       <Locus 1 copy 1>        <Locus 2 copy 1> ...
#               <Locus 1 copy 2>        <Locus 2 copy 2> ...
# }
#               SampleName="PopID"
#               SampleSize=#
#               SampleData={
# Ind1  1       <Locus 1 copy 1>        <Locus 2 copy 1> ...
#               <Locus 1 copy 2>        <Locus 2 copy 2> ...
# Ind1  2       <Locus 1 copy 1>        <??????????????> ...
#               <Locus 1 copy 2>        <??????????????> ...
# ...   ...     ...                     ...
# }
#
#
# DIYabc file format:
#
# <Project Name> <NF=NF>
# Locus 1 <A>
# Locus 2 <A>
# Pop
# Ind1  ,       <[Locus 1 copy 1][Locus 1 copy 2]>      <[Locus 2 copy 1][Locus 2 copy 2]> ...
# Ind2  ,       <[Locus 1 copy 1][Locus 1 copy 2]>      <[Locus 2 copy 1][Locus 2 copy 2]> ...
# Pop
# Ind1  ,       <[Locus 1 copy 1][Locus 1 copy 2]>      <[Locus 2 copy 1][Locus 2 copy 2]> ...
# Ind2  ,       <[Locus 1 copy 1][Locus 1 copy 2]>      <[??????????????][??????????????]> ...
# ...           ...              ...                    ...              ...
#
#
# LFMM format (assuming 1 snp at first locus, 2 snps at second locus):
#
# Ind1  Pop1    0       0       0       0       A       A       C       T       G       G ...
# Ind2  Pop1    0       0       0       0       T       T       C       C       G       G ...
# Ind1  Pop2    0       0       0       0       A       T       C       T       T       T ...
# Ind1  Pop2    0       0       0       0       T       T       0       0       0       0 ...
# ...   ...     ...     ...     ...     ...     ...     ...     ...     ...     ...     ...
#
#
# Phylip format (assuming 1 snp at first locus, 2 snps at second locus):
# (Only SNPs will be output, and loci will be concatenated in the order supplied in the whitelist.)
#
# <number of individuals> <number of base pairs>
# Ind1_Pop1 AYG ...
# Ind2_Pop1 TCG ...
# Ind1_Pop2 WYT ...
# Ind1_Pop2 TNN ...
# ...       ...    ...
#
#
# G-Phocs format:
#
# <Number of loci>
#
# <Locus name> <number of individuals> <number of sites for locus 1>
# Ind1 <locus 1 sequence>
# Ind2 <locus 1 sequence>
# Ind1 <locus 1 sequence>
# Ind2 <locus 1 sequence>
#
# <Locus name> <number of individuals> <number of sites for locus 2>
# Ind1 <locus 1 sequence>
# Ind2 <locus 1 sequence>
# Ind1 <locus 1 sequence>
# Ind2 <locus 1 ????????>
# ...   ...
#
#
# Treemix format:
#
# Pop1     Pop2
# 2,2      1,3 ...
# 3,1      1,1 ...
# 4,0      2,0 ...
# ...      ...
#
#
# Structure format:
#
#               Loc1    Loc2 ...
# Ind1  Pop1    1       1 ...
# Ind1  Pop1    1       2 ...
# Ind2  Pop1    2       1 ...
# Ind2  Pop1    2       1 ...
# Ind1  Pop2    1       3 ...
# Ind1  Pop2    2       4 ...
# Ind2  Pop2    2       0 ...
# Ind2  Pop2    2       0 ...
# ...   ...     ...     ...
#
#
# Genepop format (six digit):
#
# <Project Name>
# Locus 1
# Locus 2
# Pop
# Ind1 ,  001001        001002 ...
# Ind2 ,  002002        001001 ...
# Pop
# Ind1 ,  001002        003004 ...
# Ind2 ,  002002        000000 ...
#
#
# Allele Frequency:
#
#       Loc1_1          Loc1_2          Loc2_1          Loc2_2          Loc2_3          Loc2_4 ...
# Pop1  0.50000         0.50000         0.75000         0.25000         0.00000         0.00000 ...
# Pop2  0.25000         0.75000         0.00000         0.00000         0.50000         0.50000 ...
# ...   ...             ...             ...             ...             ...             ...
#
#
# SamBada format:
#
#       Loc1_1  Loc1_2  Loc2_1  Loc2_2  Loc2_3  Loc2_4 ...
# Ind1a 1       0       1       0       0       0 ...
# Ind1b 1       0       0       1       0       0 ...
# Ind2a 0       1       1       0       0       0 ...
# Ind2b 0       1       1       0       0       0 ...
# Ind1a 1       0       0       0       1       0 ...
# Ind1b 0       1       0       0       0       1 ...
# Ind2a 0       1       0       0       NaN     NaN ...
# Ind2b 0       1       0       0       NaN     NaN ...
# ...   ...     ...     ...     ...     ...     ...
#
#
# Bayescan format:
#
# [loci]=2
#
# [populations]=2
#
# [pop]=1
# 1     4       2       2       2
# 2     4       4       3       1       0       0
#
# [pop]=2
# 1     4       2       1       3
# 2     2       4       0       0       1       1
#
#
# GenAlEx format:
#
# 2     4       2       2       2
#                       Pop1    Pop2
# IndID PopID   Loc1            Loc2
# Ind1  Pop1    1       1       1       2 ...
# Ind2  Pop1    2       2       1       1 ...
# Ind1  Pop2    1       2       3       4 ...
# Ind2  Pop2    2       2       0       0 ...
# ...   ...     ...     ...     ...     ...
#
#
#
#
import sys,re,csv,collections,itertools
from decimal import *

if len(sys.argv) != 6:
        print "Error: improper number of arguments. Please see program header for instructions."
        print "fasta2genotype.py [fasta file] [whitelist file] [population file] [VCF file] [output name]"
        exit(1)

outname = str(sys.argv[5])
outfile = outname + ".out"
outfile_loci = outname + "_loci.out"
outfile_pops = outname + "_pops.out"
listfile = outname + ".list"

choice = int(raw_input("Output type? [1] Migrate [2] Arlequin [3] DIYABC [4] LFMM [5] Phylip [6] G-Phocs [7] Treemix [8] Haplotype: "))
while choice not in range(1,9):
	choice = int(raw_input("Not a valid option. Conversion type? [1] Migrate [2] Arlequin [3] DIYABC [4] LFMM [5] Phylip [6] G-Phocs [7] Treemix [8] Haplotype: "))

if choice == 5:
        haplo = int(raw_input("Use only variable regions or full sequences? [1] Variable regions [2] Full Sequences : "))
        while haplo not in range (1,3):
                       haplo = int(raw_input("Not a valid option. How many SNPs to keep per locus? [1] Only one [2] All : "))

if choice == 5:
        breakpoints = int(raw_input("Flag break points between loci with '!' symbol in first sequence? [1] Yes [2] No : "))
        while breakpoints not in range (1,3):
                       breakpoints = int(raw_input("Not a valid option. Flag break points between loci with '!' symbols? [1] Yes [2] No : "))

if choice == 5:
        haplotypes = int(raw_input("Type of sequences for alignment? [1] Haploid [2] Diploid : "))
        while haplotypes not in range (1,3):
                       haplotypes = int(raw_input("Not a valid option. Type of sequences for alignment? [1] Haploid [2] Diploid : "))

if choice == 7:
        one_snp = int(raw_input("How many SNPs to keep per locus? [1] Only one [2] All : "))
        while one_snp not in range (1,3):
                       one_snp = int(raw_input("Not a valid option. How many SNPs to keep per locus? [1] Only one [2] All : "))

if choice == 8:
        HaploChoice = int(raw_input("Specific output type? [1] Structure [2] Genepop [3] AlleleFreqency [4] SamBada [5] Bayescan [6] Arlequin [7] GenAlEx : "))
        while HaploChoice not in range (1,8):
                       HaploChoice = int(raw_input("Not a valid option. Specific output type? [1] Structure [2] Genepop [3] AlleleFreqency [4] SamBada [5] Bayescan [6] Arlequin [7] GenAlEx :"))

outtype = int(raw_input("Loci to use? [1] Variable [2] All: "))
while outtype not in range(1,3):
	outtype = int(raw_input("Not a valid option. Loci to use? [1] Variable [2] All: "))

UseCoverage = 0
CoverageCutoff = 0
if sys.argv[4] != 'NA' and sys.argv[4] != 'na' and sys.argv[4] != 'Na' and sys.argv[4] != 'nA':
        UseCoverage = 1
        CoverageCutoff = int(raw_input("Coverage Cutoff (number reads for locus)? Use '0' to ignore coverage: "))
        while str(type(CoverageCutoff)) != "<type 'int'>":
                CoverageCutoff = int(raw_input("Not a valid option. Coverage Cutoff (number reads for locus)? Use '0' to ignore coverage: "))

allele_threshold = 0; allele_pop_threshold = 0
allele_filter = int(raw_input("Filter for allele frequency? False alleles might bias data. [1] Yes [2] No: "))
while allele_filter not in range(1,3):
	allele_filter = int(raw_input("Not a valid option. Filter for allele frequency? False alleles might bias data. [1] Yes [2] No: "))
if allele_filter == 1:
        allele_threshold = float(raw_input("Allele frequency threshold for removal across all individuals? Use '0' to ignore this: "))
        while str(type(allele_threshold)) != "<type 'float'>":
                allele_threshold = float(raw_input("Not a valid option. Allele frequency threshold for removal across all individuals? Use '0' to ignore this: "))
        allele_pop_threshold = float(raw_input("Frequency of populations containing allele for removal across all individuals? Use '0' to ignore this: "))
        while str(type(allele_pop_threshold)) != "<type 'float'>":
                allele_pop_threshold = float(raw_input("Not a valid option. Frequency of populations containing allele for removal across all individuals? Use '0' to ignore this: "))

monomorphic_filter = int(raw_input("Remove monomorphic loci? [1] Yes [2] No: "))
while allele_filter not in range(1,3):
	allele_filter = int(raw_input("Not a valid option. Remove monomorphic loci? [1] Yes [2] No: "))


missing_data_filter = int(raw_input("Filter for missing genotypes? These might bias data. [1] Yes [2] No: "))
while missing_data_filter not in range(1,3):
	missing_data_filter = int(raw_input("Not a valid option. Filter for missing genotypes? These might bias data. [1] Yes [2] No: "))
if missing_data_filter == 1:
        locus_threshold = float(raw_input("Locus frequency threshold for locus removal across all individuals? Use '0' to ignore this: "))
        while str(type(locus_threshold)) != "<type 'float'>":
                locus_threshold = float(raw_input("Not a valid option. Locus frequency threshold for removal across all individuals? Use '0' to ignore this: "))
        locus_pop_threshold = float(raw_input("Population frequency threshold for locus removal across each population? Use '0' to ignore this: "))
        while str(type(locus_pop_threshold)) != "<type 'float'>":
                locus_pop_threshold = float(raw_input("Not a valid option. opulation frequency threshold for locus removal across each population? Use '0' to ignore this: "))
        ind_threshold = float(raw_input("Individual frequency threshold for individual removal across all loci? Use '0' to ignore this: "))
        while str(type(ind_threshold)) != "<type 'float'>":
                ind_threshold = float(raw_input("Not a valid option. Individual frequency threshold for individual removal across all loci? Use '0' to ignore this: "))
	
clipcutsite = int(raw_input("Clip cut sites? These may bias data. [1] Yes [2] No: "))
while clipcutsite not in range(1,3):
	outtype = int(raw_input("Not a valid option. Clip cut sites? These may bias data. [1] Yes [2] No: "))

cutsite1 = ""
cutsite2 = ""
checkbothends = 0
if clipcutsite == 1:
	cutsite1 = str(raw_input("Sequence of first (5') cut site? "))
	while  str(type(cutsite1)) != "<type 'str'>":
		cutsite1 = str(raw_input("Not a valid option. Sequence of first (5') cut site? "))
	cutsite2 = str(raw_input("Sequence of second (3') cut site? (Leave blank if none): "))
	while  str(type(cutsite2)) != "<type 'str'>":
		cutsite2 = str(raw_input("Not a valid option. Sequence of second (3') cut site? (Leave blank if none): "))
	checkbothends = int(raw_input("Cut site found on 5' end only [1] or both ends [2]: "))
	while checkbothends not in range(1,3):
		checkbothends = int(raw_input("Not a valid option. Cut site found on 5' end only [1] or both ends [2]: "))



# Output new file containing only variable loci
def KeepLoci(fasta,whiteloci,keeploci):
        print "Creating new fasta from whitelist..."
        carrot = ">"
        whitelist = []

        try:
                fin=open(whiteloci,"U")
                for line in fin:
                        whitelist.append(str(line.strip()))
                fin.close()
        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)


        try:
                fin=open(fasta,"U")
                fout=open(keeploci,"w")
                for line in fin:
                        if carrot in line:
                                locinum = re.sub(r'>CLocus_(\w+)_Sample_\w+_Locus_\w+_Allele_\w+', r'\1', line)
                                locinum = locinum.strip()
                                if str(locinum) in whitelist:
                                        fout.write(line)
                                        nextline = fin.next()
                                        fout.write(nextline)

                fout.close()
                fin.close()

        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)

        return 1




# Create dictionary of populations and individuals
def Pops(populations):
        print "Cataloging populations..."
        try:
		pops = csv.DictReader(open(populations,"U"), delimiter="\t", quotechar='"', dialect="excel-tab")
        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)

        popsdict = {} #Structure: {Population : {SampleID : IndividualID} }
        for i in pops:
                if i[pops.fieldnames[2]] in popsdict.keys():
                        popsdict[i[pops.fieldnames[2]]][i[pops.fieldnames[0]]] = i[pops.fieldnames[1]]
		else:
			popsdict[i[pops.fieldnames[2]]] = {i[pops.fieldnames[0]]:i[pops.fieldnames[1]]}
        return popsdict




# Create dictionary of Coverage by individual X locus
def LocusCoverage(VCFfile):
        print "Calculating loci coverage..."
        try:
                with open(VCFfile,"U") as f:
                        cov = csv.reader(f, delimiter="\t")
                        d = list(cov)

        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)
        covdict = {} # Structure: {IndividualID : {LocusID : Coverage} }
        rindex = 0
        for row in d:
                if rindex < 9:
                        rindex += 1
                        continue
                cindex = 0
                locus = int(d[rindex][2])
                for column in row:
                        ind = d[8][cindex]
                        if cindex < 9:
                                cindex += 1
                                continue
                        if ind in covdict.keys():
                                covnum = int(re.sub(r'\S+:(\d+):\S+', r'\1', d[rindex][cindex]))
                                covdict[ind][locus] = covnum
                        else:
                                covnum = int(re.sub(r'\S+:(\d+):\S+', r'\1', d[rindex][cindex]))
                                covdict[ind] = {locus : covnum}
                        cindex += 1
                rindex += 1
        return covdict




# Create dictionary of individuals, sequences, and alleles
def Seqs (listfile, clipcutsite, cutsite1, cutsite2, checkbothends, CoverageCutoff, covdict):
        print "Cataloging loci..."
        if CoverageCutoff > 0:
                printtext = "Removing genotypes below coverage threshold of %s..."
                printval = (CoverageCutoff)
                print (printtext % printval)
        newfasta = open(listfile,"U") #This is old fasta if outtype == 2
        seqsdict = {} #Structure: {SampleID : {LocusID : {AlleleID : DNAsequence} } }
        carrot = ">"

        try: # Make temporary dictionary to associate SampleID and IndividualID using population file
                # This is used for coverage cutoff option
                pops = csv.DictReader(open(sys.argv[3],"U"), delimiter="\t", quotechar='"', dialect="excel-tab")
        except IOError:
                print "Error: File does not appear to exist. Check the file and the directory path."
                exit(1)
        popsdicttemp = {}
        for i in pops:
                popsdicttemp[i[pops.fieldnames[0]]] = i[pops.fieldnames[1]]
        
        for line in newfasta:
                if carrot in line:
                        indnum = re.sub(r'>CLocus_\w+_Sample_(\w+)_Locus_\w+_Allele_\w+', r'\1', line); indnum=indnum.strip()
                        locusnum = re.sub(r'>CLocus_(\w+)_Sample_\w+_Locus_\w+_Allele_\w+', r'\1', line); locusnum=locusnum.strip()
                        allelenum = re.sub(r'>CLocus_\w+_Sample_\w+_Locus_\w+_Allele_(\w+)', r'\1', line); allelenum=allelenum.strip()
                        nextline = newfasta.next(); nextline = nextline.strip()
                        skip = 0
                        # Clip off cut sites
                        if clipcutsite == 1 and checkbothends == 2:
                                print "Clipping cut sites..."
                        	if nextline[0:len(cutsite1)] == cutsite1 and nextline[(len(nextline)-len(cutsite2)):] == cutsite2:
                        		nextline = nextline[len(cutsite1):(len(nextline)-len(cutsite2))]
                        		skip = 1
                        	elif nextline[0:len(cutsite2)] == cutsite2 and nextline[(len(nextline)-len(cutsite1)):] == cutsite1 and skip == 0:
                        		nextline = nextline[len(cutsite2):(len(nextline)-len(cutsite1))]
                        elif clipcutsite == 1 and checkbothends == 1:
                                print "Clipping cut sites..."
                        	if nextline[0:len(cutsite1)] == cutsite1:
                        		nextline = nextline[len(cutsite1):]
                        		skip = 1
                        	elif nextline[0:len(cutsite2)] == cutsite2 and skip == 0:
                        		nextline = nextline[len(cutsite2):]

                        # Produce seqsdict without coverage cutoff
                        if CoverageCutoff == 0:
                                if indnum in seqsdict.keys():
                                        if locusnum in seqsdict[indnum].keys():
                                                if allelenum not in seqsdict[indnum][locusnum].keys():
                                                        seqsdict[indnum][locusnum][allelenum] = nextline
                                                else:
                                                        print "Error, you should only have one copy of each allele for each individual/locus!"
                                        else:
                                                seqsdict[indnum][locusnum] = {allelenum:nextline}
                                else:
                                        seqsdict[indnum] = {locusnum:{allelenum:nextline}}
                                        
                        # Produce seqsdict with coverage cutoff
                        if CoverageCutoff > 0:
                                if indnum in popsdicttemp.keys():
                                        if int(locusnum) in covdict[popsdicttemp[indnum]].keys():
                                                if indnum in seqsdict.keys():
                                                        if locusnum in seqsdict[indnum].keys():
                                                                if allelenum not in seqsdict[indnum][locusnum].keys() and covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                        seqsdict[indnum][locusnum][allelenum] = nextline
                                                                elif allelenum not in seqsdict[indnum][locusnum].keys() and covdict[popsdicttemp[indnum]][int(locusnum)] < CoverageCutoff:
                                                                        pass
                                                                else:
                                                                        print "Error, you should only have one copy of each allele for each individual/locus!"
                                                        else:
                                                                if covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                        seqsdict[indnum][locusnum] = {allelenum:nextline}
                                                else:
                                                        if covdict[popsdicttemp[indnum]][int(locusnum)] >= CoverageCutoff:
                                                                seqsdict[indnum] = {locusnum:{allelenum:nextline}}

	newfasta.close()
        return seqsdict




covdict = {}
if UseCoverage == 1 and CoverageCutoff != 0:
        covdict = LocusCoverage(sys.argv[4])
        
if outtype == 1:
	Job1 = KeepLoci(sys.argv[1],sys.argv[2],listfile)
	seqsdict = Seqs(listfile, clipcutsite, cutsite1, cutsite2, checkbothends, CoverageCutoff, covdict)
	
popsdict = Pops(sys.argv[3])
num_pops = len(popsdict)

if outtype == 2:
	seqsdict = Seqs(sys.argv[1], clipcutsite, cutsite1, cutsite2, checkbothends, CoverageCutoff, covdict)




# Create dictionary with number gene copies per population
def GeneCopyCount():
        print "Counting gene copies..."
        gene_copies = {} #Structure: {PopulationID : DiploidGeneCopies}
        for i, j in popsdict.iteritems():
                gene_copies[i] = 2*len(j)
        return gene_copies




# Create dictionary with locus number and number of sites
def SeqSitesCount(listfile, cutsite1, cutsite2, checkbothends):
        print "Counting locus lengths..."
        num_sites = {} #Structure: {LocusID : NumberNucleotides}
        newfasta = open(listfile,"U")
        for i in newfasta:
                locusnum = re.sub(r'>CLocus_(\w+)_Sample_\w+_Locus_\w+_Allele_\w+', r'\1', i); locusnum=locusnum.strip()
                nextline = newfasta.next(); nextline = nextline.strip()
                skip = 0
                if clipcutsite == 1 and checkbothends == 2:
                        if nextline[0:len(cutsite1)] == cutsite1 and nextline[(len(nextline)-len(cutsite2)):] == cutsite2:
                                nextline = nextline[len(cutsite1):(len(nextline)-len(cutsite2))]
                                skip = 1
                        elif nextline[0:len(cutsite2)] == cutsite2 and nextline[(len(nextline)-len(cutsite1)):] == cutsite1 and skip == 0:
                                nextline = nextline[len(cutsite2):(len(nextline)-len(cutsite1))]
                elif clipcutsite == 1 and checkbothends == 1:
                        if nextline[0:len(cutsite1)] == cutsite1:
                                nextline = nextline[len(cutsite1):]
                                skip = 1
                        elif nextline[0:len(cutsite2)] == cutsite2 and skip == 0:
                                nextline = nextline[len(cutsite2):]
                if locusnum not in num_sites.keys():
                        num_sites[locusnum] = len(nextline)
        newfasta.close()
        return num_sites




gene_copies = GeneCopyCount()
if outtype == 1:
	num_sites = SeqSitesCount(listfile, cutsite1, cutsite2, checkbothends)
elif outtype == 2:
	num_sites = SeqSitesCount(sys.argv[1], cutsite1, cutsite2, checkbothends)




# Screen for false alleles
def AlleleRemoval(seqsdict, popsdict, gene_copies, num_sites, allele_threshold, allele_pop_threshold, allele_filter):
        #Create dictionary of loci containing dictionaries of allele counts
        allelecount = {} #Build structure: {LocusID : {Sequence : Count} }
        print "Counting alleles..."
        for x in sorted(seqsdict.iterkeys()): #Cycle through all individuals
                for p in sorted(seqsdict[x].iterkeys()): # Cycle through all loci
                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through both alleles
                                if p in allelecount.keys():
                                        if str(seqsdict[x][p][a]) not in allelecount[p].keys(): #If allele not added yet, add and set count to 1
                                                allelecount[p][str(seqsdict[x][p][a])] = 1
                                        else:
                                                allelecount[p][str(seqsdict[x][p][a])] += 1
                                else:
                                        allelecount[p] = {str(seqsdict[x][p][a]):1}

        if allele_filter == 1:
                print "Removing alleles using thresholds..."
                if allele_threshold != 0: printtext = "Allele must have an overall frequency of %s..."; printval = allele_threshold; print (printtext % printval)
                if allele_pop_threshold != 0: printtext = "Allele must be present in %s of populations..."; printval = allele_pop_threshold; print (printtext % printval)
                 
                num_inds = 0
                for k in popsdict.iterkeys():
                        num_inds += int(len(popsdict[k]))

                # Locus by locus in allelecount dictionary
                for p in sorted(allelecount.iterkeys()):
                        count = 0
                        remove_flag = 0
                        for b in sorted(allelecount[p].iterkeys()):
                                count = allelecount[p][b]
                                if float(count)/(2.0*float(num_inds)) < allele_threshold and allele_threshold != 0: #Remove alleles below threshold
                                        remove_flag = 1
                                        print "Removed allele"
                                        for x in sorted(seqsdict.iterkeys()):
                                                if p in seqsdict[x].keys():
                                                        for a in sorted(seqsdict[x][p].iterkeys()):
                                                        		if p in seqsdict[x]: #This locus might have been removed for individual already
                                                                		if seqsdict[x][p][a] == b:
                                                                				del seqsdict[x][p]
                                                                				if b in allelecount[p]: del allelecount[p][b]
                                if remove_flag == 0: #Remove alleles not present in enough populations
                                        count = 0
                                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                                flag = 0
                                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                                if p in seqsdict[x].keys():
                                                                        for a in sorted(seqsdict[x][p].itervalues()):
                                                                                if a == b:
                                                                                        flag = 1
                                                if flag == 1: count += 1
                                        if float(count)/(float(gene_copies[k])) < allele_pop_threshold and allele_pop_threshold != 0:
                                                print "Removed allele"
                                                for x in sorted(seqsdict.iterkeys()):
                                                        if p in seqsdict[x].keys():
                                                                for a in sorted(seqsdict[x][p].iterkeys()):
                                                        				if p in seqsdict[x]: #This locus might have been removed for individual already
                                                                				if seqsdict[x][p][a] == b:
                                                                						del seqsdict[x][p]
                                                                						if b in allelecount[p]: del allelecount[p][b]
        return seqsdict, num_sites




def MissingData(seqsdict, popsdict, gene_copies, num_sites, locus_threshold, locus_pop_threshold, ind_threshold):
        print "Removing loci using thresholds..."
        if locus_threshold != 0: printtext = "Locus must have an overall frequency of %s..."; printval = locus_threshold; print (printtext % printval)
        if locus_pop_threshold != 0: printtext = "Locus must have a frequency of %s in each population..."; printval = locus_pop_threshold; print (printtext % printval)
        if ind_threshold != 0: printtext = "Individual must have %s of total loci..."; printval = ind_threshold; print (printtext % printval)

        if locus_pop_threshold != 0: # Remove loci from pops when below pop threshold
                locus_pop_count = {}
                for k in sorted(popsdict.iterkeys()): #Look at one population
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                if p in locus_pop_count.keys():
                                                        if k in locus_pop_count[p].keys(): locus_pop_count[p][k] += 1
                                                        else: locus_pop_count[p][k] = 1
                                                else:
                                                        locus_pop_count[p] = {k : 1}
                for p in locus_pop_count.iterkeys():
                        for k in locus_pop_count[p].iterkeys():
                                if float(locus_pop_count[p][k])/(float(gene_copies[k])/2.0) < locus_pop_threshold:
                                        print "Removed locus..."
                                        for x in sorted(seqsdict.iterkeys()):
                                                if p in seqsdict[x] and x in popsdict[k].keys():
                                                        del seqsdict[x][p]
                                                        
        if locus_threshold != 0: # Remove loci below threshold                                                       
                num_inds = 0
                for k in popsdict.iterkeys():
                        num_inds += int(len(popsdict[k]))

                locus_count = {} #Structure: { LocusID : count }
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                if p in locus_count.keys(): locus_count[p] += 1
                                else: locus_count[p] = 1
                for p in locus_count.iterkeys():
                        if float(locus_count[p])/(float(num_inds)) < locus_threshold:
                                print "Removed locus..."
                                for x in sorted(seqsdict.iterkeys()):
                                        if p in seqsdict[x]:
                                                del seqsdict[x][p]
                                                if p in num_sites:
                                                        del num_sites[p]
                                                        
        if ind_threshold != 0: # Remove individuals below threshold
                locus_list = []
                for x in sorted(seqsdict.iterkeys()):
                        for p in sorted(seqsdict[x].iterkeys()):
                                if p not in locus_list: locus_list.append(p)
                new_num_loci = len(locus_list)
                if new_num_loci == 0:
                        print "All loci removed. Check data."
                        exit(1)

                locus_ind_count = {}
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        locus_ind_count[x] = len(seqsdict[x])

                for x in sorted(seqsdict.iterkeys()):
                        if float(locus_ind_count[x])/float(new_num_loci) < ind_threshold:
                                print "Removed individual..."
                                del seqsdict[x]

        return seqsdict, num_sites




# Screen for false alleles
def LocusRemoval(seqsdict, popsdict, gene_copies, num_sites, monomorphic_filter):
        #Create dictionary of loci containing dictionaries of allele counts
        allelecount = {} #Build structure: {LocusID : {Sequence : Count} }
        print "Counting alleles..."
        for x in sorted(seqsdict.iterkeys()): #Cycle through all individuals
                for p in sorted(seqsdict[x].iterkeys()): # Cycle through all loci
                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through both alleles
                                if p in allelecount.keys():
                                        if str(seqsdict[x][p][a]) not in allelecount[p].keys(): #If allele not added yet, add and set count to 1
                                                allelecount[p][str(seqsdict[x][p][a])] = 1
                                        else:
                                                allelecount[p][str(seqsdict[x][p][a])] += 1
                                else:
                                        allelecount[p] = {str(seqsdict[x][p][a]):1}

        if monomorphic_filter == 1:
                print "Removing monomorphic loci..."
                for p in sorted(allelecount.iterkeys()): #Remove monomorphic loci
                        if len(allelecount[p]) == 1:
                                print "Removed locus..."
                                for x in sorted(seqsdict.iterkeys()):
                                        if p in seqsdict[x]:
                                                del seqsdict[x][p]
                                                if p in num_sites:
                                                        del num_sites[p]

        return seqsdict, num_sites




if allele_filter == 1:
        seqsdict, num_sites = AlleleRemoval(seqsdict, popsdict, gene_copies, num_sites, allele_threshold, allele_pop_threshold, allele_filter)

if missing_data_filter == 1:
        seqsdict, num_sites = MissingData(seqsdict, popsdict, gene_copies, num_sites, locus_threshold, locus_pop_threshold, ind_threshold)

if monomorphic_filter == 1:
        seqsdict, num_sites = LocusRemoval(seqsdict, popsdict, gene_copies, num_sites, monomorphic_filter)



        
# Test functions
# print seqsdict.values(); print popsdict.values(); print num_sites; print gene_copies; print num_pops




# Output Migrate file
def Fasta2Migrate(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting migrate-n file..."
        try:
                OrderedLoci = []
                par=""
                fout=open(outfile,"w")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        
                fout.write(str(num_pops) + '\t' + str(len(OrderedLoci)) + "\n")
                for key in sorted(num_sites.iterkeys()):
                        fout.write("%s\t" % num_sites[key])
                fout.write('\n')
                
                

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write(str(gene_copies[k]) + '\t' + 'Pop_' + k + '\n')
                        for i in OrderedLoci: #Look at one locus
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                if int(len(str(popsdict[k][x])))>9: print "Error, Ind ID > 9 characters"; exit(1)
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind+'a').ljust(10)); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("?"); z+=1
                                                        fout.write('\n')
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind+'b').ljust(10)); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("?"); z+=1
                                                        fout.write('\n')

                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                ind = str(popsdict[k][x])
                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                        if int(a) == 0: par = 'a'
                                                                                        elif int(a) == 1: par = 'b'
                                                                                        fout.write((ind+par).ljust(10)+str(seqsdict[x][p][a])+'\n')
                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write((ind+'a').ljust(10)+str(seqsdict[x][p][a])+'\n'+(ind+'b').ljust(10)+str(seqsdict[x][p][a])+'\n')
                                                                                else:
                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                        if int(a) == 0: par = 'a'
                                                                                        elif int(a) == 1: par = 'b'
                                                                                        else: par = 'X'; break
                                                                                        fout.write((ind+par).ljust(10)+str(seqsdict[x][p][a])+'\n')
                fout.close()

        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




# Output Arlequin file
def Fasta2Arlequin(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        try:
                OrderedLoci = []
                title = raw_input("Title of Project: ")
                print "Outputting Arelequin file..."
                fout=open(outfile,"w")

                fout.write("[Profile]\n\n\t\"" + title + "\"\n\n\t\tNbSamples=" + str(num_pops))
                fout.write("\n\t\tGenotypicData=1\n\t\tGameticPhase=0\n\t\tDataType=DNA\n\t\t")
                fout.write("LocusSeparator=TAB\n\t\tMissingData=\"?\"\n\n\n[Data]\n\n\t[[Samples]]\n\n")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write("\t\tSampleName= \"Pop_" + str(k) + "\"\n\t\tSampleSize=" + str(int(gene_copies[k])/2) + "\n\t\tSampleData={\n")
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        count = 0
                                        while count < 2:
                                                ind = str(popsdict[k][x])
                                                if count == 0: fout.write(ind+'\t1\t')
                                                if count == 1: fout.write('\t\t')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                                z=0
                                                                while z < int(num_sites[i]):
                                                                        fout.write("?"); z+=1
                                                                fout.write('\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))>2:
                                                                                                print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                count += 1
                                                fout.write('\n')
                        fout.write("}\n")
                fout.close()
        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




# Output DIYABC file
def Fasta2DIYABC(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting DIYabc file..."
        try:
                OrderedLoci = []
                title = raw_input("Title of Project: ")
                fout=open(outfile,"w")

                fout.write(title + " <NM=NF>\n")

                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        fout.write(key + "\t<A>\n")

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        fout.write("Pop\n")
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        ind = str(popsdict[k][x])
                                        fout.write(ind+'\t,\t')
                                        for i in OrderedLoci: #Look at one locus
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write empty brackets
                                                        fout.write("<[][]>\t")

                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                        if (int(a)==0):
                                                                                                fout.write("<[" + str(seqsdict[x][p][a]) + "]")
                                                                                        elif (int(a)==1):
                                                                                                fout.write("[" + str(seqsdict[x][p][a]) + "]>\t")

                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write("<[" + str(seqsdict[x][p][a]) + "][" + str(seqsdict[x][p][a]) + "]>\t")
                                                                                elif int(len(seqsdict[x][p]))>2:
                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                        if (int(a)==0):
                                                                                                fout.write("<[" + str(seqsdict[x][p][a]) + "]")
                                                                                        elif (int(a)==1):
                                                                                                fout.write("[" + str(seqsdict[x][p][a]) + "]>\t")
                                        fout.write('\n')
                fout.close()
        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




# Output LFMM file
def Fasta2LFMM(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting LFMM (ped) file..."

        try:

                whitelist = []

                with open(sys.argv[2],"U") as f:
                        for line in f:
                                whitelist.append(str(line.strip()))

                with open(sys.argv[4],"U") as f:
                        vcf = csv.reader(f, delimiter="\t")
                        d = list(vcf)
                snpname = str(outname+".snp")
                fout=open(snpname,"w")

                print "Building dictionary of SNPs..."
                snpdict = {} # Structure: {IndividualID : {SnpID : [Allele1, Allele2] } }
                rindex = 0
                lastlocus = 0
                for row in d:
                        if rindex < 9:
                                rindex += 1
                                continue
                        cindex = 0
                        locus = int(d[rindex][2])
                        if str(locus) not in whitelist:
                                rindex += 1
                                continue
                        if locus == lastlocus:allele += 1
                        else: allele = 1
                        lastlocus = locus
                        fout.write(str(d[rindex][0])+'\t'+str(d[rindex][2])+'_snp'+str(allele)+'\t0\t'+str(d[rindex][1])+'\t'+str(d[rindex][3])+'\t'+str(d[rindex][4])+'\n')
                        for column in row:
                                ind = d[8][cindex]
                                if cindex < 9:
                                        cindex += 1
                                        continue
                                if ind in snpdict.keys():
                                        a1 = str(re.sub(r'(\S)/\S:\d+:\S+', r'\1', d[rindex][cindex]))
                                        a2 = str(re.sub(r'\S/(\S):\d+:\S+', r'\1', d[rindex][cindex]))
                                        if a1 == '.': allele1 = 0
                                        else:
                                                if int(a1) == 1:
                                                        allele1 = str(d[rindex][4])
                                                elif int(a1) == 0:
                                                        allele1 = str(d[rindex][3])
                                        if a2 == '.': allele2 = 0
                                        else:
                                                if int(a2) == 1:
                                                        allele2 = str(d[rindex][4])
                                                elif int(a2) == 0:
                                                        allele2 = str(d[rindex][3])
                                        snpdict[ind][str(str(locus)+"_snp"+str(allele))] = [allele1, allele2]
                                else:
                                        a1 = str(re.sub(r'(\S)/\S:\d+:\S+', r'\1', d[rindex][cindex]))
                                        a2 = str(re.sub(r'\S/(\S):\d+:\S+', r'\1', d[rindex][cindex]))
                                        if a1 == '.': allele1 = 0
                                        else: allele1 = int(a1)
                                        if a2 == '.': allele2 = 0
                                        else: allele2 = int(a2)
                                        snpdict[ind] = {str(str(locus)+"_snp"+str(allele)) : [allele1, allele2]}
                                cindex += 1
                        rindex += 1
                fout.close()

                print "Writing file..."
                fout=open(outfile,"w")

                for k in sorted(popsdict.iterkeys()): #Look at one population
                        for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                if x in popsdict[k].keys(): #If that individual's in the population
                                        ind = str(popsdict[k][x])
                                        fout.write(str(ind)+'\t'+str(k)+'\t0\t0\t0\t0')
                                        for i in sorted(snpdict[ind].iterkeys()): #Look at one locus
                                                for a in sorted(snpdict[ind][i]): #Cycle through 1 or 2 alleles
                                                        fout.write('\t'+str(a))
                                        fout.write("\n")
                fout.close()


        except IOError:
                print "Error: Problems outputting file. Make sure VCF file is included. Check the directory path."
                return 0

        return 1




# Output Phylip file
def Fasta2Phylip(num_pops, popsdict, seqsdict, gene_copies, num_sites, haplo, breakpoints, haplotypes):
        print "Outputting Phylip file..."

        try:
                

                if haplo == 1:

                        nind = sum(len(v) for v in popsdict.itervalues())
                        
                        print "   Cataloging unique sequences..."
                        globalseqsdict = {} #Structure: {LocusID : {AlleleID : DNAsequence} }
                        for ind in seqsdict.keys():
                                for locus in seqsdict[ind].keys():
                                        for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                                if locus not in globalseqsdict.keys():
                                                        globalseqsdict[locus] = { 0 : seqsdict[ind][locus][allele_id] }
                                                        
                                                else:
                                                        if allele_seq not in globalseqsdict[locus].values():
                                                                new_allele = max(globalseqsdict[locus].keys()) + 1
                                                                globalseqsdict[locus][new_allele] = allele_seq

                        print "   Building dictionary of haplotypes..."
                        haplodict = {} #Structure: {LocusID : {AlleleID : Haplotype} }
                        for locus in globalseqsdict.keys():
                                SNP_positions = []
                                for letter, letter_str in enumerate(globalseqsdict[locus][0]):
                                        default = globalseqsdict[locus][0][letter]
                                        for allele in globalseqsdict[locus].keys():
                                                if globalseqsdict[locus][allele][letter] != default:
                                                        SNP_positions.append(letter)

                                for allele in globalseqsdict[locus].keys():
                                        hap = ""
                                        for i,j in enumerate(SNP_positions):
                                                hap += globalseqsdict[locus][allele][SNP_positions[i]]

                                        if locus not in haplodict.keys():
                                                haplodict[locus] = { 0 : hap }
                                        else:
                                                haplodict[locus][allele] = hap

                        print "   Cataloging haplotypes..."
                        newseqsdict = seqsdict #Structure: {SampleID : {LocusID : {AlleleID : Haplotype} } }
                        for ind in seqsdict.keys():
                                for locus in seqsdict[ind].keys():
                                        for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                                newseqsdict[ind][locus][allele_id] = haplodict[locus][globalseqsdict[locus].keys()[globalseqsdict[locus].values().index(allele_seq)]]

                        print "   Counting haplotype SNP sites..."
                        new_num_sites = {} #Structure: {LocusID : NumberNucleotides}
                        for locus in haplodict.keys():
                                new_num_sites[locus] = len(haplodict[locus][0])

                if haplo == 2:
                        newseqsdict = seqsdict
                        new_num_sites = num_sites
                        nind = sum(len(v) for v in popsdict.itervalues()) * 2

                OrderedLoci = []
                fin=open(str(sys.argv[2]),"U")
                for line in fin:
                        OrderedLoci.append(str(line.strip()))
                fin.close()
                fout=open(outfile,"w")

                nbp = sum(new_num_sites.itervalues())
                fout.write(str(nind) + '\t' + str(nbp) + '\n')
                        
        	print "   Writing file..."
        	if haplotypes == 2:
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(newseqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write((ind+' ').ljust(10)); z=0
                                                for i in OrderedLoci: #Look at one locus
                                                        if int(len(str(popsdict[k][x])))>9: print "Error, Ind ID > 9 characters"; exit(1)
                                                        if i not in newseqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                                z=0
                                                                while z < int(new_num_sites[i]):
                                                                        fout.write("N"); z+=1
                                                                        if breakpoints == 1:
                                                                                fout.write("!")

                                                        else:
                                                                for p in sorted(newseqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                if int(len(newseqsdict[x][p]))==2:#If heterozygote
                                                                                        seq = ""
                                                                                        for n, m in enumerate(str(newseqsdict[x][p]['0'])):
                                                                                                if str(newseqsdict[x][p]['0'][n]) == str(newseqsdict[x][p]['1'][n]):
                                                                                                        seq += str(newseqsdict[x][p]['0'][n])
                                                                                                elif str(newseqsdict[x][p]['0'][n]) not in ['A','T','C','G']:
                                                                                                        ind = str(popsdict[k][x])
                                                                                                        fout.write((ind+' ').ljust(10)); z=0
                                                                                                        while z < int(new_num_sites[i]):
                                                                                                                fout.write("N"); z+=1
                                                                                                else:
                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "C":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "S"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                        seq += "M"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "Y"
                                                                                                                
                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "G":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "S"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                        seq += "R"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "K"

                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "A":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "M"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "R"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "W"

                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "T":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "Y"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "K"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                        seq += "W"
                                                                                                
                                                                                        fout.write(seq)
                                                                                        if breakpoints == 1:
                                                                                                fout.write("!")
                                                                                elif int(len(newseqsdict[x][p]))==1:#If homozygote
                                                                                        fout.write(str(newseqsdict[x][p]['0']))
                                                                                        if breakpoints == 1:
                                                                                                fout.write("!")
                                                                                else:
                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                        seq = ""
                                                                                        for n, m in enumerate(str(newseqsdict[x][p]['0'])):
                                                                                                if str(newseqsdict[x][p]['0'][n]) == str(newseqsdict[x][p]['1'][n]):
                                                                                                        seq += str(newseqsdict[x][p]['0'][n])
                                                                                                elif str(newseqsdict[x][p]['0'][n]) not in ['A','T','C','G']:
                                                                                                        ind = str(popsdict[k][x])
                                                                                                        fout.write((ind+' ').ljust(10)); z=0
                                                                                                        while z < int(new_num_sites[i]):
                                                                                                                fout.write("N"); z+=1
                                                                                                else:
                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "C":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "S"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                        seq += "M"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "Y"
                                                                                                                
                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "G":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "S"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                        seq += "R"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "K"

                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "A":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "M"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "R"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "T":
                                                                                                                        seq += "W"

                                                                                                        if str(newseqsdict[x][p]['0'][n]) == "T":
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "C":
                                                                                                                        seq += "Y"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "G":
                                                                                                                        seq += "K"
                                                                                                                if str(newseqsdict[x][p]['1'][n]) == "A":
                                                                                                                                seq += "W"                                                                                                                                                          
                                                                                        fout.write(seq)
                                                                                        if breakpoints == 1:
                                                                                                fout.write("!")
                                                fout.write('\n')
                if haplotypes == 1:
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(newseqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                gene_copies = ['0','1']
                                                for gene_copy in gene_copies:
                                                        if gene_copy == '0':
                                                                ind = str(popsdict[k][x]) + 'a'
                                                                fout.write((ind+' ').ljust(10)); z=0
                                                        if gene_copy == '1':
                                                                ind = str(popsdict[k][x]) + 'b'
                                                                fout.write((ind+' ').ljust(10)); z=0
                                                        for i in OrderedLoci: #Look at one locus
                                                                if int(len(str(popsdict[k][x])))>9: print "Error, Ind ID > 9 characters"; exit(1)
                                                                if i not in newseqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                                        z=0
                                                                        while z < int(new_num_sites[i]):
                                                                                fout.write("N"); z+=1
                                                                        if breakpoints == 1:
                                                                                fout.write("!")

                                                                else:
                                                                        for p in sorted(newseqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        if int(len(newseqsdict[x][p]))>=2:#If heterozygote
                                                                                                if int(len(newseqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                fout.write(str(newseqsdict[x][p][gene_copy]))
                                                                                                if breakpoints == 1:
                                                                                                        fout.write("!")
                                                                                        elif int(len(newseqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(newseqsdict[x][p]['0']))
                                                                                                if breakpoints == 1:
                                                                                                        fout.write("!")
                                                        fout.write('\n')
                fout.close()


        except IOError:
                print "Error: Problems outputting file. Make sure VCF file is included. Check the directory path."
                return 0

        return 1




# Output G-Phocs file
def Fasta2GPhocs(num_pops, popsdict, seqsdict, gene_copies, num_sites):
        print "Outputting G-Phocs file..."
        try:
                OrderedLoci = []
                fout=open(outfile,"w")
                
                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        
                numinds = sum(len(q) for q in popsdict.itervalues())
                
                fout.write(str(len(OrderedLoci)) + '\n\n')
				
                for i in OrderedLoci: #Look at one locus
                        fout.write(str(i) + '\t' + str(numinds) + '\t' + str(num_sites[i]) + '\n')
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write Ns
                                                        ind = str(popsdict[k][x])
                                                        fout.write((ind) + '\t'); z=0
                                                        while z < int(num_sites[i]):
                                                                fout.write("N"); z+=1


                                                else:
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                if p == i:      #If this is the right locus
                                                                        ind = str(popsdict[k][x])
                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                seq = ""
                                                                                for n, m in enumerate(str(seqsdict[x][p]['0'])):
                                                                                        if str(seqsdict[x][p]['0'][n]) == str(seqsdict[x][p]['1'][n]):
                                                                                                seq += str(seqsdict[x][p]['0'][n])
                                                                                        elif str(seqsdict[x][p]['0'][n]) not in ['A','T','C','G']:
                                                                                                ind = str(popsdict[k][x])
                                                                                        	fout.write((ind) + '\t'); z=0
                                                                                                while z < int(num_sites[i]):
                                                                                                        fout.write("N"); z+=1
                                                        					fout.write('\n')
                                                                                        else:
                                                                                                if str(seqsdict[x][p]['0'][n]) == "C":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "S"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                                seq += "M"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "Y"
                                                                                                        
                                                                                                if str(seqsdict[x][p]['0'][n]) == "G":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "S"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                                seq += "R"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "K"

                                                                                                if str(seqsdict[x][p]['0'][n]) == "A":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "M"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "R"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "W"

                                                                                                if str(seqsdict[x][p]['0'][n]) == "T":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "Y"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "K"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                                seq += "W"
                                                                                        
                                                                                fout.write(ind + '\t' + seq + '\n')
                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                fout.write(ind + '\t' + str(seqsdict[x][p]['0'])+'\n')
                                                                        else:
                                                                                print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                seq = ""
                                                                                for n, m in enumerate(str(seqsdict[x][p]['0'])):
                                                                                        if str(seqsdict[x][p]['0'][n]) == str(seqsdict[x][p]['1'][n]):
                                                                                                seq += str(seqsdict[x][p]['0'][n])
                                                                                        elif str(seqsdict[x][p]['0'][n]) not in ['A','T','C','G']:
                                                                                                ind = str(popsdict[k][x])
                                                                                        	fout.write((ind) + '\t'); z=0
                                                        					while z < int(num_sites[i]):
                                                                					fout.write("N"); z+=1
                                                        					fout.write('\n')
                                                                                        else:
                                                                                                if str(seqsdict[x][p]['0'][n]) == "C":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "S"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                                seq += "M"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "Y"
                                                                                                        
                                                                                                if str(seqsdict[x][p]['0'][n]) == "G":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "S"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                                seq += "R"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "K"

                                                                                                if str(seqsdict[x][p]['0'][n]) == "A":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "M"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "R"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "T":
                                                                                                                seq += "W"

                                                                                                if str(seqsdict[x][p]['0'][n]) == "T":
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "C":
                                                                                                                seq += "Y"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "G":
                                                                                                                seq += "K"
                                                                                                        if str(seqsdict[x][p]['1'][n]) == "A":
                                                                                                        		seq += "W"
                                                                                                                                                   
                                                                                fout.write(ind + '\t' + seq + '\n')
                        fout.write('\n')
                fout.close()

        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




# Output Treemix file
def Fasta2Treemix(num_pops, popsdict, seqsdict, gene_copies, num_sites, one_snp):
        print "Outputting Treemix file..."

        try:
        	print "   Finding unique sequences..."
                # Build dictionary of unique sequences
        	globalseqsdict = {} #Structure: {LocusID : {AlleleID : DNAsequence} }
        	for ind in seqsdict.keys():
        		for locus in seqsdict[ind].keys():
        			for allele_id, allele_seq in seqsdict[ind][locus].iteritems():
                                        if locus not in globalseqsdict.keys():
                                                globalseqsdict[locus] = { 0 : seqsdict[ind][locus][allele_id] }
                                                
                                        else:
                                                if allele_seq not in globalseqsdict[locus].values():
                                                        new_allele = max(globalseqsdict[locus].keys()) + 1
                                                        globalseqsdict[locus][new_allele] = allele_seq

        	print "   Finding SNP positions..."
                # Find sequence positions that are SNPs
                SNPpositions = {} #Structure {LocusID : [Positions] }
        	for locus in globalseqsdict.keys():
                        SNP_positions = []
                        for letter, letter_str in enumerate(globalseqsdict[locus][0]):
                                default = globalseqsdict[locus][0][letter]
                                for allele in globalseqsdict[locus].keys():
                                        if globalseqsdict[locus][allele][letter] != default and letter not in SNP_positions:
                                                SNP_positions.append(letter)
                        SNPpositions[locus] = SNP_positions

        	print "   Cataloging SNPs..."
                # Build dictionary of SNPs
                SNPdict = {} #Structure: {SampleID : {LocusID : {PositionID : {AlleleID : SNP} } } }
                for ind in seqsdict.keys():                               
                        for locus in seqsdict[ind].keys():
                                for position in SNPpositions[locus]:
                                        for allele_id in seqsdict[ind][locus]:
                                                if ind not in SNPdict.iterkeys():
                                                        SNPdict[ind] = {locus : {position: {allele_id : seqsdict[ind][locus][allele_id][position]} } }
                                                else:
                                                        if locus not in SNPdict[ind].iterkeys():
                                                                SNPdict[ind][locus] = {position: {allele_id : seqsdict[ind][locus][allele_id][position]} }
                                                        else:
                                                                if position not in SNPdict[ind][locus].iterkeys():
                                                                        SNPdict[ind][locus][position] = {allele_id : seqsdict[ind][locus][allele_id][position]}
                                                                else:
                                                                        if allele_id not in SNPdict[ind][locus][position].iterkeys():
                                                                                SNPdict[ind][locus][position][allele_id] = seqsdict[ind][locus][allele_id][position]

        	print "   Doing more SNP cataloging..."
                # Dictionary with both SNP letters present at each SNP locus
                SNPletters = {} #Structure: {LocusID : {PositionID : [Letters] } }
                
        	print "   Flagging SNPs that are not biallelic..."
                # Flag SNPs with >2 alleles
                globalSNPalleleCounter = {} #Structure: {LocusID : {PositionID : NumAlleles } }
                for ind in SNPdict.keys():                               
                        for locus in SNPdict[ind].keys():
                                for position in SNPdict[ind][locus].keys():
                                        for allele_id, allele_snp in SNPdict[ind][locus][position].iteritems():
                                                if locus not in globalSNPalleleCounter.iterkeys():
                                                        globalSNPalleleCounter[locus] = {position : 1 }
                                                        SNPletters[locus] = {position : [allele_snp]}
                                                else:
                                                        if position not in globalSNPalleleCounter[locus].iterkeys():
                                                                globalSNPalleleCounter[locus][position] = 1
                                                                SNPletters[locus][position] = [allele_snp]
                                                        else:
                                                                if allele_snp not in SNPletters[locus][position]:
                                                                        SNPletters[locus][position].append(allele_snp)
                                                                        globalSNPalleleCounter[locus][position] += 1


        	print "   Counting population SNP alleles..."
                # Dictionary with count for each SNP allele
                SNPpopCount = {} #Structure: {PopID : {LocusID : {PositionID : {SNP : Count} } } }
                for pop in popsdict.iterkeys(): #Look at one population
                        for ind in SNPdict.iterkeys(): #look at one individual
                                if ind in popsdict[pop].iterkeys():
                                        for locus in SNPdict[ind].iterkeys(): #Look at one locus
                                                for position in SNPdict[ind][locus].iterkeys(): #Look at one SNP position in locus
                                                        for allele_id, allele_snp in SNPdict[ind][locus][position].iteritems(): #Look at one SNP allele
                                                                if pop not in SNPpopCount.iterkeys():
                                                                        SNPpopCount[pop] = {locus : {position: {allele_snp : 1} } }
                                                                else:
                                                                        if locus not in SNPpopCount[pop].iterkeys():
                                                                                SNPpopCount[pop][locus] = {position: {allele_snp : 1} }
                                                                        else:
                                                                                if position not in SNPpopCount[pop][locus].iterkeys():
                                                                                        SNPpopCount[pop][locus][position] = {allele_snp : 1}
                                                                                else:
                                                                                        if allele_snp not in SNPpopCount[pop][locus][position].iterkeys():
                                                                                                SNPpopCount[pop][locus][position][allele_snp] = 1
                                                                                        else:
                                                                                                SNPpopCount[pop][locus][position][allele_snp] += 1
        
                OrderedLoci = []
                fin=open(str(sys.argv[2]),"U")
                for line in fin:
                        for pop in SNPpopCount.iterkeys():
                                for locus in SNPpopCount[pop].iterkeys():
                                        if pop in SNPpopCount.iterkeys():
                                                if locus in SNPpopCount[pop].iterkeys():
                                                        cleanedline = str(line.strip())
                                                        if cleanedline not in OrderedLoci:
                                                                OrderedLoci.append(cleanedline)

                fin.close()
                fout=open(outfile,"w")

                pops = []
                for pop in sorted(SNPpopCount.iterkeys()): #Look at one population
                        fout.write(str(pop)+' ') #Print each population as column header
                        pops.append(pop)
                fout.write('\n')

                numsnps = 0
                for locus in globalSNPalleleCounter:
                        for position in globalSNPalleleCounter[locus]:
                                if globalSNPalleleCounter[locus][position] <= 2:
                                        numsnps += 1

                #Print out counts of allele by population
                print "   Writing file..."
                locus = 0
                while locus < len(OrderedLoci):
                        locus_index = OrderedLoci[locus]
                        position = 0
                        if one_snp == 1:
                                num_SNPs_per_locus = 1
                        elif one_snp == 2: 
                                num_SNPs_per_locus = len(SNPpositions[locus_index])
                        while position < num_SNPs_per_locus:
                                position_index = SNPpositions[locus_index][position]
                                if globalSNPalleleCounter[locus_index][position_index] <= 2:
                                        for pop in pops:
                                                if locus_index in SNPpopCount[pop]:
                                                        for SNP in sorted(SNPpopCount[pop][locus_index][position_index].iterkeys()):
                                                                num_alleles = len(SNPpopCount[pop][locus_index][position_index])
                                                                this_SNP_first = 2
                                                                other = 2
                                                                if SNPletters[locus_index][position_index][0] == SNP and SNPletters[locus_index][position_index][0] < SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 1
                                                                        other = 1
                                                                if SNPletters[locus_index][position_index][1] == SNP and SNPletters[locus_index][position_index][0] < SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 0
                                                                        other = 0
                                                                if SNPletters[locus_index][position_index][0] == SNP and SNPletters[locus_index][position_index][0] > SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 0
                                                                        other = 1
                                                                if SNPletters[locus_index][position_index][1] == SNP and SNPletters[locus_index][position_index][0] > SNPletters[locus_index][position_index][1]:
                                                                        this_SNP_first = 1
                                                                        other = 0
                                                                if num_alleles == 1:#Only one SNP allele present in population                                                                       
                                                                        if this_SNP_first == 1:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNP]) + ",0")
                                                                        if this_SNP_first == 0:
                                                                                fout.write("0," + str(SNPpopCount[pop][locus_index][position_index][SNP]))
                                                                        fout.write(' ')
                                                                if num_alleles == 2:#Both SNP alleles in population
                                                                        if this_SNP_first == 1:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNP]) + "," + str(SNPpopCount[pop][locus_index][position_index][SNPletters[locus_index][position_index][other]]))
                                                                        if this_SNP_first == 0:
                                                                                fout.write(str(SNPpopCount[pop][locus_index][position_index][SNPletters[locus_index][position_index][other]]) + "," + str(SNPpopCount[pop][locus_index][position_index][SNP]))
                                                                        fout.write(' ')
                                                                        break
                                                else:
                                                        fout.write("0,0 ")
                                        fout.write("\n")
                                position += 1
                        locus += 1                                      


                fout.close()




        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




# Output file of haplotypes
def Fasta2Haplotype(num_pops, popsdict, seqsdict, gene_copies, num_sites, HaploChoice):
        try:
                OrderedLoci = []
                LociOrdered = {}
                OrderedPops = []
                PopsOrdered = {}

                fout=open(outfile,"w")

                count = 1
                for key in sorted(num_sites.iterkeys()):
                        OrderedLoci.append(key)
                        LociOrdered[key] = count
                        count += 1

                count = 1
                for key in sorted(popsdict.iterkeys()):
                        OrderedPops.append(key)
                        PopsOrdered[key] = count
                        count += 1

        	print "   Cataloging unique sequences..."
                #Create dictionary of loci containing dictionaries of unique sequences/integers
                newseqdict = {} #Build structure: {LocusID : {Sequence : UniqueInteger} }
                for x in sorted(seqsdict.iterkeys()):
                	for p in sorted(seqsdict[x].iterkeys()):
                		for a in sorted(seqsdict[x][p].iterkeys()):
                			if p in newseqdict.keys():
                				if str(seqsdict[x][p][a]) not in newseqdict[p].keys() and "?" not in seqsdict[x][p][a]:
                					locusints = []
                					for i in sorted(newseqdict[p].itervalues()):
                						locusints.append(i)
                						newseqdict[p][str(seqsdict[x][p][a])] = max(locusints) + 1
                                                elif "?" in seqsdict[x][p][a]:
                                                        newseqdict[p][str(seqsdict[x][p][a])] = 0

                			else:
                                                if "?" not in seqsdict[x][p][a]:
                                                        newseqdict[p] = {str(seqsdict[x][p][a]):1}
                                                else:
                                                        newseqdict[p] = {str(seqsdict[x][p][a]):0}

        	print "   Creating dictionary of unique haplotype integers..."
                #Convert sequences into unique integers by locus
                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through all loci
                		for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                			seqsdict[x][p][a] = newseqdict[str(p)][str(seqsdict[x][p][a])]

                if HaploChoice == 1: # Structure format
                        print "Outputting Structure file..."
                        fout.write('\t')
                        for i in OrderedLoci:
                                fout.write('\t'+str(i))
                        fout.write('\n')
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                count = 0
                                                while count < 2:
                                                        ind = str(popsdict[k][x])
                                                        fout.write(ind+'\t'+k+'\t')
                                                        for i in OrderedLoci: #Look at one locus
                                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                        fout.write('0\t')

                                                                else:
                                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                        count += 1
                                                        fout.write("\n")
                        fout.close()

                if HaploChoice == 2: #Genepop format
                        FourOrSix = int(raw_input("Genepop in four [1] or six [2] digit format? "))
                        while FourOrSix not in range (1,3):
                                FourOrSix = int(raw_input("Not a valid option. Genepop in four [1] or six [2] digit format? "))

                        title = raw_input("Title of Project: ")
                        print "Outputting Genepop file..."
                        fout.write(title+'\n')
                        for i in OrderedLoci:
                                fout.write(str(i)+'\n')

                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                fout.write('Pop\n')
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write(ind+' ,  ')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                if FourOrSix == 1:
                                                                        fout.write('0000\t')
                                                                if FourOrSix == 2:
                                                                        fout.write('000000\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        count = 0
                                                                        while count < 2:
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2))
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3))
                                                                                                        if (count==1 and int(a)==1):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2)+'\t')
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3)+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        if FourOrSix == 1:
                                                                                                                fout.write(str(seqsdict[x][p][a]).zfill(2)+str(seqsdict[x][p][a]).zfill(2)+'\t')
                                                                                                        if FourOrSix == 2:
                                                                                                                fout.write(str(seqsdict[x][p][a]).zfill(3)+str(seqsdict[x][p][a]).zfill(3)+'\t')
                                                                                                        count += 1
                                                                                                elif int(len(seqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                        if (count==0 and int(a)==0):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2))
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3))
                                                                                                        if (count==1 and int(a)==1):
                                                                                                                if FourOrSix == 1:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(2)+'\t')
                                                                                                                if FourOrSix == 2:
                                                                                                                        fout.write(str(seqsdict[x][p][a]).zfill(3)+'\t')
                                                                                count += 1
                                                fout.write("\n")
                        fout.close()

                        fout=open(outfile_pops,"w")
                        for i in sorted(PopsOrdered.iterkeys()):
                                fout.write(str(PopsOrdered[i])+'\t'+str(i)+'\n')
                        fout.close()

                if HaploChoice == 3: #Allele frequency by locus X population
                        print "Outputting allele frequency X population matrix..."
                        fout.write('\t')
                        for p in sorted(newseqdict.iterkeys()):
                                for s in sorted(newseqdict[p].itervalues()):
                                        fout.write(str(p)+'_'+str(s)+'\t') #Print each locus/allele combo as column header
                        fout.write('\n')

                        popfreq = {} #Build Structure: {PopulationID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "   Creating dictionary of population allele frequencies..."
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                if k in sorted(popfreq.keys()):
                                        for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                else:
                                        popfreq[k] = {}
                                        for p in sorted(newseqdict.iterkeys()):
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass

                        print "   Tabulating population allele frequencies..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                count = 0
                                                                while count < 2:
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))>=2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        popfreq[k][p][seqsdict[x][p][a]] += 1
                                                                                                if int(len(seqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."

                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                popfreq[k][p][seqsdict[x][p][a]] += 2
                                                                                                count += 1
                                                                        count += 1

                        print "   Writing file..."
                        #Print out frequencies of allele by population
                        for k in sorted(popfreq.iterkeys()):
                                fout.write(str(k)+'\t')
                                for p in sorted(popfreq[k].iterkeys()):
                                        total = 0
                                        for n in sorted(popfreq[k][p].itervalues()):
                                                total += n
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                if total > 0:
                                                        x = str(Decimal(str(popfreq[k][p][n])).quantize(Decimal('0.00001'))/Decimal(str(total)).quantize(Decimal('0.00001')))
                                                        fout.write(str(Decimal(str(x)).quantize(Decimal('0.00001')))+'\t')
                                                else:
                                                        fout.write(str(Decimal(str(popfreq[k][p][n])).quantize(Decimal('0.00001')))+'\t')
                                fout.write('\n')

                        fout.close()

                if HaploChoice == 4: #Sambada format
                        print "Outputting SamBada file..."
                        fout.write('\t')
                        for p in sorted(newseqdict.iterkeys()):
                                for s in sorted(newseqdict[p].itervalues()):
                                        fout.write(str(p)+'_'+str(s)+'\t') #Print each locus/allele combo as column header
                        fout.write('\n')

                        allelecount = {} #Build Structure: {IndividualID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "   Creating dictionary of allele counts..."
                        for i in sorted(popsdict.iterkeys()): #Look at one population
                                for k in sorted(popsdict[i].itervalues()): #Look at one individual
                                        if k in sorted(allelecount.keys()):
                                                for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                        if p in sorted(allelecount[k].iterkeys()):
                                                                for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                                        else:
                                                                allelecount[k][p] = {}
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                        else:
                                                allelecount[k] = {}
                                                for p in sorted(newseqdict.iterkeys()):
                                                        if p in sorted(allelecount[k].iterkeys()):
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass
                                                        else:
                                                                allelecount[k][p] = {}
                                                                for n in sorted(newseqdict[p].itervalues()):
                                                                        if n not in sorted(allelecount[k][p].values()):
                                                                                allelecount[k][p][n] = 0
                                                                        else:
                                                                                pass

                        print "   Tabulating allele counts..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one individual
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, input -1s
                                                                for a in sorted(allelecount[popsdict[k][x]][i].iterkeys()): #Cycle through 1 or 2 alleles
                                                                        allelecount[popsdict[k][x]][i][a] = -1

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        count = 0
                                                                        while count < 2:
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))>=2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                allelecount[popsdict[k][x]][p][seqsdict[x][p][a]] += 1
                                                                                                        if int(len(seqsdict[x][p]))>2:
                                                                                                                print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."

                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        allelecount[popsdict[k][x]][p][seqsdict[x][p][a]] += 2
                                                                                                        count += 1
                                                                                count += 1

                        print "   Writing file..."
                        #Print out counts of allele by individual
                        for k in sorted(allelecount.iterkeys()):
                                count = 0
                                while count < 2:
                                        if count == 0: fout.write(str(k)+'a\t')
                                        if count == 1: fout.write('\n'+str(k)+'b\t')
                                        for p in sorted(allelecount[k].iterkeys()):
                                                for n in sorted(allelecount[k][p].iterkeys()):
                                                        if count == 0:
                                                                if allelecount[k][p][n] == 0:
                                                                        fout.write('0\t')
                                                                if allelecount[k][p][n] == 1:
                                                                        allelecount[k][p][n] = 0
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == 2:
                                                                        allelecount[k][p][n] = 1
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == -1:
                                                                        fout.write('NaN\t')
                                                                if allelecount[k][p][n] > 2:
                                                                        print "Error! Program retained more than 2 of a particular allele for one individual/locus. Check code."
                                                                        exit(1)
                                                                        # This shouldn't ever happen but the code hasn't been thoroughly tested.
                                                        elif count == 1:
                                                                if allelecount[k][p][n] == 0:
                                                                        fout.write('0\t')
                                                                if allelecount[k][p][n] == 1:
                                                                        allelecount[k][p][n] = 0
                                                                        fout.write('1\t')
                                                                if allelecount[k][p][n] == -1:
                                                                        fout.write('NaN\t')
                                                                if allelecount[k][p][n] > 1:
                                                                        print "Error! Program retained more than 2 of a particular allele for one individual/locus. Check code."
                                                                        exit(1)
                                                                        # This shouldn't ever happen but the code hasn't been thoroughly tested.
                                        count += 1

                                fout.write('\n')

                        fout.close()

                if HaploChoice == 5: #Bayescan format
                        print "Outputting Bayescan file..."
                        fout.write('[loci]='+str(len(OrderedLoci))+'\n\n[populations]='+str(num_pops)+'\n')


                        popfreq = {} #Build Structure: {PopulationID : {LocusID : {AlleleInteger : Count} } }
                                     #First build empty dictionary structure with 0 count for each allele

                        print "   Creating dictionary of allele counts..."
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                if k in sorted(popfreq.keys()):
                                        for p in sorted(newseqdict.iterkeys()): #Look at one locus
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()): #Look at one allele integer
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                else:
                                        popfreq[k] = {}
                                        for p in sorted(newseqdict.iterkeys()):
                                                if p in sorted(popfreq[k].iterkeys()):
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass
                                                else:
                                                        popfreq[k][p] = {}
                                                        for n in sorted(newseqdict[p].itervalues()):
                                                                if n not in sorted(popfreq[k][p].values()):
                                                                        popfreq[k][p][n] = 0
                                                                else:
                                                                        pass

                        print "   Tabulating allele counts..."
                        #Add counts of each allele
                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                for i in OrderedLoci: #Look at one locus
                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                count = 0
                                                                while count < 2:
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))>=2:#If heterozygote
                                                                                                if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                        popfreq[k][p][seqsdict[x][p][a]] += 1
                                                                                                if int(len(seqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."

                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                popfreq[k][p][seqsdict[x][p][a]] += 2
                                                                                                count += 1
                                                                        count += 1

                        print "   Writing file..."
                        #Print out Bayescan format (Locus  \t  NumGeneCopiesPop  \t  NumAllelesAtLocus  \t  Counts  \t  For  \t  Each  \t  Allele...)
                        for k in sorted(PopsOrdered.iterkeys()):
                                fout.write('\n[pop]='+str(PopsOrdered[k])+'\n')
                                for p in sorted(popfreq[k].iterkeys()):
                                        genecopies = 0
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                genecopies += popfreq[k][p][n]
                                        fout.write(str(LociOrdered[p])+'\t'+str(genecopies)+'\t'+str(len(newseqdict[p])))
                                        for n in sorted(popfreq[k][p].iterkeys()):
                                                fout.write('\t'+str(popfreq[k][p][n]))
                                        fout.write('\n')

                        fout.close()

                        fout=open(outfile_loci,"w")
                        for i in sorted(LociOrdered.iterkeys()):
                                fout.write(str(LociOrdered[i])+'\t'+str(i)+'\n')
                        fout.close()

                        fout=open(outfile_pops,"w")
                        for i in sorted(PopsOrdered.iterkeys()):
                                fout.write(str(PopsOrdered[i])+'\t'+str(i)+'\n')
                        fout.close()

                if HaploChoice == 6: #Arlequin format
                        OrderedLoci = []
                        title = raw_input("Title of Project: ")
                        fout=open(outfile,"w")
                        print "Outputting Arlequin file..."

                        fout.write("[Profile]\n\n\t\"" + title + "\"\n\n\t\tNbSamples=" + str(num_pops))
                        fout.write("\n\t\tGenotypicData=1\n\t\tGameticPhase=0\n\t\tDataType=STANDARD\n\t\t")
                        fout.write("LocusSeparator=TAB\n\t\tMissingData=\"?\"\n\n\n[Data]\n\n\t[[Samples]]\n\n")

                        for key in sorted(num_sites.iterkeys()):
                                OrderedLoci.append(key)

                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                fout.write("\t\tSampleName= \"Pop_" + str(k) + "\"\n\t\tSampleSize=" + str(int(gene_copies[k])/2) + "\n\t\tSampleData={\n")
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                count = 0
                                                while count < 2:
                                                        ind = str(popsdict[k][x])
                                                        if count == 0: fout.write(ind+'\t1\t')
                                                        if count == 1: fout.write('\t\t')
                                                        for i in OrderedLoci: #Look at one locus
                                                                if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write ?s
                                                                        fout.write("?\t")

                                                                else:
                                                                        for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                                if p == i:      #If this is the right locus
                                                                                        for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                                if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                                elif int(len(seqsdict[x][p]))>2:
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                        if (count==0 and int(a)==0) or (count==1 and int(a)==1):
                                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                        count += 1
                                                        fout.write('\n')
                                fout.write("}\n")
                        fout.close()

                if HaploChoice == 7: # GenAlEx format
                        print "Outputting GenAlEx file..."

                        fout=open(outfile,"w")

                        num_inds = 0
                        for k in popsdict.iterkeys():
                                num_inds += int(len(popsdict[k]))

                        fout.write(str(len(OrderedLoci))+'\t'+str(num_inds)+'\t'+str(num_pops))
                        for k in sorted(popsdict.iterkeys()):
                                fout.write('\t'+str(len(popsdict[k])))

                        fout.write('\n\t\t')
                        for k in sorted(popsdict.iterkeys()):
                                fout.write('\t'+str(k))
                        fout.write('\nIndID\tPopID\t')

                        count=0
                        for i in OrderedLoci:
                                if count > 0:
                                        fout.write('\t\t')
                                fout.write(str(i))
                                count += 1
                        fout.write('\n')



                        for k in sorted(popsdict.iterkeys()): #Look at one population
                                for x in sorted(seqsdict.iterkeys()): #Cycle through all seqs by individual
                                        if x in popsdict[k].keys(): #If that individual's in the population
                                                ind = str(popsdict[k][x])
                                                fout.write(str(ind)+'\t'+str(k)+'\t')
                                                for i in OrderedLoci: #Look at one locus
                                                        if i not in seqsdict[x].keys(): #If individual doesn't have this locus, write 0s
                                                                fout.write('0\t0\t')

                                                        else:
                                                                for p in sorted(seqsdict[x].iterkeys()): #Cycle through its loci
                                                                        if p == i:      #If this is the right locus
                                                                                for a in sorted(seqsdict[x][p].iterkeys()): #Cycle through 1 or 2 alleles
                                                                                        if int(len(seqsdict[x][p]))==2:#If heterozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))==1:#If homozygote
                                                                                                fout.write(str(seqsdict[x][p][a])+'\t'+str(seqsdict[x][p][a])+'\t')
                                                                                        elif int(len(seqsdict[x][p]))>2:
                                                                                                if (int(a)==0 or int(a)==1):
                                                                                                        print "Warning, Number of individual's alleles found greater than 2! Retaining first 2 only."
                                                                                                        fout.write(str(seqsdict[x][p][a])+'\t')
                                                fout.write("\n")
                        fout.close()


        except IOError:
                print "Error: Problems outputting file. Check the directory path."
                return 0

        return 1




if choice == 1:
	Job2 = Fasta2Migrate(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 2:
	Job2 = Fasta2Arlequin(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 3:
	Job2 = Fasta2DIYABC(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 4:
	Job2 = Fasta2LFMM(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 5:
	Job2 = Fasta2Phylip(num_pops, popsdict, seqsdict, gene_copies, num_sites, haplo, breakpoints, haplotypes)
if choice == 6:
	Job2 = Fasta2GPhocs(num_pops, popsdict, seqsdict, gene_copies, num_sites)
if choice == 7:
	Job2 = Fasta2Treemix(num_pops, popsdict, seqsdict, gene_copies, num_sites, one_snp)
if choice == 8:
	Job2 = Fasta2Haplotype(num_pops, popsdict, seqsdict, gene_copies, num_sites, HaploChoice)
