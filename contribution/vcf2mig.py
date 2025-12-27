#!/usr/bin/env python3
#
# translate a vcf file to migrate
# this takes a file that either contains diploid or haploid individuals
# with a reference sequence it generates full DNA migrate file
# without a reference sequence it generages a SNP migrate file
# for more help see the help() function
#
# part of the migrate distribution
#
# Created: 2020
# Considerably modified and improved December 2025
#
# MIT licence
# (c) Peter Beerli, Tallahassee 2020-2025
import sys
import gzip
import datetime as dt
from collections import Counter

IUPACTRANS = {'A': ['A'],'C': ['C'],'G': ['G'],'T': ['T'],
              'U': ['T'],'R': ['A','G'],'Y': ['C','T'],
              'M': ['A','C'],'K': ['G','T'], 'S': ['G','C'],
              'W': ['A','T'],'H': ['A','C','T'], 'B': ['C','G','T'],
              'V': ['A','C','G'], 'D': ['A','G','T'],
              'N': ['A','C','G', 'T'], 'X': ['A','C','G', 'T'],
              '?': ['A','C','G', 'T']}
IUPAC =    'URYMKSWHBVDNX'
NONIUPAC = '?????????????'
CHROMLINK = -1
DELIM = "|/"
NMLEN = 10 

def help(args):
    print("syntax: vcf2mig --vcf vcffile.vcf")
    print("               <<--ref|--abbrevref> ref1.fasta,ref2.fasta,... | --linksnp number >  ")
    print("               <--popspec numpop ind1 ind2 .... | --pop populationfile.txt>")
    print("               <--chrom chr1,chr2,...>")
    print("                --out migrateinfile\n\nDetails:")
    print("  --vcf vcffile : a VCF file that is uncompressed or .gz, currently only")
    print("                  few VCF options are allowed, simple reference")
    print("                  and alternative allele, diploid and haploid data")
    print("                  can be used")
    print("  --abbrevref ref1.fasta,ref2.fasta,... : reference in fasta format")
    print("                  for more info see next option, returns snps + invariant counts")
    print("  --ref ref1.fasta,ref2.fasta,... : reference in fasta format")
    print("                  several references can be given, for example for")
    print("                  each chromosome, if this option is NOT present then")
    print("                  the migrate dataset will contain only the SNPs")
    print("  --allowindel   if there are indels or deletions they will be used and not deleted")
    print("  --linksnp <number|chrom>: cannot not be used with --ref; defines linkage groups of snps")
    print("                  the keyword 'chrom' will link all snps within one chromosome (the VCF tag CHROM") 
    print("                  the 'number' specifies the distance among snps that are linked")
    print("                  read from first to last snp, so if number=1000 and the first snp is at position x")
    print("                  then all snps within the x+1000 will belong to the linkage group, is done for each chrom")
    print("                  If this option and the --ref are are missing, then the resulting dataset")
    print("                  will contain single, unlinked snps")
    print("  --popspec numpop ind1,ind2,... : specify the population structure, number of populations")
    print("                  with the number of individuals for each population")
    print("                  This option excludes the option --pop; if the numbers do not match the VCF file")
    print("                  then the options takes precedence and distributes according to --popspec")
    print("  --pop popfile:  specify a file that contains a single line with (use spaces!)")
    print("                  numpop ind1 ind2 ... ")
    print("                  This option exlcudes the option --popspec")
    print("  --chrom chr1,chr2,... specify subset of chromosomes in vcf file")
    print("                  if all chromosomes are used ignore this option")
    print("  --out migratedatafile:  specify a name for the converted dataset in migrate format")
    print("  --strict: replaces all characters that are not ACGTN? with ?")    
    print("")
    print("Example:")
    print("vcf2mig.py --vcf vcffile.vcf.gz --ref ref.fasta --popspec 2 10,10 --out migratefile")
    print("vcf2mig.py --vcf vcffile.vcf --popspec 3 10,10,10 --out migratefile")
    print("vcf2mig.py --vcf vcffile.vcf --linksnp 10000 --popspec 2 5,10 --out migratefile") 
    print("")
    print(f"\n\nYou specified:{args}")


def parse_args(args):
    '''parse the commandline arguments'''
    popset        = False
    use_chrom     = None
    linkedsnps    = None
    referencefile = None
    numind = []
    numloc = []
    migratefile = None
    strictvcf = False
    refabbrev= False
    argstring = " ".join(args)
    if "--help" in argstring or "-h" in argstring or "-help" in argstring:
        help(args)
        sys.exit(-1)
        
    try:
        # search for vcffile
        key = '--vcf'
        vcffile = args[args.index(key)+1]
        # search for referencefile
        key = '--ref'
        if key in argstring:
            referencefile = args[args.index(key)+1]
        key = '--abbrevref'
        if key in argstring:
            referencefile = args[args.index(key)+1]
            refabbrev=True
        else:
            refabbrev=False
        # search for linked snps
        key = '--linksnp'
        if key in argstring:
            linkedsnps = args[args.index(key)+1]
            if linkedsnps[0] == 'C' or linkedsnps[0] == 'c':
                linkedsnps = CHROMLINK
            else:
                linkedsnps = int(linkedsnps)            
        # allow indels and deletions
        key = '--allowindel'
        if key in argstring:
            allowindel = True
        else:
            allowindel = False
        # search for populationspec
        key = '--popspec'
        if key in argstring:
            numpop = int(args[args.index(key)+1])
            numind = args[args.index(key)+2]
            numind = [int(x) for x in numind.split(',')]
            populationfile = None
            popset=True
        # search for chromosome specification     
        key = '--chrom'
        if key in argstring:
            use_chrom = args[args.index(key)+1]
            use_chrom = use_chrom.strip().split(',')
        # search for populationfile
        key = '--pop'
        if key in argstring and not popset:
            populationfile = args[args.index(key)+1]
            numind,numloc = read_populations(populationfile)
            popset=True
        if not popset:
            raise(NameError)
        # search for migratefile
        key = '--out'
        if key in argstring:
            migratefile = args[args.index(key)+1]
        key = '--strict'
        if key in argstring:
            strictvcf = True
        else:
            strictvcf = False
    except:
        print(key)
        print(args)
        help(args)
        sys.exit(-1)
    return vcffile, referencefile, linkedsnps, numind, numloc, migratefile,use_chrom,allowindel,strictvcf, refabbrev

# parses the header of the VCF file
def read_vcf_header(vcffile):
    '''
    reads vcf file header material without processing, it 
    also reads the first non-header line
    '''
    if '.gz' in vcffile:
        f = gzip.open(vcffile,'r')
        mygzip=True
    else:
        f = open(vcffile,'r')
        mygzip = False
    lines=[]
    for line in f:
        if mygzip:
            line = line.decode('ascii')
        #print("@@",line, line[0], line[0]=='#')
        if line[0]!='#':
            break
        lines.append(line.strip())
    f.close()    
    return lines

# parses VCF file header
def find_header(header,key):
    values = [h for h in header if key in h]
    return values

# parses the data part of the VCF file
def read_body(vcffile, header):
    if '.gz' in vcffile:
        f = gzip.open(vcffile,'rb')
        mygzip = True
    else:
        f = open(vcffile,'r')
        mygzip = False
    variables = header[-1].split()
    minimal=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    d = {}
    for v in variables:
        d[v.replace('#','')] = find_header(header[:-1],v.replace('#',''))
    if "FORMAT" in " ".join(variables):
        minimal.extend(["FORMAT"])
        indstart = variables.index('FORMAT')+1
    else:
        indstart = variables.index("INFO")+1
    names = variables[indstart:]
    data = []
    for line in f:
        if mygzip:
            line = line.decode('ascii')
        if line[0]=='#':
            continue
        a = line.strip().split()
        #print(line[:50])
        chrom = a[0]
        pos   = int(a[1])
        id    = a[2]
        ref   = a[3]
        alt   = a[4]
        if not allowindel:
            if len(ref)>1:
                continue
            if max(list(map(len,alt.split(','))))>1:
                continue
        qual  = a[5]
        filter = a[6]
        info = a[7]
        if "FORMAT" in minimal:
            format = a[8]
            individuals = a[9:]
        else:
            format = '.'
            individuals = a[8:]
        data.append([chrom,pos,id,ref,alt,qual,filter,info,format,individuals])     
    f.close()
    chroms = list(set([str(d[0]) for d in data]))
    return data,names,chroms


# parses vcf file, this is minimal, lots of the details will be ignored
# main goal is parsing location and ref and alt and samples
def read_vcf(vcffile):
    header = read_vcf_header(vcffile)
    data, names, chroms = read_body(vcffile,header)
    s = data[0][-1][0]
    count = sum(s.count(d) for d in DELIM)
    if count:
        ploidy = count + 1
    else:
        ploidy = 1
    return data, names, chroms, ploidy


def read_populations(populationfilename):
    f = open(populationfilename,'r')
    x = f.read().split()
    f.close()
    numpop=int(x[0])    
    numind = [int(xi) for xi in x[1:]]
    return numind,numpop

def calculate_freq(sequence, head):
    counts = Counter(sequence)
    #print(counts)
    counts = sorted(counts.items())
    counts.insert(0,("Loc", head[1:].split()[0]))
    return dict(counts)

def read_reference(file):
    references=[]
    allfiles = file.split(',')        
    for fi in allfiles:
        f = open(fi,'r')
        head = f.readline()
        mysequence = f.read().strip()
        if ">" in mysequence:
            while ">" in mysequence:
                h = mysequence.index(">")
                newseq = "".join(mysequence[:h].split())
                newsites = len(newseq)
                if strictvcf:
                    trans = newseq.maketrans(IUPAC,NONIUPAC)
                    newseq = newseq.translate(trans)
                myfreqs = calculate_freq(newseq.upper(),head)
                references.append([head,newseq,newsites,myfreqs])
                mysequence = mysequence[h:]
                head = mysequence.split('\n',1)
                mysequence = head[1]
                head = head[0]
                sites = len(mysequence)
        mysequence = "".join(mysequence.split())
        if strictvcf:
            trans = mysequence.maketrans(IUPAC,NONIUPAC)
            mysequence.translate(trans)
        myfreqs = calculate_freq(mysequence.upper(), head)
        #print("myfreqs", myfreqs)
        sites = len(mysequence)
        references.append([head, mysequence, sites, myfreqs])
        f.close()

    #print("lreflast", len(references[-1]), len(references))
    #print(references[0][0],references[0][1][:10],references[0][2],references[0][3])
    return references

def harmonize_use_chroms(use_chrom, chroms):
    
    if use_chrom == None:
        if references == None:
            use_chrom = sorted(chroms)
            return use_chrom
        if len(chroms) == len(references):
            use_chrom = sorted(chroms)
        else:
            print("Number of chromosomes and number of references mismatch")
            print(use_chrom,len(references))
            print(chroms)
            sys.exit(-10)
    else:
        lu = len(use_chrom)            
        use_chrom = list(set(use_chrom).intersection(chroms))
        if use_chrom == None:
            #check for : or = and remove those
            print("Mismatch with names in --chrom and reference file")
            chroms = [ ch.replace(':','').replace('=','') for ch in chroms ]
            use_chrom = list(set(use_chrom).intersection(chroms))
            if use_chrom == None:
                print("Warning, no match with vcf CHROM tag")
        if len(use_chrom) != lu:
            print("the --chrom list did only partly match with the vcf, used these:")
            print(f"{use_chrom}")
    return use_chrom

# takes vcf data and parses its content into a population structure
# containint all individuals and all chromosomes, snps
def create_pop_snps(vcf,begin, stop, use_chrom, ploidy):
    
    idata = []
    positions = []
    for vi,v in enumerate(vcf):
        chrom = v[0]
        if chrom not in use_chrom:
            continue
        chrom = use_chrom.index(chrom)
        pos = v[1]
        positions.append((chrom,pos))
        ref = v[3]
        alt = v[4].split(',')
        if allowindel:
            lr = len(ref)
            if lr>1:
                #deletion?
                for ai in range(len(alt)):
                    if len(alt[ai]) < lr:
                        alt[ai] = alt[ai].ljust(lr, '-')
            else: # insertion?
                la = []
                for ai in range(len(alt)):
                    la.append(len(alt[ai]))
                maxla = max(la)
                if(maxla>1):
                    ref = '@'+ref.ljust(maxla,'-')
                    for ai in range(len(alt)):
                        alt[ai] = '@'+alt[ai].ljust(maxla,'-')
                else:
                    pass #standard point mutation
        individuals = v[9]
        popindividuals = individuals[begin:stop]
        sdata = []
        for ri, rawsample in enumerate(popindividuals):
            s = rawsample.split(':')[0] #removes the quality scores
            sep = next((c for c in DELIM if c in s), None)
            sample = s.split(sep) if sep else s
            if type(sample) == list:
                for pi in range(ploidy):
                    if sample[pi] == '0':
                        sdata.append(ref)
                    elif sample[pi] != '.':
                        sdata.append(alt[int(sample[pi])-1])
                    else:
                        sdata.append('?')
        idata.append(sdata)
    return idata,positions

def convert_chrompos(positions):
    chromset = dict()
    for pi in positions:
        if  pi[0] in chromset:
            chromset[pi[0]].append(pi[1])
        else:
            chromset[pi[0]] = [pi[1]]
    return [items for items in chromset.items()]
        
        

def create_pop_references(references, snps, positions, use_chrom, ploidy):
    chrompos = convert_chrompos(positions)
    allpop = []
    for si in snps:  #population block
        #print("si", si)
        poploci = []
        #print("chrompos",chrompos)
        z = -1
        for chrom, pos in chrompos:  #locus block
            haplotypes=[references[chrom] for _ in si[chrom]] #seq for ind in locus
            #print("in ref: haplot", haplotypes)
            insertionmuts=[]
            for po in pos:   # fill in all the variants
                z += 1
                #print(f'si{z} {si[z]}')
                for hi, mut in enumerate(si[z]):
                    #print("@mut",mut, hi, po, pos, chrom, z)
                    if '@' in mut: #insertion
                        haplotypes[hi] = haplotypes[hi][:po] + '@' + haplotypes[hi][po+1:]
                        insertionmuts.append((hi,mut[1:]))
                    else:
                        haplotypes[hi] = haplotypes[hi][:po] + mut + haplotypes[hi][po+len(mut):]
            #print("in ref2: haplot", haplotypes)
            fix_haplotypes(haplotypes,insertionmuts)            
            poploci.append(haplotypes)
        allpop.append((poploci,positions))
        #sys.exit()
    return allpop
            
def fix_haplotypes(haplotypes,insertionmuts):
    search = '@'
    for hi, mut in insertionmuts:
        count = 0
        while True:
            index = haplotypes[hi].find(search)
            if index == -1:
                break
            haplotypes[hi] = (haplotypes[hi][:index] + mut + haplotypes[hi][index + 1:])
            count += 1
    return haplotypes

# Final result is stored in current_string


# writer for migrate modern format:
# data is either
#   - augmented refsequence+VCF data
#   - or snps from the VCF data 
def write_migrate(migratefile, data, freqs, positions, sites, references, names, comment):
    f = open(migratefile,'w')
    numpop = len(data)
    # header section
    #loci = len(list(set(list(zip(*positions))[0])))
    sites = list(map(len,[di[0] for di in data[0]]))
    loci = len(sites)
    chrompos = convert_chrompos(positions)
    if references != None and refabbrev == False:
        f.write(f'{numpop} {loci} {vcffile}\n')
        f.write(f"# VCF file used:      {vcffile}\n")
        f.write(f"# Translated from VCF {dt.date.today()}\n")
        f.write(f"# Reference file: {referencefile}\n")
        f.write(f"# Migrate input file: {migratefile}\n")
        f.write(f"# References augmented with VCF data file!\n")
        f.write(f"# {comment}\n")
        sitestr = " ".join([f'(s{si})' for si in sites])
        f.write(f"{sitestr}\n")
    else:
        if linkedsnps == None:
            unlinked_loci = len(sites)
            f.write(f'{numpop} {unlinked_loci} {vcffile}\n')
            f.write(f"# VCF file used:      {vcffile}\n")
            f.write(f"# Translated from VCF {dt.date.today()}\n")
            f.write(f"# Migrate input file: {migratefile}\n")
            f.write(f"# SNP data file!\n")
            f.write(f"# {comment}\n")
            sitestr = " ".join([f'(n1)' for _ in sites])
            f.write(f"{sitestr}\n")
        elif linkedsnps == CHROMLINK:
            chrompos = convert_chrompos(positions)
            sitestr = " ".join([f'(n{len(si[1])})' for si in chrompos])
            linkedloci = len(chrompos)
            f.write(f'{numpop} {linkedloci} {vcffile}\n')
            f.write(f"# VCF file used:      {vcffile}\n")
            f.write(f"# Translated from VCF {dt.date.today()}\n")
            f.write(f"# Migrate input file: {migratefile}\n")
            f.write(f"# SNP data file!\n")
            f.write(f"# {comment}\n")
            f.write(f"{sitestr}\n")
        else:
            nucs=[]
            delta = linkedsnps 
            chrompos = convert_chrompos(positions)
            for chrom in chrompos:
                nuc = 1
                x = 0
                for xi, pos in enumerate(chrom[1][1:]):
                    if pos < chrom[1][x] + delta:
                        nuc += 1
                    else:
                        nucs.append(nuc)
                        x = xi
                        nuc=1
                else:
                    nucs.append(nuc)
            loci = len(nucs)
            sitestr = " ".join([f'(n{si})' for si in nucs])
            f.write(f'{numpop} {loci} {vcffile}\n')
            f.write(f"# VCF file used:      {vcffile}\n")
            f.write(f"# Translated from VCF {dt.date.today()}\n")
            f.write(f"# Migrate input file: {migratefile}\n")
            f.write(f"# SNP data file!\n")
            f.write(f"# {comment}\n")
            f.write(f"{sitestr}\n")
    # write out the frequencies if present
    if freqs != None:
        fistr = check_freqsout(freqs)
        for fi, fis in enumerate(fistr):
            f.write(f"#*freq: {fi+1} {fis}\n")
    
    #individual name adjustments
    newnames=[]
    maxnamelen = 0
    for namepop in names:
        newnamepop = []
        for ni in namepop:
            if ploidy !=1:
                newnamepop1 = [ni+f":{plo+1}" for plo in range(ploidy)]
            else:
                newnamepop1 = [ni]
            newnamepop.extend(newnamepop1)
            maxlen = max(map(len,newnamepop1))
            if maxlen>maxnamelen:
                maxnamelen = maxlen
        newnames.append(newnamepop)
    names = newnames
    namelen = maxnamelen if maxnamelen > 10 else 10
    f.write(f"# individual name length is {namelen}!\n")
    # data section handles snps and reference augmented sequences
    for indx, (di,ni) in enumerate(zip(data,names)):
        dii = list(zip(*di))
        #print(f" {len(ni)} Pop{indx+1}")
        f.write(f" {len(ni)} Pop{indx+1}\n")
        for idx, d in enumerate(dii):
            #print(idx, end=' ')
            #print(f"{ni[idx]:<{namelen}}","".join(d))
            f.write(f'{ni[idx]:<{namelen}} {"".join(d)}\n')
    f.close()

def check_freqsout(freqs):
    fis = []
    for fi in freqs:
        if 'A' not in fi:
            fi['A']=0
        if 'C' not in fi:
            fi['C']=0
        if 'G' not in fi:
            fi['G']=0
        if 'T' not in fi:
            fi['T']=0
        if '?' not in fi:
            fi['?']=0
        
        alltotal = sum(value for key, value in fi.items() if key != 'Loc')
        total = sum(value for key, value in fi.items() if key in list('ACGT'))
        fistr = f"{fi['Loc']} ACGT={total} All={alltotal} A={fi['A']} C:{fi['C']} G:{fi['G']} T:{fi['T']} ?:{fi['?']}" 
        fis.append(fistr)
    return fis

    
if __name__ == "__main__":
    numpop = -1
    vcffile, referencefile, linkedsnps, numind, numloc, migratefile, use_chrom, allowindel, strictvcf, refabbrev = parse_args(sys.argv)
    print("Parsed options:")

    vcf,names, chroms, ploidy  = read_vcf(vcffile)
    print(f"VCF file used: {vcffile}")

    start = 0
    populations=[]
    data = [] 
    
    if referencefile != None:
        references = read_reference(referencefile)
        for ref in referencefile.split(','):
            print(f"Reference file: {ref}")
        refheaders, references, sites, freqs = list(zip(*references))
        bsnps=False
        
    else: # we only report snps and if linkedsnps !=None then the snps are linked 
        references = None
        refheaders = None
        sites = linkedsnps
        bsnps=True

    use_chrom = harmonize_use_chroms(use_chrom,chroms)

    if refabbrev:
        for fi in freqs:
            item = fi['Loc']
            fi['Loc'] = item.replace(':','').replace('=','')
        freqs = [fi for fi in freqs if fi['Loc'] in use_chrom]
    else:
        freqs = None

    for nii in numind:
        ni = nii
        populations.append(names[start:ni+start])
        data.append(create_pop_snps(vcf, start, ni+start,use_chrom,ploidy))
        start += ni
    snps, positions = zip(*data)
    if references != None and refabbrev==False:
        data = create_pop_references(references, snps,positions[0], use_chrom, ploidy)
        refdata, positions = zip(*data)
        referencefiles = referencefile.split(',')
    else:
        positions = positions[0]
        sites = len(positions)
    if references!=None and refabbrev==False:
        comment = 'Using references augmented VCF data'
        print(comment)
        #for re in refdata:
        #    for r in enumerate(re):
        #        print(r)
        #print(positions[0])
        write_migrate(migratefile, refdata, freqs, positions, sites,references, populations, comment)
    else:
        if linkedsnps == CHROMLINK:
            comment = "Using all SNPs, linked by chromosome"
        elif linkedsnps == None:
            comment = "Using unlinked SNPS"
        else:
            comment = f"Using all SNPs, linked every {linkedsnps} for each chromosome"
        print(comment)
        write_migrate(migratefile, snps,  freqs, positions, sites, references, populations, comment)
    print(f"Migrate input file: {migratefile}")
