#!/usr/bin/env python3
#
# translate a vcf file to migrate
# this takes a file that either contains diploid or haploid individuals
# with a reference sequence it generates full DNA migrate file
# without a reference sequence it generages a SNP migrate file
# for more help see the help() function
#
# part of the migrate distribution (and of git repository ....[soon to be added]
#
# MIT licence
# (c) Peter Beerli, Tallahassee 2020 
import sys
import gzip
import datetime as dt

DELIM = "|/"
NMLEN = 10

def help(args):
    print("syntax: vcf2mig --vcf vcffile.vcf <--ref ref1.fasta,ref2.fasta,... | --linksnp number >  <--popspec numpop ind1 ind2 .... | --pop populationfile.txt> <--chrom chr1,chr2,...> --out migrateinfile")
    print("  --vcf vcffile : a VCF file that is uncompressed or .gz, currently only")
    print("                  few VCF options are allowed, simple reference")
    print("                  and alternative allele, diploid and haploid data")
    print("                  can be used")
    print("  --ref ref1.fasta,ref2.fasta,... : reference in fasta format")
    print("                  several references can be given, for example for")
    print("                  each chromosome, if this option is NOT present then")
    print("                  the migrate dataset will contain only the SNPs")
    print("  --linksnp number : cannot not be used with --ref; defines linkage groups of snps")
    print("                  for example in a VCF file covering 10**9 sites, a value of 100000")
    print("                  will lead to 10 linked snp loci, if this option and the --ref are")
    print("                  are missing, then the resulting dataset will contain single, unlinked snps")
    print("  --popspec numpop ind1,ind2,... : specify the population structure, number of populations")
    print("                  with the number of individuals for each population")
    print("                  This option exlcudes the option --pop")
    print("  --pop popfile:  specify a file that contains a single line with")
    print("                  numpop ind1,ind2 in it")
    print("                  This option exlcudes the option --popspec")
    print("  --chrom chr1,chr2,... sepcify subset of chromosomes in vcf file")
    print("                  if all chromosomes are used ignore this option")
    print("  --out migratedatafile:  specify a name for the converted dataset in migrate format")    
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

        # search for linked snps
        key = '--linksnps'
        if key in argstring:
            linkedsnps = int(args[args.index(key)+1])
            print("linked snps",linkedsnps)
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
    except:
        print(key)
        print(args)
        help(args)
        sys.exit(-1)
    return vcffile, referencefile, linkedsnps, numind, numloc, migratefile,use_chrom

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
    
    #for line in lines:
    print(f"In read_header@{lines}")
    return lines

def find_header(header,key):
    values = [h for h in header if key in h]
    return values

def read_body(vcffile, header):
    if '.gz' in vcffile:
        f = gzip.open(vcffile,'rb')
        mygzip = True
    else:
        f = open(vcffile,'r')
        mygzip = False
    #print("@",header)
    variables = header[-1].split()
    #print(f"@{variables}")
    #numcol = len(variables)
    minimal=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    d = {}
    for v in variables:
        d[v.replace('#','')] = find_header(header[:-1],v.replace('#',''))
    #for di in d.keys():
    #    print(di,d[di])
    #print(variables)
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
        print(line[:50])
        chrom = a[0]
        pos   = int(a[1])
        id    = a[2]
        ref   = a[3]
        alt   = a[4]        
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


#def read_vcf2(vcffile):


def read_vcf(vcffile):
    '''
    input: a vcffile name
    ref is a single nuc
    alt can be a comma delimited list of alternatives
    output: a list with [locus, pos, ref, alt, ploidy, list with either 0,1,2 or [a1,a2]]
    exactly as [[chrom,pos,id,ref,alt,qual,filter,info,format,individuals],names,chroms,ploidy]
    '''
    header = read_vcf_header(vcffile)
    data, names, chroms = read_body(vcffile,header)
    #print(data)
    #print("first:",data[0][-1][0])
    #if "|" in data[0][-1][0]:
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
    numpop=x[0]
    numind = [int(xi) for xi in x[1:]]
    return numpop,numind

def read_reference(file):
    x = int(file)
    if x>1:
        return [[">gugus",'A'*x,x]]
    references=[]
    allfiles = file.split(',')        
    for fi in allfiles:
        f = open(fi,'r')
        head = f.readline()
        sequence = f.read().strip()
        if ">" in sequence:
            while ">" in sequence:
                h = sequence.index(">")
                newseq = "".join(sequence[:h].split())
                newsites = len(newseq)
                references.append([head,newseq,newsites])
                #print(h,[head,newseq,newsites])
                sequence = sequence[h:]
                head = sequence.split('\n',1)
                sequence = head[1]
                head = head[0]
                sites = len(sequence)
        sequence = "".join(sequence.split())
        sites = len(sequence)
        #print(h,[head,newseq,newsites])
        references.append([head,sequence, sites])
        f.close()

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
        print(use_chrom)
        print(chroms)
        for c in use_chrom:
            if c not in " ".join(chroms):
                print("mismatch between chromomose in {use_crom} and in VCF file {chroms}")
                sys.exit(-9)
    return use_chrom

def create_pop(references,vcf,begin, stop, use_chrom, ploidy):
    global snpcount
    ind = []
    snps = False
    #print(ind)
    if references != None:
        for k,r in enumerate(references):
            ref = r
            newind=[]
            for i in range(begin,stop):
                for p in range(ploidy):
                    newind.append(list(ref))
            ind.append(newind)
    else:
        snps = True
        count = 0
        for v in vcf:
            chrom = v[0]
            if chrom not in use_chrom:
                continue
            count +=1
        newind=[]
        snpcount = count
        print("variable sites in VCF:",count)
        for i in range(begin,stop):
            for p in range(ploidy):
                newind.append(list('.'*count))
        ind.append(newind)
        
    indix = range(len(use_chrom))
    count = 0
    positions = []
    for v in vcf:
        chrom = v[0]
        if chrom not in use_chrom:
            continue
        chrom = use_chrom.index(chrom)
        if snps:
            pos = count
            count += 1
            positions.append([chrom,v[1]])
        else:
            pos = v[1]
        ref = v[3]
        alt = v[4].split(',')
        individuals = v[9]
        #print(individuals)
        #print(v)
        if chrom >= len(ind):
            print(f"Problem with #loci>={chrom} and # of references={len(ind)}")
            sys.exit(-1)
        #print(">>>>>>>",begin,stop)
        for i,j in enumerate(range(begin,stop)):
            s = individuals[j].split(':')[0]
            #print(f'{s=}')
            sep = next((c for c in DELIM if c in s), None)
            # split on it (or wrap in a list if none)
            individuals2 = s.split(sep) if sep else s
            #if stop==20:
            #print("@", 1, i, j, ref,alt,individuals2, pos) #, ind[chrom][i][pos])
            #sys.exit()
            if type(individuals2)==list:
                ip = i*ploidy
                cadd=-1
                for iind in individuals2:
                    cadd += 1
                    try:
                        if iind == '0':
                            ind[chrom][ip+cadd][pos]  = ref
                        elif iind == '1':
                            #print("@@@@@", chrom,i,ip,cadd,pos,alt,ref)
                            ind[chrom][ip+cadd][pos] = alt[0]
                        else:
                            #print("@2    ", chrom, ip, alt, ref, pos, ref, individuals)
                            ind[chrom][ip+cadd][pos] = alt[1]
                            #ind[chrom][i] == "".join(xx)
                    except:
                        print("EXCEPT",ip,cadd,pos,ref,alt[0])
    
            else:                
                if individuals2 =='0':
                    #print(pos, ref, individuals)
                    ind[chrom][i][pos]  = ref
                elif individuals2 == '1':
                    ind[chrom][i][pos] = alt[0]
                else:
                    #print(chrom, i, alt, ref, pos, ref, individuals2)
                    #print("###", individuals2)
                    ind[chrom][i][pos] = alt[1]
                    #ind[chrom][i] == "".join(xx)
                #if stop==20:
                #    print(2, i, j, ref, alt, individuals2, pos, ind[chrom][i][pos])
            #if individuals[i]=='1':
            #    print(ind[chrom][i][pos], ref,alt)
            #    sys.exit(-1)
    ni=[]        
    for i in ind:
        newind=[]
        for j in i:
            newind.append("".join(j))
            #print(f'{"".join(j):6.6s}')
        ni.append(newind)
    if snps:
        return ni,positions
    else:
        return ni,None
    
#    for v in vcf:
#        chrom = v[0]-1
#        ref = references[chrom][2]
        
def write_migrate(migratefile, data, sites, references, names):
    f = open(migratefile,'w')
    if references==None:
        if linkedsnps==None:
            f.write(f" {len(data)} {sites} Translated from VCF {dt.date.today()}\n")
        else:
            f.write(f" {len(data)} {len(sites)} Translated from VCF {dt.date.today()}\n")
    else:
        print(f" {len(data[0])} {len(data[0][0])} Translated from VCF {dt.date.today()}\n")
        f.write(f" {len(data[0])} {len(data[0][0])} Translated from VCF {dt.date.today()}\n")
    f.write(    f"# VCF file used:      {vcffile}\n")
    if references==None:
        f.write(f"# SNP data file!\n")
    else:
        for ref in references:
            f.write(f"# Reference file:     {ref}\n")
    f.write(    f"# Migrate input file: {migratefile}\n")
    positions = data[0][1]
    if positions==None:
        for s in sites:
            f.write(f"(s{s}) ")
        f.write("\n")
    else:
        if linkedsnps==None:
            for pi in positions:
                f.write("(n1) ")
            f.write("\n")
        else:
            for si in sites:
                    f.write(f"(n{si}) ")
            f.write("\n")
    for i,pop in enumerate(data):
        pop1,pop2 = pop
        count = 0
        newpop = list(zip(*pop1))
        f.write(f" {len(newpop)} Pop{i+1}\n")
        #print(f" {len(newpop)} Pop{i+1}\n")
        for z,ind in enumerate(newpop):
            #print(z,len(ind))
            #for ip in range(ploidy):
            ploi = f":{((z % ploidy)+1)}"
            thename = f"{names[i][int(count/ploidy)]:{NMLEN}.{NMLEN}}".rstrip()
            #print(ploi,count, ploidy, i, len(newpop))
            nname = f"{z}{thename}{ploi}"
            print(nname)
            f.write(f"{nname:20s} ")
            print(f"{nname:20s} ")
            for locus in ind:
                f.write(f"{locus}")
            f.write("\n")
            count += 1
    f.close()

def link_snps(pos,size):
    s=[]
    h=0
    count=0
    oldpi0 = pos[0][0]
    for pi in pos:
        if pi[0]!=oldpi0:
            oldpi0=pi[0]
            h = 0
            count=0
            
        if pi[1] < h + size:
            count += 1
        else:            
            h += size
            s.append(count)
            count = 0
    s.append(count)
    return s

    
if __name__ == "__main__":

    numpop = -1
    vcffile, referencefile, linkedsnps, numind, numloc, migratefile, use_chrom = parse_args(sys.argv)
    print("parsed options")
    vcf,names, chroms, ploidy  = read_vcf(vcffile)
    print("read VCF file")
    #if numpop == -1:
    #    numpop, numind = read_populations(populationfile)
    #print("read populations")
    #else: already done in parse_args using --popspec
    print(f"VCF file used: {vcffile}")
    if referencefile != None:
        references = read_reference(referencefile)
        print("finished reading references")
        for ref in referencefile.split(','):
            print(f"Reference file: {ref}")
        refheaders, references, sites = list(zip(*references))
        snps=False
    else:
        # we only report snps and if linkedsnps !=None then the snps are linked within blocks of that size
        references = None
        refheaders = None
        sites = linkedsnps
        snps=True
        
    start = 0
    #print("@", sites)
    populations=[]
    data = [] #this contains all data for all population popxlocixindividuals
    use_chrom = harmonize_use_chroms(use_chrom,chroms)
    for ni in numind:
        #print("numind",ni)
        populations.append(names[start:ni+start])
        data.append(create_pop(references,vcf, start, ni+start,use_chrom,ploidy))
        start += ni
    if snps:
        positions = data[0][1]
        if linkedsnps!=None:
            sites = link_snps(positions,linkedsnps)
        else:
            sites = len(positions)
        referencefiles = None
    else:
        referencefiles = referencefile.split(',')
        
    write_migrate(migratefile, data, sites,referencefiles, populations)
    print(f"Migrate input file: {migratefile}")
