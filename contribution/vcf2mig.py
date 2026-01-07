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
    print("               <--bound start,stop>")
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
    print("  --bound start,stop specifies the left and right bound of the reference sequence,")
    print("                  snps are only reported within the bound")
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
    bound         = None
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
        # search for bounds
        key = '--bound'
        if key in argstring:
            b = args[args.index(key)+1]
            b0, b1 = b.split(',')
            b0 = int(b0.replace("_",""))
            b1 = int(b1.replace("_",""))
            bound = (b0,b1)

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
    return vcffile, referencefile, linkedsnps, numind, numloc, migratefile, use_chrom, allowindel, strictvcf, refabbrev, bound

# parses the header of the VCF file
def read_vcf_header(vcffile):
    opener = gzip.open if vcffile.endswith('.gz') else open
    mode = 'rt' if vcffile.endswith('.gz') else 'r'
    lines = []
    with opener(vcffile, mode, encoding='ascii', newline='') as f:
        for line in f:
            if line and line[0] != '#':
                break
            lines.append(line.strip())
    return lines

# parses VCF file header
def find_header(header,key):
    values = [h for h in header if key in h]
    return values

# parses the data part of the VCF file
def read_body(vcffile, header):
    variables = header[-1].split()
    if "FORMAT" in " ".join(variables):
        indstart = variables.index('FORMAT') + 1
    else:
        indstart = variables.index("INFO") + 1
    names = variables[indstart:]

    opener = gzip.open if vcffile.endswith('.gz') else open
    mode = 'rt' if vcffile.endswith('.gz') else 'r'

    data = []
    with opener(vcffile, mode, encoding='ascii', newline='') as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            a = line.strip().split()
            chrom = a[0]
            pos   = int(a[1])
            ref   = a[3]
            alt   = a[4]

            if not allowindel:
                if len(ref) > 1:
                    continue
                if max(map(len, alt.split(','))) > 1:
                    continue

            qual   = a[5]
            filt   = a[6]
            info   = a[7]
            if len(a) > 8:
                fmt = a[8]
                inds = a[9:]
            else:
                fmt = '.'
                inds = a[8:]

            data.append([chrom, pos, a[2], ref, alt, qual, filt, info, fmt, inds])

    chroms = list(set(d[0] for d in data))
    return data, names, chroms


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
    counts = sorted(counts.items())

    # head can be a FASTA header (">chr22 ...") OR just a contig name ("chr22"/"22"/"1")
    if head is None:
        loc = "Loc"
    else:
        h = str(head).strip()
        if h.startswith(">"):
            parts = h[1:].split()
            loc = parts[0] if parts else "NA"
        else:
            # plain contig name
            loc = h.split()[0] if h.split() else "NA"

    counts.insert(0, ("Loc", loc))
    return dict(counts)



#def calculate_freq(sequence, head):
#    counts = Counter(sequence)
#    #print(counts)
#    counts = sorted(counts.items())
#    counts.insert(0,("Loc", head[1:].split()[0]))
#    return dict(counts)

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


def _score_ref_mapping_one_chrom(refseq, vcf_rows, max_check=500):
    """
    Score two candidate POS->index mappings:
      offset = -1  => idx = POS-1  (standard VCF)
      offset =  0  => idx = POS    (shifted)
    Only uses simple SNPs with len(REF)==1.
    Returns (best_offset, score_minus1, score_0, n_used).
    """
    if not refseq or not vcf_rows:
        return (-1, 0, 0, 0)

    L = len(refseq)
    s_m1 = 0
    s_0 = 0
    n = 0

    for row in vcf_rows:
        if n >= max_check:
            break
        pos = row[1]
        ref = row[3]
        if not ref or len(ref) != 1:
            continue
        r = ref.upper()
        if r not in "ACGTN":
            continue

        i_m1 = pos - 1
        i_0 = pos
        if 0 <= i_m1 < L and refseq[i_m1].upper() == r:
            s_m1 += 1
        if 0 <= i_0 < L and refseq[i_0].upper() == r:
            s_0 += 1
        n += 1

    best = -1 if s_m1 >= s_0 else 0
    return (best, s_m1, s_0, n)


def detect_pos_offsets_by_chrom(references_full, vcf, use_chrom, max_check=500):
    """
    Determine best POS mapping per chromosome in use_chrom order.
    Returns offsets list length len(use_chrom), values in {-1,0}.
    Prints warnings to stderr if shifted mapping wins or if match rate is suspiciously low.
    """
    rows_by_name = {}
    for row in vcf:
        rows_by_name.setdefault(row[0], []).append(row)

    offsets = []
    for ci, cname in enumerate(use_chrom):
        refseq = references_full[ci]
        rows = rows_by_name.get(cname, [])
        best, s_m1, s_0, n = _score_ref_mapping_one_chrom(refseq, rows, max_check=max_check)

        if n > 0:
            r_m1 = s_m1 / n
            r_0 = s_0 / n
            best_rate = max(r_m1, r_0)
        else:
            best_rate = 0.0

        if best == 0 and n >= 5:
            print(
                f"Warning: {cname}: REF matches better with idx=POS (offset=0) "
                f"than idx=POS-1 (offset=-1): score0={s_0}/{n} score-1={s_m1}/{n}",
                file=sys.stderr
            )

        if n >= 20 and best_rate < 0.60:
            print(
                f"Warning: {cname}: low REF match rate (best {best_rate:.2f}). "
                f"FASTA/VCF may be mismatched (assembly/contig/coordinate issue).",
                file=sys.stderr
            )

        offsets.append(best)
    return offsets


def filter_vcf_by_bound(vcf, use_chrom, bound):
    """Filter VCF rows by chromosome subset and inclusive bounds (start,stop)."""
    use_set = set(use_chrom)
    if bound is None:
        return [row for row in vcf if row[0] in use_set]
    start, stop = bound
    if start > stop:
        start, stop = stop, start
    out = []
    for row in vcf:
        if row[0] not in use_set:
            continue
        if start <= row[1] <= stop:
            out.append(row)
    return out


def crop_references_with_offsets(references_full, use_chrom, bound, offsets):
    """
    Crop reference sequences according to --bound using per-chrom offset.
    bound is inclusive [start,stop] in POS space.
    Cropping slice uses idx = POS + offset, so slice is [start+off : stop+off+1).
    If bound is None, returns references_full unchanged.
    """
    if bound is None:
        return references_full[:]

    b0, b1 = bound
    if b0 > b1:
        b0, b1 = b1, b0

    refs2 = []
    for ci, r in enumerate(references_full):
        off = offsets[ci]  # -1 or 0
        start_idx = b0 + off
        stop_excl = b1 + off + 1
        if start_idx < 0:
            start_idx = 0
        if stop_excl > len(r):
            stop_excl = len(r)
        if stop_excl < start_idx:
            stop_excl = start_idx
        refs2.append(r[start_idx:stop_excl])
    return refs2


def make_positions_rel_from_abs(positions_abs, base0_by_chrom):
    """positions_rel = (chrom_index, POS - base0) used only for reference augmentation."""
    return [(ci, pos - base0_by_chrom[ci]) for (ci, pos) in positions_abs]


def prepare_inputs(vcf, use_chrom, references, freqs, bound, refabbrev, strictvcf):
    """
    Central place for:
      - filtering VCF by chrom/bound (crops positions implicitly)
      - validating POS mapping vs reference (offsets -1 or 0 per chrom)
      - cropping references to bound consistently
      - recomputing freqs for --abbrevref on the cropped segments (if bound is set)

    Returns: vcf2, references2 (or None), freqs2 (or None), offsets (or None)
    """
    vcf2 = filter_vcf_by_bound(vcf, use_chrom, bound)

    if references is None:
        return vcf2, None, None, None

    references_full = list(references)  # full-length per chrom in use_chrom order
    offsets = detect_pos_offsets_by_chrom(references_full, vcf2, use_chrom, max_check=500)

    references2 = crop_references_with_offsets(references_full, use_chrom, bound, offsets)

    # freqs handling:
    # - If --abbrevref and bound is set: recompute on cropped segments (makes sense)
    # - If --abbrevref and no bound: keep existing freqs behavior (your prior code)
    freqs2 = None
    if refabbrev:
        if bound is not None:
            freqs2 = [calculate_freq(references2[i].upper(), use_chrom[i]) for i in range(len(use_chrom))]
        else:
            # Preserve existing freqs behavior as much as possible:
            # caller's old code cleaned fi['Loc'] and filtered by use_chrom.
            # We'll just pass freqs through; main can keep its existing cleanup/filter if you want.
            freqs2 = freqs

    return vcf2, references2, freqs2, offsets


def build_dataset(vcf, names, numind, use_chrom, ploidy,
                  references, refabbrev, bound, offsets):
    """
    Builds populations + data for write_migrate.

    Returns: populations, data_out, positions_out, sites
      - positions_out is ALWAYS absolute (chrom_index, POS)
      - data_out is snps (SNP-only/abbrevref) or refdata (full --ref)
    """
    start = 0
    populations = []
    data_blocks = []

    for nii in numind:
        ni = nii
        populations.append(names[start:ni+start])
        data_blocks.append(create_pop_snps(vcf, start, ni+start, use_chrom, ploidy))
        start += ni

    snps, pos_per_pop = zip(*data_blocks)
    positions_abs = pos_per_pop[0]  # absolute POS, cropped already by vcf filtering

    if references is not None and refabbrev is False:
        # build base0_by_chrom so positions_rel matches the CROPPED reference
        if bound is not None:
            b0, b1 = bound
            if b0 > b1:
                b0, b1 = b1, b0
            base0_by_chrom = [b0] * len(use_chrom)  # positions_rel = POS - bound_start
        else:
            # no bound: idx = POS + offset, so base0 = -offset
            # offset=-1 -> base0=1 (POS-1), offset=0 -> base0=0 (POS)
            base0_by_chrom = [(-off) for off in offsets]

        positions_rel = make_positions_rel_from_abs(positions_abs, base0_by_chrom)

        ref_pairs = create_pop_references(references, snps, positions_rel, use_chrom, ploidy)
        refdata, _ignore = zip(*ref_pairs)

        data_out = refdata
        positions_out = positions_abs  # crucial: keep absolute for linkage/output
        sites = len(positions_out)
        return populations, data_out, positions_out, sites

    # SNP-only or abbrevref: just return snps
    data_out = snps
    positions_out = positions_abs
    sites = len(positions_out)
    return populations, data_out, positions_out, sites



def filter_vcf_by_bound(vcf, use_chrom, bound):
    if bound is None:
        return vcf
    start, stop = bound
    if start > stop:
        start, stop = stop, start
    use_set = set(use_chrom)
    out = []
    for row in vcf:
        if row[0] not in use_set:
            continue
        if start <= row[1] <= stop:
            out.append(row)
    return out


def make_positions_rel(positions_abs, base0_by_chrom):
    return [(ci, pos - base0_by_chrom[ci]) for ci, pos in positions_abs]

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
def create_pop_snps(vcf, begin, stop, use_chrom, ploidy):
    idata = []
    positions = []

    chrom_to_idx = {c: i for i, c in enumerate(use_chrom)}
    DELIM0, DELIM1 = DELIM
    local_allowindel = allowindel

    for v in vcf:
        ci = chrom_to_idx.get(v[0])
        if ci is None:
            continue

        pos = v[1]
        positions.append((ci, pos))

        ref = v[3]
        alt = v[4].split(',')

        if local_allowindel:
            lr = len(ref)
            if lr > 1:
                for i in range(len(alt)):
                    alt[i] = alt[i].ljust(lr, '-')
            else:
                maxla = max(len(a) for a in alt)
                if maxla > 1:
                    ref = '@' + ref.ljust(maxla, '-')
                    for i in range(len(alt)):
                        alt[i] = '@' + alt[i].ljust(maxla, '-')

        popinds = v[9][begin:stop]
        sdata = []
        sdata_append = sdata.append
        alt_get = alt.__getitem__

        for raw in popinds:
            gt, _, _ = raw.partition(':')
            if DELIM0 in gt:
                parts = gt.split(DELIM0)
            elif DELIM1 in gt:
                parts = gt.split(DELIM1)
            else:
                parts = (gt,)

            for pi in range(ploidy):
                x = parts[pi] if pi < len(parts) else '.'
                if x == '0':
                    sdata_append(ref)
                elif x != '.':
                    sdata_append(alt_get(int(x) - 1))
                else:
                    sdata_append('?')

        idata.append(sdata)

    return idata, positions
                                                                                    
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
    ref_bytes = [r.encode('ascii') for r in references]

    allpop = []
    for pop_snps in snps:
        poploci = []
        z = -1

        for chrom, pos_list in chrompos:
            haplotypes = None
            insertionmuts = []

            for po in pos_list:
                z += 1
                row = pop_snps[z]

                if haplotypes is None:
                    base = ref_bytes[chrom]
                    haplotypes = [bytearray(base) for _ in range(len(row))]

                for hi, mut in enumerate(row):
                    if not mut:
                        continue
                    if mut[0] == '@':
                        haplotypes[hi][po:po+1] = b'@'
                        insertionmuts.append((hi, mut[1:].encode('ascii')))
                    else:
                        mb = mut.encode('ascii')
                        haplotypes[hi][po:po+len(mb)] = mb

            fix_haplotypes_bytes(haplotypes, insertionmuts)
            poploci.append([h.decode('ascii') for h in haplotypes])

        allpop.append((poploci, positions))
    return allpop

def fix_haplotypes_bytes(haplotypes, insertionmuts):
    for hi, ins in insertionmuts:
        h = haplotypes[hi]
        while True:
            idx = h.find(b'@')
            if idx == -1:
                break
            h[idx:idx+1] = ins

            
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
def write_migrate(migratefile, data, freqs, positions, sites, references, names, comment, bound=None):
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
            nucs = []
            delta = linkedsnps
            chrompos = convert_chrompos(positions)  # (chrom_index, [POS...])

            # bound-fixed origin (inclusive)
            if bound is not None:
                b0, b1 = bound
                if b0 > b1:
                    b0, b1 = b1, b0
            else:
                b0 = None
                b1 = None

            for chrom_index, poslist in chrompos:
                if not poslist:
                    continue
                # unique SNP columns
                posu = sorted(set(poslist))

                groups = group_snp_positions(posu, delta=linkedsnps, bound=bound, chromlink=False)
                nucs.extend([len(g[2]) for g in groups])

                # restrict to bound if present
                #if b0 is not None:
                #    posu = [p for p in posu if b0 <= p <= b1]
                #    span = (b1 - b0 + 1)
                #    origin = b0
                #else:
                 #   # no bound: define bins from first position (consistent within chrom)
                 #   origin = posu[0]
                 #   span = posu[-1] - origin + 1

                #if span <= 0:
                #    continue

                #k = (span + delta - 1) // delta  # number of bins
                #counts = [0] * k

                #for p in posu:
                #    gi = (p - origin) // delta
                #    if 0 <= gi < k:
                #        counts[gi] += 1

                #nucs.extend(counts)

            loci = len(nucs)
            sitestr = " ".join([f'(n{si})' for si in nucs])

            #@@@            
            f.write(f'{numpop} {loci} {vcffile}\n')
            f.write(f"# VCF file used:      {vcffile}\n")
            f.write(f"# Translated from VCF {dt.date.today()}\n")
            f.write(f"# Migrate input file: {migratefile}\n")
            f.write(f"# SNP data file!\n")
            f.write(f"# {comment}\n")
            f.write(f"{sitestr}\n")
    # write out the frequencies if present
    if freqs is not None:
        freqs_out = freqs

        # Only do the old estimate-based expansion if we do NOT already have exact per-group freqs
        if (refabbrev and (linkedsnps not in (None, CHROMLINK))
            and not _freqs_already_grouped(freqs)):
        # ... your existing split_freq_by_lengths expansion that produced 1 0.1, 1 0.2, ...
            freqs_out = expanded

        fistr = check_freqsout(freqs_out)
        for fi, fis in enumerate(fistr):
            f.write(f"#*freq: {fi+1} {fis}\n")
        #print("DEBUG write_migrate freqs len [freqs]", len(freqs) if freqs else None,
        #      "Locs", [fi["Loc"] for fi in freqs[:3]],
        #      "A", [fi.get("A") for fi in freqs[:3]],
        #      file=sys.stderr)
        

    #if freqs != None:
    #    fistr = check_freqsout(freqs)
    #    for fi, fis in enumerate(fistr):
    #        f.write(f"#*freq: {fi+1} {fis}\n")
    
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

def split_freq_by_lengths(freq, lengths, loc_prefix):
    """
    Create one freq dict per segment length, using the A/C/G/T proportions
    from `freq` and allocating expected counts by segment length.
    `lengths` are in BASES (not SNP counts).
    """
    total = sum(freq.get(b, 0) for b in "ACGT")
    if total <= 0:
        p = {b: 0.25 for b in "ACGT"}
    else:
        p = {b: freq.get(b, 0) / total for b in "ACGT"}

    out = []
    for si, L in enumerate(lengths, start=1):
        # expected counts with rounding; then fix drift to sum exactly to L
        counts = {b: int(round(p[b] * L)) for b in "ACGT"}
        diff = L - sum(counts.values())

        # deterministic drift fix: add/subtract following highest-prob bases
        bases = sorted("ACGT", key=lambda b: (-p[b], b))
        j = 0
        while diff != 0:
            b = bases[j % 4]
            if diff > 0:
                counts[b] += 1
                diff -= 1
            else:
                if counts[b] > 0:
                    counts[b] -= 1
                    diff += 1
            j += 1

        fi = {"Loc": f"{loc_prefix}.{si}"}
        fi.update(counts)
        fi["?"] = 0
        out.append(fi)
    return out

    
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
        fistr = f"{fi['Loc']} ACGT={total} All={alltotal} A={fi['A']} C={fi['C']} G={fi['G']} T={fi['T']} ?={fi['?']}" 
        fis.append(fistr)
    return fis

def compute_invariant_group_freqs_bound_fixed(
    references,          # list of reference strings in use_chrom order (or mapped order you pass in)
    positions,           # list of (chrom_index, POS_abs) for VCF records (VCF POS is 1-based)
    use_chrom,           # list of chrom names in the same order as references
    linkedsnps,          # None, CHROMLINK, or integer delta
    bound=None,          # (start, stop) inclusive, POS coordinates
):
    """
    Emits freq dicts for migrate:

      - one 'chrom/window total' line: Loc="<chrom> 0"
      - then one 'invariant-by-group' line per locus: Loc="<chrom> <g>"

    Semantics:
      * If bound is not None and linkedsnps is int delta:
            groups are base-fixed bins starting at bound_start.
      * If bound is None and linkedsnps is int delta:
            groups use the SAME rolling/SNP-anchored rule as your (n...) construction:
              start at first SNP; group covers [start, start+delta-1];
              next group starts at first SNP > previous_end.
      * If linkedsnps == CHROMLINK:
            one group per contig (within bound if present, else whole contig).
      * If linkedsnps is None:
            emit only Loc "<chrom> 0" totals (no per-group lines).

    Invariant counts exclude variable sites (unique POS).
    """

    # Normalize bound
    if bound is not None:
        b0, b1 = bound
        if b0 > b1:
            b0, b1 = b1, b0
        bound_len = b1 - b0 + 1
    else:
        b0, b1 = None, None
        bound_len = None

    # Build unique SNP positions per chromosome index
    varpos_by_chrom = {}
    for ci, pos in positions:
        varpos_by_chrom.setdefault(ci, set()).add(pos)

    out = []

    for ci, refseq in enumerate(references):
        chromname = use_chrom[ci] if use_chrom else str(ci)
        varpos = varpos_by_chrom.get(ci, set())
        if not varpos:
            # no SNPs => no loci in (n...) either; skip for consistency
            continue

        ref_len = len(refseq)

        # Decide contig window for CHROMLINK or totals when bound given
        if bound is not None:
            window_start_pos = b0
            window_end_pos   = b1
            # if already cropped to bound, use directly
            if bound_len is not None and ref_len == bound_len:
                base_ref = refseq
                base_ref_start_pos = b0  # POS of base_ref[0]
            else:
                # slice [b0..b1] from full contig
                s = max(0, b0 - 1)
                e = min(ref_len, b1)  # exclusive
                base_ref = refseq[s:e]
                base_ref_start_pos = b0
        else:
            # no bound: for CHROMLINK we use full contig
            base_ref = refseq
            base_ref_start_pos = 1
            window_start_pos = 1
            window_end_pos   = ref_len

        # Helper to get a reference slice for POS interval [p0..p1] inclusive
        def _slice_ref(p0, p1):
            # clamp to contig bounds
            if p0 < window_start_pos:
                p0 = window_start_pos
            if p1 > window_end_pos:
                p1 = window_end_pos
            if p1 < p0:
                return ""
            # map POS to indices in base_ref
            i0 = p0 - base_ref_start_pos
            i1 = p1 - base_ref_start_pos + 1  # exclusive
            if i0 < 0:
                i0 = 0
            if i1 > len(base_ref):
                i1 = len(base_ref)
            return base_ref[i0:i1]

        # Build groups (list of (g_start_pos, g_end_pos, snp_positions_in_group_set))
        groups = []

        if linkedsnps is None:
            # only totals later
            pass

        elif linkedsnps == CHROMLINK:
            # one group per contig window
            # include only SNPs within window (bound if present)
            snps_in = {p for p in varpos if window_start_pos <= p <= window_end_pos}
            groups = [(window_start_pos, window_end_pos, snps_in)]

        elif isinstance(linkedsnps, int) and linkedsnps > 0:
            delta = linkedsnps

            if bound is not None:
                # base-fixed bins starting at b0
                span = window_end_pos - window_start_pos + 1
                k = (span + delta - 1) // delta
                for gi in range(k):
                    g0 = window_start_pos + gi * delta
                    g1 = min(window_end_pos, g0 + delta - 1)
                    snps_in = {p for p in varpos if g0 <= p <= g1}
                    # IMPORTANT: include empty groups? In migrate header you include them, yes (n0).
                    # But your current (n...) rolling code does NOT. For bound-fixed you DO include all bins.
                    groups.append((g0, g1, snps_in))

            else:
                # NO bound: use rolling/SNP-anchored groups to MATCH your (n...) logic
                posu = sorted(varpos)
                i = 0
                while i < len(posu):
                    g0 = posu[i]
                    g1 = min(ref_len, g0 + delta - 1)
                    # collect SNPs within [g0..g1]
                    j = i
                    snps = set()
                    while j < len(posu) and posu[j] <= g1:
                        snps.add(posu[j])
                        j += 1
                    groups.append((g0, g1, snps))
                    i = j

        else:
            # invalid linkedsnps
            continue

        # Decide "chrom total" window for Loc "<chrom> 0"
        # If we have groups, total should cover the union span from first group start to last group end.
        if groups:
            tot0 = groups[0][0]
            tot1 = groups[-1][1]
        else:
            tot0 = window_start_pos
            tot1 = window_end_pos

        total_seq = _slice_ref(tot0, tot1)
        tot = Counter(total_seq.upper())
        chrom_tot = {
            "Loc": f"{chromname} 0",
            "A": tot.get("A", 0),
            "C": tot.get("C", 0),
            "G": tot.get("G", 0),
            "T": tot.get("T", 0),
            "?": tot.get("?", 0) + tot.get("N", 0) + tot.get("X", 0),
        }
        out.append(chrom_tot)

        # If unlinked, stop after totals
        if linkedsnps is None:
            continue

        # Now compute per-group invariant freqs (ref composition in group window minus SNP sites)
        for gi, (g0, g1, snps_in) in enumerate(groups, start=1):
            seq = _slice_ref(g0, g1)
            c = Counter(seq.upper())
            gfreq = {
                "A": c.get("A", 0),
                "C": c.get("C", 0),
                "G": c.get("G", 0),
                "T": c.get("T", 0),
                "?": c.get("?", 0) + c.get("N", 0) + c.get("X", 0),
            }

            # subtract each SNP position once (unique POS) within group
            for p in snps_in:
                if p < g0 or p > g1:
                    continue
                idx = p - g0  # 0-based within seq
                if 0 <= idx < len(seq):
                    b = seq[idx].upper()
                    if b in "ACGT":
                        gfreq[b] -= 1
                    else:
                        gfreq["?"] -= 1

            # clamp
            for b in "ACGT?":
                if gfreq[b] < 0:
                    gfreq[b] = 0

            out.append({
                "Loc": f"{chromname} {gi}",
                "A": gfreq["A"],
                "C": gfreq["C"],
                "G": gfreq["G"],
                "T": gfreq["T"],
                "?": gfreq["?"],
            })

    return out

def compute_invariant_group_freqs_bound_fixed_outdated(
    references,          # list of reference strings in use_chrom order (or mapped order you pass in)
    positions,           # list of (chrom_index, POS_abs) for VCF records (VCF POS is 1-based)
    use_chrom,           # list of chrom names in the same order as references
    linkedsnps,          # None, CHROMLINK, or integer delta
    bound=None,          # (start, stop) inclusive, POS coordinates
):
    """
    Returns a list of freq dicts suitable for check_freqsout(), containing:
      - one 'chrom total' line per chromosome window: Loc="<chrom> 0"
      - then one 'invariant-by-group' line per linkage group: Loc="<chrom> <g>"

    Grouping semantics:
      - If bound is given: groups are bound-fixed bins starting at bound_start.
      - If bound is None: groups are SNP-span-fixed bins starting at min SNP POS on that contig.
      - If linkedsnps == CHROMLINK: one group per contig window.
      - If linkedsnps is None: only emits chrom totals (Loc "<chrom> 0") and no per-group lines.

    Invariant counts exclude variable sites (unique POS per chromosome).
    """

    # Normalize bound
    if bound is not None:
        b0, b1 = bound
        if b0 > b1:
            b0, b1 = b1, b0
        bound_len = b1 - b0 + 1
    else:
        b0, b1 = None, None
        bound_len = None

    # Build unique SNP position sets per chromosome index
    varpos_by_chrom = {}
    for ci, pos in positions:
        varpos_by_chrom.setdefault(ci, set()).add(pos)

    out = []

    for ci, refseq in enumerate(references):
        chromname = use_chrom[ci] if use_chrom else str(ci)

        varpos = varpos_by_chrom.get(ci, set())

        # Skip contigs with no SNPs in positions[] (keeps freqs consistent with loci list)
        # If you prefer to still emit "<chrom> 0" totals for all references, remove this continue,
        # but then your output can contain freqs for contigs that never appear as loci.
        if not varpos:
            continue

        # Decide window and origin
        if bound is not None:
            # bound-fixed window
            window_start_pos = b0
            window_end_pos   = b1
            origin = b0

            # If references are already cropped to bound, use them directly.
            # Otherwise slice the full contig.
            if bound_len is not None and len(refseq) == bound_len:
                window = refseq
            else:
                start_idx = window_start_pos - 1
                stop_idx  = window_end_pos      # exclusive
                if start_idx < 0:
                    start_idx = 0
                if stop_idx > len(refseq):
                    stop_idx = len(refseq)
                window = refseq[start_idx:stop_idx]
        else:
            # NO bound: use SNP-span window so group counts match loci implied by SNP span.
            window_start_pos = min(varpos)
            window_end_pos   = max(varpos)
            origin = window_start_pos

            start_idx = window_start_pos - 1
            stop_idx  = window_end_pos      # exclusive
            if start_idx < 0:
                start_idx = 0
            if stop_idx > len(refseq):
                stop_idx = len(refseq)
            window = refseq[start_idx:stop_idx]

        win_len = len(window)
        if win_len <= 0:
            continue

        # ---- (A) Chromosome total composition for the window
        tot = Counter(window.upper())
        chrom_tot = {
            "Loc": f"{chromname} 0",
            "A": tot.get("A", 0),
            "C": tot.get("C", 0),
            "G": tot.get("G", 0),
            "T": tot.get("T", 0),
            "?": tot.get("?", 0) + tot.get("N", 0) + tot.get("X", 0),
        }
        out.append(chrom_tot)

        # If unlinked SNPs, we emit only chrom total (sanity) and stop here
        if linkedsnps is None:
            continue

        # Decide group lengths in this window
        if linkedsnps == CHROMLINK:
            group_lengths = [win_len]
        elif isinstance(linkedsnps, int) and linkedsnps > 0:
            delta = linkedsnps
            span = win_len
            k = (span + delta - 1) // delta
            if k <= 0:
                continue
            group_lengths = [delta] * (k - 1) + [span - delta * (k - 1)]
        else:
            # unknown/invalid linkedsnps: no groups
            continue

        # Precompute per-group totals by scanning window in group-sized chunks
        group_totals = []
        offset = 0
        for L in group_lengths:
            chunk = window[offset:offset+L]
            c = Counter(chunk.upper())
            group_totals.append({
                "A": c.get("A", 0),
                "C": c.get("C", 0),
                "G": c.get("G", 0),
                "T": c.get("T", 0),
                "?": c.get("?", 0) + c.get("N", 0) + c.get("X", 0),
                "_len": L,
            })
            offset += L

        # Optional debug (safe even for single-group)
        #if group_totals:
        #    print("DEBUG", chromname, "win_len", win_len,
        #          "groups", len(group_totals),
        #          "first group ACGT?", {b: group_totals[0][b] for b in "ACGT?"},
        #          "last group ACGT?", {b: group_totals[-1][b] for b in "ACGT?"},
        #          file=sys.stderr)

        # Subtract variable sites: remove the reference base at each unique SNP position
        if linkedsnps == CHROMLINK:
            # all SNPs fall into group 1 (index 0) within this window
            for pos in varpos:
                if pos < window_start_pos or pos > window_end_pos:
                    continue
                idx = pos - origin  # 0-based into window
                if 0 <= idx < win_len:
                    b = window[idx].upper()
                    if b in "ACGT":
                        group_totals[0][b] -= 1
                    else:
                        group_totals[0]["?"] -= 1
        else:
            # delta-linked
            delta = linkedsnps
            for pos in varpos:
                if pos < window_start_pos or pos > window_end_pos:
                    continue
                gi = (pos - origin) // delta
                idx = pos - origin
                if 0 <= gi < len(group_totals) and 0 <= idx < win_len:
                    b = window[idx].upper()
                    if b in "ACGT":
                        group_totals[gi][b] -= 1
                    else:
                        group_totals[gi]["?"] -= 1

        # Clamp after subtraction
        for g in group_totals:
            for b in "ACGT?":
                if g[b] < 0:
                    g[b] = 0

        #if group_totals:
        #    print("DEBUG emit first/last",
        #          group_totals[0]["A"], group_totals[-1]["A"],
        #          file=sys.stderr)

        # Emit per-group invariant freqs
        for gi, g in enumerate(group_totals, start=1):
            out.append({
                "Loc": f"{chromname} {gi}",
                "A": g["A"],
                "C": g["C"],
                "G": g["G"],
                "T": g["T"],
                "?": g["?"],
            })
    return out


def _refseq_list_from_references(references):
    if references is None:
        return None
    if isinstance(references, (list, tuple)) and references and isinstance(references[0], str):
        return list(references)
    if isinstance(references, list) and references and isinstance(references[0], (list, tuple)):
        if len(references[0]) >= 2 and isinstance(references[0][1], str):
            return [r[1] for r in references]
    if isinstance(references, (list, tuple)) and len(references) >= 2:
        seqs = references[1]
        if isinstance(seqs, (list, tuple)) and seqs and isinstance(seqs[0], str):
            return list(seqs)
    raise TypeError(f"Unrecognized references structure: {type(references)}")

def _freqs_already_grouped(freqs):
    # True if we have entries like "chrom group" (e.g., "1 1") besides the "chrom 0" total
    # i.e. second token is an integer
    try:
        for fi in freqs:
            loc = str(fi.get("Loc", "")).strip().split()
            if len(loc) >= 2 and loc[1].isdigit() and int(loc[1]) >= 1:
                return True
    except Exception:
        pass
    return False

def _reference_dict_from_headers_and_seqs(refheaders, refseqs):
    """
    Build a dict {chromname: seq} using the first token after '>' in each FASTA header.
    Also stores a cleaned variant with ':' and '=' removed to mirror your earlier cleanup.
    """
    d = {}
    for h, s in zip(refheaders, refseqs):
        name = h.strip()
        if name.startswith(">"):
            name = name[1:].strip()
        name = name.split()[0] if name else ""
        if not name:
            continue
        d[name] = s
        d[name.replace(":", "").replace("=", "")] = s
    return d

def group_snp_positions(posu, *, delta, bound=None, chromlink=False):
    """
    posu: sorted unique POS (1-based)
    Returns list of groups, each group is (g0, g1, snps_set)
    """
    if not posu:
        return []

    if chromlink:
        if bound is not None:
            b0, b1 = bound
            if b0 > b1: b0, b1 = b1, b0
            snps = {p for p in posu if b0 <= p <= b1}
            return [(b0, b1, snps)]
        else:
            return [(posu[0], posu[-1], set(posu))]

    if bound is not None:
        b0, b1 = bound
        if b0 > b1: b0, b1 = b1, b0
        span = b1 - b0 + 1
        k = (span + delta - 1) // delta
        groups = []
        for gi in range(k):
            g0 = b0 + gi * delta
            g1 = min(b1, g0 + delta - 1)
            snps = {p for p in posu if g0 <= p <= g1}
            groups.append((g0, g1, snps))
        return groups

    # NO bound: rolling groups
    groups = []
    i = 0
    while i < len(posu):
        g0 = posu[i]
        g1 = g0 + delta - 1
        j = i
        snps = set()
        while j < len(posu) and posu[j] <= g1:
            snps.add(posu[j])
            j += 1
        groups.append((g0, g1, snps))
        i = j
    return groups


if __name__ == "__main__":
    numpop = -1
    vcffile, referencefile, linkedsnps, numind, numloc, migratefile, use_chrom, allowindel, strictvcf, refabbrev, bound = parse_args(sys.argv)
    print("Parsed options:")

    vcf,names, chroms, ploidy  = read_vcf(vcffile)
    print(f"VCF file used: {vcffile}")

    start = 0
    populations=[]
    data = [] 
    freqs = None
    
    if referencefile != None:
        references = read_reference(referencefile)
        for ref in referencefile.split(','):
            print(f"Reference file: {ref}")
        refheaders, references, sites, freqs = list(zip(*references))
        bsnps=False
        #print("@@", freqs)
    else: # we only report snps and if linkedsnps !=None then the snps are linked 
        references = None
        refheaders = None
        sites = linkedsnps
        bsnps=True
        #if freqs!=None:
        #freqs=None

    # after reading vcf/names/chroms/ploidy and (maybe) references/freqs:

    use_chrom = harmonize_use_chroms(use_chrom, chroms)

    # Prepare inputs (bound filtering + validation + reference cropping)
    vcf, references, freqs2, offsets = prepare_inputs(
        vcf=vcf,
        use_chrom=use_chrom,
        references=references,
        freqs=freqs,
        bound=bound,
        refabbrev=refabbrev,
        strictvcf=strictvcf
    )
    
    # Preserve your original abbrevref freq cleanup if bound is None (unchanged behavior)
    if references is not None and refabbrev:
        if bound is None:
            # your old behavior:
            freqs = list(freqs2)
            for fi in freqs:
                item = fi['Loc']
                fi['Loc'] = item.replace(':','').replace('=','')
            freqs = [fi for fi in freqs if fi['Loc'] in use_chrom]
        else:
            # recomputed freqs are already aligned with use_chrom
            freqs = freqs2
    else:
        # full --ref mode or SNP-only mode
        if refabbrev:
            freqs = freqs2
        else:
            freqs = None

    # Build dataset (snps or refdata) + absolute positions for linkage/writeout
    populations, data_out, positions_out, sites = build_dataset(
        vcf=vcf,
        names=names,
        numind=numind,
        use_chrom=use_chrom,
        ploidy=ploidy,
        references=references,
        refabbrev=refabbrev,
        bound=bound,
        offsets=offsets
    )

    # Decide comment exactly as before, but you can append bound info if you want
    if references is not None and refabbrev is False:
        comment = "Using references augmented VCF data"
    else:
        if linkedsnps == CHROMLINK:
            comment = "Using all SNPs, linked by chromosome"
        elif linkedsnps is None:
            comment = "Using unlinked SNPS"
        else:
            comment = f"Using all SNPs, linked every {linkedsnps} for each chromosome"

    if bound is not None:
        b0, b1 = bound
        if b0 > b1:
            b0, b1 = b1, b0
        comment += f" with bound [{b0},{b1}]"


    if refabbrev:
        # references currently is a sequence list (or tuple) of contigs, not necessarily filtered to use_chrom
        refseqs_all = _refseq_list_from_references(references)  # returns list of seq strings
        refdict = _reference_dict_from_headers_and_seqs(refheaders, refseqs_all)

        # now build list in use_chrom order; raise a helpful error if missing
        refseqs = []
        missing = []
        for ch in use_chrom:
            if ch in refdict:
                refseqs.append(refdict[ch])
            else:
                missing.append(ch)
        if missing:
            raise ValueError(f"Missing reference contigs for: {missing[:5]}{'...' if len(missing)>5 else ''}")
        
        freqs = compute_invariant_group_freqs_bound_fixed(
            references=refseqs,          # list of strings in use_chrom order
            positions=positions_out,            # list of (chrom_index, POS_abs)
            use_chrom=use_chrom,
            linkedsnps=linkedsnps,
            bound=bound,                    # (start, stop) or None
        )
        #print("DEBUG freqs lines:", len(freqs), "first Loc:", freqs[0]["Loc"], "last Loc:", freqs[-1]["Loc"], file=sys.stderr)

        
    print(comment)

    # Final write: always pass absolute positions_out
    write_migrate(migratefile, data_out, freqs, positions_out, sites, references, populations, comment, bound=bound)
    print(f"Migrate input file: {migratefile}")

