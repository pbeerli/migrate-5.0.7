#!/usr/bin/env python
#
# converts an aligned stacks .fa file to the migrate format
# reads should be labeled like this
# >CLocus_12706_Sample_1_Locus_34105_Allele_0 [BayOfIslands_s088.fq; groupI, 125578, +]
#    CLocus_#1                      Allele_#2  Location_#3 Sample_#4
#
# options are
# -i inputfile
# -o outputfile
# -r     remove trailing N
# Default will check list for keywords Clocus and then pick next element
#                                      Allele and then pick next element split and take first
#                                      assume that the next element with [ is a location
#                                      element after allele 0/1 label with fastq filenname with sample name
#
# We assume that the fasta file has a standard string from Stacks with the addition that the
# filename must contain the location and sample identification, for example Location_Sample 
# like this BayOfIslands_s088
# the script will parse a string like this:
# >CLocus_12706_Sample_1_Locus_34105_Allele_0 [BayOfIslands_s088.fq; groupI, 125578, +]  
# (1) read a the line as string, discard the newline at the end
# (2) split the string on the character '_', this leaves a list:
# h=['>CLocus','12706','Sample','1','Locus','34105','Allele','0 [BayOfIslands', 's088.fq; groupI, 125578, +]']  
# (3) using index() we find Clocus and Allele an pick the next element to fill locus
# (4) the string after allele is then split on the empty () returning
#  temp = ['0', '[BayOfIslands'], these will be stored as allele and location
# (5) the sample will be extracted from the list element allele_index+2 and split on the '.' 
# (6) the next line is the complete sequence for the header=[locus,location,sample,allele]
#
# Distributed under the MIT open source license
# (c) 2015 Peter Beerli, Tallahassee
# (c) 2025 Peter Beerli, Tallahassee: Updated to work with Python 3, and fix for slight format variants
#                                     in the fasta file
###################################################################################################### 
import sys

def read_lines(f,headers,sequences):
    while True:
        header = f.readline()
        if not header: 
            break
        sequence = f.readline()
        if not sequence: 
            break
        h = header.rstrip().split('_')
        index = h.index('>CLocus')
        locus = h[index+1]
        index = h.index('Allele')
        temp = h[index+1].split()
        allele, location = temp[0],temp[1][1:]
        temp = h[index+2].split('.')
        sample = temp[0].replace("]","")
        headers.append([locus,location,sample,allele])
        sequences.append(sequence.rstrip())
        #print(locus,location,sample,allele,end=' ')
        #print(sequence)

 
def read_stacks(file):
    headers = []
    sequences = []
    if file != sys.stdin:
        with open(file,'r') as f:
            read_lines(f,headers,sequences)
    else:
        read_lines(sys.stdin,headers,sequences)
    return headers,sequences

def column(matrix, i):
    return [row[i] for row in matrix]

def find_unique_locations(locations):
    return list(set(locations))

def group(headers,sequences, col1, locations, col2, loci):
    data = zip(headers,sequences)
    populations = []
    for l in locations:
        pop = list(filter(lambda x: x[0][col1] == l , data))
        poploci = []
        for lo in loci:
            ll = list(filter(lambda x: x[0][col2] == lo , pop))
            poploci.append(ll)
            #print "Pop:", l, "Locus:",lo, ll[0][0][0],len(ll)
        populations.append(poploci)
    return populations

def findN(x):
    try:
        return list(x[1]).index('N')
    except:
        return len(x[1])
        

def find_trailing_Ns(populations):
    trailingNs = [0]*len(populations[0])
    for pop in populations:
        for locus,i in zip(pop,range(len(pop))):
            a = max(list(map(findN,locus)))
            if a>trailingNs[i]:
                trailingNs[i] = a
    return trailingNs

def to_migrate(populations,locations,filename):
    loc = iter(locations)
    output =[]
    output.append("%li %li data: %s" % (len(populations), len(populations[0]),filename))
    output.append("# Automatically translated by stacks2mig.py")
    output.append("# Converter can be found: http://www.peterbeerli.com/software")
    output.append("# MIT license, 2015, Peter Beerli")
    if trailNtrue:
        lastNs = find_trailing_Ns(populations)
        sites = list(map(str,lastNs))
        output.append("%s" % " ".join(sites))
    else:
        siteslist = populations[0]
        sites = []
        for sit in siteslist:
            print(sit)
            sites.append(str(len(sit[0][1])))
        output.append("%s" % " ".join(sites))

    for pop in populations:
        samples = list(map(len,pop))
        samples = list(map(str,samples))
        output.append("%s %s" % ( " ".join(samples),  next(loc))) 
        if trailNtrue:
            tN = iter(lastNs)
        for locus in pop:
            if trailNtrue:
                bmax= next(tN)
            #output.append("# %s" % locus[0][0][0])
            #output.append("# %s" % locus)
            for si in locus:
                sss = ("%-.8s" %si[0][-2])+":"+si[0][-1]
                if trailNtrue:
                    output.append("%-10.10s %s" % (sss, si[1][:bmax]))
                else:
                    output.append("%-10.10s %s" % (sss,si[1]))
    return output

if __name__ == '__main__':
    options = sys.argv[1:]
    infilename = sys.stdin
    outfilename = sys.stdout
    trailNtrue = False
    reorder = False
    if "-h" in options:
        print ("stacks2mig.py")
        print ("syntax: python stacks2mig.py < -i stacksfile.fa >  <-o migratefile> <-r>")
        print ("<> means options")
        print (" -i filename # expects a stacks fasta file, if ommited reads from standard input")
        print (" -o filename # outputfilename to write the migrate infile to, if omitted writes to standard out")
        print (" -r          # removes trailing \'N\', if omitted does not remove trailing \'N\'")
        print (" -l          # interactively asks about reordering the locations in the outputfile")
        print (" The title for each sequence must look like this:")
        print (" >CLocus_12706_Sample_1_Locus_34105_Allele_0 [BayOfIslands_s088.fq; groupI, 125578, +]")
        print ("         xxxxx                             x  xxxxxxxxxxxx xxxx")
        print (" where the xxx are extracted for the conversion")
        print (" The sequence must be on one line!")
        print ("")
        print (" MIT license, Peter Beerli 2025")
        print()
        sys.exit(-1)
    if "-i" in options:
        index = options.index("-i") 
        infilename = options[index+1]
    if "-o" in options:
        index = options.index("-o") 
        outfilename = options[index+1]
    if "-r" in options:
        trailNtrue = True
    if "-l" in options:
        reorder = True
    headers,sequences = read_stacks(infilename)
    locations = find_unique_locations(column(headers,1))
    loci = find_unique_locations(column(headers,0))
    if reorder:
        for i in range(1,len(locations)+1):
            print (i, locations[i-1])
        newlocations = [int(i) for i in input("Enter the numbers in the order you like:\n").split()]
        if min(newlocations)==0:
            locations = [locations[i] for i in newlocations]
        if min(newlocations)==1:
            locations = [locations[i-1] for i in newlocations]
    print (locations)
    populations = group(headers,sequences, 1, locations, 0, loci)
    output = to_migrate(populations,locations,infilename)
    if outfilename != sys.stdout:
        of = open(outfilename,'w')
        for oi in output:
            of.write(oi+'\n')
    else:
        for oi in output:
            print( oi)

