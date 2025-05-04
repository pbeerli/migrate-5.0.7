#!/usr/bin/env python
#
import sys
import os
infile = sys.argv[1]
def read_migrate_infile(filename):
    data = []
    f = open(filename,'r')
    first = True
    second = False
    while True:
        fi = f.readline()
        print fi, first,second 
        if fi[0]=='#':
            continue
        if first:
            first=False
            second=True
            line = fi.rstrip().split()
            numpop = int(line[0])
            numloci = int(line[1])
            data = [[] for locus in range(numloci)]
            for locus in range(numloci):
                popblock = [[] for pop in range(numpop)]
                data[locus] = popblock[:]
            try:
                title = line[2:]
            except:
                title=""
            continue
        if second:
            second = False
            sites = [int(xi) for xi in fi.rstrip().split()]
            break
        if second == False and first == False:
            break
    #print data
    #sys.exit()
    for pop in range(numpop):
        while True:
            try:
                fi = f.readline()
            except:
                break
            #if pop==numpop-1 and locusdone:
            #    break
            if fi[0]=='#':
                continue
            line = fi.rstrip().split()
            print line
            if len(line)>numloci:
                indarray = [int(xi) for xi in line[:numloci]]
            else:
                indarray = [int(line[0]) for xi in range(numpop)]
            print indarray
            print numloci
            print zip(indarray,range(numloci))
            locusdone=False
            for indpop,locus in zip(indarray,range(numloci)):
                
                #s = data[locus][pop]
                for ind in range(indpop):
                    while True:
                        fi = f.readline()
                        if fi[0] != '#':
                            break 
                    line = fi.rstrip()
                    #print line
                    indname = line[:10]
                    if '\t' in indname:
                        indname = indname.split('\t')
                    print "@",indname,"@",locus
                    data[locus][pop].append(line[10:])
            break
    return data

data = read_migrate_infile(infile)
locus = 0
for di in data:
    # locus
    locus += 1
    # find number of alleles:
    alleles = []
    for d in di:
        alleles.extend(list(set(d)))
    alleles = list(set(alleles))
    print locus, len(alleles)

