#!/usr/bin/env python
# convert structure data into migrate data
# reads
#(1) Individual code number --- for HGDP-CEPH individuals, these code
#    numbers match the numbers previously used.
#(2) Population code number assigned by us.
#(3) Population name.
#(4) Geographic information about the population.
#(5) Geographic region.
import sys
#
def read(file):
    f = open(file,'r')
    locinames = f.readline().split()
    data=[]
    for i in f:
        data.append(i.split())
    return [locinames,data]


if __name__ == '__main__':
    
    if len(sys.argv)<2:
            print "Syntax: python structure2mig.py filename codenum"
            print "filename is the file containing data formatted for structure[? aka Wang and Rosenberg 2007]"
            print "codenum 1 population code"
            print "        2 population name"
            print "        3 geographic region name"
            print "        4 geographic code"
            sys.exit(-1)
    else:
        filename = sys.argv[1]
        codenum = int(sys.argv[2])
        
    [locinames, data] = read(filename)
    indcode = [x[0] for x in data]
    popcode = [x[1] for x in data]
    popname = [x[2] for x in data]
    geoname = [x[3] for x in data]
    geocode = [x[4] for x in data]
    raw     = [x[5:] for x in data]
    if codenum==1:
        popname2 = [popcode[i] for i in xrange(0,len(popcode),2)]
    if codenum==2:
        popname2 = [popname[i] for i in xrange(0,len(popname),2)]
    if codenum==3:
        popname2 = [geoname[i] for i in xrange(0,len(geoname),2)]
    if codenum==4:
        popname2 = [geocode[i] for i in xrange(0,len(geocode),2)]
        
    indcode2 = [indcode[i] for i in xrange(0,len(indcode),2)]
    left = [raw[i] for i in xrange(0,len(raw),2)]
    right = [raw[i] for i in xrange(1,len(raw),2)]
    #
    sites=[]
    for i in range(0,len(left[0])):
        a = [x[i] for x in left]
        b = [x[i] for x in right]
        all = [float(x) for x in list(set(a+b))]
        all2 = sorted(all)
        all2.pop(0)
        oldj = all2[0]
        sum = 0.0
        for j in all2:
            sum += j-oldj
            oldj = j
        site = str(int(round(sum/(len(all2)-1))))+' '
        sites.append(site)
    thesites =  "".join(sites)
    #
    #    
    #
    a = zip(popname2,left,right)
    a.sort()
    [popname2,left,right] = zip(*a)
    #
    print len(list(set(popname2))), len(locinames), "/", "wang and rosenberg"
    print '#@M',thesites
    for ind in range(0,len(left)):
        if popname2[ind]!=popname2[ind-1]:
            print popname2.count(popname2[ind]), popname2[ind]
        print "%-10.10s" % indcode2[ind],
        for loc in range(0,len(left[0])):
            if left[ind][loc]=='-9':
                left[ind][loc]='?'
            if right[ind][loc]=='-9':
                right[ind][loc]='?'
            print "%s/%s " % (left[ind][loc],right[ind][loc]),
        print
    
