# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 22:08:17 2019

@author: huang
"""

import pyBigWig
# import numpy as np

infile = "MH63.seedling.H3K4me3.RNA.fwd.RPKM.bw"
# outfile = "test.bw"
binsize = 5
chrom = 'all'

bw = pyBigWig.open(infile)
# outbw = pyBigWig.open(outfile,'w')

chrom_dict = bw.chroms()

for i in chrom_dict.keys():
    if chrom == 'all' or i == chrom:
#        outbw.addHeader([(i,chrom_dict[i])],maxZooms=0)
        flag = None
        count = 0
        for j in range(0,chrom_dict[i],binsize):
            value = bw.values(i,j,j+1)[0]
            if flag != None and flag != value:
                print("{0}\t{1}\t{2}\t{3}".format(i,j-binsize*count,j,flag))
                count = 0
            flag = value
            count += 1
        print("{0}\t{1}\t{2}\t{3}".format(i,j-binsize*(count-1),chrom_dict[i],flag))
        
bw.close()