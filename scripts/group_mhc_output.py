import sys
import os
import pandas as pd
import numpy as np

first = True
out = open(snakemake.output[0], "w")
for e in snakemake.input:
    if os.path.getsize(e) > 0:
        mhcout = open(e, 'r')
        alleles = next(mhcout).split('\t')
        header = next(mhcout).rstrip().split('\t')
        if first:
            allele = ''
            for i in range(0, len(header)):
                #print(header[i].rstrip())
                if i < len(alleles):
                    #print(alleles[i])
                    if alleles[i] != '':
                        allele = alleles[i].rstrip() + '_'
                header[i] = allele + header[i].rstrip()
            header[len(header) -1] = "NB"
            first = False
            #print(header)
            out.write('\t'.join(header) + '\n')
        for line in mhcout:
            out.write(line)
out.close()
