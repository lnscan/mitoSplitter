##20230102

import os
import sys 
import pysam
import pandas as pd

bulkl  = sys.argv[1]
mitofa = sys.argv[2]
resdir = sys.argv[3]

bulkfs = []
with open(bulkl, 'r') as bulkr:
    for line in bulkr.readlines():
        bulkfs.append(line.rstrip())

mitofasta = pysam.FastaFile(mitofa)
maxBP = mitofasta.lengths[0]
mitoseq = mitofasta.fetch(mitofasta.references[0], 0, maxBP)
reflist = []
for i in range(maxBP):
    reflist.append(str(i+1) + "_" + mitoseq[i].upper())

whole_af_var = pd.DataFrame()
for eachbulk in bulkfs:
    af_var  = pd.read_csv(eachbulk, sep = "\t", header = 0, index_col = 0)
    if eachbulk == bulkfs[0]:
        whole_af_var = af_var.copy(deep = True)
    else:
        whole_af_var = pd.merge(whole_af_var, af_var, how = 'outer', left_index = True, right_index = True) 

whole_af_var = whole_af_var[~whole_af_var.index.isin(reflist)]
whole_af_var = whole_af_var.fillna(0)
whole_af_var.to_csv(resdir + '/bulk_all.af.txt', sep = "\t")

