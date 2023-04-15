##20221212

import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

start = time.time()
readstart = time.gmtime()
print("start at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readstart))

prefix = sys.argv[1]

def pivot_mtx(mtxdir):
	scdata = sc.read_10x_mtx(mtxdir)
	scvmat = scv.DataFrame(scdata.X)
	scvmat.index = scdata.obs.index
	scvmat.columns = scdata.var.gene_ids
	return(scvmat.T)

def cal_ATCG_af(prefix, base):
	basemat = pivot_mtx(prefix + "." + base + ".alt")
	covmat = pivot_mtx(prefix + ".coverage." + base + ".alt")
	afmat = basemat / (covmat + 0.00000001)
	afmat.index = afmat.index + base
	return(afmat)

afA = cal_ATCG_af(prefix, 'A')
afC = cal_ATCG_af(prefix, 'C')
afG = cal_ATCG_af(prefix, 'G')
afT = cal_ATCG_af(prefix, 'T')

af = pd.concat([afA, afC, afG, afT], axis = 0)
af = af.fillna(0)
#af = af.loc[(af!=0).any(1), (af!=0).any(0)]

af.to_csv(prefix + ".var_af.txt", sep = "\t")

end = time.time()
readend = time.gmtime()
print("end at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readend))
print("running time " + str(end - start) + " seconds.")
