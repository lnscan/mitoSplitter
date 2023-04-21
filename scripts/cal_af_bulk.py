##20221213

import os
import sys
import pandas as pd

prefix = sys.argv[1]

sample = os.path.basename(prefix)

pdcov = pd.read_csv(prefix + ".coverage.txt", header = None, index_col = 0, sep = " ")
pdcov.columns = ['Counts']

def cal_af_ATCG(prefix, base, pdcov):
	df = pd.read_csv(prefix + "." + base + ".txt", header = None, index_col = 0, sep = " ")
	df.columns = ['Counts', 'mq', 'refbase']
	df = df[df['refbase'].str.contains(base) == False]
	df = df.drop(['mq'], axis = 1)
	df['Af'] = df['Counts'] / pdcov['Counts']
	df = df.dropna(axis = 0)
	df.index = df.index.astype("str") + '_' + base
	#df.index = df.index.astype("str") + "_r" + df.refbase + '_' + base
	df = df.drop(['Counts', 'refbase'], axis = 1)
	return(df)

afA = cal_af_ATCG(prefix, 'A', pdcov)
afC = cal_af_ATCG(prefix, 'C', pdcov)
afG = cal_af_ATCG(prefix, 'G', pdcov)
afT = cal_af_ATCG(prefix, 'T', pdcov)

af = pd.concat([afA, afC, afG, afT], axis = 0)
af.columns = [sample]
af.to_csv(prefix + ".af.txt", sep = "\t")

