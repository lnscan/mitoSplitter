## check correlation between bulk references

import sys 
import numpy as np
import pandas as pd
from itertools import combinations

bulk_af = sys.argv[1] ## bulk RNA-seq mtRNA matrix
cor_war = sys.argv[2] ## correlation warning value

outfile = bulk_af.replace("txt", "cor.txt")
outp = open(outfile, 'w')

bulk = pd.read_csv(bulk_af, sep = "\t", header = 0, index_col = 0)
for combins in combinations(np.unique(bulk.columns), 2): 
    bulk_com = bulk[list(combins)]
    spear = bulk_com.corr().loc[combins[0],combins[1]]
    sspear = str(spear)
    outp.write(f"{combins[0]}\t{combins[1]}\t{sspear}\n")
    if spear > float(cor_war):
        print(f"Warning: The correlation between {combins[0]} and {combins[1]} is {sspear}, which is greater than {cor_war}. This suggests that the genetic differences between the two samples are relatively small.")
      
