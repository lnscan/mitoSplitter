##20221209

import sys
import scipy
import pandas as pd

rdir = sys.argv[1]
altbase = sys.argv[2]
rawpd = pd.read_csv(rdir + "." + altbase + ".alt.txt", sep = " ", header = None)

rawpd.columns = ['Pos', 'Barcodes', 'Count']

codpos, uniqpos = rawpd["Pos"].factorize(sort = True)
rawpd["codPos"] = codpos + 1
codbc, uniqbc = rawpd['Barcodes'].factorize(sort = True)
rawpd['codBC'] = codbc + 1

h1 = '''%%MatrixMarket matrix coordinate integer general
%metadata_json: {"format_version": 2, "software_version": "3.1.0"}'''
h2 = str(len(uniqpos)) + " " + str(len(uniqbc)) + " " + str(rawpd.shape[0])
wmtx = open(rdir + "." + altbase + ".alt/matrix.mtx", 'w')
wmtx.write("%s\n%s\n" % (h1, h2))
wmtx.close()
rawpd[["codPos", "codBC", "Count"]].to_csv(rdir + "." + altbase + ".alt/matrix.mtx", mode = "a", index = None, header = None, sep = " ")

uniqpos = uniqpos.astype(str) + "_"
uniqpospd = pd.DataFrame(uniqpos)
uniqpospd['gs'] = uniqpos
uniqpospd['ex'] = "Gene Expression"
uniqpospd.to_csv(rdir + "." + altbase + ".alt/features.tsv.gz", compression = 'gzip', index = None, header = None, sep = "\t")
uniqbcpd = pd.DataFrame(uniqbc)
uniqbcpd.to_csv(rdir + "." + altbase + ".alt/barcodes.tsv.gz", compression = 'gzip', index = None, header = None)
