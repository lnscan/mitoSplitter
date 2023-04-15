# -*- coding: utf-8 -*-
## used for plot ROC curve

import os
import re
import sys
import time
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import seaborn as sns
import scrublet as src
import matplotlib.pyplot as plt
from collections import Counter as cc
from sklearn import metrics
from sklearn import preprocessing
from sklearn.manifold import TSNE
from matplotlib.ticker import MultipleLocator
from sklearn.semi_supervised import LabelSpreading

#bulk_af = sys.argv[1]
#mix_af = sys.argv[2]
#mixnum = sys.argv[1]
#squal = sys.argv[2] ## sc qual
#bqual = sys.argv[3] ## bulk qual

#bulk_af = "MIX5/outs/bulk_q" + bqual +".af.txt"
#mix_af = "MIX" + mixnum + "/outs/whole_q" + squal +"_all.var_af.txt"
#prefix = "MIX" + mixnum + "/outs/sc" + squal +"_bulk" + bqual

bulk_af = sys.argv[1]
mix_af  = sys.argv[2]
prefix  = sys.argv[3]

start = time.time()
readstart = time.gmtime()
print(prefix, "start at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readstart))

if not os.path.exists(prefix):
    os.mkdir(prefix)

#mix_dge = 'dge.txt'
bulk = sc.read_csv(bulk_af, delimiter = "\t", first_column_names = True).T
#bulk = sc.read(bulk_af,cache=True).T
mix = sc.read_csv(mix_af, delimiter = "\t", first_column_names = True).T

sc.pp.highly_variable_genes(bulk, flavor='seurat_v3', span=0.3, n_top_genes=2000)
hvf_info = bulk.var
#hvf_info
hvf_info.to_csv(prefix + '/mitoSpliter.hvf_meta.scanpy.txt', sep='\t')
hvf = pd.DataFrame(bulk.var.index[bulk.var.highly_variable])    # True/False
hvf.to_csv(prefix + '/mitoSpliter.hvf.scanpy.txt', index=False, header=False)
plt.figure(figsize=(5,5))
x1 = list(bulk.var.means)
y1 = list(bulk.var.variances_norm)
x2 = list(bulk.var.means[bulk.var.highly_variable])
y2 = list(bulk.var.variances_norm[bulk.var.highly_variable])
plt.scatter(x1, y1, color='grey', s=40)
plt.scatter(x2, y2, color='red', s=40)
plt.xlabel("Mean frequency of variants", fontsize=12)
plt.ylabel("Dispersions of variants", fontsize=12)    # norm
plt.savefig(prefix + '/mitoSpliter.variant_selection.pdf')
plt.close()

#singlet_bcs = pd.read_csv('221122_mitospliter/result/pbmc_3/scrublet/single_bcs.list',header=None)
singlet_bcs = mix.obs_names
hvf2 = np.intersect1d(hvf, mix.var.index)
bcs2 = np.intersect1d(singlet_bcs, mix.obs_names)
mix_tmp = mix[bcs2, hvf2]
hvf2 = hvf2[np.sum(mix_tmp.X, axis=0)>0]
bulk2 = bulk[:, hvf2]
mix2 = mix[bcs2, hvf2]
hvf2df = pd.DataFrame({'hvf2':hvf2})
hvf2df.to_csv(prefix + '/mitoSpliter.hvf2.scanpy.txt', index=False, header=False)

##### Correlation-clusters_tmp
bulk_mix = ad.concat([bulk2,mix2], merge = "same")
corr = pd.DataFrame(np.corrcoef(bulk_mix.X))
corr.columns = bulk_mix.obs_names
corr.index = bulk_mix.obs_names
corr = corr.loc[mix2.obs_names,bulk2.obs_names]
corr['max']=corr.max(axis=1)
corr['clusters']=corr.idxmax(axis=1)
corr.to_csv(prefix + '/mitoSpliter.clusters_corr.txt',index=True,header=True,sep='\t')

plt.figure(figsize=(5,5))
plt.hist(corr['max'], bins = 50)
plt.xlabel("Max correlation", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
#plt.show()
plt.savefig(prefix + '/mitoSpliter.corr_hist.pdf')
plt.close()

##### LPA-clusters_final
data = mix[bcs2,:]
cells = list(data.obs_names)
# cluster labels processing
clu = corr
clu = clu[-(np.isnan(clu['max']))]
label_coding = {}
label_decoding = {-1 : "Unassigned"}
coding = 1
for bs in bulk.obs_names:
    label_coding[bs] = coding
    label_decoding[coding] = bs
    coding +=1
dict_clu = {}
for bc in list(clu.index.values):
    dict_clu[bc] = int( label_coding[ clu.loc[bc,'clusters'] ] )

labels = [dict_clu.get(c, -1) for c in cells]
labels = np.array(labels) 
print('Unlabels Number', list(labels).count(-1))

# LPA
label_prop_model = LabelSpreading(kernel='rbf', alpha=0.2, max_iter=100)
# fitting
label_prop_model.fit(data.X, labels)
# predicting
#Y_pred = label_prop_model.predict(data.X)
#predicted_clusters = [label_decoding[i] for i in Y_pred]
Y_pred_prob = label_prop_model.predict_proba(data.X)
Y_pred_prob = pd.DataFrame(Y_pred_prob)
Y_pred_prob[pd.isna(Y_pred_prob)] = 0
Y_pred_prob.index = cells
Y_pred_prob.columnsumns = label_coding.keys()
Y_pred_prob.to_csv(prefix + '/mitoSpliter.clusters_prob.txt',index=True,header=True,sep='\t')

Y_pred_prob_maxid = Y_pred_prob.idxmax(axis = 1)
Y_pred_prob_max = Y_pred_prob.max(axis = 1)
Y_pred_prob_cluster = [label_decoding[i+1] for i in Y_pred_prob_maxid]
Y_pred = pd.DataFrame({'cluster_id' : Y_pred_prob_maxid,
                       'cluster' : Y_pred_prob_cluster,
                       'prob' : Y_pred_prob_max})
Y_pred.index = cells
Y_pred.to_csv(prefix + '/mitoSpliter.clusters.txt',index=True,header=True,sep='\t')

end = time.time()
readend = time.gmtime()
print(prefix, "end at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readend))
print(prefix, "running time " + str(end - start) + " seconds.")
