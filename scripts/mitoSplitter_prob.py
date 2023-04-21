# -*- coding: utf-8 -*-
## used for multiplexed samples using mitoSplitter

import os
import re
import sys
import time
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import seaborn as sns
#import scrublet as src
import matplotlib.pyplot as plt
from collections import Counter as cc
from sklearn import metrics
from sklearn import preprocessing
from sklearn.manifold import TSNE
from matplotlib.ticker import MultipleLocator
from sklearn.semi_supervised import LabelSpreading

## read files
bulk_af = sys.argv[1] ## bulk RNA-seq mtRNA matrix
mix_af  = sys.argv[2] ## single-cell RNA-seq mtRNA matrix
prefix  = sys.argv[3] ## output dir

start = time.time()
readstart = time.gmtime()
print(prefix, "start at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readstart))

if not os.path.exists(prefix):
    os.mkdir(prefix)

bulk = sc.read_csv(bulk_af, delimiter = "\t", first_column_names = True).T
mix = sc.read_csv(mix_af, delimiter = "\t", first_column_names = True).T

## identification of highly variable features(sites) from bulk matrix
sc.pp.highly_variable_genes(bulk, flavor='seurat_v3', span=0.3, n_top_genes=2000)
hvf_info = bulk.var
hvf_info.to_csv(prefix + '/mitoSplitter.hvf_meta.scanpy.txt', sep='\t')
hvf = pd.DataFrame(bulk.var.index[bulk.var.highly_variable])    # True/False
hvf.to_csv(prefix + '/mitoSplitter.hvf.scanpy.txt', index=False, header=False)
plt.figure(figsize=(5,5))
x1 = list(bulk.var.means)
y1 = list(bulk.var.variances_norm)
x2 = list(bulk.var.means[bulk.var.highly_variable])
y2 = list(bulk.var.variances_norm[bulk.var.highly_variable])
plt.scatter(x1, y1, color='grey', s=40)
plt.scatter(x2, y2, color='red', s=40)
plt.xlabel("Mean frequency of variants", fontsize=12)
plt.ylabel("Dispersions of variants", fontsize=12)    # norm
plt.savefig(prefix + '/mitoSplitter.variant_selection.pdf')
plt.close()

## Overlap of bulk hvf and single-cell variants

#Remove singlets using Scrublet 
#mix_dge = sys.argv[4] ## single-cell RNA-seq gene expression matrix
#dge = sc.read(mix_dge).T
#expected_dbl_rate = {0:[0.000,0.008],
#                     1:[0.008,0.008],
#                     2:[0.016,0.007],
#                     3:[0.023,0.008],
#                     4:[0.031,0.008],
#                     5:[0.039,0.007],
#                     6:[0.046,0.008],
#                     7:[0.054,0.007],
#                     8:[0.061,0.008],
#                     9:[0.069,0.007],
#                     10:[0.076,0.008]}
#bc_num_raw = dge.n_obs
## cell/var filtering
#sc.pp.filter_cells(dge, min_genes=200)
#sc.pp.filter_genes(dge, min_cells=3)
## expected_doublet_rate calculating
#bc_num = dge.n_obs
#nn = int(bc_num/1000)
#nn2 = expected_dbl_rate[nn]
#nn3 = nn2[0]+nn2[1]*(bc_num-nn*1000)/1000
#nn3 = round(nn3,3)
## Scrublet object initializing
#scrub = scr.Scrublet(dge.X, expected_doublet_rate=nn3)
## predict doubelts
#doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,min_cells=3)
## output
#doublets_info = pd.DataFrame( {'Barcode':dge.obs_names,
#                               'Predicted_doublets':predicted_doublets,
#                               'Doublet_scores':doublet_scores} )
#doublets_info.to_csv('mitoSplitter.doublets_info.scrublet.txt',index=False,header=True,sep='\t')
## singlet/doublet barcodes
#singlet_bcs = doublets_info.Barcode[( bool(i) for i in (1-doublets_info.Predicted_doublets) )]
#doublet_bcs = doublets_info.Barcode[doublets_info.Predicted_doublets]

#singlet_bcs = pd.read_csv('xx/scrublet/single_bcs.list',header=None) ## read singlet barcode list here if known
singlet_bcs = mix.obs_names
hvf2 = np.intersect1d(hvf, mix.var.index)
bcs2 = np.intersect1d(singlet_bcs, mix.obs_names)
mix_tmp = mix[bcs2, hvf2]
hvf2 = hvf2[np.sum(mix_tmp.X, axis=0)>0]
bulk2 = bulk[:, hvf2]
mix2 = mix[bcs2, hvf2]
hvf2df = pd.DataFrame({'hvf2':hvf2})
hvf2df.to_csv(prefix + '/mitoSplitter.hvf2.scanpy.txt', index=False, header=False)

## Correlation-clusters
bulk_mix = ad.concat([bulk2,mix2], merge = "same")
corr = pd.DataFrame(np.corrcoef(bulk_mix.X))
corr.columns = bulk_mix.obs_names
corr.index = bulk_mix.obs_names
corr = corr.loc[mix2.obs_names,bulk2.obs_names]
corr['max']=corr.max(axis=1)
corr['Clusters']=corr.idxmax(axis=1)
corr.to_csv(prefix + '/mitoSplitter.clusters_corr.txt',index=True,header=True,sep='\t')

plt.figure(figsize=(5,5))
plt.hist(corr['max'], bins = 50)
plt.xlabel("Max correlation", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
#plt.show()
plt.savefig(prefix + '/mitoSplitter.corr_hist.pdf')
plt.close()

## LPA-clusters
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
    dict_clu[bc] = int( label_coding[ clu.loc[bc,'Clusters'] ] )

labels = [dict_clu.get(c, -1) for c in cells]
labels = np.array(labels) 
print('Unlabels Number', list(labels).count(-1))

# LPA
label_prop_model = LabelSpreading(kernel='rbf', alpha=0.2, max_iter=100)
# fitting
label_prop_model.fit(data.X, labels)
# predicting
Y_pred_prob = label_prop_model.predict_proba(data.X)
Y_pred_prob = pd.DataFrame(Y_pred_prob)
Y_pred_prob[pd.isna(Y_pred_prob)] = 0
Y_pred_prob.index = cells
Y_pred_prob.columnsumns = label_coding.keys()
Y_pred_prob.to_csv(prefix + '/mitoSplitter.clusters_prob.txt',index=True,header=True,sep='\t')

Y_pred_prob_maxid = Y_pred_prob.idxmax(axis = 1)
Y_pred_prob_max = Y_pred_prob.max(axis = 1)
Y_pred_prob_cluster = [label_decoding[i+1] for i in Y_pred_prob_maxid]
Y_pred = pd.DataFrame({'Cluster_id' : Y_pred_prob_maxid,
                       'Cluster' : Y_pred_prob_cluster,
                       'prob' : Y_pred_prob_max})
Y_pred.index = cells
Y_pred.to_csv(prefix + '/mitoSplitter.clusters.txt',index=True,header=True,sep='\t')

## TSNE plot
colors = ["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
           "#7E6148FF", "#B09C85FF", "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
           "#FFDC91FF", "#EE4C97FF", "#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", "#BA6338FF",
           "#5DB1DDFF", "#802268FF"]
mix3=mix2
mix3=mix3[np.sum(mix3.X, axis=1)>0,:]
mix3_result=Y_pred.loc[mix3.obs_names,'Cluster']
X_tsne = TSNE(n_components=2).fit_transform(mix3.X)
X_tsne_data = np.vstack((X_tsne.T, mix3_result)).T
df_tsne = pd.DataFrame(X_tsne_data, columns=['tSNE1','tSNE2','Cluster'])
df_tsne = df_tsne.sort_values(by='Cluster', ascending=True)
bulknum = df_tsne['Cluster'].nunique()
plt.figure(figsize=(8,8))
sns.scatterplot(data=df_tsne, x='tSNE1', y='tSNE2', hue='Cluster', palette=colors[0:bulknum], size=0.1, linewidth=0)
ax=plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
plt.title('tSNE visualization')
plt.savefig(prefix + '/mitoSplitter.tSNE.pdf')
plt.close()

end = time.time()
readend = time.gmtime()
print(prefix, "end at " + time.strftime("%Y-%m-%d %H:%M:%S" ,readend))
print(prefix, "running time " + str(end - start) + " seconds.")
