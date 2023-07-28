## remove doublets using scrublet in mtx format

import os
import sys 
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import scrublet as scr 
import matplotlib.pyplot as plt 

prefixid = sys.argv[1] ## output dir
removebc = sys.argv[2] ## remove barcodes list before doublet detection
countdir = sys.argv[3] ## /path/to/filtered_feature_bc_matrix

## read remove bc list
rmbc = []
with open(removebc, 'r') as rr: 
    for line in rr.readlines():
        rmbc.append(line.rstrip())

## read cellranger output
data = sc.read_10x_mtx(countdir, cache = True)

data = data[~data.obs.index.isin(rmbc), :]
sc.pp.highly_variable_genes(data, flavor='seurat_v3', span=0.3, n_top_genes=2000) ## calculate hvf info for each gene in data

## run scrublet
bc_num_raw = data.n_obs
# cell/var filtering
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)
# expected_doublet_rate calculating
bc_num = data.n_obs
nn = 0.008 * (bc_num / 1000)
nn = round(nn,3)
# Scrublet object initializing
scrub = scr.Scrublet(data.X, expected_doublet_rate=nn)
# predict doubelts
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,min_cells=3)
# output
doublets_info = pd.DataFrame( { 'Barcode':data.obs_names,
                                'Predicted_doublets':predicted_doublets,
                                'Doublet_scores':doublet_scores} )
doublets_info.to_csv(prefixid + '/doublets_info.scrublet.txt',index=False,header=True,sep='\t')
# singlet/doublet barcodes
singlet_bcs = doublets_info.Barcode[( bool(i) for i in (1-doublets_info.Predicted_doublets) )]
doublet_bcs = doublets_info.Barcode[doublets_info.Predicted_doublets]
## if scrublet cannot auto find threshold, make nn the threshold
#singlet_bcs = doublets_info.Barcode[( bool(i) for i in (doublets_info.Doublet_scores < nn) )]
#doublet_bcs = doublets_info.Barcode[( bool(i) for i in (doublets_info.Doublet_scores >= nn) )]

singlet_bcs.to_csv(prefixid + '/scrublet_singlet.list',header = False, index = False)
doublet_bcs.to_csv(prefixid + '/scrublet_doublet.list',header = False, index = False)

