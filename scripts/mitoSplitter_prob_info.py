# -*- coding: utf-8 -*-
## used for plot ROC curve

import os
import re
import sys
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

bulk_clu = sys.argv[1]
prefix = sys.argv[2]

Y_pred = pd.read_csv(prefix + '/mitoSplitter.clusters.txt', sep = "\t", index_col = 0)
true_labels = pd.read_csv(bulk_clu, header=0, index_col=0, sep='\t', encoding='utf-8')
count_labels = cc(true_labels['Cluster'].tolist())

subindex = Y_pred.index.intersection(true_labels.index)
Y_pred = Y_pred.loc[subindex, ]
true_labels = true_labels.loc[subindex, ]

Y_pred['True_cluster'] = true_labels.loc[Y_pred.index, 'Cluster']
Y_pred.to_csv(prefix + '/mitoSplitter.clusters_with_truelabels.txt',index=True,header=True,sep='\t')

colors = ["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
           "#7E6148FF", "#B09C85FF", "#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#6F99ADFF",
           "#FFDC91FF", "#EE4C97FF", "#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", "#BA6338FF",
           "#5DB1DDFF", "#802268FF"]

bc_left = []
accuracy = []
recall = []
precision = []
f1 = []
fpr = []
tpr = []
threshold = []
auc = []
filter_probs = [0.9999999999]

for filter_prob in filter_probs:
    dd_i = Y_pred[Y_pred['Prob'] >= filter_prob]
    bc_left.append(dd_i.shape[0])
    true_l = dd_i['True_cluster']
    pred_l = dd_i['Cluster']
    indvs = np.unique(pred_l).tolist()
    acc = metrics.accuracy_score(true_l, pred_l)
    rec = metrics.recall_score(true_l, pred_l, average = "macro")
    eachrec = metrics.recall_score(true_l, pred_l, average = None)
    prec = metrics.precision_score(true_l, pred_l, average = "macro")
    eachprec = metrics.precision_score(true_l, pred_l, average = None)
    f1sc = metrics.f1_score(true_l, pred_l, average = "macro")
    eachf1sc = metrics.f1_score(true_l, pred_l, average = None)
    accuracy.append(acc)
    recall.append(rec)
    precision.append(prec)
    f1.append(f1sc)
    
    perform = pd.DataFrame({'Accuracy' : [acc],
                           'Recall' : [rec],
                           'Precision' : [prec],
                           'F1' : [f1sc],
                           'State' : ['Combine'],
                           'Sample' : ['MIX']})
    
    ## draw ROC curve for each sample
    eachacc = []
    eachauc = []

    plt.figure(figsize=(5,5))
    plt.title("Receiver Operating Characteristic")
    for eachlabel in set(true_l):
        true_l0 = [1 if x == eachlabel else 0 for x in true_l]
        pred_l0 = [1 if x == eachlabel else 0 for x in pred_l]
        accuracy0 = metrics.accuracy_score(true_l0, pred_l0)
        recall0 = metrics.recall_score(true_l, pred_l, average = "macro")
        precision0 = metrics.precision_score(true_l, pred_l, average = "macro")
        f10 = metrics.f1_score(true_l, pred_l, average = "macro")
        perform.loc[len(perform)] = [accuracy0, recall0, precision0, f10, 'Individual', eachlabel]

        #tmpdf = pd.DataFrame({'true_l' : true_l, 'true_l0' : true_l0, 
        #                      'pred_l' : pred_l, 'pred_l0' : pred_l0})
        fpr0, tpr0, threshold0 = metrics.roc_curve(true_l0, pred_l0)
        auc0 = metrics.auc(fpr0, tpr0)
        eachacc.append(accuracy0)
        eachauc.append(auc0)
        plt.plot(fpr0, tpr0, color=colors[indvs.index(eachlabel)], label = u'%s AUC = %0.3f' % (eachlabel, auc0))
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig(prefix + '/mitoSplitter_filter_prob' + str(filter_prob) +'.sample_ROC.pdf')
    plt.close()
    eachperform = pd.DataFrame({'Sample' : indvs,
                                'Accuracy' : eachacc,
                                'Recall' : eachrec,
                                'Precision' : eachprec,
                                'F1' : eachf1sc,
                                'AUC' : eachauc})
    eachperform.to_csv(prefix + '/mitoSplitter_filter_prob' + str(filter_prob) + '_each_perform.txt', sep = '\t')

    ## draw ROC curve for each prob
    true_l0 = preprocessing.label_binarize(true_l, classes = np.unique(true_l))
    pred_l0 = preprocessing.label_binarize(pred_l, classes = np.unique(pred_l))
    fpr0, tpr0, threshold0 = metrics.roc_curve(true_l0.ravel(), pred_l0.ravel())
    auc0 = metrics.auc(fpr0, tpr0)
    auc.append(auc0)
    
dd_line = pd.DataFrame({'Filter_prob' : filter_probs,
                        'Barcodes_left' : bc_left,
                        'Accuracy' : accuracy,
                        'Recall' : recall,
                        'Precision' : precision,
                        'F1' : f1,
                        #'fpr' : fpr,
                        #'tpr' : tpr,
                        #'thresholds' : threshold,
                        'AUC' : auc})
dd_line.to_csv(prefix + '/mitoSplitter_filter_prob_lineplot.txt', sep = '\t')
