#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 07:42:57 2020

@author: meyerct6
"""
import pandas as pd
import numpy as np
import glob
import pickle

# =============================================================================
# Find all the sc-Impute datasets 
# after running R code in /data/isolated_pop_data_for_imputation-SCIMPUTE
# =============================================================================
tf_list = pd.read_table('master_tf_list.txt')
#Fit to the transcription factors in list
genes =list(tf_list['gene_symbol']) + ['INS','GCG']
df=pd.DataFrame()
for fil in glob.glob('./data/isolated_pop_data_for_imputation-SCIMPUTE/*/*.csv'):
    tmp = pd.read_csv(fil,index_col=0).T
    tmp = tmp.iloc[:,np.in1d(list(tmp),genes)]
    df = df.append(tmp)

df.fillna(0,inplace=True)
    
with open('./data/test_train_indices.p','rb') as f:
    T = pickle.load(f)

df.loc[T['test_cellID']].to_csv('./data/scimpute_genes/scimpute_test.csv')
df.loc[T['train_cellID']].to_csv('./data/scimpute_genes/scimpute_train.csv')
