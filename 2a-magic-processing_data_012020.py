#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 07:36:02 2020

@author: meyerct6
"""
import magic
import pandas as pd
import scprep
import numpy as np
import pickle
# =============================================================================
# Read in data from 2 protocols
# =============================================================================
for e in [1,2]:
    df = pd.read_csv(f'./data/filtered_dataset_x{e}.csv',index_col=0)
    df.fillna(0,inplace=True)

# =============================================================================
# Run Magic
# =============================================================================
    #Start MAGIC
    magic_op = magic.MAGIC(n_jobs=24,knn=5,decay=1,t=3,knn_dist='euclidean',n_pca=100)
    #Run data preprocessing
    df = scprep.normalize.library_size_normalize(df)    
    df = scprep.transform.sqrt(df)
    #What transcription factors to include
    tf_list = pd.read_table('./data/master_tf_list.txt')
    #Fit to the transcription factors in list
    genes =list(tf_list.loc[np.in1d(tf_list['gene_symbol'],list(df)),'gene_symbol'].values) + ['INS','GCG']
    df_magic = magic_op.fit_transform(df, genes=genes)
    df_magic.to_csv(f'./data/magic_imputed_genes_completeTF_x{e}.csv')

df = pd.DataFrame()
look_up_table = pd.DataFrame()
for e in [1,2]:
    df = df.append(pd.read_csv(f'./data/magic_imputed_genes_completeTF_x{e}.csv',index_col=0))
    look_up_table = look_up_table.append(pd.read_csv(f'./data/cellID-lookuptable_x{e}.csv',index_col=0))
    df.fillna(0,inplace=True)
    
with open('test_train_indices.p','rb') as f:
    T = pickle.load(f)

df.loc[T['test_cellID']].to_csv('./data/magic_imputed_genes/magicTF_test.csv')
df.loc[T['train_cellID']].to_csv('./data/magic_imputed_genes/magicTF_train.csv')
look_up_table.to_csv('./data/cellID-lookuptable.csv')