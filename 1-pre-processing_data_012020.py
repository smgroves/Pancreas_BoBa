#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:53:47 2020

@author: meyerct6
"""

import loompy
import numpy as np
import pandas as pd
from collections import OrderedDict
import scprep
import sklearn.model_selection as ms
import pickle
import magic
#References:
#comparision: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8388285
#scImpute https://github.com/Vivianstats/scImpute
#MAGIC:   https://nbviewer.jupyter.org/github/KrishnaswamyLab/MAGIC/blob/master/python/tutorial_notebooks/emt_tutorial.ipynb

# =============================================================================
# Read through each sample and add to loompy
# =============================================================================
#For each protocol
for e in [1,2]:
    #Samples for a protocol
    sampls = ['x'+str(e)+'_S'+str(i) + 'c' for i in [3,4,5,6]]
    tds = OrderedDict()
    #Read in the loompy data
    for tp in sampls:
        loom_fn = f'./data/scbeta_indrops/01_Stages_3_to_6/data/complete_processing/{tp}.processed.loom'
        tds[tp] = loompy.connect(loom_fn)
    
    #Create a dataframe
    df = pd.DataFrame()
    #create data set
    for tp in sampls:
        df = df.append(pd.DataFrame(tds[tp][:,:],columns=tds[tp].ca.CellID,index=tds[tp].ra.Gene).T)

    #Create look-up table of cellular features
    look_up_table = pd.DataFrame([])
    for tp in sampls:
        tmp = pd.DataFrame({'CellStage':tds[tp].ca.CellStage,
                            'CellID':tds[tp].ca.CellID,
                            'CellLabel':tds[tp].ca.DetailedLabels,
                            'CellProtocol':tds[tp].ca.CellProtocol})
        look_up_table = look_up_table.append(tmp,ignore_index=True)
    look_up_table.set_index('CellID',inplace=True)


    #Now run scprep filtering
    df.fillna(0,inplace=True)

    #scprep.plot.plot_library_size(master_df,cutoff=2000)
    #Filter out cells with small library size
    df = scprep.filter.filter_library_size(df,cutoff=2000)
        
    #Filter out rare genes
    #Keep transcription factors no matter what
    tf_list = pd.read_table('./data/master_tf_list.txt')
    genes = tf_list['gene_symbol']
    to_keep = df.loc[:,genes[np.in1d(genes,list(df))]]
    df = df.loc[:,df.columns[~np.in1d(list(df),genes)]]
    df = scprep.filter.filter_rare_genes(df,cutoff=1,min_cells=100)
    #Add back in TFs
    df = pd.concat((to_keep,df),axis=1)    
    df.to_csv('./data/filtered_dataset_x'+str(e)+'.csv',index=True)
    look_up_table.to_csv(f'./data/cellID-lookuptable_x{e}.csv',index=True)    
    del df

# =============================================================================
# Select training and testing datasets
# =============================================================================
#Now combine the datasets and select the training and testing subsets
df = [];look_up_table=pd.DataFrame()
for e in [1,2]:
    df = df + list(pd.read_csv(f'./data/filtered_dataset_x{e}.csv',usecols=[0]).iloc[:,0])
    look_up_table = look_up_table.append(pd.read_csv(f'./data/cellID-lookuptable_x{e}.csv',index_col=0))
    
#Split testing and training dataset
kf = ms.StratifiedKFold(n_splits=5,random_state=1234)
train_index, test_index = next(kf.split(df, look_up_table.loc[df,'CellLabel']))
    
T = {'test_cellID':[df[i] for i in test_index],'test_index':test_index,'train_index':train_index,'train_cellID':[df[i] for i in train_index]}
with open('./data/test_train_indices.p','wb') as f:
    pickle.dump(T,f)
    
# =============================================================================
# Prepare data for scimpute
# =============================================================================
df = pd.DataFrame()
look_up_table = pd.DataFrame()
for e in [1,2]:
    df = df.append(pd.read_csv(f'./data/filtered_dataset_x{e}.csv',index_col=0))
    look_up_table = look_up_table.append(pd.read_csv(f'./data/cellID-lookuptable_x{e}.csv',index_col=0))
df.fillna(0,inplace=True)

for cl in look_up_table['CellLabel'].unique():
    print(cl)
    df.loc[df.index[np.in1d(df.index,look_up_table[look_up_table['CellLabel']==cl].index)]].T.to_csv(f'./data/isolated_pop_data_for_imputation-SCIMPUTE/{cl}_scimpute.csv',index=True)


