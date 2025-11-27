# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 18:52:54 2025

@author: dnt
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import Levenshtein as lv

path_data = '../data/'
file_data = 'df_concat_cdsRdRp_labels_143654.csv'
df_data = pd.read_csv (path_data + file_data, index_col = 'seqID')
print (df_data.shape)

df_index = pd.DataFrame()
path = '../supervised_training3_control/'

df_con = pd.DataFrame()
for file in os.listdir(path):
    if file.endswith('.csv')&file.startswith('df_dcrcp'):
        print (file)
        
        df = pd.read_csv (path + file, index_col = 'seqID')
        
        indexlst = df.index.tolist()
        
        df_data_indexed = df_data.loc[indexlst]
        print(df_data_indexed.shape)
        df_data_indexed = df_data_indexed.drop_duplicates(keep = 'first')
        
        df_con = pd.concat([df_con, df_data_indexed])


        
df_data_indexedH = df_con[df_con['host_Standardization'] == 'Homo sapiens']
df_data_indexednonH = df_con[df_con['host_Standardization'] != 'Homo sapiens']


seqRef_h, reqRef_nonh = df_data_indexedH.iloc[0, 0], df_data_indexednonH.iloc[0, 0]

distance_htoh, distance_htononh, distance_nonhtoh, distance_nonhtononh = [], [], [], []

for i in range(df_data_indexedH.shape[0]):
    seq_h = df_data_indexedH.iloc[i, 0]
    distance0 = lv.ratio(seq_h, seqRef_h)
    distance_htoh.append(distance0)

    distance1 = lv.ratio(seq_h, reqRef_nonh)
    distance_htononh.append(distance1)
    
    
for j in range(df_data_indexednonH.shape[0]):
    seq_nonh = df_data_indexednonH.iloc[j, 0]
    distance2 = lv.ratio(seq_nonh, seqRef_h)
    distance_nonhtoh.append(distance2)

    distance3 = lv.ratio(seq_nonh, reqRef_nonh)
    distance_nonhtononh.append(distance3)
    
df_dist_h, df_dist_nonh = pd.DataFrame(), pd.DataFrame()

df_dist_h['similarity_htoh'] = distance_htoh
df_dist_h['similarity_htononh'] = distance_htononh

df_dist_nonh['similarity_nonhtoh'] = distance_nonhtoh
df_dist_nonh['similarity_nonhtononh'] = distance_nonhtononh

df_dist_h.to_csv ('df_h_similarity_all' + '.csv')
df_dist_nonh.to_csv ('df_nonh_similarity_all' + '.csv')
        
        