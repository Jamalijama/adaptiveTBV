# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 08:51:36 2024

@author: Jing Li, Small steps make change. dnt_seq@163.com
"""


from sklearn.manifold import TSNE
# from sklearn.datasets import load_iris
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import umap
import math
from math import cos,sin,radians
from sklearn.linear_model import LinearRegression
import os
from sklearn.utils import shuffle



amino_table = ['I', 'D', 'M', 'H', 'E', 'W', 'R', 'L', 'Y', 'Q', 'G', 'A', 'S', 'P', 'C', 'T', 'V', 'F', 'N', 'K']
nt_table = ['t','c', 'a', 'g']
dnt_table = [nt1+nt2 for nt1 in nt_table for nt2 in nt_table]
dnts_table = [dnt1+dnt2 for dnt1 in dnt_table for dnt2 in dnt_table]
codon_table = [nt1+nt2+nt3 for nt1 in nt_table for nt2 in nt_table for nt3 in nt_table]
codon_table1 = codon_table.copy()
codon_table1.remove('taa')
codon_table1.remove('tag')
codon_table1.remove('tga')
codonpair_table = [codon0 + codon1 for codon0 in codon_table1 for codon1 in codon_table1]
dnt_category = ['n12', 'n23','n31']
dntpair_category = ['n12m12', 'n23m23','n31m31','n12n31','n23m12','n31m23']
dnts_cols_list = ['Freq_'+ dnt + '_' + dnt_cat for dnt_cat in dnt_category for dnt in dnt_table]
dntpair_cols_list = ['Freq_'+ dnts + '_' + dnts_cat for dnts_cat in dntpair_category for dnts in dnts_table]
codon_cols_list = ['Freq_'+ codon for codon in codon_table]
codonpair_cols_list = ['Freq_'+ codonpair for codonpair in codonpair_table]
amino_cols_list = ['Freq_'+ amino for amino in amino_table]
# full_cols_list0 = dnts_cols_list + dntpair_cols_list + codon_cols_list + codonpair_cols_list + amino_cols_list

dcr_set_name_list = ['dnts', 'codons', 'aminos', 'DCRcodonpair']
dcr_set_list = [dnts_cols_list,codon_cols_list, amino_cols_list,  dntpair_cols_list+codonpair_cols_list]

cnames = {
'deepskyblue':          '#00BFFF',
'purple':               '#800080',
'orchid':               '#DA70D6',
'lightgreen':           '#90EE90',
'darkgreen':            '#006400',
'lightblue':            '#ADD8E6',
'yellow':               '#FFFF00',
'yellowgreen':          '#9ACD32',
'deeppink':             '#FF1493',
'burlywood':            '#DEB887',
'red':                  '#FF0000',
'indianred':            '#CD5C5C',
'darkred':              '#8B0000',
    }
color_num_list = list (range(1,16,1))
color_dict = dict(zip(color_num_list,cnames.values()))
color_list0 = list(color_dict.values())


color_list_all = []
for i in color_list0:
    color_list_all.append(i)
# print (color_list_all)

path_label = '../data/'
file_label = 'virusesformodeling0827.csv'
df_labels0 = pd.read_csv (path_label + file_label, index_col = 'ID')
label_lst0 = df_labels0['host'].tolist()
order_lst0 = df_labels0['order'].tolist()

print (df_labels0['host'].value_counts())
# print (df_labels0['order'].value_counts())

hosts = ['Carnivora', 'Tick', 'Human', 'Rodentia', 'Artiodactyla', 'Chiroptera', 'Primate']
hosts_new = ['Animal', 'Tick', 'Human', 'Animal', 'Animal', 'Animal', 'Animal']

dicth = dict(zip(hosts, hosts_new))
print (dicth)

hostTypes = hosts_new

dict_host = dict(zip(hosts_new, range(1,4,1)))
print (dict_host)
dict_host2 = dict(zip(range(1,4,1), hostTypes))
print (dict_host2)

label_lst00 = [dicth[i] for i in label_lst0]


label_lst = [dict_host[i] for i in label_lst00]

y_types = list(set(label_lst))
print (y_types)

hue_order = ['Tick', 'Mammal', 'Human']

path = '../counting/'
path_us = './'
for file in os.listdir(path):
    if file.endswith ('.csv') & file.startswith('df_full_DCR_counting_virusesformodeling0827'):
        
        gene = 'polymerase'
        print (file, gene)
        
        df_DCR0 = pd.read_csv (path + file, index_col = 'seqID')    
        df_all = pd.concat([df_DCR0, df_labels0], axis = 1)        
        df_all['label'] = label_lst
        
        for kk in range(10):
        
            df_all_s = shuffle(df_all, random_state = kk)
            # print (df_all_s['host'].value_counts())
            
            df_all_h = df_all_s[df_all_s['host'] == 'Human']
            df_all_t = df_all_s[df_all_s['host'] == 'Tick']
            
            df_all_others = df_all_s[(df_all_s['host'] != 'Human')&(df_all_s['host'] != 'Tick')]
            # print (df_all_others['host'].value_counts())
            
            df_all_hs = df_all_h.sample(n = 200)
            df_all_ts = df_all_t.sample(n = 200)
            
            df_con = pd.concat([df_all_hs, df_all_ts, df_all_others], axis = 0)
            
            num = df_con.shape[0]
            
            print (df_all.shape, df_con.shape, df_labels0.shape)
            
            print (df_con['host'].value_counts())
    
    
            
            index_lst = df_con.index.tolist()
            label_lst1 = df_con['label'].tolist()
            label_lst2 = df_con['host'].tolist()
            # print (label_lst2)
            
            
            label_lst3 = [dicth[jj] for jj in label_lst2]
            # print (label_lst3)
            
            df_con['host'] = label_lst3
            
            print (df_con['host'])
            
            y_types = list(set(label_lst3))
            # print (y_types)
            label_lst3 = ['Mammal' if i == 'Animal' else i for i in label_lst3]
    
            y_num = len(y_types)
            label_cate_list = y_types
            label_c_list = list(range(1,8,1))
            # print (label_c_list)
        
            dict1 = dict(zip(label_cate_list,label_c_list))
            dict2 = dict(zip(label_cate_list,color_dict))
            label_list1 = [dict1[i] for i in label_cate_list]
            # print (label_list1)
            color_list = [color_dict[i] for i in label_list1]
            # print(color_list)
            
            # print (df_DCR.value_counts('subtype')[:50])
        
            for set_i in range(3,4,1):
                set_name = dcr_set_name_list[set_i]
                dcr_set = dcr_set_list[set_i]
                # print (set_name, len(dcr_set))
                data = np.array (df_con[dcr_set])
                # print (data.shape)
            
                X_tsne = TSNE(learning_rate = 100).fit_transform(data)
                X_pca = PCA(n_components = 2).fit_transform(data)
                X_umap = umap.UMAP(n_components = 2).fit_transform(data, label_lst1) #, min_dist=0.001
                # print (X_tsne.shape)
                # print ('X_pca', X_pca.shape)
                print (X_umap.shape)
        
        
                df_tsne = pd.DataFrame (X_tsne,columns = ['t_SNE1','t_SNE2'])
                df_tsne = (df_tsne - df_tsne.min()) / (df_tsne.max() - df_tsne.min())
                df_tsne ['host'] = label_lst3
                
                df_umap = pd.DataFrame (X_umap, columns = ['UMAP1','UMAP2'])
                df_umap = (df_umap - df_umap.min()) / (df_umap.max() - df_umap.min())
                df_umap ['host'] = label_lst3    
                print (df_umap.shape)
        
                df_pca = pd.DataFrame (X_pca,columns = ['PCA1','PCA2'])
                df_pca = (df_pca - df_pca.min()) / (df_pca.max() - df_pca.min())    
                df_pca ['host'] = label_lst3
                
                file_name0 = 'sns_scatterplot_tSNE_PCA_UMAP_tickV_sampled_' + '_' + set_name + '_' + str(kk) + '_3labels_20250118.svg'
                # file_name1 = 'df_PCA_tickV_sampled' + '_' + set_name + '_' + str(kk) + '.csv'
                # file_name2 = 'df_tSNE_tickV_sampled' + set_name + '_' + str(kk) + '.csv'
                # file_name3 = 'df_UMPA_tickV_sampled' + set_name + '_' + str(kk) + '.csv'
        
            
                plt.figure(figsize=(3,12), facecolor='white')
                plt.subplot(311)
                sns.scatterplot(data = df_tsne, x = 't_SNE2', y = 't_SNE1',\
                                markers=True, hue = 'host',\
                                palette= {'Tick': '#3478b1',\
                                         'Human': '#c9403a',\
                                         'Mammal': '#e49d39'},\
                                hue_order = hue_order)
                plt.legend(scatterpoints=1)
    
                plt.subplot(312)
                sns.scatterplot(data = df_umap, x = 'UMAP2', y = 'UMAP1',\
                                markers=True, hue = 'host',\
                                palette={'Tick': '#3478b1',\
                                         'Human': '#c9403a',\
                                         'Mammal': '#e49d39'},\
                                hue_order = hue_order)
                    
                plt.legend(scatterpoints=1)
                
                plt.subplot(313)
                # for y_i in range(y_num):
                #     y_ = y_types[y_i]
                    # color = color_list0[y_i]
                #     plt.xlim([-0.1,1.1])
                #     plt.ylim([-0.1,1.1])
        
                    # df_X_pca_label = df_pca[df_pca['host'] == y_]
                #     # print (df_X_pca_label.shape)
                sns.scatterplot(data = df_pca, x = 'PCA2', y ='PCA1',\
                                markers=True, hue = 'host',\
                                palette={'Tick': '#3478b1',\
                                         'Human': '#c9403a',\
                                         'Mammal': '#e49d39'},\
                                hue_order = hue_order) # ,x_estimator=np.mean
        
                plt.legend(scatterpoints=1)
                plt.savefig(path_us +'11' + file_name0, dpi = 300, bbox_inches = 'tight')
                plt.close()

                plt.figure (figsize = (4,4))
                sns.set(font_scale=2)
                sns.set(style='white')
                sns.pairplot(data=df_pca,hue = 'host', kind="scatter", palette={'Tick': '#3478b1', 'Human': '#c9403a', 'Mammal': '#e49d39'}, hue_order = hue_order)
                sns.set(rc={'figure.figsize':(15,5)})
                sns.set_theme(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=True, rc=None)
                plt.savefig('11Pairplot_PCA_' + set_name + '_' + str(kk) + '_3labels_20250118.svg', dpi = 600,bbox_inches = 'tight')
                plt.close()


              
                df_pca['accession'] = index_lst
                
                print (df_con.shape, df_pca.shape, df_umap.shape, df_tsne.shape)
                df_con['PCA1_' + set_name] = df_pca['PCA1'].tolist()
                df_con['PCA2_' + set_name] = df_pca['PCA2'].tolist()
                
                df_con['t_SNE1_' + set_name] = df_tsne['t_SNE1'].tolist()
                df_con['t_SNE2_' + set_name] = df_tsne['t_SNE2'].tolist()
                
                df_con['UMAP1_' + set_name] = df_umap['UMAP1'].tolist()
                df_con['UMAP2_' + set_name] = df_umap['UMAP2'].tolist()
                
                
            df_con.to_csv (path_us + 'df_DCRcodonpair_polymerase_CDS_tickV_sampled_PCAtSNEUMAP_added' + '_' + str(kk) +'.csv', encoding = 'gbk')            
            
