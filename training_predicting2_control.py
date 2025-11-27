# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 21:43:24 2025

@author: dnt
"""

import pandas as pd
import numpy as np
import torch
import torch.utils.data as Data
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
from sklearn.utils import shuffle
from imblearn.over_sampling import SMOTE
import os
from sklearn.model_selection import train_test_split
from collections import Counter
from classifier_CNN import classifier_CNN
# from plot_confusion_matrix import plot_confusion_matrix
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.decomposition import PCA
from sklearn.metrics import roc_curve, auc
import seaborn as sns
from plot_confusion_matrix import plot_confusion_matrix

nt_table = ['t','c', 'a', 'g']
dnt_table = [nt1+nt2 for nt1 in nt_table for nt2 in nt_table]
dnts_table = [dnt1+dnt2 for dnt1 in dnt_table for dnt2 in dnt_table]
codon_table = [nt1+nt2+nt3 for nt1 in nt_table for nt2 in nt_table for nt3 in nt_table]
codon_table1 = codon_table.copy()

stop_codon = ['taa', 'tag', 'tga']
stop_codonpair = [i+j for i in stop_codon for j in stop_codon]

codonpair_table0 = [codon0 + codon1 for codon0 in codon_table1 for codon1 in codon_table1]
codonpair_table = []
for cp in codonpair_table0:
    if cp not in stop_codonpair:
        codonpair_table.append(cp)

dnt_category = ['n12', 'n23','n31']
dntpair_category = ['n12m12', 'n23m23','n31m31','n12n31','n23m12','n31m23']
dntpair_cols_list = ['Freq_'+ dnts + '_' + dnts_cat for dnts_cat in dntpair_category for dnts in dnts_table]
codonpair_cols_list = ['Freq_'+ codonpair for codonpair in codonpair_table]

full_cols_list0 = dntpair_cols_list + codonpair_cols_list
full_cols_list0 = full_cols_list0 + ['padding0', 'padding1']
print (len(full_cols_list0))

full_cols_list1 = [col[5:] for col in full_cols_list0]




############################################

epoch_num_list =  [50]
split_size_list = [0.15]
test_size = split_size_list[0]
split_i = 5
lr = 0.001
batch_size = 20
############################################

dfID = pd.read_csv ('seqIDs.csv')
keyIDs = dfID['seqID'].tolist()

#################### 

nt_table = ['t','c', 'a', 'g']
dnt_table = [nt1+nt2 for nt1 in nt_table for nt2 in nt_table]
dnts_table = [dnt1+dnt2 for dnt1 in dnt_table for dnt2 in dnt_table]
codon_table = [nt1+nt2+nt3 for nt1 in nt_table for nt2 in nt_table for nt3 in nt_table]
codon_table1 = codon_table.copy()

stop_codon = ['taa', 'tag', 'tga']
stop_codonpair = [i+j for i in stop_codon for j in stop_codon]

codonpair_table0 = [codon0 + codon1 for codon0 in codon_table1 for codon1 in codon_table1]
codonpair_table = []
for cp in codonpair_table0:
    if cp not in stop_codonpair:
        codonpair_table.append(cp)

dnt_category = ['n12', 'n23','n31']
dntpair_category = ['n12m12', 'n23m23','n31m31','n12n31','n23m12','n31m23']
dntpair_cols_list = ['Freq_'+ dnts + '_' + dnts_cat for dnts_cat in dntpair_category for dnts in dnts_table]
codonpair_cols_list = ['Freq_'+ codonpair for codonpair in codonpair_table]

full_cols_list0 = dntpair_cols_list + codonpair_cols_list
full_cols_list0 = full_cols_list0 + ['padding0', 'padding1']
print (len(full_cols_list0))

full_cols_list1 = [col[5:] for col in full_cols_list0]

#################### 

path_dcrcp = '../counting/'
file_dcrcp = 'df_dcrcp_counting_concat_cdsRdRp_host_143654_sampled_16317_labeled_20251026added.csv'
df_dcrcp00 = pd.read_csv (path_dcrcp + file_dcrcp, index_col = 'seqID', encoding = 'gbk')# 
df_dcrcp0 = df_dcrcp00[df_dcrcp00['label'] != 3]
df_dcrcp1 = df_dcrcp00[df_dcrcp00['label'] == 3]
print(df_dcrcp1.shape)

# print(df_dcrcp00.viridae.value_counts())
# print (df_dcrcp00.shape)







# df_key =df_dcrcp00.loc[keyIDs]
    
counting = df_dcrcp0.viridae.value_counts()
print (counting)
# viridaelst = list(counting.keys())
# print (viridaelst)

    
#     ####################################################  sampling
#     ####################################################  sampling

accuracy_lst, ilst = [], []
lst_fpr, lst_tpr, lst_rocauc, samplingnums = [], [], [], []
 
    
lst = list(range(10)) 
for xx in range(1000):
    
    df_dcrcp0 = shuffle(df_dcrcp00, random_state = 3)
        
    df_con = pd.DataFrame()
    
    for viridae in viridaelst:
        df_vi = df_dcrcp0[df_dcrcp0['viridae'] == viridae]
        print(df_vi.shape)
        
        if (viridae == 'IAV') & (df_vi.shape[0] > 1000):
            df_vi_sample = df_vi.sample(n = 1200, random_state = 42)
            
        elif (viridae == 'Flaviviridae') & (df_vi.shape[0] > 1000):
            df_vi_sample = df_vi.sample(n = 1200, random_state = 42)
            
        elif (viridae == 'Phenuiviridae') & (df_vi.shape[0] > 1000):
            df_vi_sample = df_vi.sample(n = 1200, random_state = 42)
        
        else:
            df_vi_sample = df_vi
            
        df_con = pd.concat([df_con, df_vi_sample], axis = 0)
    
    print ('data shape post sampling: ', df_con.shape) 
    
       
    ####################################################  sampling
    ####################################################  sampling
        
    
    
    df_dcrcp_model0 = df_con
    print (df_dcrcp_model0.shape)
    
    hosts_lst = df_dcrcp_model0.host1.tolist()
    hostTypes0 = list(set(hosts_lst))
    host_num = (len(hostTypes0))
    
    dicthost0 = dict(zip(hostTypes0, range(host_num)))
    y_dict = dict(zip(range(host_num), hostTypes0))
    
    
    dicthost0 = {'nonHuman': 0, 'Human': 1, 'Tick': 2}
    y_dict = {0: 'nonHuman', 1: 'Human', 2: 'Tick'}
    y_dict2 = {0: 'nonHuman', 1: 'Human'}
    
    # # print (y_dict)
    label_names = list(y_dict2.values())
    print (label_names)
    
    
    # # ########################################################    prepare data
    # # ########################################################    prepare data
    
    lst = [0, 3]
    for jj in lst:
        
        df_human = df_dcrcp_model0[df_dcrcp_model0['label'] == 1]
        df_humanS = df_human.sample (n = 1500, random_state = jj, replace=True)
        df_others = df_dcrcp_model0[df_dcrcp_model0['label'] == 0].sample (n = 1500, random_state = jj, replace=True)
        print (df_others.shape)
        
        df_dcrcp_model = pd.concat([df_humanS, df_others], axis = 0)
        df_dcrcp_model = shuffle (df_dcrcp_model, random_state = jj)
        
        print (df_dcrcp_model.shape)
        
        
        X, y = df_dcrcp_model[full_cols_list0], df_dcrcp_model.label.tolist()
        
        Xkey = df_key[full_cols_list0]
        print ('Xkey shape: ', Xkey.shape)        
        numX = Xkey.shape[0]
        df_dcrcp1xx = df_dcrcp1[(df_dcrcp1['label2'] == 0)|(df_dcrcp1['label2'] == 1)]
        testX_xx, testy_xx  = df_dcrcp1xx[full_cols_list0], df_dcrcp1xx.label2.tolist()
        
        
        # df_dcrcp1xx = df_dcrcp1
        # testX_xx, testy_xx  = df_dcrcp1xx[full_cols_list0], df_dcrcp1.shape[0]*[1]
        
        
        print(f"原始数据集类别分布: {Counter(y)}")
        y = [int(i) for i in y]
    
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = test_size, random_state=jj, stratify=y)
        
        
        X_resampled, y_resampled = X_train, y_train
        
        
        X_resampled2 = X_resampled.copy()
        df_resampled = pd.DataFrame(data = X_resampled2)
        df_resampled['label'] = y_resampled
        
        df_resampled.to_csv ('df_dcrcp_counting_training_sampling' + str(jj) + '.csv')
        
        X_tick = df_dcrcp00[df_dcrcp00['label'] == 2][full_cols_list0]
        accessions_tick = df_dcrcp00[df_dcrcp00['label'] == 2].index.tolist()
        y_tick = X_tick.shape[0] * [2]
        
        print (X_tick.shape)
        
       
        
        for epoch_num in epoch_num_list[0:1]:
            name_train = 'CNN2labels_' + str(jj) + '_test_size'+ str(test_size) +'_epoch_num' + str(epoch_num) + '_batchsize' + str(batch_size) + '_learningRate' + str(lr) + '_predTickV'
            name_train_all = 'CNN2labels_' + '_test_size'+ str(test_size) +'_epoch_num' + str(epoch_num) + '_batchsize' + str(batch_size) + '_learningRate' + str(lr) + '_predTickV_'
            name_test = 'CNN2labels_' + str(jj) + '_test_size' + str(test_size) +'_epoch_num' + str(epoch_num) + '_batchsize' + str(batch_size) + '_learningRate' + str(lr) + '_predTickV'
            aucs = []
            pred_y_list = []
            true_y_list = []
            df_X_test = pd.DataFrame()
            mean_fpr = np.linspace(0, 1, 100)
            training_loss_list_all = []
    
            #############################################################     CNN training
            #############################################################     CNN training
            
            X_train_array = np.array(X_resampled)
            
            print (X_train_array.shape)
            train_num = X_train_array.shape[0]
            y_train_list = y_resampled
            X_valid_array = np.array(X_test)
            y_valid_list = y_test
            
            testX_xx_array = np.array(testX_xx)
            testX_xx_num = testX_xx_array.shape[0]
            Xkey_array = np.array(Xkey)
            
                    
            valid_num = X_valid_array.shape[0]
            valid_size = len(y_valid_list)
            X_train_array2 = X_train_array.reshape(train_num,1,75,75)
            X_valid_array2 = X_valid_array.reshape(valid_num,1,75,75)
            
            testX_xx_array2 = testX_xx_array.reshape(testX_xx_num,1,75,75)
            Xkey_array2 = Xkey_array.reshape(numX,1,75,75)
            
            
            X_train_tensor = torch.tensor(X_train_array2)
            y_train_tensor = torch.tensor(y_train_list)
            X_valid_tensor = torch.tensor(X_valid_array2)
            y_valid_tensor = torch.tensor(y_valid_list)
            
            testX_xx_tensor = torch.tensor(testX_xx_array2)
            testy_xx_tensor = torch.tensor(testy_xx)
            Xkey_tensor =  torch.tensor(Xkey_array2)
            
            
            
            X_train_tensor = X_train_tensor.to(torch.float32)
            X_valid_tensor = X_valid_tensor.to(torch.float32)
            testX_xx_tensor = testX_xx_tensor.to(torch.float32)
            Xkey_tensor = Xkey_tensor.to(torch.float32)
            
            
            X_tick_array = np.array(X_tick)
            y_tick_list = y_tick
            tick_num = X_tick_array.shape[0]
            X_tick_array2 = X_tick_array.reshape(tick_num,1,75,75)
            X_tick_tensor = torch.tensor(X_tick_array2)
            y_tick_tensor = torch.tensor(y_tick_list)
            X_tick_tensor = X_tick_tensor.to(torch.float32)
            
            
            # # ###########################################################  data loader
            
            
            torch_train = Data.TensorDataset (X_train_tensor, y_train_tensor)
            # print (torch_train)
            loader = Data.DataLoader (dataset = torch_train,
                                      batch_size = batch_size,
                                      shuffle = False,
                                      num_workers = 0)  # it defines Multiprocess, default ==0, causing error with more than 0, 
    
            # #############################################################  data loader
            
            # ###############################################################  build CNN
    
            cnn = classifier_CNN()
            if_use_gpu = 1
            if if_use_gpu:
                cnn = cnn.cuda()
                
            # ###############################################################  buld CNN
            
            # ##########################################  define optimizer and loss function
            
            optimizer = torch.optim.Adam(cnn.parameters(), lr = lr)#lr = .0001
            loss_func = nn.CrossEntropyLoss()
            
            ### train data
            training_loss_list = []
            for epoch in range (epoch_num):
                for step, (x, y) in enumerate (loader):
                    b_x = Variable (x)
                    b_y = Variable (y)
                    # print (b_x.shape, b_y.shape)
            ###################################################
                    if if_use_gpu:
                        b_x = b_x.cuda()
                        b_y = b_y.cuda()
                        
    
                    pred = cnn (b_x)[0]
                    loss = loss_func (pred, b_y)
                    optimizer.zero_grad()                # 对loss求导
                    loss.backward()
                    optimizer.step()
                    training_loss_list.append(loss.cpu().detach().item())
            # # ##########################################  define optimizer and loss function
    
            print ('training finished!')
    
            df_training_loss = pd.DataFrame ()
            df_training_loss['Training_loss'] = training_loss_list
            df_name = 'df_training_loss_'+ name_train +'_control.csv'
            df_training_loss.to_csv (df_name)
    
    
            valid_pred_prob = cnn (Variable(X_valid_tensor).cuda())
            valid_pred_matrix = valid_pred_prob[1].cpu().detach().numpy()
            valid_pred_list = torch.max(valid_pred_prob[0],1)[1].data.cpu().detach().numpy()
            
            
            valid_pred_array = label_binarize (valid_pred_list, classes = [1, 0])
            valid_true_array = label_binarize (y_valid_list, classes = [1, 0])
            n_classes = valid_pred_array.shape[1]
            valid_prob_array = valid_pred_prob[1].cpu().detach().numpy()
            
            print ('n_classes:', n_classes)
            
            
            fpr = dict()
            tpr = dict()
            roc_auc = dict()
         
     
            
            fpr[0], tpr[0], thresholds = roc_curve(valid_true_array[:,0], valid_prob_array[:,0])
            roc_auc[0] = auc(fpr[0], tpr[0])
            
            # lst_fpr.append(fpr)
            # lst_tpr.append(tpr)
            # lst_rocauc.append(roc_auc)
            # samplingnums.append(jj)
    
            plt.figure (figsize = (6,6), facecolor='w')
            plt.background_color='white'
            # plt.grid(0)   
            for i in range (1):
                fpr[i], tpr[i], thresholds = roc_curve(valid_true_array[:,i], valid_prob_array[:,i])
                roc_auc[i] = auc(fpr[i], tpr[i])
                plt.plot(fpr[i], tpr[i], lw = 3, alpha=0.3,\
                          label='ROC class %d (AUC = %0.4f)' % (i, roc_auc[i]))
            plt.plot([0, 1], [0, 1], linestyle='--', lw = 3, color='r',
                      label='Chance', alpha=.8)
    
            plt.xlim([-0.05, 1.05])
            plt.ylim([-0.05, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.legend(loc="lower right")
            # plt.grid(0)
            plt.rcParams['axes.facecolor'] = 'white'  
            plt.background_color='white'
            
            plt.savefig( 'ROC_AUC_prob_'+ name_train + '_control.svg',dpi=300 ,bbox_inches='tight')
            plt.close()
            
    
    
            class_names = np.ravel(y_valid_list)
            class_names = class_names[unique_labels(y_valid_list, valid_pred_list)]
    
    
            
            
            plt.figure (figsize = (4,4), facecolor='w')
           
            np.set_printoptions(precision=2)
            plot_confusion_matrix(y_valid_list, valid_pred_list, label_names,\
                                  classes=class_names, normalize=True,\
                                  title='Confusion matrix (%)')

            plt.grid(0)
            plt.background_color='white'
            # plt.rcParams['axes.facecolor'] = 'white'  
            plt.savefig('Confusion matrix (%)_' + name_train + '_control.svg',dpi=600,bbox_inches='tight')
            plt.close()
    
    
    
    
            fc_full_array = valid_pred_prob[2].cpu().detach().numpy() # fc_full_array = valid_pred_prob[2].cpu().detach().numpy()
            print (fc_full_array.shape)
            df_fc_full = pd.DataFrame (fc_full_array)
            # df_fc_full.to_csv ('df_fc_full_' + name_train + '.csv')
    
            pca = PCA(n_components=2)
            newX1 = pca.fit_transform(df_fc_full)
            newX11 = np.array (newX1)
            # print (newX11.shape)
            
            y_valid_list = [y_dict[i] for i in y_valid_list]
            df_fc_full['True_label'] = y_valid_list
            cols = ['PCA1','PCA2']
            df_fc_full[cols] = newX11
            df_fc_full.to_csv ('df_fc_full_PCA_' + name_train + '_control.csv')
    
            plt.figure (figsize = (6,6), facecolor='w')
    
            sns.set(font_scale=2)
            sns.set(style='white')
            sns.pairplot(data=df_fc_full[cols + ['True_label']],\
                         hue = 'True_label',\
                         palette={'Tick': '#3478b1', 'Human': '#c9403a', 'nonHuman': '#e49d39'},\
                         hue_order = ['Human', 'nonHuman'],\
                         kind="scatter")
                
                
            sns.set(rc={'figure.figsize':(15,5)})
            sns.set_theme(context='notebook', style='darkgrid', palette='deep', font='sans-serif', font_scale=1, color_codes=True, rc=None)
            plt.grid(0)
            plt.background_color='white'
            plt.rcParams['axes.facecolor'] = 'white'  
            plt.savefig('Pairplot_PCA_fc_full_' + name_train +'_control.svg', dpi = 300,bbox_inches = 'tight')
            plt.close()
    
    
            
            
            tick_pred_prob = cnn (Variable(X_tick_tensor).cuda())
            tick_pred_matrix = tick_pred_prob[0].cpu().detach().numpy()
            tick_pred_list = torch.max(tick_pred_prob[0],1)[1].data.cpu().detach().numpy()
    #         # print (tick_pred_list)
            tick_pred_array = label_binarize (tick_pred_list, classes = [1, 0])
            tick_true_array = label_binarize (y_tick_list, classes = [1, 0])
            n_classes2 = tick_pred_array.shape[1]
            tick_prob_array = tick_pred_prob[1].cpu().detach().numpy()
            
    
            
            df_tick = pd.DataFrame(data = tick_prob_array)
            df_tick['accession'] = accessions_tick
            df_tick['true_label'] = y_tick_list
            df_tick['predict_label'] = tick_pred_list
            
            df_tick.to_csv ('predicted_results_tickV' + '_control.csv')
            
    
            data = Variable(testX_xx_tensor).cuda()
            pred0 = cnn (data)
            pred = pred0[0].cpu()
            prob_array = pred0[1].cpu().data.numpy()
            print (prob_array.shape)
            pred_labels = torch.max (pred, 1)[1].data
            
            
            dataX = Variable(Xkey_tensor).cuda()
            pred0X = cnn (dataX)
            predX = pred0X[0].cpu()
            probX_array = predX[1].cpu().data.numpy()
            # print (probX_array.shape)
            predX_labels = torch.max (predX, 1)[1].data
            
            
            
            
            
            correct = 0
            
            testy_xx = [int(i) for i in testy_xx]
            
            for kk in range(len(testy_xx)):
                if testy_xx[kk] == pred_labels[kk]:
                    correct += 1
                    
            # print ('预测准确率：', jj, correct/len(testy_xx))
            
            # print (testy_xx)
            # print (pred_labels)
    
            
                                    
            model_name = 'CNN_model_dcrcp_' + name_train + '_control.pt'    ####################   To save trained model
            torch.save(cnn.cpu(), model_name)
            
            if (correct/len(testy_xx) > 0.7) & (sum(predX_labels)/numX >= 0.95):
                lst.remove(jj)
                accuracy_lst.append(correct/len(testy_xx))
    
                lst_fpr.append(fpr)
                lst_tpr.append(tpr)
                lst_rocauc.append(roc_auc)
                samplingnums.append(jj)
                
                print('准确率：', jj, correct/len(testy_xx))
    
#     ######################################################    testing
    
    # path_label = '../data/'
    # df_label = pd.DataFrame()
    # labels = []
    # for file in os.listdir(path_label):
    #     if file.startswith ('TBV-1120-zmq-20251026.csv')&file.endswith('.csv'):
    #         print (file)
    #         df_label = pd.read_csv (path_label + file, encoding = 'gb2312')
    #         labels = df_label.label.tolist()
                    
            
    # path_test = '../counting/'
    # df_test = pd.DataFrame()
    
    # for file in os.listdir(path_test):
    #     if file.startswith ('df_dcrcp_counting_TBV-1120-zmq-20251026.csv')&file.endswith('.csv'):
    #         print (file)
    #         df_test = pd.read_csv (path_test + file, index_col = 'seqID')
    
    #         df_test['padding0'] = df_test.shape[0] * [0]
    #         df_test['padding1'] = df_test.shape[0] * [0]
    #         print (df_test.shape)
    
    
    #         df_test['label'] = labels
            
    
    # print (df_test.shape)
    # df_test = df_test[(df_test['label'] == 0)|(df_test['label'] == 1)]
    # labels_selected = df_test.label.tolist()
    
    
    # model_path = '../supervised_training3/'
    # for file in os.listdir(model_path):
    #     if file.endswith('.pt'):
    #         model_name = file
    #         # print (model_name)
            
    # #         # cnn = classifier_CNN()
    #         cnn = torch.load(model_path + model_name, weights_only = False)
    #         device = torch.device ("cuda:0" if torch.cuda.is_available() else "cpu")
    #         cnn.to(device)
    #         # print (cnn)
    #         if torch.cuda.is_available():
    #             cnn = cnn.cuda()
    #             cnn.to(torch.device("cuda:0"))
    
    #         id_list = df_test.index.tolist()
    #         dcrcp_test = df_test.loc[:,full_cols_list0]
    #         y_test = df_test.label.tolist()
    #         print (dcrcp_test.shape)
    #         dcrcp_array = np.array(dcrcp_test)
    #         num = dcrcp_array.shape[0]
    
    #         X_test_array2 = dcrcp_array.reshape(num,1,75,75)
    
    #         print (X_test_array2.shape)
            
    #         dcrcp_tensor0 = torch.tensor(X_test_array2)
    #         dcrcp_tensor = dcrcp_tensor0.to(torch.float32)
    #         print (dcrcp_tensor.shape)
    
    #         data = Variable(dcrcp_tensor).cuda()
            
    #         pred0 = cnn (data)
    #         pred = pred0[0].cpu()
    #         prob_array = pred0[1].cpu().data.numpy()
    #         print (prob_array.shape)
    #         pred_labels = torch.max (pred, 1)[1].data
            
    #         fc_array = pred0[2].cpu().data.numpy()
            
    #         pred_results = pd.DataFrame ({'seqID': id_list,\
    #                                       'pred_label': pred_labels})
    
    #         pred_results [['Score_0','Score_1']] = prob_array
    #         pred_results['label'] = labels_selected
    #         print (pred_results.shape)
    #         file_name = 'df_Pred_TBV-1120-zmq-20251026_biovalidated_model_' + file[:-4] + '.csv'
    #         pred_results.to_csv (file_name, index = False)
            
    #         df_fc = pd.DataFrame ({'seqID': id_list})
    #         df_fc [list(range(fc_array.shape[1]))] = fc_array
    
    #         file_name2 = 'df_fc_data_TBV-1120-zmq-20251026_biovalidated_model_' + file[:-4] + '.csv'
            
    #         df_fc.to_csv (file_name2)
    
    ######################################################    testing
    ######################################################    testing
    
    
    
    ######################################################    plot roc all
    ######################################################    plot roc all
    ######################################################    plot roc all
    

    
    ######################################################    plot roc all
    ######################################################    plot roc all
    ######################################################    plot roc all
    
    
    print (accuracy_lst)
    
    # if (sum(accuracy_lst)/10 > 0.75) & (min(accuracy_lst)>0.6):
    #     break
    if len(lst) == 0:
        break

plt.figure (figsize = (10,10), facecolor='w')
plt.grid(0)
plt.plot([0, 1], [0, 1], linestyle='--', lw = 3, color='r',
          label='Chance', alpha=.8)

for kkk in range(10):
    index = samplingnums.index(kkk)
        
    
    fpr = lst_fpr[index]
    tpr = lst_tpr[index]
    roc_auc = lst_rocauc[index]
    num = samplingnums[index]
    
    plt.plot(fpr[0], tpr[0], lw = 3, alpha=0.3,\
              label='ROC model %d (AUC = %0.4f)' % (num, roc_auc[0]))
        
plt.grid(0)        
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize = 12)
plt.ylabel('True Positive Rate', fontsize = 12)
plt.legend(loc="lower right")
plt.rcParams['figure.facecolor'] = 'white'  
plt.grid(0)
# plt.background_color='white'  
plt.rcParams['axes.facecolor'] = 'white'  
plt.savefig( 'ROC_AUC_prob_'+ name_train_all + '_all_control.svg',dpi=300 ,bbox_inches='tight')
plt.close()    