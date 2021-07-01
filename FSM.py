# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 17:49:03 2019

@author: ppttssy
"""

import random
import numpy as np
from sklearn import metrics
from sklearn.metrics import auc
from sklearn import preprocessing
from sklearn.svm import OneClassSVM
from sklearn.metrics import ndcg_score
from sklearn.model_selection import KFold
from sklearn.metrics import precision_recall_curve


alls_half=[]
Gold_dict={}
Train_drug=[]
alls_half_dict={}
dict_drug_index={}
dict_index_drug={}
Filter_final_dicts={}
floder = KFold(n_splits=2,shuffle=True)
with open("data//drugs.txt",'r') as r:
        index=0    
        for i in r:
                i=i.strip('\n')
               # Filter_final_dicts[i]=1
                dict_drug_index[i]=index
                dict_index_drug[index]=i
                index+=1
with open("data//Gold_hypertension.txt",'r') as f:
        for i in f:
            i=i.strip('\n').split('\t')
            if i[0] in dict_drug_index and i[1] in dict_drug_index:
                Gold_dict[i[0]+'_'+i[1]]=1
                Gold_dict[i[1]+'_'+i[0]]=1

with open("data//target_tanimoto.txt",'r') as r:
    index=0
    for line in r:
        line=line.strip('\n').split('\t')
        if line[0] in dict_drug_index and line[1] in dict_drug_index and line[0]!=line[1]:
              alls_half.append(line[0]+'_'+line[1])
              alls_half_dict[line[0]+'_'+line[1]]=index
              alls_half_dict[line[1]+'_'+line[0]]=index
              index+=1

def ddmatrixs(position):
    ddmatrix=np.zeros((len(dict_drug_index),len(dict_drug_index)))
    with open("data//"+str(position),'r') as r:
        for line in r:
            line=line.strip('\n').split('\t')
            if line[0]+'_'+line[1] in alls_half_dict :
                ddmatrix[dict_drug_index[line[0]]][dict_drug_index[line[1]]]=float(line[2])
                ddmatrix[dict_drug_index[line[1]]][dict_drug_index[line[0]]]=float(line[2])      
    return ddmatrix

#def Friendship(position):
    

if __name__=='__main__':
    
    
    '''
    Friendship feature
    '''
    X_pathway,X_atc,X_Off,X_CHE,X_seperation=ddmatrixs("pathway_tanimoto.txt"),ddmatrixs("atc_tanimotol.txt"),ddmatrixs("offside_tanimoto.txt"),ddmatrixs("chemical_tanimoto.txt"),ddmatrixs("target_tanimoto.txt")
    temp1=np.hstack((X_pathway,X_atc))
    temp2=np.hstack((temp1,X_Off))
    temp3=np.hstack((temp2,X_seperation))
    temp4=np.hstack((temp3,X_CHE))
    temp4 = preprocessing.scale(temp4)
    Drug_combine_X=np.zeros((len(alls_half),1525))
    
    for index,value in enumerate(alls_half):
        Drug_combine_X[index]=temp4[dict_drug_index[value.split("_")[0]]]+temp4[dict_drug_index[value.split('_')[1]]]

    # for ssk in range(100):
    #     w1 = open("D://drug//result//hypertension_result_305_"+str(ssk)+".txt",'w')
    #     w1.close()
  #  a=[]
    for ssk in range(100):
        # w1 = open("D://drug//result//hypertension_result_305_"+str(ssk)+".txt",'a')
        alls_halfs=[]
        alls_halfs_dict={}
        y_s=[]
        index_1=0
        Drug_combine_X_exteart=[]
        for index in range(len(alls_half)):
            if alls_half[index] in Gold_dict:      
                y_s.append(1)
    
            else:
                y_s.append(0)
            alls_halfs.append(alls_half[index])
            alls_halfs_dict[alls_half[index]]=index_1
            Drug_combine_X_exteart.append(index)
            index_1+=1           
            
        y_s=np.array(y_s)      
        Drug_combine_X_1=Drug_combine_X[Drug_combine_X_exteart]
        for train_index,val_index in floder.split(alls_halfs,y_s):  
            G,P=[],[]
            test_solo=[]
            train_solo=[]
            test_all_gold=[]
            test_all=[]
            train_solo_fea=[]
            for value in train_index:
                    train_solo_fea.append(alls_halfs[value].split('_')[0])
                    train_solo_fea.append(alls_halfs[value].split('_')[1])  
            for value in train_index:
                if y_s[value]==1:
                    G.append(dict_drug_index[alls_halfs[value].split('_')[0]])
                    G.append(dict_drug_index[alls_halfs[value].split('_')[1]]) 
                    train_solo.append(alls_halfs[value].split('_')[0])
                    train_solo.append(alls_halfs[value].split('_')[1])  
            for value in val_index:
                    test_all.append(alls_halfs[value])
                    if dict_drug_index[alls_halfs[value].split('_')[0]] not in G and dict_drug_index[alls_halfs[value].split('_')[0]] not in P:                          
                        P.append(dict_drug_index[alls_halfs[value].split('_')[0]])
                    if dict_drug_index[alls_halfs[value].split('_')[1]] not in G and dict_drug_index[alls_halfs[value].split('_')[1]] not in P:                              
                        P.append(dict_drug_index[alls_halfs[value].split('_')[1]])
                        
            G2=[]                            
            G=list(set(G))
            P=list(set(P))
            test_solo=list(set(test_solo)) 
            train_solo=list(set(train_solo))     
            train_solo_fea=list(set(train_solo_fea))     
            y_train_solo=np.zeros(len(train_solo_fea))
            for index,value in enumerate(train_solo_fea):
                G2.append(dict_drug_index[value])
                if value in train_solo:
                    y_train_solo[index]=1
                    
            
            Filter=[]                    
            for index in range(100):
                for train_solo2,val_solo in floder.split(temp4[G2],y_train_solo):  
                  #  if sum(y_train_solo[train_solo2])>1 and sum(y_train_solo[val_solo])>1:
                      clf = OneClassSVM(gamma='auto',kernel='linear').fit(temp4[G2][train_solo2][y_train_solo[train_solo2]==1])
                      y_s_svr=clf.score_samples(temp4[G2][val_solo])
                      fpr, tpr, thresholds =  metrics.roc_curve(y_train_solo[val_solo], y_s_svr)
                      y = tpr - fpr
                      Youden_index = np.argmax(y)  # Only the first occurrence is returned.
                      optimal_threshold = thresholds[Youden_index]

                      Filter.append(optimal_threshold)
            clfs = OneClassSVM(gamma='auto',kernel='linear').fit(temp4[G])
            Essthr=np.mean(Filter)  
            
            
            allstest=[]
            y_s_svr=clfs.score_samples(temp4[P])
            for index_2 in range(len(y_s_svr)):
                      allstest.append((dict_index_drug[P[index_2]],y_s_svr[index_2]))

            allstest.sort(key=lambda x:-x[1])
            Train_all=[]
            Train_drug=[]
            for index in range(len(allstest)):  
                    if   allstest[index][1] > Essthr:       
                        Train_drug.append(allstest[index][0])


            Train_drug=list(set(train_solo)|set(Train_drug))
            for index in range(len(Train_drug)):
                          for index_2 in range(len(Train_drug)):
                            if index>index_2:
                                Train_all.append(Train_drug[index]+'_'+Train_drug[index_2])
                                Train_all.append(Train_drug[index_2]+'_'+Train_drug[index])  
            
            final_test=list(set(test_all)&set(Train_all))    
            X_train=[]
            for index in range(len(train_index)):
                if  alls_halfs[train_index[index]] in Gold_dict: 
                    X_train.append(Drug_combine_X_1[alls_halfs_dict[alls_halfs[train_index[index]]]])
            X_train=np.array(X_train)


            X_test,y_test=np.zeros(((len(test_all)),1525)),np.zeros(len(test_all))
            for index in range(len(test_all)):
                X_test[index]=Drug_combine_X_1[alls_halfs_dict[test_all[index]]]
                if  test_all[index] in Gold_dict:
                    y_test[index]=1
        
            clf = OneClassSVM(gamma='auto',kernel='rbf').fit(X_train)
            y_svr=clf.score_samples(X_test)

            for index in range(len(test_all)):
                if test_all[index] not in Train_all:
                    y_svr[index]=-100000
            one_drug_list_valid2=[]
            for index in range(len(test_all)):
                one_drug_list_valid2.append((test_all[index],y_svr[index]))
            one_drug_list_valid2.sort(key=lambda x:-x[1])     
            
            
            
            '''
            predicted drug combinations
            '''
            # for index in range(int(len(test_all)/20)):
            #     w1.write(one_drug_list_valid2[index][0]+'\t'+str(one_drug_list_valid2[index][1])+'\n')       
          
            
            '''
            evalutation
            '''
            fpr, tpr, thresholds = metrics.roc_curve(y_test, y_svr, pos_label=1)
            precision, recall, pr_thresholds = precision_recall_curve(y_test, y_svr)
            aupr_score = auc(recall, precision) 
            print(str(np.trapz(tpr,fpr)))
            print(aupr_score)    
            print(ndcg_score(np.reshape(np.array(y_test),(1, -1)),np.reshape(np.array(y_svr),(1,-1))))

    
    
        # w1.close()