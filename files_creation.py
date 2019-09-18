# -*- coding: utf-8 -*-
"""
this script should run first, it will create all needed files for the rest of the program.

"""
import numpy as np
import pandas as pd
import utiles
import pickle
import scipy
RBPs = pd.read_csv('data files\\RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])
utiles.createFolder('figures')
utiles.createFolder('temporary files')

WTDic = {}
for wt in WTs:
    WTDic[len(wt)] = wt 
#==============================================================================
#create the distance files
WT = {}
for i in range(len(proteins)):
    WT[proteins[i]] = WTs[i]
    
Ldic = {}    
for key,val in WT.items():
    Ldic[key] = len(val)
    
for key,val in WT.items():
    L = len(val)
    Data = pd.read_excel('data files\\Data.xlsx',key)
    D=[]
    for p in proteins:
        D.append(Data[p][Data['site length']==L])
    distances = pd.DataFrame(D).T
    distances.to_csv('temporary files\\dist_'+key+'.csv')
#==============================================================================
#create SNPS all proteins
for key,val in WT.items():
    utiles.Create_SNPS(val,key+' all snps')
#==============================================================================
#create di-nucs only MS2    
utiles.Creat_Di(WT['MS2'],'MS2 all di_nuc')
#create data files for Fully Connected
#============================================================================== 
for key,val in WT.items():
    
    Data = pd.read_excel('data files\\Data.xlsx',key)
    Data = Data[Data['site length']==Ldic[key]]
    Data['seq'] = Data['seq'].astype(str).str.upper()
    DataC = Data[Data['prefix']=='C']
    DataGC = Data[Data['prefix']=='GC']
    
    DataC = {'X':utiles.OneHot(list(DataC['seq'])),'Y':DataC['score']}
    DataGC ={'X':utiles.OneHot(list(DataGC['seq'])),'Y':DataGC['score']}
    
#    DataC = np.asarray(DataC)
#    DataC = np.swapaxes(DataC,1,2)
#    
#    DataGC = np.asarray(DataGC)
#    DataGC = np.swapaxes(DataGC,1,2)
    with open('temporary files\\DataForFC_'+key+'_C.txt','wb') as fp:
        pickle.dump(DataC,fp) 
    with open('temporary files\\DataForFC_'+key+'_GC.txt','wb') as fp:
        pickle.dump(DataGC,fp)            
#create data files for CNN
#============================================================================== 
L = 50      #preset value of the length of the input for the CNN


Dic={'PP7':0,'MS2':1,'Qb':2}
full_mat = scipy.io.loadmat('data files\\structure_matrix.mat')

with open('data files\\TSS.txt','r') as f:
    TSS = f.readline()
with open('data files\\mc.txt','r') as f:
    mc = f.readline()
LTSS = len(TSS)
for key,val in WT.items():
    Data = pd.read_excel('data files\\Data.xlsx',key)
    Data = Data.replace(np.nan, '', regex=True)
    
    mat = full_mat['struct_mat'][0,Dic[key]]
    FullSeqs = []
     #sequecne encoding
    for i in range(len(Data)):
        tmp = TSS+Data['prefix'][i]+Data['seq'][i]+Data['suffix'][i]+mc
        FullSeqs.append(tmp[LTSS-10:LTSS-10+L].upper())
    OHSeqs = utiles.OneHot(FullSeqs)
    
    #structure encoding
    FULLData=[]
    for i in range(len(Data)):
        SeqMat = np.asarray(OHSeqs[i])
        StructMat = mat[i,0]
        L1 = len(Data['prefix'].iloc[i])
        Lseq = StructMat.shape[1]
        LeftZeros = np.zeros([5,10+L1])
        RightZeros = np.zeros([5,L-10-L1-Lseq])
        Struc = np.concatenate((LeftZeros,StructMat,RightZeros),axis=1) 
        FullMat = np.concatenate((SeqMat.T,Struc),axis=0)
        
        FULLData.append(FullMat.T)
    FULLData = {'X':FULLData,'Y':Data['score']}
    with open('temporary files\\DataForCNN_'+key+'.txt','wb') as fp:
        pickle.dump(FULLData,fp)
    

#%% check for common sequences in the library
MS2 = pd.read_excel('data files\\Data.xlsx','MS2')
PP7 = pd.read_excel('data files\\Data.xlsx','PP7')
Qb = pd.read_excel('data files\\Data.xlsx','Qb')


MS2['seq'] = MS2['seq'].astype(str).str.lower()
MS2 = MS2.sort_values(by=['seq','prefix'])
PP7['seq'] = PP7['seq'].astype(str).str.lower()
PP7 = PP7.sort_values(by=['seq','prefix'])
Qb['seq'] = Qb['seq'].astype(str).str.lower()
Qb = Qb.sort_values(by=['seq','prefix'])


All_Data = pd.DataFrame(['AAAAA',1,1,1]).T
All_Data.columns = ['seq','MS2','PP7','Qb']

com_seq = []
MS2scores = []
PP7scores = []
Qbscores =[]
M=[]
P = []
Q = []
CGC=[]
for m in range(len(MS2)):
    flag=0
    for p in range(len(PP7)):
        if flag:
            break
        #if MS2 and PP7 has seq in common
#        if MS2['seq'].iloc[m]<PP7['seq'].iloc[p]:
#            break
        if MS2.iloc[m,0]==PP7.iloc[p,0] and MS2.iloc[m,1]==PP7.iloc[p,1]:
            for q in range(len(Qb)):
                if flag:
                    break
#            if Qb['seq'].iloc[q]>PP7['seq'].iloc[p]:
#                break                
                if Qb.iloc[q,0] == PP7.iloc[p,0] and Qb.iloc[q,1] == PP7.iloc[p,1]:
#                    print(111)
                    com_seq.append(PP7.iloc[p,0])
                    MS2scores.append(MS2['score'].iloc[m])
                    PP7scores.append(PP7['score'].iloc[p])
                    Qbscores.append(Qb['score'].iloc[q])
                    Q.append(Qb.index[q])
                    M.append(MS2.index[m])
                    P.append(PP7.index[p])
                    CGC.append(Qb['prefix'].iloc[q])
                    flag=1
                    
All_Data = pd.DataFrame(com_seq)
All_Data.columns = ['seq']
All_Data['MS2'] = MS2scores
All_Data['PP7'] = PP7scores
All_Data['Qb'] = Qbscores
All_Data['Qidx'] = Q
All_Data['Midx'] = M
All_Data['Pidx'] = P
All_Data['CGC'] = CGC
All_Data.to_csv('temporary files\\common_sequences.csv')


