# -*- coding: utf-8 -*-
"""
this script should run first, it will create all needed files for the rest 
of the program.
the files include the encoding of the data for preparing it to be used as 
an input for the nets
"""
import numpy as np
import pandas as pd
import utiles
import pickle
import scipy.io as sio
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

RBPs = pd.read_csv(dir_path+'/data/RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])
# utiles.createFolder('figures')
utiles.createFolder('temporary')

#change this to get the script to check the common sequences between the 3 libraries
check_common = False

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
    Data = pd.read_excel('data/Data.xlsx',key)
    D=[]
    for p in proteins:
        D.append(Data[p][Data['site length']==L])
    distances = pd.DataFrame(D).T
    distances.to_csv('temporary/dist_'+key+'.csv')
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
    
    Data = pd.read_excel('data/Data.xlsx',key)
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
    with open('temporary/DataForFC_'+key+'_C.txt','wb') as fp:
        pickle.dump(DataC,fp) 
    with open('temporary/DataForFC_'+key+'_GC.txt','wb') as fp:
        pickle.dump(DataGC,fp) 
           
#create data files for CNN
#============================================================================== 
L = 50      #preset value of the length of the input for the CNN


Dic={'PP7':0,'MS2':1,'Qb':2}
full_mat = sio.loadmat('data/structure_matrix.mat')

with open('data/TSS.txt','r') as f:
    TSS = f.readline()
with open('data/mc.txt','r') as f:
    mc = f.readline()
LTSS = len(TSS)
for key,val in WT.items():
    Data = pd.read_excel('data/Data.xlsx',key)
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
    with open('temporary/DataForCNN_'+key+'.txt','wb') as fp:
        pickle.dump(FULLData,fp)
    

#%% create a temporary file that will contatn the index of the Run 
tmpdf = pd.DataFrame([0])
tmpdf.columns = ['Num']
tmpdf.to_csv('temporary/RunNum.csv')
