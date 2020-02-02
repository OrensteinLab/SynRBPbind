'''
this file contain all general functions for this project
please note that the file should be in the same folder as the rest of the scripts


'''
import numpy as np
import pandas as pd
import os
import itertools
#import matplotlib.pyplot as plt
from tensorflow.keras.models import clone_model
import scipy.stats as stat

#create a folder
#==============================================================================
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
#==============================================================================
#this function is from stuck overflow - computes the distance of one string from the other
def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]
#==============================================================================
#this function is from stuck overflow - create all possible sequences with distance 
#n_mut from the sequence wt
def polymorph(wt, n_mut):
    for locs in itertools.combinations(range(len(wt)), n_mut):
        this_word = [[char] for char in wt]
        for loc in locs:
            orig_char = wt[loc]
            this_word[loc] = [l for l in "ACGU" if l != orig_char]
        for poss in itertools.product(*this_word):
            yield ''.join(poss)     
#==============================================================================
# def density_scatter( x , y, ax = None, sort = True, bins = 50, **kwargs )   :
#     """
#     Scatter plot colored by 2d histogram
#     """
#     if ax is None :
#         fig , ax = plt.subplots()
#     data , x_e, y_e = np.histogram2d( x, y, bins = bins)
#     z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False )

#     # Sort the points by density, so that the densest points are plotted last
#     if sort :
#         idx = z.argsort()
#         x, y, z = x[idx], y[idx], z[idx]

#     ax.scatter( x, y, c=z, **kwargs )
#     return ax    

#==============================================================================
def OneHot(Data):
    #function to one hot encode the data, input: list of strings ; output onehot matrix
    OneHotDic= {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],
                         'T':[0,0,0,1],'U':[0,0,0,1]}
    OHdata=[]
    for i in range(len(Data)):
        tmp = Data[i]
        oh = list(OneHotDic[char] for char in tmp)
        OHdata.append(oh)
    return OHdata 
#==============================================================================
def Creat_Di(seq,file_name):
    muts=[]
    for nuc1 in range(len(seq)):
        for nuc2 in range(len(seq)):
            if nuc1<nuc2:
                muts.append(seq[0:nuc1]+'A'+seq[nuc1+1:nuc2]+'A'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'A'+seq[nuc1+1:nuc2]+'C'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'A'+seq[nuc1+1:nuc2]+'G'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'A'+seq[nuc1+1:nuc2]+'U'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'C'+seq[nuc1+1:nuc2]+'A'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'C'+seq[nuc1+1:nuc2]+'C'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'C'+seq[nuc1+1:nuc2]+'G'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'C'+seq[nuc1+1:nuc2]+'U'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'G'+seq[nuc1+1:nuc2]+'A'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'G'+seq[nuc1+1:nuc2]+'C'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'G'+seq[nuc1+1:nuc2]+'G'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'G'+seq[nuc1+1:nuc2]+'U'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'U'+seq[nuc1+1:nuc2]+'A'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'U'+seq[nuc1+1:nuc2]+'C'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'U'+seq[nuc1+1:nuc2]+'G'+seq[nuc2+1::])
                muts.append(seq[0:nuc1]+'U'+seq[nuc1+1:nuc2]+'U'+seq[nuc2+1::])
    Muts = pd.DataFrame(muts)            
    Muts.columns = ['seq']            
    Muts.to_csv('temporary/'+file_name+'.csv',index=False)
    return Muts


#create SNPS
#==============================================================================
def Create_SNPS(seq,file_name):
    #get the seq - wildtype
    #the snps are saved in a csv file with the name that was given
    muts=[]
    for i in range(len(seq)):
        muts.append(seq[0:i]+'A'+seq[i+1::])
        muts.append(seq[0:i]+'C'+seq[i+1::])
        muts.append(seq[0:i]+'G'+seq[i+1::])
        muts.append(seq[0:i]+'U'+seq[i+1::])
    Allmuts = pd.DataFrame(muts)
    Allmuts.columns = ['seq']
    Allmuts.to_csv('temporary/'+file_name+'.csv',index=False)
    return Allmuts
#==============================================================================
def CreateROC(NScore,PScore):
    P=len(PScore)
    N=len(NScore)
    allScores = np.c_[np.concatenate((NScore,PScore)),
                      np.concatenate((np.zeros(N),np.ones(P)))]
    allScores = np.random.permutation(allScores)
    Sorted = allScores
    Sorted = np.array(sorted(Sorted, key=lambda a_entry: a_entry[0]))
    TN1=np.cumsum(Sorted[:,1]==0)
    TP1=np.cumsum(Sorted[:,1])
    AUC = round(np.trapz(np.array(TN1)/N,np.array(TP1)/P),3)

    return AUC
#==============================================================================
def N_fold_CV(X,Y,folds,model,params={'optimizer':'adam','epochs':1,'loss':'mse',
                                      'batch_size':16}):
    pos_thresh = 3.5
    neg_thresh = 0
    Rand = np.random.rand(len(Y),1).argsort(axis=0).ravel()
    #shuffle the data
    ShuffeledX = list(X[i] for i in Rand)
    ShuffeledY = list(Y[i] for i in Rand)    
    vec = np.round(np.linspace(0,len(Y)-1,folds+1)).astype(int)
    
    Pears = []
    AUC = []
    for i in range(len(vec)-1):
        curr_model = clone_model(model)
        #prepare the train and test data
        if( i!=0 and i!=(len(vec)-1)):
            X_test  = ShuffeledX[vec[i]:vec[i+1]]
            X_train = ShuffeledX[:vec[i]]+ ShuffeledX[vec[i+1]:]
            Y_test  = ShuffeledY[vec[i]:vec[i+1]]
            Y_train = ShuffeledY[:vec[i]]+ ShuffeledY[vec[i+1]:]
        else:
            if i==0:
                X_test  = ShuffeledX[vec[i]:vec[i+1]]
                X_train = ShuffeledX[vec[i+1]:]
                Y_test  = ShuffeledY[vec[i]:vec[i+1]]
                Y_train = ShuffeledY[vec[i+1]:]                       
            if i==len(vec)-1:
                X_test  = ShuffeledX[vec[i]:vec[i+1]]
                X_train = ShuffeledX[:vec[i]]
                Y_test  = ShuffeledY[vec[i]:vec[i+1]]
                Y_train = ShuffeledY[:vec[i]]        
    
        curr_model.compile(optimizer=params['optimizer'],loss=params['loss'])
        curr_model.fit(np.asarray(X_train),np.asarray(Y_train),batch_size = params['batch_size'],epochs=params['epochs'])    
        Predictions = curr_model.predict(np.asarray(X_test))
        Pearson = stat.pearsonr(Predictions.ravel(),np.array(Y_test))
        Pears.append(Pearson[0])
        #AUC
        pos_score = Predictions[np.array(Y_test)>=pos_thresh]
        Neg_score = Predictions[np.array(Y_test)<neg_thresh]
        AUC.append(CreateROC(Neg_score,pos_score))  
#        all_data['fold_{0}'.format(i)]={'Y_test':Y_test,'pred':Predictions}
    Mean_AUC = np.mean(AUC)
    Mean_Pearson = np.mean(Pears)
    return {'AUC':Mean_AUC,'pearson':Mean_Pearson}         
    
    
    
    
    
    
#==============================================================================