# -*- coding: utf-8 -*-
"""
this script will perform WT-specific models evaluation and hyper parameters selection
to change the number of sets tested in this section and the number of folds in the CV
pleas see variables:
    # Nsets
    # Folds
these variables are located before the wt-specific model selction (lines 37-38) 


"""

import os
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import pandas as pd
from tensorflow.keras.layers import Flatten,Conv1D,Dense
from tensorflow.keras.models import Sequential
# from keras.layers import Flatten,Conv1D,Dense
# from keras.models import Sequential
import scipy.stats as stat
import pickle
import numpy as np
import random
import utiles
import sys

RBPs = pd.read_csv('data/RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])
#TODO for check the code
proteins = [proteins[0]]
#%% Fully connected:
utiles.createFolder('temporary/wt_specific_model_selection')
utiles.createFolder('temporary/wt_specific_final_models')
Nsets = 10
Folds = 10
#define parameters space
#==============================================================================
Nodes = list(range(5,50,5))
Layers = [1,2]
Activations = ['tanh','linear','relu']
epochs = list(range(20,110,10))

ParamsDic = {'Nodes':Nodes,'Layers':Layers,'Activations':Activations,'epochs':epochs}
prot = sys.argv[1]
proteins = [prot]


#==============================================================================
for p in proteins:
    
    SNPS = pd.read_csv('temporary/'+p+' all snps.csv') 
    OH = utiles.OneHot(SNPS['seq'])
    OH = np.asarray(OH)


    for Prefix in ['C','GC']:
        random.seed(0)
        np.random.seed(0)
        #random sets
        Sets={}
        for i in range(Nsets):
            Set={}
            for key,val in ParamsDic.items(): 
                tmp = random.sample(val,1)
                Set[key] = tmp
            Sets['set{0}'.format(i)] = Set
        library = p+'_'+Prefix
        print(f'{library}\n')

                
        print('Loading data\n') 
        #load data 
        with open('temporary/DataForFC_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)
        #orginize data
        #assign 10% as validation set
        Ndata = len(Data['X'])
        idxs = random.sample(range(Ndata),int(Ndata/5*4))
        idxs.sort(reverse=True)
        YY = list(Data['Y'])
        XX = Data['X']
        
        val_X = []
        val_Y = []
        for i in idxs:
            val_X.append(XX[i])
            val_Y.append(YY[i])
            del YY[i]
            del XX[i]
        
        Data['X'] = XX
        Data['Y']  = pd.Series(YY)
        #Random parameters sets
        
        results = {}
        pears = []
        print('Random search:\n')
        for i in range (len(Sets)):
            print(f'set {i} out of {len(Sets)}/n')
            Set = Sets['set{0}'.format(i)]
            ##build model
            model = Sequential()
            model.add(Flatten())
            model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
            
            for l in range(Set['Layers'][0]):
                 model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
                 
            #add the output layer
            model.add(Dense(1,activation= 'linear'))
        
            ##N fold function
            tmp = utiles.train_test_model(X_train = val_X,Y_train = list(val_Y),
                                          X_test =Data['X'] ,Y_test = list(Data['Y']),
                                          model = model,
                                          params={'optimizer':'adam','epochs':Set['epochs'][0],
                                                  'loss':'mse','batch_size':16} )
            #Set['epochs'][0]
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary/wt_specific_model_selection/'+library+
                  '_results_random_parameters.txt','wb') as fb:    
            pickle.dump(results,fb)
        #find maximum point
        idx = np.argmax(pears)
        Set = Sets['set{0}'.format(idx)]
        
        #surrounding parameters sets
#==============================================================================
        Sets={}
        i=0
        
        for n in range(Set['Nodes'][0]-5,Set['Nodes'][0]+6):
            if n<1:
                continue
            for e in range(Set['epochs'][0]-15,Set['epochs'][0]+16,5):
                if e<1:
                    continue
                tmp = Set.copy()
                tmp['epochs'] = [e]
                tmp['Nodes'] = [n]
                Sets['set{0}'.format(i)] = tmp
                i+=1
        
        #10 fold for each one
        results = {}
        pears = []
        print('Grid search:/n')
        for i in range (len(Sets)):
            print(f'set {i} out of {len(Sets)}/n')
            Set = Sets['set{0}'.format(i)]
            ##build model
            model = Sequential()
            model.add(Flatten())
            model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
            
            for l in range(Set['Layers'][0]):
                 model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
                 
            #add the output layer
            model.add(Dense(1,activation= 'linear'))
        
            ##N fold function
            tmp = utiles.train_test_model(X_train = val_X,Y_train = list(val_Y),
                                  X_test =Data['X'] ,Y_test = list(Data['Y']),
                                  model = model,
                                  params={'optimizer':'adam','epochs':Set['epochs'][0],
                                          'loss':'mse','batch_size':16} )
            # tmp = utiles.N_fold_CV(X = Data['X'],Y = list(Data['Y']),
            #                        folds = Folds ,model = model,
            #                        params={'optimizer':'adam','epochs':Set['epochs'][0],
            #                                'loss':'mse','batch_size':16} )
            
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary/wt_specific_model_selection/'+library+
                  '_results_surround_parameters.txt','wb') as fb:    
            pickle.dump(results,fb)    
        #find maximum point
        idx = np.argmax(pears)
        Set = Sets['set{0}'.format(idx)]
        
        #validate the model
#==============================================================================
        print('Performance evaluation over 80%/n')
        Sets={}
        Sets['set0'] = Set
        model = Sequential()
        model.add(Flatten())
        model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
        
        for l in range(Set['Layers'][0]):
             model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
             
        #add the output layer
        model.add(Dense(1,activation= 'linear'))
        
        final_pears,final_std = utiles.N_fold_CV1(X = val_X,Y = val_Y
                                       ,folds = Folds,model = model,
                                       params={'optimizer':'adam',
                                               'epochs':Set['epochs'][0],
                                               'loss':'mse','batch_size':16}) 
        # model.compile(optimizer='adam',loss='mse')
        # model.fit(np.asarray(Data['X']),np.asarray(list(Data['Y'])),batch_size = 8,epochs=Set['epochs'][0]) 

        # Predictions = model.predict(np.asarray(val_X))
        
        print('CNN {0} got pearson: {1} +- {2}, AUC: {3} +- {4}'
              .format(library,final_pears['pearson'],final_std['pearson'],
                      final_pears['AUC'],final_std['AUC'])) 

        print('\n the selected parameters are:\n')
        for key,val in Set.items():
            if type(val)==list:
                print('{0} : {1}, '.format(key,val[0]))
            else:
                print('{0} : {1}, '.format(key,val))
        #train and save the model for future use
#=============================================================================        
        print('Train the final model')
        Sets={}
        print('train the full model - WT')
        Sets['set0'] = Set
        with open('temporary/DataForFC_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)        
        model = Sequential()
        model.add(Flatten())
        model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
        
        for l in range(Set['Layers'][0]):
             model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
             
        #add the output layer
        model.add(Dense(1,activation= 'linear'))
        model.compile(optimizer='adam',loss='mse')
        model.fit(np.asarray(Data['X']),np.asarray(list(Data['Y'])),batch_size = 8,epochs=Set['epochs'][0],verbose=0) 
       
        #save the model
        filename = library+'_model.h5'
        model.save(filename)
        # filename = library+'_model.sav'
        # pickle.dump(model, open('temporary/wt_specific_final_models/'+filename, 'wb'))
        
        #predict the SNPS binding intensity
        OH1 = np.asarray(OH)
        predictions = model.predict(OH1)
        DF = pd.DataFrame(predictions)
        DF.to_csv('temporary/'+p+Prefix+'_predictions.csv',index=False)
