# -*- coding: utf-8 -*-
"""
this script will perform all models evaluation and hyper parameters selection
"""

 
import pandas as pd
from keras.layers import Flatten,Conv1D,Dense
from keras.models import Sequential
#import scipy.stats as stat
import pickle
import numpy as np
import random
import utiles
RBPs = pd.read_csv('data files\\RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])
#%% Fully connected:
utiles.createFolder('temporary files\\wt specific model selection')
utiles.createFolder('temporary files\\wt specific final models')
Nsets = 1
Folds = 2
#define parameters space
#==============================================================================
Nodes = list(range(5,50,5))
Layers = [1,2]
Activations = ['tanh','linear','relu']
epochs = list(range(20,110,10))
ParamsDic = {'Nodes':Nodes,'Layers':Layers,'Activations':Activations,'epochs':epochs}
#==============================================================================
for p in proteins:
    SNPS = pd.read_csv('temporary files\\'+p+' all snps.csv') 
    OH = utiles.OneHot(SNPS['seq'])
    OH = np.asarray(OH)
    for Prefix in ['C','GC']:

        library = p+'_'+Prefix
        
        #random sets
        Sets={}
        for i in range(Nsets):
            Set={}
            for key,val in ParamsDic.items(): 
                tmp = random.sample(val,1)
                Set[key] = tmp
            Sets['set{0}'.format(i)] = Set
                
         
        #load data 
        with open('temporary files\\DataForFC_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)
        #orginize data
            
        
        #Random parameters sets
        results = {}
        pears = []
        for i in range (len(Sets)):
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
        #TODO change this line so the function will get the parameters from here
            tmp = utiles.N_fold_CV(X = Data['X'],Y = list(Data['Y']),
                                   folds = Folds ,model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            #Set['epochs'][0]
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary files\\wt specific model selection\\'+library+
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
            for e in range(Set['epochs'][0]-15,Set['epochs'][0]+16,5):
                tmp = Set.copy()
                tmp['epochs'] = [e]
                tmp['Nodes'] = [n]
                Sets['set{0}'.format(i)] = tmp
                i+=1
        
        #10 fold for each one
        results = {}
        pears = []
        for i in range (len(Sets)):
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
        
            tmp = utiles.N_fold_CV(X = Data['X'],Y = list(Data['Y']),
                                   folds = Folds ,model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary files\\wt specific model selection\\'+library+
                  '_results_surround_parameters.txt','wb') as fb:    
            pickle.dump(results,fb)    
        #find maximum point
        idx = np.argmax(pears)
        Set = Sets['set{0}'.format(idx)]
        
        #train model on all data
#==============================================================================
        Sets={}
        Sets['set0'] = Set
        model = Sequential()
        model.add(Flatten())
        model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
        
        for l in range(Set['Layers'][0]):
             model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
             
        #add the output layer
        model.add(Dense(1,activation= 'linear'))
        model.compile(optimizer='adam',loss='mse')
        model.fit(np.asarray(Data['X']),np.asarray(list(Data['Y'])),batch_size = 8,epochs=Set['epochs'][0]) 
        #save the model
        filename = library+'_model.sav'
        pickle.dump(model, open('temporary files\\wt specific final models\\'+filename, 'wb'))
        
        #predict the SNPS binding intensity
        OH1 = np.asarray(OH)
        predictions = model.predict(OH1)
        DF = pd.DataFrame(predictions)
        DF.to_csv('temporary files\\'+p+Prefix+'_predictions.csv',index=False)
#%% CNN
utiles.createFolder('temporary files\\whole library model selection')
utiles.createFolder('temporary files\\whole library final models')
Nsets = 2
Folds = 2

#define parameters space
#==============================================================================
Nodes = list(range(5,50,5))
Layers = [1,2]
Activations = ['tanh','linear','relu']
epochs = list(range(20,110,10))
kernels = list(range(4,35))
ker_length = list(range(4,10))
ParamsDic = {'Nodes':Nodes,'Layers':Layers,'Activations':Activations,
             'epochs':epochs,'kernels':kernels,'ker_length':ker_length}
#==============================================================================
for p in proteins:
        library = p
        
        #random sets
        Sets={}
        for i in range(Nsets):
            Set={}
            for key,val in ParamsDic.items(): 
                tmp = random.sample(val,1)
                Set[key] = tmp
            Sets['set{0}'.format(i)] = Set
                
        
        #load data 
        with open('temporary files\\DataForCNN_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)
        #orginize data
            
#==============================================================================        
        #Random parameters sets
        results = {}
        pears = []
        for i in range (len(Sets)):
            Set = Sets['set{0}'.format(i)]
            ##build model
            model = Sequential()
            model.add(Conv1D(filters = Set['kernels'][0],
                             kernel_size = Set['ker_length'][0],
                             strides = 1,padding='valid',
                             input_shape = (50,9),
                             activation='relu'))   
            model.add(Flatten())
            for l in range(Set['Layers'][0]):
                 model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
                 
            #add the output layer
            model.add(Dense(1,activation= 'linear'))
        
            ##N fold function
            tmp = utiles.N_fold_CV(X = Data['X'],Y = list(Data['Y']),
                                   folds = Folds ,model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
            with open('temporary files\\whole library model selection\\'+library+
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
            for e in range(Set['epochs'][0]-15,Set['epochs'][0]+16,5):
                for kl in range(Set['kernels']-5,Set['kernels']+5):
                    for kw in range(Set['ker_length']-3,Set['ker_length']+3):
                        tmp = Set.copy()
                        tmp['kernels'] =[ kl]
                        tmp['ker_length'] =[ kw]
                        tmp['epochs'] = [e]
                        tmp['Nodes'] = [n]
                        Sets['set{0}'.format(i)] = tmp
                        i+=1
        
        #10 fold for each one
        results = {}
        pears = []
        for i in range (len(Sets)):
            Set = Sets['set{0}'.format(i)]
            ##build model
            model = Sequential()
            model.add(Conv1D(filters = Set['kernels'][0],
                             kernel_size = Set['ker_length'][0],
                             strides = 1,padding='valid',
                             input_shape = (50,9),
                             activation='relu'))   
            model.add(Flatten())
            for l in range(Set['Layers'][0]):
                 model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
                 
            #add the output layer
            model.add(Dense(1,activation= 'linear'))
        
            tmp = utiles.N_fold_CV(X = Data['X'],Y = list(Data['Y']),
                                   folds = Folds ,model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary files\\whole library model selection\\'+library
                  +'_results_surround_parameters.txt','wb') as fb:    
            pickle.dump(results,fb)    
        #find maximum point
        idx = np.argmax(pears)
        Set = Sets['set{0}'.format(idx)]
        
        #train model on all data
#==============================================================================
        Sets={}
        Sets['set0'] = Set
        model = Sequential()
        model.add(Flatten())
        model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
        
        for l in range(Set['Layers'][0]):
             model.add(Dense(Set['Nodes'][0],activation=Set['Activations'][0]))
             
        #add the output layer
        model.add(Dense(1,activation= 'linear'))
        model.compile(optimizer='adam',loss='mse')
        model.fit(np.asarray(Data['X']),np.asarray(list(Data['Y'])),
                  batch_size = 8,epochs=Set['epochs'][0]) 
        
        #save the model
#==============================================================================
        filename = library+'_model.sav'
        pickle.dump(model, open('temporary files\\whole library final models\\'+filename, 'wb'))