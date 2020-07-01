"""
this script will perform whole-library models evaluation and hyper parameters selection
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

import scipy.stats as stat
import pickle
import numpy as np
import random
import utiles
import sys
random.seed(0)
RBPs = pd.read_csv('data/RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])
#TODO for check the code
proteins = [proteins[0]]

prot = sys.argv[1]
proteins = [prot]

#%% CNN
utiles.createFolder('temporary/whole_library_model_selection')
utiles.createFolder('temporary/whole_library_final_models')
Nsets = 10
Folds = 10

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
                
        print('Loading data\n')
        #load data 
        with open('temporary/DataForCNN_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)
        #orginize data
        #assign 20% as validation set
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
#==============================================================================        
        #Random parameters sets
        results = {}
        pears = []
        print('Random search:\n')
        for i in range (len(Sets)):
            print(f'set {i} out of {len(Sets)}/n')
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
            tmp = utiles.train_test_model(X_train = val_X,Y_train = list(val_Y),
                                          X_test =Data['X'] ,Y_test = list(Data['Y']),
                                          model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
            with open('temporary/whole_library_model_selection/'+library+
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
            for e in range(Set['epochs'][0]-10,Set['epochs'][0]+11,5):
                if e<1:
                    continue
                for kl in range(Set['kernels'][0]-5,Set['kernels'][0]+5):
                    if kl<1:
                        continue
                    for kw in range(Set['ker_length'][0]-3,Set['ker_length'][0]+3):
                        if kw<1:
                            continue
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
        print('Grid search:/n')
        for i in range (len(Sets)):
            print(f'set {i} out of {len(Sets)}/n')
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
        
            tmp = utiles.train_test_model(X_train = val_X,Y_train = list(val_Y),
                                          X_test =Data['X'] ,Y_test = list(Data['Y']),
                                          model = model,
                                   params={'optimizer':'adam','epochs':Set['epochs'][0],
                                           'loss':'mse','batch_size':16} )
            results['set{0}'.format(i)] = {'Set':Set,'meas':tmp}
            pears.append(tmp['pearson'])
        with open('temporary/whole_library_model_selection/'+library
                  +'_results_surround_parameters.txt','wb') as fb:    
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
        
        print('CNN {0} got pearson: {1} +- {2}, AUC: {3} +- {4}'
              .format(library,final_pears['pearson'],final_std['pearson'],
                      final_pears['AUC'],final_std['AUC'])) 
        print('\n the selected parameters are:\n')
        for key,val in Set.items():
            if type(val)==list:
                print('{0} : {1}, '.format(key,val[0]))
            else:
                print('{0} : {1}, '.format(key,val))
        #train model on all data
#==============================================================================
        print('Train the final model')
        with open('temporary/DataForCNN_'+library+'.txt','rb') as fb:
            Data = pickle.load(fb)
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
        model.save(library+'_model.h5')
        # filename = library+'_model.sav'
        # pickle.dump(model, open('temporary/whole_library_final_models/'+filename, 'wb'))

