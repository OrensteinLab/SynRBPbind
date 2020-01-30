
"""
this script will run the WT-specific models
"""
import pandas as pd
from tensorflow.keras.layers import Flatten,Conv1D,Dense, Dropout
from tensorflow.keras.models import Sequential
#import scipy.stats as stat
import pickle
import numpy as np
import random
import utiles


Qb = {'C':{'Nodes':39 ,'layers':1 },'GC':{'Nodes':39 ,'layers':1 }}
PP7 = {'C':{'Nodes':10 ,'layers':2 },'GC':{'Nodes':10 ,'layers':2 }}
MS2 = {'C':{'Nodes':38 ,'layers':1 },'GC':{'Nodes':32 ,'layers':1 }}
Parameters = {'Qb':Qb,'PP7':PP7,'MS2':MS2}

for prot,Params in Parameters.items():
    #load all snps for prediction
    Data = pd.read_csv('temporary/'+prot+' all snps.csv')
    Data_to_pred = utiles.OneHot(list(Data['seq']))

    for cgc,Set in Params.items():    
        print('{0}_{1}'.format(prot,cgc))
        #load the training data
        with open('temporary/DataForFC_{0}_{1}.txt'.format(prot,cgc),'rb') as fb: 
            train_data = pickle.load(fb)
             
        
        #build model
        model = Sequential()
        model.add(Flatten())
        model.add(Dropout(0.2))
        for l in range(Set['layers']):
             model.add(Dense(Set['Nodes'],activation='tanh'))
        
        model.add(Dense(1,activation= 'linear')) 
        model.compile(optimizer='adam',loss='mse')
        
        #train model
        X_train =  np.asarray(train_data['X']).astype(float)
        Y_train = np.asarray(list(train_data['Y']))
        
        model.fit(X_train, Y_train, batch_size = 8,epochs=70) 
        
        # model.fit(np.asarray(train_data['X']),np.asarray(list(train_data['Y'])),
        #           batch_size = 8,epochs=70) 
        
        #predict using model
        pred = model.predict(np.asarray(Data_to_pred))
        #save predictions
        Data[cgc] = pred

    Data.to_csv('temporary/SNPS_predictions_{0}.csv'.format(prot))        
    

                