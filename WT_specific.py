
"""
this script will run the WT-specific models
if you will choose to train it it will train it,
if you will choose to predict on sequences using it
"""
import pandas as pd
from tensorflow.keras.layers import Flatten,Conv1D,Dense, Dropout
from tensorflow.keras.models import Sequential,load_model
#import scipy.stats as stat
import pickle
import numpy as np
import random
import utiles
import optparse

Qb = {'C':{'Nodes':39 ,'layers':1 },'GC':{'Nodes':39 ,'layers':1 }}
PP7 = {'C':{'Nodes':10 ,'layers':2 },'GC':{'Nodes':10 ,'layers':2 }}
MS2 = {'C':{'Nodes':38 ,'layers':1 },'GC':{'Nodes':32 ,'layers':1 }}
Parameters = {'Qb':Qb,'PP7':PP7,'MS2':MS2}

RBPs = pd.read_csv('data/RBPs.csv')
WTs = list(RBPs['wt'])
proteins = list(RBPs['name'])

WT = {}
for i in range(len(proteins)):
    WT[proteins[i]] = WTs[i]
    
Ldic = {}    
for key,val in WT.items():
    Ldic[key] = len(val)  
    
# getOptions
#==============================================================================
def getOptions():
    parser = optparse.OptionParser()    
    
    parser.add_option("--trainfile", action="store", type="str", default="data/Data.xlsx",
                      help = "directory and filename of the .xlsx file containing the data ")
    
    parser.add_option("--train", action="store", type="int", default=1,
                      help = "train the model? yes = 1, No = 0")  
    
    parser.add_option("--predict", action="store", type="int", default=1,
                      help = "predict on SNPS using the model?  yes = 1, No = 0")   
    
    (options, args) = parser.parse_args()    
    return options  


options = getOptions()

if options.train:
    for prot,Params in Parameters.items():
        #read the train data, seperate it to the appropriate lengths and encode it
        key=prot
        Data = pd.read_excel(options.trainfile,key)
        Data = Data[Data['site length']==Ldic[key]]
        Data['seq'] = Data['seq'].astype(str).str.upper()
        DataC = Data[Data['prefix']=='C']
        DataGC = Data[Data['prefix']=='GC']
        
        DataC = {'X':utiles.OneHot(list(DataC['seq'])),'Y':DataC['score']}
        DataGC ={'X':utiles.OneHot(list(DataGC['seq'])),'Y':DataGC['score']}
    
        with open('temporary/DataForFC_'+key+'_C.txt','wb') as fp:
            pickle.dump(DataC,fp) 
        with open('temporary/DataForFC_'+key+'_GC.txt','wb') as fp:
            pickle.dump(DataGC,fp) 
            
        #perform the training process           
        for cgc,Set in Params.items():    
            print('training: {0}_{1}'.format(prot,cgc))
            #load the training data
            with open('temporary/DataForFC_{0}_{1}.txt'.format(prot,cgc),'rb') as fb: 
                train_data = pickle.load(fb)       
            
            #build model
            model = Sequential()
            model.add(Flatten(input_shape = (Ldic[prot],4)))
            model.add(Dropout(0.2))
            for l in range(Set['layers']):
                 model.add(Dense(Set['Nodes'],activation='tanh'))
            
            model.add(Dense(1,activation= 'linear')) 
            model.compile(optimizer='adam',loss='mse')
            
            #train model
            X_train =  np.asarray(train_data['X']).astype(float)
            Y_train = np.asarray(list(train_data['Y']))
            
            model.fit(X_train, Y_train, batch_size = 8,epochs=70) 
                       
            model.save('{0}_{1}.h5'.format(prot,cgc))
        
if options.predict:
    for prot,Params in Parameters.items():
        #load all snps for prediction if needed
        Data = pd.read_csv('temporary/'+prot+' all snps.csv')
        Data_to_pred = utiles.OneHot(list(Data['seq']))
    
        for cgc,Set in Params.items():    
            print('predicting: {0}_{1}'.format(prot,cgc))
            #load the model
            model = load_model('{0}_{1}.h5'.format(prot,cgc))

            #predict using model
            pred = model.predict(np.asarray(Data_to_pred).astype(float))
            #save predictions
            Data[cgc] = pred
    
        Data.to_csv('output_files/SNPS_predictions_{0}.csv'.format(prot))             
            
