'''
this file is a class that trains and predict on all sequences that in distance Distances from every WT
notice the order of proteins! 
MS2 Qb PP7 
'''
import pickle
import pandas as pd
import numpy as np
from tensorflow.keras.models import Sequential, clone_model,load_model
# from keras.models import Sequential, clone_model
from tensorflow.keras.layers import Flatten,Conv1D,Dense
import os
import scipy.io as sio
import itertools
import math
#uncomment this to avoid GPU usage while training models
#os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
#os.environ["CUDA_VISIBLE_DEVICES"]="0"
import random
import utiles
class DisSeqs:
    
        
    def __init__(self,Lengths=[19,20,25],Distances=[3,4,5,6,7],
                 SeqsInBatch=20000, struct=1, SampleSize=2_000_000):
# =============================================================================
#       Initialize parameters
# =============================================================================
        self.Files = []                                                         #file list
        self.Lengths = Lengths                                                  #lengths of wild types 
        self.Distances = Distances                                             #distances of the sequences 
        self.OneHotDic= {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],
                         'T':[0,0,0,1],'U':[0,0,0,1]}                           #one hot dictionary
        RBPs = pd.read_csv('data/RBPs.csv')
        WTs = list(RBPs['wt'])
        self.WTDic = {}
        for wt in WTs:
            self.WTDic[len(wt)] = wt        
        
        self.Proteins = list(RBPs['name'])
#        self.Proteins = ['MS2','Qb','PP7']                                      #list of proteins to use
        
        with open('data/TSS.txt','r') as f:
            self.TSS = f.readline()

        with open('data/mc.txt','r') as f:
            self.mc = f.readline()
            
        self.SeqsInBatch = SeqsInBatch                                          #maximal number of sequences at a time
        self.MS2_ModelParams = {'epochs':40,'batchsize':16,'kernels':18,
                                'kernelL':6,'nuerons':20}        
        self.PP7_ModelParams = {'epochs':40,'batchsize':16,'kernels':11,
                                'kernelL':5,'nuerons':20}  
        self.Qb_ModelParams = {'epochs':40,'batchsize':16,'kernels':18,
                                'kernelL':5,'nuerons':20}  
        
        self.ModelParams = {'Qb':self.Qb_ModelParams,'MS2':self.MS2_ModelParams,
                            'PP7':self.PP7_ModelParams}
        
        # self.ModelParams = {'epochs':40,'batchsize':16,'kernels':18,
        #                     'kernelL':5,'nuerons':20} 
        # self.ModelParams = {'epochs':40,'batchsize':16,'kernels':18,
        #                     'kernelL':5,'nuerons':20}                           #model parameters
        self.SampleSize = SampleSize                                            #number of random sequences to test, -1=ALLsequences
        self.struct = struct
        self.CGCDic = {'MS2':'C','PP7':'C', 'Qb':'GC'}
        self.ALLseqs=0
        self.RandomIndexes=[]
        self.CurrIndex=0
        self.CurrRandIndex = 0
        self.seqNums=[]
        self.SequencNumber1()
        self.CreateRandomIndexes2()
#        print(self.RandomIndexes)
#==============================================================================

    def CreateFileList(self):
        #this function will create the file names list to create
        for L in self.Lengths:
            for d in self.Distances:
                self.Files.append('w{0}_d{1}'.format(L,d))
                
#============================================================================== 
                
    def SequencNumber(self):
        for L in self.Lengths:
            for d in self.Distances:
                self.ALLseqs += self.nCr(L,d)*3**d
                
                
    def SequencNumber1(self):
        
        tmpAll=0
        for L in self.Lengths:
            for d in self.Distances:
                tmpAll+= self.nCr(L,d)*3**d
            self.ALLseqs += tmpAll
            self.seqNums.append(tmpAll)
            tmpAll=0        
#============================================================================== 

    def nCr(self,n,r):
        f = math.factorial
        return f(n) / f(r) / f(n-r)     

#==============================================================================
    def CreateRandomIndexes(self):
        self.RandomIndexes = random.sample(range(int(self.ALLseqs)),self.SampleSize)
        self.RandomIndexes = sorted(self.RandomIndexes)
        
    def CreateRandomIndexes2(self):
        tmp=0
        for i in range(len(self.seqNums)):
            if not i:
                Rand = random.sample(range(int(self.seqNums[i])),int(self.SampleSize/(len(self.seqNums))))
            else:
                Rand = random.sample(range(tmp,int(tmp+self.seqNums[i])),int(self.SampleSize/(len(self.seqNums))))
            self.RandomIndexes = self.RandomIndexes+Rand
            tmp +=int(self.seqNums[i])
        
        self.RandomIndexes = sorted(self.RandomIndexes)
                    
#==============================================================================    

    #this function is from stuck overflow
    def polymorph(self,wt, n_mut):
        for locs in itertools.combinations(range(len(wt)), n_mut):
            this_word = [[char] for char in wt]
            for loc in locs:
                orig_char = wt[loc]
                this_word[loc] = [l for l in "ACGU" if l != orig_char]
            for poss in itertools.product(*this_word):
                yield ''.join(poss)  

#==============================================================================

    def CreateSeqs(self,distances=[1],WTs=[19]):
    # this function will generate all possible sequnces of distance d from WildType wt
        for wt in WTs:
            for d in distances:
                SeqsList = self.polymorph(self.WTDic[wt],d)
                with open('temporary/w{0}_d{1}.txt'.format(wt,d),'w')as f:
                    for i in SeqsList:
                        f.write(i+'\n')
                f.close()

#==============================================================================
                
    def BuildModel(self,protein):
        #build the model according to the given parameters
        Params = self.ModelParams[protein]
        model = Sequential()
        if self.struct:
            model.add(Conv1D(filters = Params['kernels'],
                             kernel_size =  Params['kernelL'],
                             strides = 1,padding='valid',
                             input_shape = (50,9),activation='relu'))
        else:
            model.add(Conv1D(filters = Params['kernels'],
                             kernel_size =  Params['kernelL'],
                             input_shape = (50,4),activation='relu'))                        
        model.add(Flatten())
        model.add(Dense(Params['nuerons'],activation='tanh'))
        model.add(Dense(1,activation='linear'))
#        self.model = model
        return model

#==============================================================================        

    def TrainModel(self,protein,model):
        Params = self.ModelParams[protein]
        #copy the model
        Model = clone_model(model)
        
        #read training data and 
        #get the data
        with open('temporary/DataForCNN_'+protein+'.txt','rb') as fp:
            FullData = pickle.load(fp)
        
        #get the scores   
        Y = FullData['Y']
        
        #reshape the data so the model could load it
        X_train = np.asarray(FullData['X'])
        X_train = np.swapaxes(X_train,1,2)
        Y_train = np.array(Y)
        
        #there is an option to not use the structur data
        if not self.struct:
            X_train = X_train[:,0:4,:]  
            
        X_train = np.swapaxes(X_train,2,1)
        #fit the model
        Model.compile(optimizer='adam',loss='mse')
        Model.fit(X_train,Y_train,batch_size = Params['batchsize'],
                  epochs=Params['epochs'])
        
        #save the trained model
        filename = protein+'_model.h5'
        Model.save(filename)
        # filename = protein+'_model.sav'
        # pickle.dump(Model, open(filename, 'wb'))
        
        return Model
        
#==============================================================================  
        
    def TrainModels(self,protein):
        #check if the input 'protein' is a list of proteins and train all the models needed
        
        if type(protein)==str:
            model = self.BuildModel(protein)
            self.TrainModel(protein,model)
        else:
            if type(protein==list):
                for p in protein:
                    model = self.BuildModel(p)
                    self.TrainModel(p,model)
                    
#==============================================================================  

    def ReadNLines(self,file):
        #this function reads self.SeqInBatch lines from file, if it reaches
        #the end of the file it return EOF=1 else EOF=0 
        
        EOF = 0
        i = 0
        lines = []
        for line in file:
            
            if self.CurrIndex == self.RandomIndexes[self.CurrRandIndex]:
                lines.append(line)
                i += 1
                self.CurrRandIndex +=1
                if self.CurrRandIndex==len(self.RandomIndexes):
                    file.close()
                    EOF=1
                    return lines,EOF

                if i==self.SeqsInBatch:
                    return lines,EOF
            self.CurrIndex+=1
            
        EOF = 1
        file.close()
     
        return lines,EOF            
#==============================================================================     
    def SaveTempFile(self,Data):
        #save a file to send to RNAfold
        file = 'temp.txt'
        with open(file,'w') as f:
            for i in Data:
                f.write(i+'\n')
        f.close()    
#==============================================================================        
    def OneHot(self,Data):
        #function to one hot encode the data
        OHdata=[]
        for line in Data:
            tmp = (line).split()
            oh = list(self.OneHotDic[char] for char in tmp[0])
            OHdata.append(oh)
        return OHdata
#==============================================================================

#    def SelfCheck(self):
#        Check=1
#        
    
#==============================================================================  
    def Predict(self,Train_Models=0):
        
        #that is the main function of this script
        #it either loads or train model for each protein, and test it
        #the test data is being loaded from pre-set files 
        #(can be created using the CreateSeqs function) then the data is passed 
        #to RNAfold, then to matlab script to transform the paranetasis to matrix 
        #the binding intensity is being predicted for each batch and saved as
        # csv file according to the batch number
        
        #read the run number
        Run = pd.read_csv('temporary/RunNum.csv')
        Run = Run['Num'].iloc[0]
        Run+=1 
        #save the run number for the next run
        tmpdf = pd.DataFrame([Run])
        tmpdf.columns = ['Num']
        tmpdf.to_csv('temporary/RunNum.csv')
        
        utiles.createFolder('temporary/predictions{0}'.format(Run))
        if Train_Models:
            self.TrainModels(self.Proteins)
            MS2model = load_model('MS2_model.h5')
            Qbmodel = load_model('Qb_model.h5')
            PP7model = load_model('PP7_model.h5')

            # MS2model = pickle.load(open('MS2_model.sav', 'rb'))
            # Qbmodel = pickle.load(open('Qb_model.sav', 'rb'))
            # PP7model = pickle.load(open('PP7_model.sav', 'rb'))

        else:
            MS2model = load_model('MS2_model.h5')
            Qbmodel = load_model('Qb_model.h5')
            PP7model = load_model('PP7_model.h5')

        #prepare the TSS and mc to go iinto the models
        L = len(self.TSS)
        OHTSS = list(self.OneHotDic[char] for char in self.TSS)
        OHmc = list(self.OneHotDic[char] for char in self.mc)     
        EOF = 1
        self.CreateFileList()
        BatchNumber = 0
        for Fidx in range(len(self.Files)):
            print('file number {0} begins'.format(Fidx))
            if EOF:
                file1 = open('temporary/'+self.Files[Fidx]+'.txt','r')
                EOF = 0
                
            while not EOF:
                BatchNumber += 1
                print('read data')
                Data,EOF = self.ReadNLines(file1)                                #read the Data to currently work with
                print(len(Data))
                if not Data:
                    continue
                print('save batch')
                self.SaveTempFile(Data)
                print('RNAfold')
                os.system('RNAfold -i temp.txt > parFile.txt')                  #RNAfold - saves the results to file named parFile
                # os.system('C:/Program Files (x86)/ViennaRNA Package/RNAfold -i temp.txt > parFile.txt')                  #RNAfold - saves the results to file named parFile

                print('MATLAB script')
                os.system(' matlab -nodisplay -nodesktop -r "run Struct.m;quit"')        #matlab script for matrix annotation
                mat = sio.loadmat('structure_matrix_tmp.mat')                #load the matrix
                mat = mat['struct_mat'][0,0]       
                SeqToSkip =list(map(int, list(pd.read_csv('Bugs.csv'))))        #read the indexes of the sequences that couldtnt be interpert
                #all indexes are compatible with python (py_idx = mat_idx -1)
                SeqToSkip = sorted(SeqToSkip, reverse=True)[:-1]
                #loop to remove un interpeted sequences
                for i in sorted(SeqToSkip, reverse=True):
                    del Data[i]
                print('encode Data')
                Seqs = Data
                #encode sequences
                Data = self.OneHot(Data)
                #concat sequnces to TSS and C or GC according to the better performing
                DataC   = []
                DataGC  = []
                print('organize data')
                for i in range(len(Data)):
                    StructMat = mat[i,0]
                    Lseq = len(Data[i])

                    OH = Data[i]
                    # C prefix - orginize data
                    L1=1
                    Mod = len(OH)%3
                    if Mod==0:
                        tmp = OHTSS+[self.OneHotDic['C']]+OH+[self.OneHotDic['C']]+[self.OneHotDic['U']]+OHmc
                    else:
                        if Mod==1:
                            tmp = OHTSS+[self.OneHotDic['C']]+OH+[self.OneHotDic['U']]+OHmc
                        else:
                            tmp = OHTSS+[self.OneHotDic['C']]+OH+OHmc
                            
                    tmp = tmp[L-10:L+40]
                    tmp = np.asarray(tmp)
                    if self.struct:
                        LeftZeros = np.zeros([5,10+L1])
                        RightZeros = np.zeros([5,50-10-L1-Lseq]) 
                        Struc = np.concatenate((LeftZeros,StructMat,RightZeros),axis=1)
                        FullMat = np.concatenate((tmp.T,Struc),axis=0)                
                        DataC.append(FullMat)
                    else:
                        DataC.append(tmp.T)
                    
                    # GC prefix - orginize data
                    L1=2
                    if Mod==0:
                        tmp = OHTSS+[self.OneHotDic['G']]+[self.OneHotDic['C']]+OH+[self.OneHotDic['U']]+OHmc
                    else:
                        if Mod==1:
                            tmp = OHTSS+[self.OneHotDic['G']]+[self.OneHotDic['C']]+OH+OHmc
                        else:
                            tmp = OHTSS+[self.OneHotDic['G']]+[self.OneHotDic['C']]+OH+[self.OneHotDic['C']]+[self.OneHotDic['U']]+OHmc                    
#                    tmp = OHTSS+[self.OneHotDic['G']]+[self.OneHotDic['C']]+OH+OHmc
                    tmp = tmp[L-10:L+40]
                    tmp = np.asarray(tmp)
                    if self.struct:
                        LeftZeros = np.zeros([5,10+L1])
                        RightZeros = np.zeros([5,50-10-L1-Lseq]) 
                        Struc = np.concatenate((LeftZeros,StructMat,RightZeros),axis=1)
                        FullMat = np.concatenate((tmp.T,Struc),axis=0)                
                        DataGC.append(FullMat)
                    else:
                        DataGC.append(tmp.T)
                
                #convert the data to a form that the model can predict on    
                DataC = np.asarray(DataC)
                DataC = np.swapaxes(DataC,1,2)
                DataGC = np.asarray(DataGC)
                DataGC = np.swapaxes(DataGC,1,2)
                print('predict on the data')
                #predict on all proteins
                MS2predictions = MS2model.predict(DataC)
                Qbpredictions  = Qbmodel.predict(DataGC)
                PP7predictions = PP7model.predict(DataC)
                #save the predictions
                print(str(len(Data))+'_'+str(len(MS2predictions)))
                DF = np.concatenate([MS2predictions,Qbpredictions,PP7predictions],axis=1)
                DF = pd.DataFrame(DF)
                DF.columns = self.Proteins
                DF['seq'] = Seqs
                DF.to_csv('temporary/predictions{0}/predictions{1}.csv'.format(Run,BatchNumber),index=False)
