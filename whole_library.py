import pandas as pd
import os
# import seaborn as sns
# import matplotlib.pyplot as plt
import utiles
from DisSeqs_class import DisSeqs
import optparse


# getOptions
#==============================================================================
def getOptions():
    parser = optparse.OptionParser()    
    
    # parser.add_option("--trainfile", action="store", type="str", default="data\Data.xlsx",
    #                   help = "directory and filename of the .xlsx file containing the data ")
    
    parser.add_option("--train", action="store", type="int", default=0,
                      help = "train the model? yes = 1, No = 0")  
    
    parser.add_option("--predict", action="store", type="int", default=0,
                      help = "predict on SNPS using the model?  yes = 1, No = 0")   

    parser.add_option("--createseqs", action="store", type="int", default=0,
                      help = "would you like to create the sequences with the"+
                      " specified hamming distances?  yes = 1, No = 0") 

    parser.add_option("--mindist", action="store", type="int", default=3,
                      help = "minimum hamming distance from wt")   
    
    parser.add_option("--maxdist", action="store", type="int", default=4,
                      help = "maximum hamming distance from wt")      
    (options, args) = parser.parse_args()    
    return options  


options = getOptions()

Distances = list(range(options.mindist,options.maxdist+1))


#predict binding intensity for a sample of sequences
#initialize the class
dis = DisSeqs(Distances=Distances,SeqsInBatch=50_000,SampleSize=1000)

#train the models if asked to
if options.train:
    dis.TrainModels()
    
#create files of the sequences with desired distances from the wt
if options.createseqs:
    dis.CreateSeqs(distances=Distances)
    
#train or load the models and predict on all sampled sequecnes.
if options.predict:
    dis.Predict()

    #%% combine the predictions files
    
    Run = pd.read_csv('temporary/RunNum.csv')
    RunNum = Run['Num'].iloc[0]
    
    
    List = os.listdir('temporary/predictions{0}'.format(RunNum))
    i=0
    for f in List:
       if not i:
           DF = pd.read_csv('temporary/predictions{0}/'.format(RunNum)+f)
           i=1
       else:
           DF1 = pd.read_csv('temporary/predictions{0}/'.format(RunNum)+f)
           DF = pd.concat([DF,DF1],axis=0)
    DF.to_csv('temporary/Predictions_from_run_{0}.csv'.format(RunNum),index=False)  
    print('Done fitst part \n')     

