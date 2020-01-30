import pandas as pd
import os
# import seaborn as sns
# import matplotlib.pyplot as plt
import utiles
from DisSeqs_class import DisSeqs

plot = False

#predict binding intensity for a sample of sequences
#initialize the class
dis = DisSeqs(Lengths=[19,20,25],Distances=[3,4],SeqsInBatch=50_000,SampleSize=500)
#create files of the sequences with desired distances from the wt
dis.CreateSeqs(distances=[3,4],WTs=[19,20,25])
#train or load the models and predict on all sampled sequecnes.
dis.Predict(Train_Models=False)

print('Done fitst part \n')

#%% combine the predictions files

#change the RunNum according to the one tou need
RunNum = 1


List = os.listdir('temporary/predictions{0}'.format(RunNum))
i=0
for f in List:
   if not i:
       DF = pd.read_csv('temporary/predictions{0}/'.format(RunNum)+f)
       i=1
   else:
       DF1 = pd.read_csv('temporary/predictions{0}/'.format(RunNum)+f)
       DF = pd.concat([DF,DF1],axis=0)
DF.to_csv('temporary/Predict.csv',index=False)       

#%% create sequences lists for validation experiments
#the list of sequences of double and single binders at the end is for proteins 
#p1 and p2, please change the names as you like.

DF1 = pd.read_csv('temporary/Predict.csv')

p1 = 'Qb'
p2 = 'PP7'
p3= 'MS2'

#threshold for binders and non binders
pos_thresh = 3.5
neg_thresh = 0
#how many sequences you like
Nsequences = 20
#minimum edit distance between two sequences
mindist = 4
Q1=[]
Q2=[]
Q3=[]
Q4=[]
for i in range(len(DF1)):
    if DF1[p1].iloc[i]>pos_thresh and DF1[p2].iloc[i]>pos_thresh and DF1[p3].iloc[i]<neg_thresh:
        Q1.append(DF1['seq'].iloc[i])
    if DF1[p1].iloc[i]>pos_thresh and DF1[p2].iloc[i]<neg_thresh and DF1[p3].iloc[i]<neg_thresh:
        Q2.append(DF1['seq'].iloc[i])
    if DF1[p1].iloc[i]<neg_thresh and DF1[p2].iloc[i]<neg_thresh and DF1[p3].iloc[i]<neg_thresh:
        Q3.append(DF1['seq'].iloc[i])
    if DF1[p1].iloc[i]<neg_thresh and DF1[p2].iloc[i]>pos_thresh and DF1[p3].iloc[i]<neg_thresh:
        Q4.append(DF1['seq'].iloc[i])
    
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

double=[Q1[0]]
for seq in Q1:
    flag=0
    for seq1 in double:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1
    if not flag:
        double.append(seq)
    
singlep1=[]
for seq in Q2:
    flag=0
    for seq1 in double:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1    
    for seq1 in singlep1:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1    
    if not flag:
        singlep1.append(seq)
        
singlep2=[]
for seq in Q4:
    flag=0
    for seq1 in double:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1    
    for seq1 in singlep1:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1  
    for seq1 in singlep2:
        if utiles.levenshteinDistance(seq,seq1)<mindist:
            flag=1             
    if not flag:
        singlep2.append(seq)        
        
#delete the seqs that in the library
ORGp1 = pd.read_excel('data/Data.xlsx',p1)
ORGp2 = pd.read_excel('data/Data.xlsx',p2)     

ORGp1['seq'] = ORGp1['seq'].astype(str).str.upper()
ORGp2['seq'] = ORGp2['seq'].astype(str).str.upper()
ORGp1 = list(ORGp1['seq'])
ORGp2 = list(ORGp2['seq'])

Double=[]
for seq in double:
    seq = seq[:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        continue
    Double.append(seq)

Singlep1=[]
for seq in singlep1:
    seq = seq[:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        continue
    Singlep1.append(seq)

Singlep2=[]
for seq in singlep2:
    seq = seq[:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        continue
    Singlep2.append(seq)

doubleB = pd.DataFrame(Double)
singleP1 = pd.DataFrame(Singlep1)
singleP2 = pd.DataFrame(Singlep2)

utiles.createFolder('output_files')
doubleB.to_csv('output_files/'+p1+p2+'double_binders.csv',index=False)
singleP1.to_csv('output_files/'+p1+'_single_binder.csv',index=False)
singleP2.to_csv('output_files/'+p2+'_single_binder.csv',index=False)




