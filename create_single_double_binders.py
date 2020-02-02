
import pandas as pd
import os
# import seaborn as sns
# import matplotlib.pyplot as plt
import utiles
import optparse

# getOptions
#==============================================================================
def getOptions():
    parser = optparse.OptionParser()    

    parser.add_option("--RunNum", action="store", type="int", default=-1,
                      help = "the Run index of the program, if not specified "+
                      "the program will use the last run")  
 
    (options, args) = parser.parse_args()    
    return options  


options = getOptions()

if options.RunNum==-1:
    Run = pd.read_csv('temporary/RunNum.csv')
    RunNum = Run['Num'].iloc[0]   
    
#%% create sequences lists for validation experiments
#the list of sequences of double and single binders at the end is for proteins 
#p1 and p2, please change the names as you like.

DF1 = pd.read_csv('temporary/Predictions_from_run_{0}.csv'.format(RunNum))

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
    

if (not len(Q1)) or (not len(Q2)) or (not len(Q3)) or (not len(Q4))   :
    print('Please run the prediction process on larger sample size')

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

# utiles.createFolder('output_files')
doubleB.to_csv('output_files/'+p1+p2+'double_binders.csv',index=False)
singleP1.to_csv('output_files/'+p1+'_single_binder.csv',index=False)
singleP2.to_csv('output_files/'+p2+'_single_binder.csv',index=False)




