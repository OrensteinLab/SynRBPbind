import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import utiles
from DisSeqs_class import DisSeqs

#predict binding intensity for a sample of sequences
dis = DisSeqs(Lengths=[19,20,25],Distances=[3,4],SeqsInBatch=50_000,SampleSize=1_000_000)
dis.CreateSeqs(distances=[1,2],WTs=[19,20,25])
dis.Predict(0)

#%% combine the predictions files


RunNum = 5
List = os.listdir('predictions{0}'.format(RunNum))
i=0
for f in List:
   if not i:
       DF = pd.read_csv('predictions{0}\\'.format(RunNum)+f)
       i=1
   else:
       DF1 = pd.read_csv('predictions{0}\\'.format(RunNum)+f)
       DF = pd.concat([DF,DF1],axis=0)
       
       
#%% 
common = pd.read_csv('temporary files\\common_sequences.csv')

DF1=DF
p1 = 'Qb'
p2 = 'PP7'
p3= 'MS2'
thresh = 3.5
sns.set_style('white')
fig,ax = plt.subplots(dpi=300,figsize=(10,10))

m1 = min(DF1[p1])-1
m2 = min(DF1[p2])-1
M1 = max(DF1[p1])+1
M2 = max(DF1[p2])+1
#

cax = fig.add_axes([1.02, 0.13, 0.05, 0.75]) 
sns.kdeplot(DF1[p2],DF1[p1],cmap='Reds',shade=True,bw=0.5,ax=ax,cbar_ax=cax,cbar=True)

sns.scatterplot(common[p2],common[p1],alpha=0.2,ax=ax, size=1,legend=False)

ax.set_xlim([m2,M2])
ax.set_ylim([m1,M1])

ax.plot([m2,M2],[thresh,thresh],linestyle='dashed',color='olive')
ax.plot([thresh,thresh],[m1,M1],linestyle='dashed',color='olive')

plt.savefig(p1+p2+'map{0}.png'.format(RunNum),bbox_inches='tight')
plt.savefig(p1+p2+'map{0}.eps'.format(RunNum), format='eps', dpi=1000 ,bbox_inches='tight')
#%% create sequences lists for validation experiments
pos_thresh = 3.5
neg_thresh = 0
Nsequences = 20
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
ORGp1 = pd.read_excel('data files\\Data.xlsx',p1)
ORGp2 = pd.read_excel('data files\\Data.xlsx',p2)     

ORGp1['seq'] = ORGp1['seq'].astype(str).str.upper()
ORGp2['seq'] = ORGp2['seq'].astype(str).str.upper()
ORGp1 = list(ORGp1['seq'])
ORGp2 = list(ORGp2['seq'])

for i in range(len(double)):
    seq=double[i][:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        del double[i]
    double[i] = seq 
    
for i in range(len(singlep1)):
    seq=singlep1[i][:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        del singlep1[i]
    singlep2[i] = seq   
     
for i in range(len(singlep2)):
    seq=singlep2[i][:-1]
    if (seq in ORGp1) or (seq in ORGp2):
        del singlep2[i]
    singlep1[i] = seq

doubleB = pd.DataFrame(double)
singleP1 = pd.DataFrame(singlep1)
singleP2 = pd.DataFrame(singlep2)

utiles.createFolder('output files')
doubleB.to_csv('output files\\'+p1+p2+'double_binders.csv',index=False)
singleP1.to_csv('output files\\'+p1+'_single_binder.csv',index=False)
singleP2.to_csv('output files\\'+p2+'_single_binder.csv',index=False)






