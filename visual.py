
"""
this file create all the visualizations of the papaer
1. histograms of distances
1. RNA plots:
    a. all plots in the same scale
    b. all plots in different scales
    
2. bar plots for comparison between different models:
    a. pearson correlation
    b. AUC

3. logo and mutations maps????
"""
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as stat
import utiles
utiles.createFolder('figures')
#
#import seaborn as sns
#from sklearn.linear_model import LogisticRegression
#from Gtools import letterAt
#from sklearn.feature_selection import chi2
#from keras.models import clone_model
#from keras import Sequential
#from keras.layers import Dense, Dropout, TimeDistributed, Conv1D, MaxPooling1D,Flatten,Conv2D
#from keras.callbacks import EarlyStopping
#import scipy.stats as stat
#from sklearn.model_selection import train_test_split


# normalization of the matrix:
###############################################################################     
def Norm_Mat(Predicted1):
    tmp = Predicted1>0
    if tmp.any():
        P = Predicted1.copy()    
        P[P<0]=0
        Mx = np.max(np.max(P))
        Mn = np.min(np.min(P[P!=0]))
        P = (P-Mn)/Mx/2
        P[P<0]=0
        P = P+0.5
        P[P==0.5]=0
    else:
        P = np.zeros(Predicted1.shape)
    tmp = Predicted1<0
    if tmp.any():
        Ne=Predicted1.copy() 
        Ne[Ne>0]=0
        Ne = np.abs(Ne)
        Mx = np.max(np.max(Ne))
        Mn = np.min(np.min(Ne[Ne!=0]))
        Ne = (Ne-Mn)/Mx/2
        Ne[Ne<0]=0
        Ne = 0.5-Ne
        Ne[Ne==0.5]=0
    else:
        Ne = np.zeros(Predicted1.shape) 
        
    mat = np.array(P)+np.array(Ne)

    return mat
#this function plots the RNA 
#==============================================================================  
def plotRNA(x,y,Scores,Fname,Wild):
    l = 0.6
    FSZ=20
#    Scores1 = Scores
    OrgScores = Scores.copy()
    Min = np.min(np.min(np.min(OrgScores)))
    Max = np.max(np.max(np.max(OrgScores)))
    Scores = Norm_Mat(Scores)
    im = plt.imshow(OrgScores,cmap = 'bwr_r')
    plt.show()
    colormap=matplotlib.cm.get_cmap('bwr_r')
    fig = plt.figure(num=None, figsize=(20,20), dpi=200, facecolor='w', edgecolor='k')

    for i in range(len(x)):
        #nucleotide A
        plt.fill([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color=colormap(Scores[0,i]))
        plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color='k')
        #nucleotide C
        plt.fill([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color=colormap(Scores[1,i]))
        plt.plot([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')        
        #nucleotide G
        plt.fill([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color=colormap(Scores[2,i]))
        plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color='k')
        #nucleotide U
        plt.fill([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color=colormap(Scores[3,i]))
        plt.plot([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')   

        if Wild[i]=='A':
            if Scores[0,i]>0.7 or Scores[0,i]<0.2 :
                plt.text(x[i],y[i]+l/4,r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i],y[i]+l/4,r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')
                
        if Wild[i]=='C':
            if Scores[1,i]>0.7 or Scores[1,i]<0.2 :
                plt.text(x[i]+l/4,y[i],r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i]+l/4,y[i],r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold') 
                
        if Wild[i]=='G':
            if Scores[2,i]>0.7 or Scores[2,i]<0.2 :
                plt.text(x[i],y[i]-l/4,r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i],y[i]-l/4,r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')
                
        if Wild[i]=='U':
            if Scores[3,i]>0.7 or Scores[3,i]<0.2 :
                plt.text(x[i]-l/4,y[i],r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')    
            else:
                plt.text(x[i]-l/4,y[i],r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')    
    #legend    
    i=0
    x=[10]
    y=[7.5]
    FSZ1=12
    #boxes
    plt.plot([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color='k')
    #letters
    plt.text(x[i],y[i]+l/4,'A',fontsize=FSZ1,ha = 'center',va='center')
    plt.text(x[i]+l/4,y[i],'C',fontsize=FSZ1,ha = 'center',va='center')
    plt.text(x[i],y[i]-l/4,'G',fontsize=FSZ1,ha = 'center',va='center')
    plt.text(x[i]-l/4,y[i],'U',fontsize=FSZ1,ha = 'center',va='center')
    
    plt.xlim([0,max(x)+5])
    plt.ylim([0,max(x)+5])
    plt.axis('off')
    plt.title(Fname,fontsize = 25)
    cax = fig.add_axes([0.85, 0.5, 0.02, 0.2]) 
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_ticks([Min,(Min+Max)/2,Max])
    cbar.set_ticklabels([float(format(Min,'.1g')),0,float(format(Max,'.4g'))])
    plt.savefig(Fname+'.png')
    plt.show()


def MiniPlotRNA(x,y,Scores,Wild):
    l = 0.6
    FSZ=18  
    colormap=matplotlib.cm.get_cmap('bwr_r')
    #ha = 'center',va='center'
    for i in range(len(x)):
        #nucleotide A
        plt.fill([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color=colormap(Scores[0,i]))
        plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color='k')
        #nucleotide C
        plt.fill([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color=colormap(Scores[1,i]))
        plt.plot([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')        
        #nucleotide G
        plt.fill([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color=colormap(Scores[2,i]))
        plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color='k')
        #nucleotide U
        plt.fill([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color=colormap(Scores[3,i]))
        plt.plot([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')   

              
        if Wild[i]=='A':
            if Scores[0,i]>0.7 or Scores[0,i]<0.2 :
                plt.text(x[i],y[i]+l/4,r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i],y[i]+l/4,r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')
                
        if Wild[i]=='C':
            if Scores[1,i]>0.7 or Scores[1,i]<0.2 :
                plt.text(x[i]+l/4,y[i],r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i]+l/4,y[i],r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold') 
                
        if Wild[i]=='G':
            if Scores[2,i]>0.7 or Scores[2,i]<0.2 :
                plt.text(x[i],y[i]-l/4,r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')
            else:
                plt.text(x[i],y[i]-l/4,r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')
                
        if Wild[i]=='U':
            if Scores[3,i]>0.7 or Scores[3,i]<0.2 :
                plt.text(x[i]-l/4,y[i],r'$\bullet$',fontsize=FSZ,color='white',ha = 'center',va='center',weight = 'bold')    
            else:
                plt.text(x[i]-l/4,y[i],r'$\bullet$',fontsize=FSZ,color='k',ha = 'center',va='center',weight = 'bold')    
                
#%% plot both C and GC results on the same scale
dic = {'MS2':19, 'Qb':20, 'PP7':25}
proteins = ['MS2','Qb','PP7']
CGC = ['C','GC']
Lengths = [19,20,25]
x_pos_dic = {'MS2':0.55, 'Qb':0.63, 'PP7':0.7}
wt19 = 'ACAUGAGGAUCACCCAUGU'
wt25 = 'UAAGGAGUUUAUAUGGAAACCCUUA'
wt20 = 'AUGCAUGUCUAAGACAGCAU'
# predict on all snps
import pickle

#for p in proteins:
#    SNPS = pd.read_csv('temporary files\\'+p+' all snps.csv') 
#    for cgc in ['C','GC']:
#        Model = pickle.load(open('temporary files\wt specific final models\\'+p+'_'+cgc+'_model.sav','rb'))
#        OH = utiles.OneHot(SNPS['seq'])
#        OH = np.asarray(OH)
#        predictions = Model.predict(OH)
#        DF = pd.DataFrame(predictions)
#        DF.to_csv('temporary files\\'+p+cgc+'_predictions.csv',index=False)
    
    
n=7
X_25 = [1,2,3,4,5,6,7,8,9,10,11,12,13.25,13.25,12,11,10,9,8,7,5,4,3,2,1]
Y_25 = [3+n,3+n,3+n,3+n,3+n,3.5+n,3+n,3+n,3+n,3+n,3.75+n,3.25+n,3+n,2+n,1.75+n,
        1.25+n,2+n,2+n,2+n,2+n,2+n,2+n,2+n,2+n,2+n]

X_20 = [2,3,4,5,6,7,8,9,10,11,12,11,10,9,8,7,5,4,3,2]
Y_20 = [3+n,3+n,3+n,3+n,3.5+n,3+n,3+n,3+n,3+n,3.5+n,2.5+n,1.5+n,2+n,2+n,2+n,2+n
        ,2+n,2+n,2+n,2+n]

X_19 = [1,2,3,4,5,6,7,8,9,10,10,9,8,7,5,4,3,2,1]
Y_19 = [3+n,3+n,3+n,3+n,3+n,3.5+n,3+n,3+n,3.5+n,3+n,2+n,1.5+n,2+n,2+n,2+n,2+n,
        2+n,2+n,2+n]
l = 0.6
FSZ=7   



tmp = []            
AllData = {}           
for p in proteins:
    for cgc in CGC:  
        L = dic[p] 
        Data = np.array(pd.read_csv('temporary files\\'+p+cgc+'_predictions.csv'))
        tmp.append(Data) 
        AllData [p+cgc] = Data   

for p in proteins:
    tmp=[]
    
    for cgc in CGC:
        L = dic[p]
        Data = np.array(pd.read_csv('temporary files\\'+p+cgc+'_predictions.csv'))
        tmp.append(Data)         
    
    TMP = np.concatenate(tmp,axis=0)
    P = TMP.copy()
    P[P<0]=0
    Mx = np.max(np.max(P))
    Mn = np.min(np.min(P[P!=0]))
    P = (P-Mn)/Mx/2
    P[P<0]=0
    P = P+0.5
    P[P==0.5]=0
    MX = Mx
    Ne=TMP.copy() 
    Ne[Ne>0]=0
    Ne = np.abs(Ne)
    Mx = np.max(np.max(Ne))
    Mn = np.min(np.min(Ne[Ne!=0]))
    Ne = (Ne-Mn)/Mx/2
    Ne[Ne<0]=0
    Ne = 0.5-Ne
    Ne[Ne==0.5]=0            
    
    mat = np.array(P)+np.array(Ne)
    
    
    NormAllData = {}    
    first = 0   
    for cgc in CGC:    
        L = dic[p] 
        NormAllData[p+cgc] = mat[first:first+len(AllData[p+cgc])]
        first = first+len(AllData[p+cgc])
            
    l=1
    FSZ=20   
    
    im = plt.imshow(TMP,cmap = 'bwr_r')
    plt.show()
    colormap=matplotlib.cm.get_cmap('bwr_r')
    fig = plt.figure(num=None, figsize=(20,20), dpi=300,
                     facecolor='w', edgecolor='k')
    disty = 0
    for cgc in CGC:
        L = dic[p]
        distx = 0
        
        #prepare the matrix for plot:
        MatToPlot = NormAllData[p+cgc]
        MatToPlot = np.asmatrix(MatToPlot[:,0]) 
        MatToPlot = np.reshape(MatToPlot,[-1,4]).T 
        if L==19:
            MiniPlotRNA(np.array(X_19)+distx,np.array(Y_19)+disty,MatToPlot,wt19)
        if L==20:
            MiniPlotRNA(np.array(X_20)+distx,np.array(Y_20)+disty,MatToPlot,wt20)            
        if L==25:
            MiniPlotRNA(np.array(X_25)+distx,np.array(Y_25)+disty,MatToPlot,wt25)            
                
#            distx = distx+14
        disty = disty -4
      
    #at the end of all plots       
    #legend
    l=2
    Mn = TMP.min()    
    i=0
    x=[12]
    y=[-2]
    plt.plot([x[i],x[i]-l/2,x[i]-l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]-l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]+l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]-l/2,y[i]],color='k')
    plt.plot([x[i],x[i]-l/2,x[i]+l/2,x[i]],[y[i],y[i]+l/2,y[i]+l/2,y[i]],color='k')
    plt.text(x[i]-0.05,y[i]+l/4,'A',fontsize=FSZ)
    plt.text(x[i]+l/4,y[i]-0.05,'C',fontsize=FSZ)
    plt.text(x[i]-0.1,y[i]-l/4-0.1,'G',fontsize=FSZ)
    plt.text(x[i]-l/4-0.1,y[i]-0.05,'U',fontsize=FSZ)
    plt.xlim([0,20])
    plt.ylim([-5,15])
    #plt.axis('off')
    plt.yticks([5.33,9.33],['GC','C'],fontsize=20)

#    plt.yticks([1.33,5.33,9.33],['PP7',r'Q$\beta$','MS2'],fontsize=20)
#    plt.xticks([5,20,35],['Wild Type 19','Wild Type 20','Wild Type 25'],fontsize=20)
    plt.title(p+'C and GC',fontsize = 25)
    cax = fig.add_axes([x_pos_dic[p], 0.5, 0.01, 0.2]) 
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_ticks([Mn,(Mn+MX)/2,MX])
    cbar.set_ticklabels([float(format(Mn,'.3g')),0,float(format(MX,'.3g'))])
    plt.savefig('figures\\'+p+'_library.png')
    plt.show()            

