list of files:
data files:
RBPs.csv -
	contain a list of N proteins and thie known binding site
Data.xlsx - 
	contains the experimental data that the models are based on,
	one tab for each protein, the tab contain columns of the following:
	seq,prefix,	suffix,site length,score,site structure,min free energy,
	and one column for each of the proteins that specify the edit distance of the sequence from the WT of the protein. 
	(see for example the supplumentry Data.xlsx file in the "data files" folder. 

structure_matrix.mat -
	mat file of array of N cells each of them contain another array of cells that
	each of them holds a matrix that encodes the structure of the sequence.
mc.txt  -
	the end of transcript site, string of nucleotides
TSS.txt
	the transcript start site, string of nucleotides

Examples for all files formation can be seen in the "data files" folder.	

scripts:
user that wish to use this program should run the files in the next order:
1. file_creation.py - create some temporary files of the encoded data.
2. model selection.py - Tests models and hyper parameters to find the optimal one.
3. SeqsPredict.py - this file runs a class that contained in "DisSeqs_class.py" 
				  it will create the validation sequences for the whole library model
				  and predict based their binding intensity based on the models that 
				  were choosen in "model selection.py" it will create a csv file that
				  containes all sequences and predictions
4. visual.py - 	create some of the graphs and visualize the data, in this file there is a usage of X,Y coordinates,
				this is specific for the WTs we used, please change that to fit your data

