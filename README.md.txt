The list of files necessary to run the program:
please note that the run time of this program depends on many variables
and can be between 30 minutes to a few days.
the fianl output of this program should be the output of seqs_predict.py,
this will be the 3 files with single and double binders of the proteins.
one more file type that is an output for this program is the model predictionon
all SNPS with the same length as the wt of the protein.

installation:
there is no need to install the software but the dependencies specified below.

Under data folder:

RBPs.csv -
	contains a list of N proteins and their known wild type binding sites
	(one site for each protein) including their known RNA secondary structure

Data.xlsx - 
	contains the experimental data that the models are based on,
	one tab for each protein, the tab contain columns of the following:
		seq,prefix,suffix,site length,score,site structure,min free energy,
		and one column for each of the proteins that specify the edit distance
		of the sequence from the WT of the protein. 
		(see for example Data.xlsx file in the data folder.) 

structure_matrix.mat -
	(cell = type of variable in MATLAB)
	mat file with array of cells (length as the number of proteins),
	each of this cells contain array of cells - every one of holds 
	a matrix that encodes the structure of the sequence.

mc.txt  -
	the 5' flank of the binding site transcript.

TSS.txt
	the 3' flank of the binding site transcript.

Examples for all file formats are in the data folder.	

Code (python scripts):
dependencies: to run the program sucsefully you will need the following
programs installed on your computer:
	RNAfold (Vienna package)
	MATLAB
	python3.6 including tensorflow 2, keras, scipy, numpy, pandas, itertools,
	random, pickle
	
This is the order in which we recommend running the scripts:

1. file_creation.py - creates temporary files of the encoded data.
					  it uses the data in the data.xls file to create:
					  files of the edit distances (for plot of the histograms)
					  files the contain all SNPS of each WT in the RBPs.csv file
					  file with all Di-nucleotide mutations only for the MS2 protein.
					  this script create the file for the wt-specific models by:
							seperate the C and GC library
							one-hot encodes the data
							save it as a dictionary with the data and the binding scores.
					  the script also creates the data for the whole-library models by:
							adding to the sequence the mc and TSS, 
							cuting the sequence to length of 50
							adding the secondary structure in the file "structure_matrix.mat"
							padding with zeros the structure matrix of the added flanks.
					  it also creates a file containing the number zero (for later use and change
					  in the "DisSeq_class.py" file)

					  the input for this script is the data files in the data folder.
					  the output will be in the temporary folder. the outputs are files
					  with one hot encoded data for bothe the "whole-library" and "wt-specific" models.
					  to run this:
					  

Optional: 
2. model_selection.py - tests models and hyper parameters to find the optimal model for each protein
			This script takes a lot of time to complete.
			This can be skipped, a pre-trained model is loaded or a model is trained with a specified set of
			hyper parameters.
			To use our pre-trained CNN models you can use the program as is.
			if you wish to train the CNN models yourself you can change the "Train_Models" parameter in 
			"SeqsPredict.py" line 15 to "True"
						
3. seqs_predict.py - runs a class contained in DisSeqs_class.py
                    It creates the validation sequences for the whole library model
                    and predicts based their binding intensity based on the models that 
                    were choosen under model selection.py
                    It generates a csv file that containes all given sequences and predictions for them.
                    When runing the main method in this class ("Predict"), there will be an option
	                to train a model or load an existing model.
                    The training of all 3 models takes about 5 minutes.

4. Run_WT_specific.py - trainsthe WT-specific models (6 of them) and uses them to predict on all 1-hamming
						distance from WT sequecnes. 
						Each model provides prediction only for sequences with the same length as the original WT




How to run the program?

all scripts do not need any input from the command lind, just hit "play"
for example:
from the command line
1. "python file_creation.py"
2. "python seqs_predict.py"
3. "python Run_WT_specific.py"

