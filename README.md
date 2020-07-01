The list of files necessary to run the program:
please note that the run time of this program depends on many variables
and can be between 30 minutes to a few days.
the fianl output of this program should be the output of seqs_predict.py,
this will be the 3 files with single and double binders of the proteins.
one more file type that is an output for this program is the model prediction
all SNPS with the same length as the wt of the protein.

installation:
------------
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
programs installed on your computer (unix):
	RNAfold (Vienna package)
	MATLAB
	python3.6 including tensorflow 2, keras, scipy, numpy, pandas, itertools,
	random, pickle

This is the order in which we recommend running the scripts:

1. files_creation.py - creates temporary files of the encoded data.
					  it uses the data in the data.xlsx file to create:
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
					  "python files_creation.py"
					  

Optional: 
2. whole_library_model_selection.py - tests models and hyper parameters for the whole-library models to find the optimal model for each protein
			This script takes a lot of time to complete.
			This can be skipped, a pre-trained model is loaded or a model is trained with a specified set of
			hyper parameters. This script takes as an input the protein name.

Optional: 			
3. WT_specific_model_selection.py - tests models and hyper parameters for the wt-specific models to find the optimal model for each protein
			This script takes a lot of time to complete.
			This can be skipped, a pre-trained model is loaded or a model is trained with a specified set of
			hyper parameters. This script takes as an input the protein name.
			To train the model with specified parameters or load our pre-trained model please read
			the section of the "whole_library" and "WT_specific"
						
4. whole_library.py - runs a class contained in DisSeqs_class.py
					it can train the models/ or predict using a pre-trained models.
                    It generates a csv file that containes all given sequences and predictions for them.
                    The training of all 3 models takes about 5 minutes.
					argumants:
					--train : 0 - don't tain the models, just load them, 1: train them using the
								   data in the file above
								   defult: 0
					--predict: 0 dont use the model for predictions, 1: predict on all 1-hamming
								  distance from WT sequecnes. 
								  defult: 0				
					--createseqs: create the sequnces that are within mind and maxd distance from thw WT
								  defult: 0
					--mind : 	   minimal hamming distance between the sequences that are tested and the WT.
									defult : 3
					--maxd : 	   maximum hamming distance between the sequences that are tested and the WT.
									defult:4

5. WT_specific.py - trains the WT-specific models (6 of them) and uses them to predict on all 1-hamming
						distance from WT sequecnes. 
						Each model provides prediction only for sequences with the same length as the
						original WT
						argumants:
						--trainfile : specify the directory+file name of the xlsx file containing the 
									  train data. 
									  defult: 'data/Data.xlsx'
						--train : 0 - don't tain the models, just load them, 1: train them using the
									   data in the file above
									   defult: 1
						--predict: 0 dont use the model for predictions, 1: predict on all 1-hamming
									  distance from WT sequecnes. 
									  defult: 1
					
6. create_single_double_binders.py - this will use the files created by the "whole_library" script after
									 runnind in prediction mode to seperate the sequences to double and 
									 single binders, it will also delete from the sequences list any
									 sequence that appears on the original data set.
									 argumants:
									 --RuunNum: the run index of the prediction process, so the program 
												will know what data to use.
												defult: last run.


How to run the program?

All scripts expect specific file names in the data folder.

1. "python file_creation.py"
   The script accept all files in the folder data
   and generates files in folder: temporary

2. optional:
	"python WT_specific_model_selection.py protein"
	The script select hyper parameters for WT_specific models of 'protein' 
	used in this project (for C and GC libraries).

3. optional:
	"python whole_library_model_selection.py protein"
	The script select hyper parameters for whole library model of 'protein' 
	used in this project.
	
4. "python whole_library.py --train 1 --predict 1 --createseqs 1"
   For train set the variable to ..., for predict run this way...
   The train data should be in folder data under the name Data.xlsx
   the prediction will be on synethetic sequecnes that will create during the program.


5. "python WT_specific.py"
	The script accept "trainfile", "train" and "predict" argumants.
	Creates files in the "output_files" folder containing all SNPS predictions
	(if --predict 0 is not specified)

6. "python create_single_double_binders.py"
	Creates files in the output_files folder containing single and double binders
