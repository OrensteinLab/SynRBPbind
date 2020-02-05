Please note that the run time of this program depends on many variables 
and can be between 30 minutes to a few days.

Installation
------------
There is no need to install the software but the dependencies are:
Programs installed on your Unix computer:
- RNAfold 2.4.13(Vienna package)
- MATLAB 2018
- python3.6.6 with packages: tensorflow 2.0.0, pandas 0.25.3.

Necessary files
---------------
Under 'data' folder (contains example of all files):

RBPs.csv -
	contains a list of N proteins and their known wild type binding sites
	(one site for each protein) including their known RNA secondary structure.

Data.xlsx - 
	contains the experimental data the models are based on,
	a single worksheet for each protein. The sheet contains the columns:
		seq,prefix,suffix,site length,score,site structure,min free energy,
		and one column for each of the proteins that specify the edit distance
		of the sequence from the WT of the protein. 

structure_matrix.mat -
	(cell = type of variable in MATLAB)
	a mat file with array of cells (length as the number of proteins),
	each of this cells contain array of cells - every one of holds 
	a matrix that encodes the structure of the sequence.
	[can be generated by Matlab script Struct.m]

mc.txt  -
	the 5' flank of the binding site transcript.

TSS.txt
	the 3' flank of the binding site transcript.


Sscript and arguments
---------------------
We recommend running the scripts in the following order.

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
2. model_selection.py - tests models and hyper parameters to find the optimal model for each protein.
			This script takes a lot of time to complete.
			This can be skipped, a pre-trained model is loaded or a model is trained with a specified set of
			hyper parameters.
			To train the model with specified parameters or load our pre-trained model please read the
			the section of the "whole_library" and "WT_specific"

						
3. whole_library.py - runs a class contained in DisSeqs_class.py
					it can train the models/ or predict using a pre-trained models.
                    It generates a csv file that containes all given sequences and predictions for them.
                    The training of all 3 models takes about 5 minutes.
					arguments:
					--train : 0 - don't tain the models, just load them, 1: train them using the
								   data in the file above
								   defult: 0
					--predict: 0 dont use the model for predictions, 1: predict on all 1-hamming
								  distance from WT sequences. 
								  defult: 0				
					--createseqs: create the sequnces that are within mind and maxd distance from thw WT
								  defult: 0
					--mind : 	   minimal hamming distance between the sequences that are tested and the WT.
									defult : 3
					--maxd : 	   maximum hamming distance between the sequences that are tested and the WT.
									defult:4

4. WT_specific.py - trains the WT-specific models (6 of them) and uses them to predict on all 1-hamming
						distance from WT sequences. 
						Each model provides prediction only for sequences with the same length as the
						original WT
						arguments:
						--trainfile : specify the directory+file name of the xlsx file containing the 
									  train data. 
									  defult: 'data/Data.xlsx'
						--train : 0 - don't tain the models, just load them, 1: train them using the
									   data in the file above
									   defult: 1
						--predict: 0 don't use the model for predictions, 1: predict on all 1-hamming
									  distance from WT sequences. 
									  defult: 1
					
5. create_single_double_binders.py - this will use the files created by the "whole_library" script after
									 runnind in prediction mode to seperate the sequences to double and 
									 single binders, it will also delete from the sequences list any
									 sequence that appears on the original data set.
									 arguments:
									 --RuunNum: the run index of the prediction process, so the program 
												will know what data to use.
												defult: last run.


How to run the program?
-----------------------
All scripts expect specific file names in the data folder.

1. "python files_creation.py"
   The script expects specific files in 'data' folder (see above)
   and generates training and testing files in folder 'temporary'

2. Optional:
	"python model_selection.py"
	The script selects hyper-parameters for all models used in this project.
	It expects specific files in 'temporary' folder (see above)
	and generates trained model files (in the same folder) and training output in 'temporary' folder.
	
3. "python whole_library.py --train 1 --predict 1 --createseqs 1"
   The script trains a model based on data files created in 'temporary' folder (see above)
   and predicts binding to sequences in specific files created in 'temporary' folder (see above).
   It ouputs predicted intensities to temporary/Predictions_from_run_<num>.csv

4. "python WT_specific.py --train 1 --predict 1 --trainfile data/Data.xlsx"
   The script trains a model based on WT-specific data files in 'temporary' folder (see above)
   and predict binding to all single nucleotide polymorphisms (SNPs).
   It creates files in the 'output_files' folder containing all SNPS predictions.

5. "python create_single_double_binders.py"
   It creates files in the 'output_files' folder containing single and double binders.
