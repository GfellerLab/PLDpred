#############

PLDpred is a predictor of peptides length distribution (PLD) of HLA-I based on HLA sequences using linear regression model.

###############

The prediction of ligand length distribution is computed with the function PLDpred, taking as input a vector or dataframe of HLA sequences. Important to mention that predictions are only available for HLA-gene A, B or C. HLA-alleles and HLA-genes are defined separately using arguments "allele" and "gene".

#### INSTALLATION ####

install.packages("devtools")
devtools::install_github("GfellerLab/PLDpred")
library(PLDpred)

#### USAGE #####

PLDpred(sequences, allele='ListOfAlleles', gene='ListOfGenes', output='path/filename')

* sequences : Vector or dataframe with HLA sequences where each element correspond to one amino acid. For dataframe, lines are HLA-allele and column are sequence position.  
* allele : list of characters (i.e.'HLA-A01:01).
* gene : list of characters, if allele unknown. ('A','B','C')

#### OUTPUT ####

PLDpred returns a vector or dataframe of peptides length distribution of HLA class I with sequence index and corresponding HLA-I allele, HLA-I gene or gene predicted (if allele=NULL and gene=NULL). Output could be saved as in csv format with filename defines in the output argument. 

#### Example ####

data(sequences.test)
PLDpred::PLDpred(sequences = sequences.test[,-1], allele=sequences.test[,1], gene=NULL, output="./test_PLDpred")


#### HELP ####

All details are in description of the function: 
?PLDpred::PLDpred


