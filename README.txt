#############

PLDpred is a predictor of peptides length distribution (PLD) of HLA-I based on HLA sequences using linear regression model.

###############

The prediction of ligand length distribution is computed with the function PLDpred, taking as input a vector or dataframe of HLA sequences. Important to mention that predictions are only available for HLA-gene A, B or C. HLA-alleles and HLA-genes are defined separately using arguments "allele" and "gene". If HLA-gene is unknown, it will be predict based on gene-specific position weight matrix score. 

#### INSTALLATION ####

install.packages("devetools")
devtools::install_github("GfellerLab/PLDpred")

#### USAGE #####

sequences <- read.table('HLA/sequence/file')
PLDpred(sequences, allele='ListOfAlleles', gene='ListOfGenes', output='path/filename')

* sequences : Vector or dataframe with the HLA sequences were each element correspond to one amino acid. For dataframe, lines are HLA and column are sequence position.  
* allele : list of characters (i.e.'HLA-A01:01).
* gene : list of characters, if allele unknown. ('A','B','C')

#### OUTPUT ####

PLDpred returns a vector or dataframe of peptides length distribution of HLA class I with sequence index and corresponding HLA-I allele, HLA-I gene or gene predicted (if allele=NULL and gene=NULL). Output could be saved as in csv format with filename defines in the output argument. 

#### HELP ####

All details are in description of the function: 
?PLDpred::PLDpred


