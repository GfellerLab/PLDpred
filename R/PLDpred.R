
#' Prediction of HLA class I ligand length distribution
#'
#' Prediction of ligand length distributions based on specific positions of HLA sequences. A Linear regression was trained independently for each HLA-gene (A,B and C) with HLA-peptidomics from mass spectrometry dataset of multi-allelic samples.
#' For that reason, HLA-gene must be specify for an accurate prediction. If HLA-gene is unknown, it will be predicted based on scores calculated from Gene-specific Position weigth matrix (PWM). PLDpred also restricts the HLA sequence to alpha1
#' and alpha2 domaines based on 10 highly conserved positions.
#'
#' @param sequences vector or dataframe (row = HLA and column = position) of HLA sequences.
#' @param allele list of character. Allele of HLA-I (e.g.'HLA-A01:01').
#' @param gene list of character. HLA-gene A, B or C.
#' @param output pathway and name for the final output (csv file).
#'
#'
#' @import foreach
#' @import Matrix
#' @import glmnet
#' @importFrom stats setNames
#' @importFrom utils write.csv
#' @importFrom utils data
#'
#' @return vector or dataframe with index sequence, allele/gene/predicted gene with predicted ligand length distribution  (Possibility to save prediction in csv format, see output argument).
#' @export


PLDpred <- function(sequences, allele=NULL, gene=NULL, output){

  allele <- as.character(allele)
  gene <- as.character(gene)

  # Setting
  aa = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  len <- sapply(8:14, function(x) paste('l',x,sep=''))

  # Function
  pred_glmnet_LOO <-function(fit_mul,test_sparse){
    pred<-list()
    for (i in 1:length(len)){
      pred[[i]]<- predict.cv.glmnet(fit_mul, t(as.matrix(test_sparse)), type="response", s="lambda.min")[,,1][i]
    }
    p=do.call(cbind.data.frame,pred)
    #p=setNames(cbind.data.frame(al,p), c('allele',len))
    return(p)
  }

  # Sequences definition
  data("PosCons")
  data("AACons")

  sequences <- sequence_format(sequences)

  if (is.null(dim(sequences))){
    sequences <- sequence_definition(sequences)
  }else{
    sequences <- lapply(1:nrow(sequences), function(i){
    d <- sequence_definition(sequences[i,])
    return(d)
    })
    empty.seq<-sapply(1:length(sequences), function(w){
      if (is.null(sequences[[w]])){
        return(w)
      }else{
        return(0)
      }
    }) # HLA with no a1 and a2 sequences found.

    if (sum(empty.seq)>0){
      print(paste("Warning: Could not find a1 and a2 domaine in sequences",empty.seq,sep=" "))
      sequences<-do.call('rbind', sequences[[-empty.seq]])
    }else{
      sequences<-do.call('rbind', sequences)
      empty.seq <- NULL
    }
  }

  # Gene definition

  HLAgene<- c('A','B','C','G')

  if(is.null(allele) & is.null(gene)){
    data("PWMgenes")
    if (is.null(dim(sequences))){
      score <-PLDpred::ScorePWM(sequences,PosWM)
      gene <- HLAgene[which.max(score)]
    }else{
      score <- lapply(1:nrow(sequences), function(x){
       sc <- ScorePWM(sequences[x,],PosWM)
      })
      gene <- sapply(score, function(d){
        g <- NULL
        g <- HLAgene[which.max(d)]
        return(g)
      })
    }
  }
  if (!is.null(allele) & is.null(gene)){
    gene <- substr(allele,5,5)
  }

  # Format sequence in 20 dimensions
  sequences <- Sequence_vector(sequences, pHF = NULL)

  #### Prediction ####
  prediction <- Prediction_LD(empty.seq, sequences, pred_sparse, allele, gene)

  # save prediction in csv file
  if (!is.null(output)){
    write.csv(prediction, file=paste(output,".csv",sep=""),row.names = F)
  }
  return(prediction)
}