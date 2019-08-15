
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


PLDpred <- function(sequences, allele=NULL, gene=NULL, output=NULL){

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
  data(PosCons)
  data(AACons)
  data(PWMgenes)

  sequences <- sequence_format(sequences)

  if (is.null(dim(sequences))){
    sequences <- sequence_definition(sequences)
  }else{
    sequences <- lapply(1:nrow(sequences), function(i){
    d <- sequence_definition(sequences[i,])
    return(d)
    })

    # HLA with no a1 and a2 sequences found.
    empty.seq<-sapply(1:length(sequences), function(w){
      if (is.null(sequences[[w]])){
        return(w)
      }else{
        return(0)
      }
    })

    if (sum(empty.seq)>0){
      warning(paste("Could not find a1 and a2 domaine in sequences",empty.seq,sep=" "))
      sequences<-do.call('rbind', sequences[[-empty.seq]])
    }else{
      sequences<-do.call('rbind', sequences)
      empty.seq <- NULL
    }
  }


  # Prediction gene from a1 and a2 domain
    data("PWMgenes")
    HLAgene<- c('A','B','C','G')

    if (is.null(dim(sequences))){
      score <-PLDpred::ScorePWM(sequences,PosWM)
      gene_pred <- HLAgene[which.max(score)]
    }else{
      score <- lapply(1:nrow(sequences), function(x){
        sc <- ScorePWM(sequences[x,],PosWM)
      })
      gene_pred<- sapply(score, function(d){
        g <- NULL
        g <- HLAgene[which.max(d)]
        return(g)
      })
    }

  # Gene definition from allele or gene (argument PLDpred function)

    if (!is.null(gene)){
      gene <- as.character(gene)
    }
    if (!is.null(allele) & is.null(gene)){
      allele <- as.character(allele)
      gene <- substr(allele,5,5)
    }
  if (!is.null(gene) | !is.null(allele)){
    if (all(gene==gene_pred)==F){
      warning(paste("Genes or Alleles gived are not corresponding to our gene prediction, possible mismatch between gene/allele arguments and sequences (PLD prediction are still based on PLDpred gene/allele arguments)"), sep="")
    }else{
      gene_pred <- NULL
    }
  }else{
    gene <- gene_pred
    gene_pred <- NULL
  }

  # Format sequence in 20 dimensions
  sequences <- Sequence_vector(sequences)

  #### Prediction ####
  prediction <- Prediction_LD(empty.seq, sequences, allele, gene, gene_pred)

  # save prediction in csv file
  if (!is.null(output)){
    write.csv(prediction, file=paste(output,".csv",sep=""),row.names = F)
  }
  return(prediction)
}
