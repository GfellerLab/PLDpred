
#' Prediction of HLA class I ligand length distribution
#'
#' Prediction_LD allows the prediction of HLA-gene if it is unknown based on position weight matrix.
#'
#' @param empty.seq HLA without a1 and a2 domains found in sequence.
#' @param sequences HLA sequences (vector or dataframe)
#' @param allele Allele of HLA (e.g.'HLA-A01:01')
#' @param gene Gene of HLA
#'
#' @import foreach
#' @import Matrix
#' @import glmnet
#' @importFrom stats setNames
#'
#' @return dataframe with gene/allele and the corresponding ligand length distribution predicted with linear regression model (Save in csv format).
#' @export


Prediction_LD <- function(empty.seq, sequences, allele, gene, gene_pred){

  len <- sapply(8:14, function(x) paste("l",x,sep=''))
  HLAgene <- c('A','B','C')

  # transform sequence in sparse matrix
  pred_sparse<- sparse.model.matrix(~., sequences)

  #### Prediction ####

  data(model_lm)
  data(model_pos)

  prediction_LR_gene <- function(hg){
    if (gene[hg] == 'G'){
      h <- 3
    }else{
      h <- which(HLAgene==gene[hg])
    }

    pos_cor <- unlist(model_pos[[h]])

    prediction <- t(data.frame(predict.cv.glmnet(model_lm[[h]], data.matrix(sequences[hg,model_pos[[h]]]), type="response", s="lambda.min")[,,1]))

    prediction<- t(apply(prediction,1,norm_dis)) # Normalisation to 1

    rownames(prediction)<-c()

    return(prediction)
  }

  if (is.null(dim(sequences))){
    prediction<-prediction_LR_gene(1)
  }else{
    prediction <- t(sapply(1:nrow(sequences), prediction_LR_gene))
  }
  prediction <- setNames(as.data.frame(prediction), len)


  if(is.null(empty.seq)){
    index = 1:nrow(sequences)
  }else{
    index = 1:nrow(sequences)[-empty.seq]
  }

  if (is.null(allele)){
    prediction <- cbind.data.frame(index=index,Gene=gene,prediction)
  }else{
    prediction <- cbind.data.frame(index=index,allele=allele,prediction)
  }

  if (!is.null(gene_pred)){
    prediction <- cbind.data.frame(index=prediction[,1], GenePred = gene_pred, prediction[,-1])
  }

  return(prediction)
}
