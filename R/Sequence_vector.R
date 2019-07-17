
#' Sequence in 20 dimension vectors
#'
#' Create a one-hot encoding scheme of HLA sequences (per 20 amino-acids).
#'
#' @param Sequence vector or dataframe of HLA sequences.
#'
#' @import stringr
#'
#' @return dataframe where each line represent sequences coded as 20 dimensions vector per position and ligand length distribution for each allele if pHF not null
#' @export
Sequence_vector <- function(Sequence){
  colnames_Freq <- c()
  Frequency_vector <- NULL
  vector_vector_allele <- NULL
  aa <-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")


  # Sequence Vector dataframe

  if (is.null(dim(Sequence))){
    FreqV <- NULL
    lenseq <- length(Sequence)
    for (i in 1:lenseq){
      freq_pos<- rep(0, times=length(aa))
      freq_pos[which(Sequence[i]== aa)] <- freq_pos[which(Sequence[i]==aa)]+1
      FreqV <-  c(FreqV, freq_pos)
    }
    Frequency_vector <- rbind.data.frame(FreqV,FreqV)[1,]
  }else{
    for (j in 1:nrow(Sequence)){
      FreqV <- NULL
      lenseq <- ncol(Sequence)
      for (i in 1:ncol(Sequence)){
        freq_pos<- rep(0, times=length(aa))
        freq_pos[which(Sequence[j,as.numeric(i)]== aa)] <- freq_pos[which(Sequence[j,as.numeric(i)]==aa)]+1
        FreqV <-  c(FreqV, freq_pos)
      }
      Frequency_vector <- rbind.data.frame(Frequency_vector,FreqV)
    }
  }

  # Colnames of the matrix containing the position and the aa
  c=0
  for (i in 1:lenseq){
    for (a in 1:20){
      c=c+1
      tmp <- c(aa[a],i)
      colnames_Freq[c]<- str_c(tmp, collapse = "_")
    }
  }

  colnames(Frequency_vector) <- colnames_Freq

  # Control :
  # sequence length is respected:
  if (all(apply(Frequency_vector, 1, sum)!= lenseq)){
    stop('Sequence created does not have sequence length respected')
  }
  # each sequence is equal to vector of length position * 20
  if (ncol(Frequency_vector) != lenseq*20){
    stop('new vector size not correct')
  }

  # Reverse manipulation to see if the sequence are similar
  control <- c()
  for (j in 1:nrow(Frequency_vector)){
    s<- NULL
    for (i in seq(1,(lenseq)*20,20)){
      if (length(which(Frequency_vector[j,i:(i+19)] == 1))>0){
        s <- c(s, aa[which(Frequency_vector[j,i:(i+19)] == 1)])
      }else{
        stop('position is empty')
      }
    }
    if (is.null(dim(Sequence))){
      control <- all(Sequence==s)
    }else{
      control[j] <- all(Sequence[j,]==s)
    }
  }

  if (all(control!=T)){
    stop('Mistake in sequence, data not reversible')
  }
  #print(paste("Reversible data:", all(control)==TRUE, sep=" "))

  return(Frequency_vector)
}

