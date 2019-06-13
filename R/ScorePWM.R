#' Score sequence similarity per HLA-gene
#'
#' Sum of amino acids frequencies at specific position in the gene-specific PWMs corresponding to sequences.
#'
#' @param sequences HLA-I amino acids sequences
#' @param PWMgene Position Weight Matrix for each HLA-gene
#'
#'
#' @return Score of sequences similarity with HLA-gene
#' @export
#'
#'
ScorePWM<-function(sequences, PWMgene){

  aa = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

  final_score <- sapply(PWMgene, function(p){
    if (is.null(dim(sequences))){
      all_score<-sum(sapply(1:length(sequences), function(x){
        score <- p[which(aa==sequences[x]),x]
        return(score)
      }))

    }else{
      all_score <-NULL
      for (i in 1:nrow(sequences)){
        score <-0
        for (j in 1:ncol(sequences)){
          score <- score + p[which(aa==sequences[i,j]),j]
        }
        all_score <- c(all_score, score)
      }
    }
  return(all_score)})
  return(final_score)
}
