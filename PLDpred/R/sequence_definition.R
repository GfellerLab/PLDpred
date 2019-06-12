
#' Delimitation of HLA sequences for prediction
#'
#' Define a1 and a2 domains in sequence using 10 highly conserved positions.
#'
#' @param seq protein sequence of HLAs
#'
#' @return
#' @export

sequence_definition <- function(seq){
  data('PosCons')
  data('AACons')
  s <- NULL
  score=1
  for (i in 1:(length(seq)-180)){
    if (seq[i]=="S"){
      for (j in 1:length(PosCons)){
        if (seq[i+PosCons[j]] == AACons[j]){
          score = score+1
        }
      }
      if ((score/10) >= 0.7){
        s<- seq[i:(i+180)]
        break
      }
      score=1
    }
  }
  return(s)
}
