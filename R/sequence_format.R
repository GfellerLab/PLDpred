#' Sequence format
#'
#' Correct the automatic conversion of "T" to TRUE or "F" to FALSE. It also change NA and empty position to X.
#'
#' @param sequences list of amino acids
#'
#' @return format HLA sequences
#' @export
#'
sequence_format<-function(sequences){
  sequences[sequences=="TRUE"] <- "T"
  sequences[sequences=="FALSE"] <- "F"
  sequences[is.na(sequences)==T] <- "X"
  sequences[is.null(sequences)==T] <- "X"
  #colnames(sequences)[1] <- "allele"
  return(sequences)
}
