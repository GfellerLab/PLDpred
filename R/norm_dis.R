#' normalisation ligand length distribution
#'
#' Normalisation of ligand length distribution when prediction of one of ligang length is negative.
#'
#' @param x ligand length distribution to normalise
#'
#' @return Vector with ligand length distribution normalised.
#'
#' @export

norm_dis <- function(x){
  if (min(x)<0){
    x[which(x<0)] <- 0
    x <- (x/sum(x))
  }else{
    x<- (x/sum(x))
  }
  return(x)
}
