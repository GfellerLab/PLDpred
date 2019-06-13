% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLDpred.R
\name{PLDpred}
\alias{PLDpred}
\title{Prediction of HLA class I ligand length distribution}
\usage{
PLDpred(sequences, allele = NULL, gene = NULL, output)
}
\arguments{
\item{sequences}{vector of amino acids or dataframe.}

\item{allele}{list of character. Allele of HLA-I (e.g.'HLA-A01:01').}

\item{gene}{list of character. HLA-gene A, B or C.}

\item{output}{pathway and name for the final output (csv file).}
}
\value{
vector or dataframe with index sequence, allele/gene/predicted gene with predicted ligand length distribution  (Possibility to save prediction in csv format, see output argument).
}
\description{
Prediction of ligand length distributions based on specific position of HLA sequences. A Linear regression was trained independently for each HLA-gene (A,B and C) with mass spectrometry dataset of multi-allelic samples.
For that reason, HLA-gene must be specify for an accurate prediction. If the HLA-gene is unknown, it will be predicted based on score calculated from Gene-specific Position weigth matrix (PWM). The score correspond to
the sum of PWM value corresponding to each element of the HLA sequences. The gene is assigned to the PWM highest score. PLDpred also restrict the HLA sequence to alpha1 and alpha2 domaines based on 10 highly conserved positions.
}