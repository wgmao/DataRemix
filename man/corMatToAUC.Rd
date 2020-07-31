\name{CorMatToAUC}
\alias{CorMatToAUC}
\title{Evaluate meanAUC/meanAUPR/medianAUC/medianAUPR of Recovering Known Pathways/Network}
\usage{
corMatToAUC(data, GS, objective = "mean.AUC")
}                      
\arguments{
\item{data}{Correlation matrix with the dimension gene-by-gene}

\item{GS}{GS stands for gold standard which is a binary matrix with the dimension gene-by-pathway. The number of rows should match that of data and the genes should be in the same order. This pathway matrix can be obtained from pathway databased, e.g., MSigDB}

\item{objective}{\itemize{\item "mean.AUC" (default option): Output is c(meanAUPR, meanAUC) and the objective is meanAUC (average AUC across all pathways)
\item "mean.AUPR"" : Output is c(meanAUC, meanAUPR) and the objective is meanAUPR
\item "median.AUC" : Output is c(medianAUPR, medianAUC) and the objective is medianAUC (median AUC across all pathways)
\item "median.AUPR" : Output is c(medianAUC, medianAUPR) and the objective is medianAUPR
}}
}



\value{
\item{}{A vector c(meanAUPR, meanAUC). By default, the biological objective is the average AUC acorss all pathways. The output will change based on objective}
}
\description{
Given the task of recovering known pathways/network, calculate the average AUC value and average AUPR value based on gene correlation matrix and the gold stanard, e.g., the known pathway matrix.
}
