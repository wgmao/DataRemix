\name{CorMatToAUC}
\alias{CorMatToAUC}
\title{Evaluate meanAUC/meanAUPR of Recovering Known Pathways/Network}
\usage{
corMatToAUC(data, GS)
}                      
\arguments{
\item{data}{Correlation matrix with the dimension gene-by-gene}

\item{GS}{GS stands for gold standard which is a binary matrix with the dimension gene-by-pathway. The number of rows should match that of data and the genes should be in the same order.}
}

\value{
\item{}{A vector c(meanAUPR, meanAUC). By default, the biological objective is the average AUC acorss all pathways}
}
\description{
Given the task of recovering known pathways/network, calculate the average AUC value and average AUPR value based on gene correlation matrix and the gold stanard, e.g., the known pathway matrix.
}
