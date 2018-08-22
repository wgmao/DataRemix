\name{CorMatToAUC}
\alias{CorMatToAUC}
\title{Network Improvement}
\usage{
corMatToAUC(data, GS)
}                      
\arguments{
\item{data}{Correlation matrix with the dimension gene-by-gene}

\item{GS}{GS stands for gold standard which is a binary matrix with the dimension gene-by-pathway. The number of rows should match that of data.}
}

\value{
\item{}{A vector c(meanAUPR, meanAUC). By default, the biological objective is the average AUC acorss all pathways}
}
\description{
Network Improvement
}
