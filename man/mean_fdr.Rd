% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_rates.R
\name{mean_fdr}
\alias{mean_fdr}
\title{Calculates average FDR when ground truth is available}
\usage{
mean_fdr(zm, id)
}
\arguments{
\item{zm}{Matrix with posterior samples of cluster assignments}

\item{id}{Vector of true cluster assignments (ground truth)}
}
\value{
Average FDR over posterior samples
}
\description{
Calculates average FDR when ground truth is available
}
\examples{
truePartition <- c(50,50,50,50)
maxPartitionSize<- length(truePartition)
uniqueNumberRecords <- sum(truePartition)
id <- rep(1:uniqueNumberRecords, times=rep(1:maxPartitionSize, times=truePartition))
numberFields <- 5
numberCategories <- rep(10,5)
trueBeta <- 0.01
simulatedData <- SimData(truePartition, numberFields, numberCategories, trueBeta)
posteriorESCD <- SampleCluster(data=simulatedData, Prior="ESCD", burn=0, nsamples=10)
mean_fdr(zm = posteriorESCD$Z, id)
}
