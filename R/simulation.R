## function to simulate data ##

#' Generates a simulated dataset based on a true partition
#'
#' @param true_L Vector of size max cluster size with number of clusters of each size   
#' @param nfields Number of fields
#' @param ncat Vector with number of categories per field
#' @param true_beta Distortion probability for the fields
#' @return Simulated data set
#' @export
#' @examples
#' truePartition <- c(2,2,2,2)
#' numberFields <- 2
#' numberCategories <- rep(5,2)
#' trueBeta <- 0.01
#' SimData(truePartition, numberFields, numberCategories, trueBeta)
SimData <- function(true_L, nfields, ncat, true_beta) {
	true_M <- length(true_L)
	true_K <- sum(true_L)
	N <- sum(c(1:true_M) * true_L) #100
	id <- rep(1:true_K, times = rep(1:true_M, times = true_L))
	## generate true entities
	y <- cbind()
	for (l in 1:nfields) {
		y <- cbind(y, sample.int(n = ncat[l], size = true_K, 
			replace = TRUE))
	}
	data <- matrix(NA, nrow = N, ncol = nfields)
	for (i in 1:N) {
		for (l in 1:nfields) {
			if (runif(1) < true_beta) {
				data[i, l] <- sample.int(n = ncat[l], 
					size = 1)
			} else {
				data[i, l] <- y[id[i], l]
			}
		}
	}
	x <- DataRemap(data)
	return(x)
}
