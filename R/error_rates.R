## functions to compute FNR and FDR ##

#' @title Calculates false discovery rate (FDR) when the ground truth is available
#' 
#' @description \loadmathjax False discovery rate (FDR) of the estimated record linkage (partition) based on the ground truth is defined as (Steorts, 2015)
#' \mjsdeqn{FDR = \frac{FP}{CL + FP}}
#' where FP is the number of false positives (not linked under the truth but linked under the estimate) and CL is the number of correct links (true positives). 
#' If both FP=0 and CL=0, define FDR = 0. 
#' 
#' FDR can be also defined as \mjseqn{FDR = 1 - Precision}, where \mjseqn{Precision = CL/(CL+FP)}.
#' 
#' \code{fdr_fun} calculates FDR for an estimated partition, and \code{mean_fdr} calculates average FDR based on posterior samples of partition.
#' 
#' @param z Vector of cluster assignments
#' @param id Vector of true cluster assignments (ground truth)
#' @param use_apply Logical (default F), whether to use \code{apply()} to calculate the rate. 
#' Setting \code{use_apply = T} may be slower but memory efficient when \code{length(z)} is large.
#' @return False discovery rate (FDR)
#' 
#' @references Steorts, R. C. (2015). Entity resolution with empirically motivated priors. Bayesian Analysis, 10(4), 849-875.
#' @seealso \code{\link{fnr_fun}}, \code{\link{mean_fnr}}
#' @export
#' @examples
#' nclusters_per_size <- c(50,50,50,50)
#' numberFields <- 5
#' numberCategories <- rep(10,5)
#' trueBeta <- 0.01
#' # generate simulated data
#' simulatedData <- SimData(nclusters_per_size, nfields = numberFields, ncat = numberCategories, true_beta = trueBeta)
#' # Fit ESCD model
#' posteriorESCD <- SampleCluster(data=simulatedData, Prior="ESCD", burn=5, nsamples=10)
#' # true number of clusters
#' trueK = sum(nclusters_per_size)
#' # true cluster membership vector
#' trueid = rep(1:trueK, times=rep(1:length(nclusters_per_size), times=nclusters_per_size))
#' # FDR calculation for a single estimate
#' fdr_fun(posteriorESCD$Z[1,], trueid)
#' # average FDR calculation
#' mean_fdr(posteriorESCD$Z, trueid)
#' 
fdr_fun <- function(z, id, use_apply = F) {
  if(!use_apply){ # use_apply = F, faster, but memory intensive
    estimatelinks = outer_equal(z)
    diag(estimatelinks) = F
    truelinks = outer_equal(id)
    diag(truelinks) = F
    denom = sum(estimatelinks)
    if(denom == 0){
      return(0)
    }else{
      return(sum(estimatelinks & !truelinks)/denom)
    }
  }else{ # use_apply = T, slower, but not memory intensive
    if (n_matches_fun(z) == 0) {
      return(0)
    } else {
      return(sum(vapply(X = c(1:length(z)), FUN = i_false_det_fun, 
                        FUN.VALUE = 1, z. = z, id. = id))/(2 * n_matches_fun(z)))
    }
  }
}


#' @title Calculates false negative rate (FNR) when the ground truth is available
#' 
#' @description \loadmathjax False negative rate (FNR) of the estimated record linkage (partition) based on the ground truth is defined as (Steorts, 2015)
#' \mjsdeqn{FNR = \frac{FN}{CL + FN}}
#' where FN is the number of false negatives (linked under the truth but not linked under the estimate) and CL is the number of correct links (true positives). 
#' If both FN=0 and CL=0, define FNR = 0. 
#' 
#' FNR can be also defined as \mjseqn{FNR = 1 - Recall}, where \mjseqn{Recall = CL/(CL+FN)}.
#'
#' \code{fnr_fun} calculates FNR for an estimated partition, and \code{mean_fnr} calculates average FNR based on posterior samples of partition.
#'
#' @param z Integer vector of cluster assignments
#' @param id Integer vector of true cluster assignments (ground truth)
#' @param use_apply Logical (default F), whether to use \code{apply()} to calculate the rate. 
#' Setting \code{use_apply = T} may be slower but memory efficient when \code{length(z)} is large.
#' @return False negative rate (FNR)
#' 
#' @references Steorts, R. C. (2015). Entity resolution with empirically motivated priors. Bayesian Analysis, 10(4), 849-875.
#' @seealso \code{\link{fdr_fun}}, \code{\link{mean_fdr}}
#' @export
#' @examples
#' nclusters_per_size <- c(50,50,50,50)
#' numberFields <- 5
#' numberCategories <- rep(10,5)
#' trueBeta <- 0.01
#' # generate simulated data
#' simulatedData <- SimData(nclusters_per_size, nfields = numberFields, ncat = numberCategories, true_beta = trueBeta)
#' # Fit ESCD model
#' posteriorESCD <- SampleCluster(data=simulatedData, Prior="ESCD", burn=5, nsamples=10)
#' # true number of clusters
#' trueK = sum(nclusters_per_size)
#' # true cluster membership vector
#' trueid = rep(1:trueK, times=rep(1:length(nclusters_per_size), times=nclusters_per_size))
#' # FNR calculation for a single estimate
#' fnr_fun(posteriorESCD$Z[1,], trueid)
#' # average FNR calculation
#' mean_fnr(posteriorESCD$Z, trueid)
#' 
fnr_fun <- function(z, id, use_apply = F) {
  if(!use_apply){ # use_apply = F, faster, but memory intensive
    estimatelinks = outer_equal(z)
    diag(estimatelinks) = F
    truelinks = outer_equal(id)
    diag(truelinks) = F
    denom = sum(truelinks)
    if(denom == 0){
      return(0)
    }else{
      return(sum(!estimatelinks & truelinks)/denom)
    }
  }else{ # use_apply = T, slower, but not memory intensive
    if (n_matches_fun(id) == 0) {
      return(0)
    } else {
      return(sum(vapply(X = c(1:length(z)), FUN = i_false_neg_fun, 
                        FUN.VALUE = 1, z. = z, id. = id))/(2 * n_matches_fun(id)))
    }  
  }
}


#' @rdname fdr_fun 
#' @param zm Matrix with posterior samples of cluster assignments, where each row corresponds to one sample from the posterior
#' @export
mean_fdr <- function(zm, id, use_apply = F) {
	fdr_vec <- apply(X = zm, MARGIN = 1, FUN = fdr_fun, 
		id = id, use_apply = use_apply)
	return(mean(fdr_vec))
}

#' @rdname fnr_fun 
#' @param zm Matrix with posterior samples of cluster assignments, where each row corresponds to one sample from the posterior
#' @export
mean_fnr <- function(zm, id, use_apply = F) {
  fnr_vec <- apply(X = zm, MARGIN = 1, FUN = fnr_fun, 
                   id = id, use_apply = use_apply)
  return(mean(fnr_vec))
}

# auxiliary functions to compute fnr and fdr

# faster version of base::outer(z, z, "==")
outer_equal<- function(z){ 
  z = as.integer(z)
  n = length(z)
  Y = rep.int(z, rep.int(n, n))
  robj <- z==Y
  dim(robj) = c(n,n)
  robj
}


i_false_det_fun <- function(i, z., id.) {
	i_match <- which(z.[-i] == z.[i])
	i_match_true <- which(id.[-i] == id.[i])
	i_common <- length(intersect(i_match, i_match_true))
	i_false_det <- length(i_match) - i_common
	return(i_false_det)
}

i_false_neg_fun <- function(i, z., id.) {
	i_match <- which(z.[-i] == z.[i])
	i_match_true <- which(id.[-i] == id.[i])
	i_common <- length(intersect(i_match, i_match_true))
	i_false_neg <- length(i_match_true) - i_common
	return(i_false_neg)
}

# calculates all possible (undirected) links in the given partition.
n_matches_fun <- function(z) {
	sizes <- table(z)
	n_match <- sum(vapply(X = sizes, FUN = function(s) {
		return(choose(s, 2))
	}, FUN.VALUE = 1))
	return(sum(n_match))
}
