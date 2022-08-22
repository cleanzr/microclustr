#' @useDynLib microclustr
#' @importFrom Rcpp sourceCpp
#' @importFrom stats rbinom rgamma runif

NULL

## Wrapper function ##
# Inputs: N - the sample size based on the dataset           #
#         Prior - Specify Prior Type: ESCNB, DP,PY, ESCP, ESCB, ESCBshift          #
#         param - Parameter for prior type specified         #
# Output: Parameters for probability of existing cluster (A) #
#         and new cluster (B) for reseating algorithms       #
SetParam <- function(Prior, param, N) {
	if (Prior == "DP") {
		theta <- param[1]
		A <- seq(N)
		B <- rep(theta, N)
	}
	if (Prior == "PY") {
		theta <- param[1]
		delta <- param[2]
		A <- seq(N) - delta
		B <- rep(theta, N) + seq(N) * delta
	}
	if (Prior == "ESCNB") {
		rnb <- param[1]
		pnb <- param[2]
		beta <- ((1 - pnb)^rnb)/(1 - ((1 - pnb)^rnb))
		A <- (seq(N) + rnb)
		B <- (seq(N) + 1) * beta * rnb
	}
  if (Prior == "ESCP"){
    lambda <- param[1]
    A <- rep(lambda,N)
    B <- (seq(N) + 1) * lambda*exp(-lambda)/(1-exp(-lambda)) # can replace lambda with 1, for both A and B
  }
  if (Prior == "ESCB"){
    Nbinom <- param[1]
    pbinom <- param[2]
    temp = (1-pbinom)^Nbinom
    A <- c(Nbinom - seq(Nbinom), rep(0,N-Nbinom))
    B <- (seq(N) + 1) * Nbinom * temp/(1-temp)
  }
  if (Prior == "ESCBshift"){
    Nbinom <- param[1]
    pbinom <- param[2]
    temp = (1-pbinom)^Nbinom
    temp2 = (seq(Nbinom)+1)/(seq(Nbinom))*(Nbinom - seq(Nbinom)+1)*pbinom/(1-pbinom)
    A <- c(temp2, rep(0,N-Nbinom))
    B <- (seq(N) + 1) * temp
  }
	return(list(A = A, B = B))
}

SetParamESCD <- function(mus, N, M) {
	A <- (seq(N) + 1) * (c((mus[-1])/mus[-(M + 1)],
		rep(0, N - M)))
	B <- (seq(N) + 1) * mus[1]
	return(list(A = A, B = B))
}


## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD, ESCP, ESCB, ESCBshift      #
# Output: Lower and Upper bound on support of the distribution #
#         of Prior parameters                                  #

lowupSlice <- function(Prior) {
	if (Prior == "DP") {
		lo <- 0
		up <- Inf
	}
	if (Prior == "PY") {
		lo <- c(0, 0)
		up <- c(Inf, 1)
	}
	if (Prior == "ESCNB") {
		lo <- c(0, 0)
		up <- c(Inf, 1)
	}
	if (Prior == "ESCD") {
		lo <- c(0, 0, 0)
		up <- c(Inf, Inf, 1)
	}
  if (Prior == "ESCP"){
    lo <- 0
    up <- Inf
  }
  if (Prior == "ESCB"||Prior == "ESCBshift"){
    lo <- c(0,0)
    up <- c(Inf,1)
  }
	return(list(lo = lo, up = up))
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD, ESCP, ESCB, ESCBshift      #
#         N - number of record in data                          #
#         Nbinom: For ESCB, the maximum size of cluster, for ESCBshift, the initial value of Nbinom.
# Output: Initial values for hyperpameters                      #

SetInit <- function(Prior, N, Nbinom = NULL) {
  if (Prior == "DP") {
    #concentration parameter
    theta <- N/2
    x1 <- as.vector(theta)
  }
  if (Prior == "PY") {
    #concentration parameter
    theta <- N/2
    # discount parameter
    beta <- 0.5
    x1 <- c(theta, beta)
  }
  if (Prior == "ESCNB") {
    rnb <- 1
    pnb <- 0.5
    x1 <- c(rnb, pnb)
  }
  if (Prior == "ESCD") {
    alpha <- 1
    rnb <- 1
    pnb <- 0.5
    x1 <- c(alpha, rnb, pnb)
  }
  if (Prior == "ESCP"){
    lambda <- 1
    x1 <- lambda
  }
  if (Prior == "ESCB"||Prior == "ESCBshift"){
    Nbinom <- Nbinom
    pbinom <- 0.1
    x1 <- c(Nbinom,pbinom)
  }
  return(x1)
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, , ESCNB, ESCD, ESCP, ESCB, ESCBshift #
# Output: Indicator for hyperparameter sampling             #

SamInd <- function(Prior) {
	if (Prior == "DP") {
		#concentration parameter
		samind <- c(1)
	}
	if (Prior == "PY") {
		# concentration parameter
		# discount parameter
    samind <- c(1, 1)
	}
	if (Prior == "ESCNB") {
		# NB parameters r and p
		samind <- c(1, 1)
	}
	if (Prior == "ESCD") {
		# NB parameters r and p, alpha fixed
		samind <- c(0, 1, 1)
	}
  if (Prior == "ESCP"){
    # poisson
    samind <- c(1)
  }
  if (Prior == "ESCB"||Prior == "ESCBshift"){
    # binomial
    samind <- c(0,1) # Nbinom fixed (ESCBshift will not use slice sampling; ignore)
  }
	return(samind)
}


## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD, ESCP, ESCB, ESCBshift       #
#         N - number of record in data                          #
# Output: Vector of hyperpameter values for prior on parameters #
#         of specified partition prior                          #

hpriorpar <- function(Prior, N) {
	if (Prior == "DP") {
		# Gamma(ag1, bg1) prior over concentration parameter
		ag <- 1
		bg <- 1/(N/2)
		hpriorpar <- c(ag, bg)
	}
	if (Prior == "PY") {
		# Gamma(ag, bg) prior over concentration parameter
		bg <- 1/(N/2)
		ag <- 1
		# Beta(ab,bb) prior over discount parameter
		ab <- 1
		bb <- 1
		hpriorpar <- c(ag, bg, ab, bb)
	}
	if (Prior == "ESCNB" || Prior == "ESCD") {
		# Gamma(ag, bg) prior over r
		bg <- 1
		ag <- 1
		# Beta(ab,bb) prior over  p
		ab <- 2
		bb <- 2
		hpriorpar <- c(ag, bg, ab, bb)
	}
  if (Prior == "ESCP" ){
    # Gamma(ag1, bg1) prior over lambda
    ag <- 1
    bg <- 1
    hpriorpar <- c(ag, bg)
  }
  if (Prior == "ESCB"||Prior == "ESCBshift"){
    # Beta(0.5,0.5) prior over pbinom
    ab <- 0.5
    bb <- 0.5
    hpriorpar <- c(ab, bb)
  }
	return(hpriorpar)
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: DP, PY, ESCNB, ESCD, ESCP, ESCB, ESCBshift      #
# Output: Parameter names used to save posterior samples     #

NameParam <- function(Prior, nfields) {
	if (Prior == "DP") {
		cnames <- c("theta")
		for (i in 1:nfields) {
			cnames <- c(cnames, sprintf("beta_%d",
				i))
		}
	}
	if (Prior == "PY") {
		cnames <- c("theta", "delta")
		for (i in 1:nfields) {
			cnames <- c(cnames, sprintf("beta_%d",
				i))
		}
	}
	if (Prior == "ESCNB") {
		cnames <- c("r", "p")
		for (i in 1:nfields) {
			cnames <- c(cnames, sprintf("beta_%d",
				i))
		}
	}
	if (Prior == "ESCD") {
		cnames <- c("alpha", "r", "p")
		for (i in 1:nfields) {
			cnames <- c(cnames, sprintf("beta_%d",
				i))
		}
	}
  if (Prior=="ESCP"){
    cnames <- c("lambda")
    for(i in 1:nfields){
      cnames <- c(cnames, sprintf("beta_%d",i))
    }
  }
  if(Prior == "ESCB"||Prior == "ESCBshift"){
    cnames <- c("N","p")
    for(i in 1:nfields){
      cnames <- c(cnames, sprintf("beta_%d",i))
    }
  }
	return(cnames)
}

# Computes proportion of each values in a field
# Used to compute empirical distribution of data
calcProp <- function(i, data) {
	return(table(data[, i])/sum(table(data[, i])))
}

# For each field, remap values to a consecutive
# list of integers starting with 1
remap <- function(i, data) {
	return(as.numeric(factor(data[, i], labels = seq(1,
		length(unique(data[, i]))))))
}

# Sample from a Dirichlet
rdirichletf <- function(params) {
	p <- sapply(params, function(x) rgamma(1, x))
	return(p/sum(p))
}

# Find Beta distribution parameters for distortion
# probailities, given the mean and sd
abeta <- function(bmean, bsd) {
	return(bmean * (1 - bmean) * (bmean/bsd^2) -
		bmean)
}


# selects a subset of points (called "block") that agree on some
# randomly choosen features. Chaperones are chosen within the block
# for a few consecutive iterations.
select_block <- function(x, size_block) {
	N <- dim(x)[1] #  number of datapoints
	L <- dim(x)[2] #  number of fields
	n_fields <- rbinom(1, L, prob = 0.75) # select how many fields the chaperones will agree on
	block <- c(1:N)
	if (n_fields > 0) {
		fields_to_agree <- sample(c(1:L), size = n_fields,
			replace = FALSE) # select which fields the chaperones will agree on
		x0 <- x[sample(c(1:N), 1), ] #select a point at random (to be used as \pivot\" point)"
		for (l in 1:n_fields) {
			if (length(block) > size_block) {
				block <- block[which(x[block, fields_to_agree[l]] ==
					x0[fields_to_agree[l]])]
			}
		}
	}
	return(block)
}

#' Remap data to a list of consecutive integers per field
#'
#' @param data Data frame containing only categorical variables
#' @return Data frame of remapped values to consecutive
#'         list of integers
#' @export
#' @examples
#' nclusters_per_size <- c(10,10,10,10)
#' numberFields <- 5
#' numberCategories <- rep(10,5)
#' trueBeta <- 0.01
#' data <- SimData(nclusters_per_size, numberFields, numberCategories, trueBeta)
#' DataRemap(data)
DataRemap <- function(data) {
	nfields <- ncol(data)
	x <- sapply(1:nfields, remap, data = data)
	return(x)
}

## Function Details ##
# Main sampler function: chaperones algorithm for cluster    #
# assignments and Slice sampler for priors' parameters and   #
# distortion probablities of the fields.                     #
#                                                            #
# Inputs: Data - e.g "Italy"                                 #
#         Prior - specify Prior Type: DP,PY,ESCNB, ESCD, ESCP, ESCB, ESCBshift      #
#         burn - MCMC burn-in period                         #
#         nsamples - MCMC iterations                         #
#         truncds - truncation value for field distortions   #
#         samds - indicator to sample field distortions      #
#         spacing - chaperones restricted Gibbs moves        #
#         w - step size for slice sampler                    #
#         m - number of steps                                #
#         hpriorpar - vector of prior hyperparameters        #
#         block_flag - TRUE for non-uniform chaperones       #
# Output: Posterior samples for prior parameters and         #
#         cluster assingments                                #
# Summary: This is a wrapper function that requires the user #
#          specify the data, prior and MCMC parameters       #
#          and returns  posterior samples of cluster         #
#          assignments and hyperparameters.                  #

#' @title Perform entity resolution with random partition priors.
#' 
#' @description \code{SampleCluster} runs Bayesian entity resolution model for categorical dataset under various possible random partition priors.
#'  The possible choices of random partition prior for the \code{Prior} argument are:
#'  \itemize{
#'    \item \code{"DP"}: Random partition induced by Dirichlet process (DP), also known as Chinese restaurant process prior.
#'    \item \code{"PY"}: Random partition induced by Pitman-Yor process (PY).
#'    \item \code{"ESCNB"}: Exchangeable sequence of clusters (ESC) model with zero-truncated negative binomial distribution.
#'    \item \code{"ESCD"}: ESC model with Dirichlet distribution.
#'    \item \code{"ESCP"}: ESC model with zero-truncated Poisson distribution.
#'    \item \code{"ESCB"}: ESC model with zero-truncated binomial distribution with fixed maximum cluster size (given by \code{Nbinom}).
#'    \item \code{"ESCBshift"}: ESC model with shifted binomial distribution, non-fixed maximum cluster size. 
#'  }
#' The \strong{Details} section below describes the properties of these priors. Especially for the details of ESC models, see Betancourt, Zanella and Steorts (2020) and Lee and Sang (2022).
#' 
#' The independent fields model of Steorts, Hall, and Fienberg (2016) has been used to model the noisy categorical data.
#' Specifically, the density vector \eqn{\theta_l} within each field is fixed as the empirical distribution of the data; See Section 4 of Betancourt et al. (2020) for detailed specification.
#'    
#' 
#' @details The choice of random partition prior significantly affects the Bayesian entity resolution results.
#' There are two major properties of random partition priors that affects the entity resolution result: \strong{microclustering} and \strong{balancedness}. 
#' 
#' \strong{Microclustering.} Traditional exchangeable random partition models such as the one induced by Dirichlet process (DP) or Pitman-Yor process (PY), assumes that 
#' the number of data points in each cluster grows linearly with the total number of data points. 
#' This growth assumption is not desirable in entity resolution tasks, where the cluster represents duplicates of the record so that
#' cluster sizes remain small, even for large data sets. 
#' The \emph{microclustering} property (Zanella, Betancourt and others, 2016) holds when cluster sizes grow sublinearly with the total number of data points, 
#' and the exchangeable sequence of clusters (ESC) model (Betancourt et al, 2020) is a very general class of random partition models that possess microclustering property. 
#' ESC models can directly control the prior distribution of the cluster sizes (\eqn{\mu}), which can be either zero-truncated negative binomial distribution, zero-truncated Poisson, and many others.
#' 
#' \strong{Balancedness.} Traditional exchangeable random partition models such as Chinese restaurant process induced by DP, often possesses the "rich-get-richer" property. 
#' This gives tendency that some few clusters become large and dominates the whole dataset a priori, leading to an unbalanced partition.
#' Lee and Sang (2022) studied the balancedness of random partition models, and characterized \emph{balance-averse} and \emph{balance-seeking} properties
#' which corresponds to favoring unbalanced and balanced partition in terms of probability, respectively.  
#' While the microclustering property has a similar rationale by limiting the growth rate of the largest cluster, 
#' the balancedness property analyzes how it assigns probabilities to partitions with different levels of balancedness in non-asymptotic point of view and they complement each other.
#' 
#' The table below summarizes the properties with different choice of random partition priors:
#' 
#' | \code{Prior} | Microclustering | Balancedness | 
#' | ----- | ----- | ----- | 
#' | \code{"DP"} | No | balance-averse | 
#' | \code{"PY"} | No | balance-averse | 
#' | \code{"ESCNB"} | Yes | balance-averse | 
#' | \code{"ESCD"} | Yes | N/A |  
#' | \code{"ESCP"} | Yes | balance-neutral | 
#' | \code{"ESCB"} | Yes, bounded microclustering | balance-seeking |  
#' | \code{"ESCBshift"} | Yes | balance-seeking |
#' 
#' Here \code{Prior = "ESCD"} is neither balance-averse nor balance-seeking, 
#' and \code{Prior = "ESCB"} has \emph{bounded microclustering property} (Betancourt et al, 2022), where the maximum size of the cluster is upper bounded by the fixed hyperparameter \code{Nbinom}.
#' To let the binomial parameter (number of trials) also be random, use \code{Prior = "ESCBShift"}.
#' Using balance-seeking prior leads to assigning less prior probability to partition with many singleton clusters, thus regularizing the emergence of singleton clusters.
#'
#' For the posterior sampling of cluster assignments (partition), a tailored MCMC scheme named \emph{chaperones algorithm} (Zanella, Betancourt and others, 2016) has been used.
#' 
#' @md
#' 
#' @param data Integer matrix, noisy categorical data where row corresponds to records and column corresponds to fields. 
#' @param Prior Specify partition prior: "DP", "PY", "ESCNB", "ESCD", "ESCP", "ESCB", "ESCBshift". See \strong{Details} section for their properties and differences.
#' @param burn MCMC burn-in period
#' @param nsamples MCMC iterations after burn-in
#' @param spacing Thinning for chaperones algorithm (default 1000)
#' @param block_flag TRUE for non-uniform chaperones (default)
#' @param beta_known default NULL, specify if true distortion probabilities (with length same as \code{ncol(data)}) is assumed to be known.
#' @param Nbinom default NULL, specify if Prior is "ESCB" which represents the maximum
#'
#' @return List with posterior samples for cluster assignments (Z),
#'         Prior parameters and distortion probabilities (Params)
#' @references   
#' Steorts, R. C., Hall, R., & Fienberg, S. E. (2016). A Bayesian approach to graphical record linkage and deduplication. Journal of the American Statistical Association, 111(516), 1660-1672.
#' 
#' Zanella, G., Betancourt, B., Miller, J. W., Wallach, H., Zaidi, A., & Steorts, R. C. (2016). Flexible models for microclustering with application to entity resolution. Advances in neural information processing systems, 29.
#' 
#' Betancourt, B., Zanella, G., & Steorts, R. C. (2020). Random partition models for microclustering tasks. Journal of the American Statistical Association, 1-13.
#' 
#' Betancourt, B., Sosa, J., & Rodríguez, A. (2022). A prior for record linkage based on allelic partitions. Computational Statistics & Data Analysis, 172, 107474.
#'
#' Lee, C. J., & Sang, H. (2022). Why the Rich Get Richer? On the Balancedness of Random Partition Models. Proceedings of the 39th International Conference on Machine Learning (ICML), PMLR 162:12521 - 12541.
#' @export
#' @examples
#' nclusters_per_size <- c(50,50,50,50)
#' numberFields <- 5
#' numberCategories <- rep(10,5)
#' trueBeta <- 0.01
#' # generate simulated data
#' simulatedData <- SimData(nclusters_per_size, numberFields, numberCategories, trueBeta)
#' # Fit ESCD model
#' posteriorESCD <- SampleCluster(data=simulatedData, Prior="ESCD", burn=0, nsamples=10)
#' 
SampleCluster <- function(data, Prior, burn, nsamples,
	spacing = 1000, block_flag = TRUE, beta_known = NULL, Nbinom = NULL){
  # sanity check
  if( !(Prior == "DP" || Prior == "PY" || Prior == "ESCNB" || Prior == "ESCD" || Prior == "ESCP" || Prior == "ESCB" || Prior == "ESCBshift")){
    stop("invaild Prior argument, see ?Samplecluster for possible random partition priors")
  }
  if( burn < 0 || nsamples < 1) stop("burn should be nonnegative integer and nsamples should be positive integer") 
  if(is.null(Nbinom) & (Prior== "ESCB")) stop("specify Nbinom value for Prior= ESCB")
  
	x <- DataRemap(data)
	N <- nrow(x)
	nfields <- ncol(x)
	proportions <- lapply(1:nfields, calcProp, data = x)
	# Initial clustering
	z0 <- sample(c(1:N), N, replace = TRUE) 
	# For ESCB, there is an upper bound of cluster size. Start with singleton partition
	if(Prior == "ESCB") z0 <- 1:N
	# No gaps
	z <- as.numeric(factor(z0,
	labels = seq(1,  length(unique(z0)))))
	# For ESCBshift, specify initial value of Nbinom
	if(Prior == "ESCBshift") Nbinom = max(tabulate(z))
	# Initial Prior parameter values
	x1 <- SetInit(Prior, N, Nbinom)
	samind <- SamInd(Prior)
	# Beta prior for distortion parameters
	bmean <- 0.005
	bsd <- 0.01
	cb <- abeta(bmean, bsd)
	db <- cb/bmean - cb
	hpriords <- c(cb, db)
	pdist <- cb/(cb + db)
	# Initial value for distortion parameters
	betas <- rep(pdist, nfields)
	truncds <- 0.1
  # If beta_known is known..
	if(!is.null(beta_known)){
	  if(length(beta_known)!=nfields) stop("wrong beta_known length")
	  betas <- beta_known
	}
	# hyperparameter values for priors on Prior parameters
	hpriorpar <- hpriorpar(Prior, N)

	# lo - lower limit distribution support (slice)
	# up - upper limit distribution support (slice)
    loup <- lowupSlice(Prior)
	lo <- loup$lo
	up <- loup$up

	s <- matrix(0, nsamples, length(x1) + nfields)
	zM <- matrix(0, nsamples, N)

	lx <- loglikxSP(betas, x, z, proportions)
	ll <- 1
	# number of iterations for each random block of chaperons
	itbl <- 50

	if (Prior == "DP" | Prior == "PY" | Prior ==
		"ESCNB"|Prior == "ESCP"|Prior == "ESCB"||Prior == "ESCBshift") {
		for (i in 1:(burn + nsamples)) {
			if ((i%%10) == 0) {
				print(paste("iter=", i))
			}
			Khat <- length(unique(z))
			Nk <- as.vector(table(z))
			# univariate slice sampler for all parameters
			w <- 1 # step size for slice sampler
			m <- 10 # limit on steps for slice sampler
			if(Prior!="ESCBshift"){
			  x1 <- unislicem(x1, N, Khat,lx, Nk, hpriorpar, w, m, lo, up, Prior, samind)
			}else{ # if prior is ESCBshift

			  # 1. first sample Nbinom (x1[1])
			  # In this implementation we trucate the support of Nbinom; TODO: full draw
			  lengthNbinom = 20 # maximum Nmax samples.
			  # support of Nbinom (assuming unimodal)
			  Ncandidates = (max(Nk)-1):(max(Nk)-1+lengthNbinom-1)
			  tempmat = matrix(sequence(from = Ncandidates[1]-Nk +1, by= 1, nvec = rep(lengthNbinom,Khat)), ncol = Khat)
			  temp = rowSums(lfactorial(tempmat)) # product of denominator term
			  logprobs = lbeta(N-Khat + hpriorpar[1], Ncandidates*Khat - N + Khat + hpriorpar[2]) - log(Ncandidates) + Khat*lfactorial(Ncandidates) - temp

			  # draw Nbinom from this log-prob. Gumbel trick
			  gumbels = -log(-log(runif(lengthNbinom)))
			  x1[1] = Ncandidates[which.max(logprobs + gumbels)]
			  # 2. then sample p.
			  x1[2] = stats::rbeta(1, N-Khat + hpriorpar[1], x1[1]*Khat - N + Khat + hpriorpar[2])
			}
			lods <- 0
			upds <- truncds
			if(is.null(beta_known)){
			  betas <- unislicespb1(betas, x, z, proportions,
				  hpriords, w, m, lods, upds, x1, N,
				  Khat, Nk, hpriorpar, Prior)
      }
			# parameters of reallocation algortihm #
			AB <- SetParam(Prior, x1, N)
			A <- AB$A
			B <- AB$B
			if (block_flag) {
				for (ss in 1:(spacing/itbl)) {
					bl <- select_block(x, 200)
					z1 <- Web_SamplerSP_fbl(x, z, bl,
						A, B, betas, proportions, 1, itbl)
					# no gaps in cluster assigment #
					z <- as.numeric(factor(z1[1, ], labels = seq(1,
						length(unique(z1[1, ])))))
				}
			} else {
				z1 <- Web_SamplerSP(x, z, A, B, betas,
					proportions, 1, spacing)
				# no gaps in cluster assigment #
				z <- as.numeric(factor(z1[1, ], labels = seq(1,
					length(unique(z1[1, ])))))
			}
			lx <- loglikxSP(betas, x, z, proportions)
			#loglik[i] <- lx
			if ((i > burn)) {
				zM[ll, ] <- z
				s[ll, ] <- c(x1, betas)
				ll <- ll + 1
			}
		}
	}

	if (Prior == "ESCD") {
		M <- 100 # max cluster size
		Khat <- length(unique(z))
		tz <- table(z)
		Nk <- as.vector(tz)
		if (max(tz) >= M) {
			M <- max(tz)
		}
		Lm <- table(factor(tz, levels = 1:M))
		alpha0 = x1[1]
		r0 = x1[2]
		p0 = x1[3]
		lmu0 = lgamma(seq(M) + r0) + (seq(M)) * log(p0) +
			r0 * log(1 - p0) - log(1 - (1 - p0)^r0) -
			lgamma(r0) - lgamma(seq(M) + 1)
		mu0 = exp(lmu0)
		for (i in 1:(burn + nsamples)) {
			if ((i%%10) == 0) {
				print(paste("iter=", i))
			}
			# univariate slice sampler for all parameters
			w <- 1 # step size for slice sampler
			m <- 10 # limit on steps for slice sampler
			x1 <- unislicemESCD(x1, lx, Lm, mu0,
				hpriorpar, w, m, lo, up, samind)
			lods <- 0
			upds <- truncds
			betas <- unislicespb1(betas, x, z, proportions,
				hpriords, w, m, lods, upds, x1, N,
				Khat, Nk, hpriorpar, Prior)

			# parameters of reseating algortihm
			alpha0 <- x1[1]
			r0 <- x1[2]
			p0 <- x1[3]

			# Sampling mu from dirichlet of size max cluster size
			lmu0 <- lgamma(seq(M) + r0) + (seq(M)) *
				log(p0) + r0 * log(1 - p0) - log(1 -
				(1 - p0)^r0) - lgamma(r0) - lgamma(seq(M) +
				1)
			mu0 <- exp(lmu0)
			alphasDM <- c(alpha0 * mu0 + Lm, abs(alpha0 *
				(1 - sum(mu0))))
			mus <- rdirichletf(alphasDM)
			AB <- SetParamESCD(mus, N, M)
			A <- AB$A
			B <- AB$B
			if (block_flag) {
				for (ss in 1:(spacing/itbl)) {
					bl <- select_block(x, 200)
					z1 <- Web_SamplerSP_fbl(x, z, bl,
						A, B, betas, proportions, 1, itbl)
					# no gaps in cluster assigment #
					z <- as.numeric(factor(z1[1, ],
					labels = seq(1, length(unique(z1[1, ])))))
				}
			} else {
				z1 <- Web_SamplerSP(x, z, A, B, betas,
					proportions, 1, spacing)
				# no gaps in cluster assigment #
				z <- as.numeric(factor(z1[1, ], labels = seq(1,
					length(unique(z1[1, ])))))
			}
			Khat <- length(unique(z))
			tz <- table(z)
			Nk <- as.vector(tz)
			if (max(tz) >= M) {
				M <- max(tz)
			}
			Lm <- table(factor(tz, levels = 1:M))
			lx <- loglikxSP(betas, x, z, proportions)
			if ((i > burn)) {
				zM[ll, ] <- z # cluster samples
				s[ll, ] <- c(x1, betas) # parameters samples
				ll = ll + 1
			}
		}
	}
	colnames(s) <- NameParam(Prior, nfields)
	return(list(Z = zM, Params = s))
}
