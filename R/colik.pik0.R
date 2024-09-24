# sourceCpp( 'dPikL0.cpp' )

################################################################
# helper function for updating ancestral node and likelihood terms
.update.alpha.cpp <- function(u,v, tree, fgy, A)
{
	pacorate <- update_alpha0(tree$mstates[,u] 
	  , tree$mstates[,v]
	  , fgy$.F
	  , fgy$.Y
	  , A
	)
	pacorate[[1]] <- as.vector( pacorate[[1]] )
	pacorate
}




################################################################################
#' Compute the coalescent log-likelihood.
#'
#' Computes the log-likelihood using a coalescent (or structured coalescent)
#' genealogical model based on a user-supplied demographic process.
#'
#' @param tree A DatedTree object
#' @param theta A named numeric vector or named list of parameter values used by
#'   the demographic model
#' @param demographic.process.model See \code{\link[phydynR]{build.demographic.process}}
#' @param x0 A named vector of initial conditions required by the model. This
#'   includes demes and any other dynamic variables.
#' @param t0 The time of origin of the process. Should predate the root of the tree.
#' @param res Integer number of time steps to use when simulating model.
#' @param integrationMethod If simulating an ODE (ordinary differential equation)
#'    model, this provides the integration routine corresponding to options in
#'    deSolve.
#' @param timeOfOriginBoundaryCondition If TRUE, will return -Inf if the root of
#'   the tree precedes the time of origin.
#' @param maxHeight It will only count internode intervals in the likelihood that
#'   occur after maxHeight years before present. Useful for large trees and when
#'   you do not want to model the entire demographic history.
#' @param forgiveAgtY If number of extant lineages exceeds simulated population
#'   size, return -Inf if this value is zero, or forgive the discrepancy if zero.
#'   If between zero and one, only forgive the discrepancy if this proportion of
#'   lineages is less than the given value.
#' @param AgtY_penalty If number of extant lineages exceeds simulated population
#'   size, penalize likelihood with value L*AgtY_penalty where L is the cumulative
#'   coalescent rate within the given internode interval. 0<= AgtY_penalty <= Inf.
#' @param returnTree If TRUE, a copy of the tree is also returned, which includes
#'   the inferred states of lineages and likelihood terms at each internal node.
#' @param step_size_res Parameter for the ODE solver; it is the default number 
#'   of timesteps to use when solving coalescent equations in each internode 
#'   interval
#' @param likelihood Toggle likelihood approximation to be used 
#'   (QL fast/approximate, PL1 faster better approximation, PL2 slow/good 
#'   approximation. See \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006546}{Volz & Siveroni 2018} 
#'   for details.
#'   
#'
#' @return The coalescent (or structured coalescent) log likelihood (numeric).
#' @author Erik Volz
#' 
#' @export
#'
#' @examples
#' # A simple exponential growth model with birth rates beta and death rates gamma:
#' dm <- build.demographic.process(births = c(I = 'parms$beta * I'),
#'                                 deaths = c(I = 'parms$gamma * I'),
#'                                 parameterNames = c('beta', 'gamma'),
#'                                 rcpp = FALSE,
#'                                 sde = FALSE)
#'
#' # simulate a tree based on the model:
#' tre <- sim.co.tree(list(beta = 1.5, gamma = 1),
#'                         dm,
#'                         x0  = c(I = 1),
#'                         t0 = 0,
#'                         sampleTimes = seq(10, 15, length.out=50),
#'                         res = 1000)
#'
#' # Compute a likelihood
#' colik(tre,
#'       list(beta = 1.5, gamma = 1),
#'       dm,
#'       x0 = c(I = 1),
#'       t0 = -1,
#'       res = 1e3)
colik = colik.pik <- function(tree, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , step_size_res = 10 # for adaptive ode solver, set to value < 1
  , likelihood = c( 'PL2', 'PL1', 'QL' )
  , PL2 = FALSE
) {
 
	if ( tree$maxHeight >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	likelihood = likelihood[1] 
	if ( likelihood == 'PL2' ){
		l = colik.pik.fgy(tree 
		  , tfgy
		  , timeOfOriginBoundaryCondition
		  , maxHeight
		  , forgiveAgtY 
		  , AgtY_penalty
		  , returnTree
		  , step_size_res
		  , PL2 = TRUE 
		)
	} else if (likelihood== 'PL1'){
		l = colik.pik.fgy(tree 
		  , tfgy
		  , timeOfOriginBoundaryCondition
		  , maxHeight
		  , forgiveAgtY 
		  , AgtY_penalty
		  , returnTree
		  , step_size_res
		  , PL2 =FALSE 
		)
	} else if ( likelihood == 'QL'){
		l = colik.fgy1(tfgy 
		  , tree
		  , integrationMethod = integrationMethod 
		  , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
		  , maxHeight = maxHeight 
		  , forgiveAgtY  = forgiveAgtY
		  , AgtY_penalty = AgtY_penalty
		  , returnTree = returnTree
		)
	} else{
		stop('Specify a valid likelihood method such as QL or PL2')
	}
	l
}

#' Compute a structured coalescent likelihood given a dated genealogy and a demographic history in FGY format
#' @inheritParams colik
#' @param tfgy TO DO
#' 
#' @return The coalescent (or structured coalescent) log likelihood (numeric).
#' @author Erik Volz
#'
#' @export 
colik.pik.fgy <- function(tree
  , tfgy
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , step_size_res = 10
  , PL2 = FALSE
) {
	if (tfgy[[1]][1] < tfgy[[1]][2] ){
		or <- order( tfgy[[1]] , decreasing=TRUE )
		tfgy[[1]] <- tfgy[[1]][or]
		tfgy[[2]] <- tfgy[[2]][or]
		tfgy[[3]] <- tfgy[[3]][or]
		tfgy[[4]] <- tfgy[[4]][or]
	}
	t0 <- tail( tfgy[[1]], 1)
	if ( min(tree$maxHeight, maxHeight) >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	
	# validate input TODO may be slow
	tfgy[[2]] <- lapply( tfgy[[2]] , function(x) pmax( x, 0 ) )
	tfgy[[3]] <- lapply( tfgy[[3]] , function(x) pmax( x, 0 ) )
	tfgy[[4]] <- lapply( tfgy[[4]] , function(x) setNames(pmax(0, x ), names(x) ) )
	
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		#ih <- 1+floor( length(tfgy[[1]]) * h / (tree$maxSampleTime- tail(tfgy[[1]],1))  )
		ih <- min( length(tfgy[[1]]), 1+floor( length(tfgy[[1]]) * h / (tree$maxSampleTime- tail(tfgy[[1]],1))  ) )
		list( .F = tfgy[[2]][[ih]]
		 , .G = tfgy[[3]][[ih]]
		 , .Y = tfgy[[4]][[ih]]
		 )
	}
	.newNodes.order <- function(nn, extantLines)
	{
		# gives order to traverse new nodes if heights are concurrent
		if (length(nn) <=1 ){
			return (nn)
		}
		o <- c()
		nrep <- 0
		while (!(all(nn %in% o) ) ){
			for (a in setdiff( nn, o) ){
				d <- tree$daughters[a, ] 
				if (all( d%in%union( union(1:tree$n,extantLines), o) )){ # note include samples here, since sample with zero edge length will not be extant
					o <- c( o, a )
				}
			}
			nrep <- nrep + 1
			if (nrep > 1 + length(nn)){
				o <- sample(nn, replace=F)
				break
			}
		}
		o
	}
	
	if (is.null( tree$n ) ) tree$n <- length( tree$sampleTimes)
	if (is.null(tree$m)) tree$m <- length( tfgy[[4]][[1]] )
	if (ncol(tree$sampleStates)==1 & tree$m == 2){ # make exception for unstructured models
		tree$sampleStates <- cbind( tree$sampleStates, rep( 0 , tree$n) )
	}
	#if (is.null( tree$lstates)) 
	{
		tree$lstates <- matrix(NA, ncol = tree$n + tree$Nnode, nrow = tree$m)
		#~ tree$lstates[1,] <- 1
		tree$lstates[,1:tree$n ] <- t( tree$sampleStates )
	}
	tree$mstates <- tree$lstates + 1 - 1
	if (is.null( tree$ustates)) tree$ustates <- matrix(0, ncol = tree$n + tree$Nnode, nrow = tree$m)
	
	eventTimes <- unique( sort(tree$heights) )
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]
	}
	S <- 1
	L <- 0
	
	# NOTE extant refers to node indices with parent lineage extant 
	extantAtEvent_nodesAtHeight <- eventTimes2extant( eventTimes, tree$heights, tree$parentheight ) #1-2 millisec
	extantAtEvent_list <- extantAtEvent_nodesAtHeight[[1]]
	nodesAtHeight <- extantAtEvent_nodesAtHeight[[2]]
	#variables to track progress at each node
	tree$A <- c() #h->A
	tree$lnS <- c()
	tree$lnr <- c()
	tree$ih <- c()
	loglik <- 0
	for (ih in 1:(length(eventTimes)-1))
	{
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
	
	  #get A0, process new samples, calculate state of new lines
	  extantLines <- extantAtEvent_list[[ih]]
	  if (length(extantLines) > 1 ){
	    A0 <- rowSums(tree$mstates[,extantLines]) #TODO faster compute this on fly
          } else if (length(extantLines)==1){ A0 <- tree$mstates[,extantLines] }
        
		## update mstates, L 
		# <new code here>
		pik0 <- cbind( tree$mstates[, extantLines ]  ,rep(0, tree$m ))
		if (PL2){
			pik1L1 <- solvePikL1(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
			 ,h0
			 ,h1
			 ,pik0
			 ,step_size_res
			) 
		} else{
			pik1L1 <- solvePikL0(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
			 ,h0
			 ,h1
			 ,pik0
			 ,step_size_res
			) 
		}
		pik1 <- pmax( pik1L1[[1]] , 0)
		pik1 <- pik1[ , 1:(ncol(pik1)-1)]
		if (is.matrix(pik1)){
			pik1 <- t( t(pik1) / colSums( pik1 ) )
		} else{
			pik1 <- pik1 / sum(pik1)
		}
		tree$mstates[, extantLines] <- pik1
		L <- L + pik1L1[[2]] 
		# </new code>
		
		#if applicable: update ustate & calculate lstate of new line
		newNodes <- nodesAtHeight[[ih+1]]
		newNodes <- newNodes[newNodes > tree$n] 
		# NOTE re-evaluation of A here appears to be quite important 
		if (length(extantLines)>1){
			A <- rowSums( tree$mstates[, extantLines] ) #TODO speedup
		} else{
			A <- tree$mstates[, extantLines]
		}
		{
			YmA <- (sum(fgy$.Y) - length(extantLines))
			if (YmA < 0){
				if (length(extantLines)/length(tree$tip.label)  > forgiveAgtY){
					L <- Inf
				} else{
					L <- L + L * abs(YmA) * AgtY_penalty
				}
			}
			
			for (alpha in .newNodes.order(newNodes, extantLines) )
			{
				u <- tree$daughters[alpha,1]
				v <- tree$daughters[alpha,2]
				#pa_corate <- .update.alpha(u,v, tree, fgy, A) # 
				pa_corate <- .update.alpha.cpp(u,v, tree, fgy, A) 
				tree$coalescentRates[alpha] <- pa_corate[[2]] 
				if (is.nan(L))
				{
					warning('is.nan(L)')
					#L <- (h1 - h0) * pa_corate[[2]]
					L <- Inf
				}
				tree$lstates[,alpha] <- pa_corate[[1]]
				tree$mstates[,alpha] <- pa_corate[[1]]
				
				# trace
				tree$A <- rbind( tree$A, A)
				tree$lnS <- c( tree$lnS, -L )
				tree$lnr <- c( tree$lnr, log(pa_corate[[2]]) )
				tree$ih <- c( tree$ih, h1)
				
				# update lik
				loglik <- loglik + log( pa_corate[[2]] ) - L ;
				if (is.infinite( loglik) | is.na(loglik) ){
					if (returnTree){
						return(list( loglik = -Inf
						 , tree = tree ))
					} else{
						return(-Inf)
					}
				}
				
			} #else if (length(newNodes)>1) {
			#	stop("Likelihood for concurrent internal nodes not yet implemented")
			#} 
		}
		
		# if coalescent occurred, reset cumulative hazard function
		if (length(newNodes) > 0) {
			L <- 0 
		}
	}
	if (returnTree){
		return(list( loglik = loglik
		 , tree = tree ))
	}
	return( loglik)
} 

