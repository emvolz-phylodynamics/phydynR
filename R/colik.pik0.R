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
colik = colik.pik <- function(tree, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , step_size_res = 10 # for adaptive ode solver, set to value < 1
) {
 
	if ( tree$maxHeight >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	colik.pik.fgy(tree 
	  , tfgy
	  , timeOfOriginBoundaryCondition
	  , maxHeight
	  , forgiveAgtY 
	  , AgtY_penalty
	  , returnTree
	  , step_size_res
	)
}



colik.pik.fgy = colik.pik0.fgy <- function(tree
  , tfgy
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , step_size_res = 10
) {
	if (tfgy[[1]][1] < tfgy[[1]][2] ) stop('tfgy must be in order of decreasing time.')
	t0 <- tail( tfgy[[1]], 1)
	if ( min(tree$maxHeight, maxHeight) >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	
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
		pik0 <- cbind( tree$mstates[, extantLines ]  ,rep(0, m ))
		pik1L1 <- solvePikL0(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
		 ,h0
		 ,h1
		 ,pik0
		 ,step_size_res
		) 
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

