# This version breaks problem in to pieces to facilitate unit testing, but won't be as fast as pure cpp version
# used for direct comparison of old rcolgem methods and new phydynr methods


##############################################################

.solve.Q.A.L.boost <- function(h0, h1, A0, L0, tree, tfgy)
{ # uses C implementation
	out <- solveQALboost0(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
	 , h0
	 , h1
	 , L0
	 , A0
	 #, tree$maxHeight
	)
	out
}

.solve.Q.A.L.boost1 <- function(h0, h1, A0, L0, tree, tfgy)
{ # solves eqns for logit transform of Q
	out <- solveQALboost1(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
	 , h0
	 , h1
	 , L0
	 , A0
	)
	out
}


.solve.Q.A.L.deSolve <- function(h0, h1, A0, L0, tree, tfgy, integrationMethod='lsoda')
{
#~ dQLcpp
#~ function (x, FF, GG, YY, A0)
	m <- length(A0)
	.dQL <-  function( t, y, parms, ...){
		# t : h
		ih <- min( length(tfgy[[1]]), 1+floor( length(tfgy[[1]]) * t / (tree$maxSampleTime- tail(tfgy[[1]],1))  ) )
		list( as.vector( dQLcpp( y, tfgy[[2]][[ih]]
		 , tfgy[[3]][[ih]]
		 , tfgy[[4]][[ih]]
		 , A0 )))
	}
	times <- seq(h0, h1, length.out = 2)
	y0 <- c( as.vector( diag( length(A0))), L0)
	tt <- tryCatch(o <- ode(y0, times,  func = .dQL, parms=list() , method = integrationMethod),warning=function(w) w)
	if (is(tt,"warning")) return(tt)
	y1 <- o[nrow(o), -1]
	L1 <- unname( y1[length(y1)] )
	QQ <- matrix( y1[-length(y1)], nrow = m, ncol = m, byrow=TRUE)
	QQ <- QQ  / rowSums( QQ)
	##AA <- QQ %*% A0
	# Igor: i think this is the right way but this AA is not used anyway
	#AA <- t(QQ) %*% A0
	#AA <- sum(A0) * AA / sum(AA)
	list( QQ , NA , L1)
}
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

colik.modular0 <- function(tree, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , expmat = FALSE # not yet implemented
  , finiteSizeCorrection=F
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , returnTree = FALSE
  , solveODE = 0
) 
{
 
	if ( tree$maxHeight >  (tree$maxSampleTime- t0) ){
		warning('t0 occurs after root of tree. Results may be innacurate.')
	}
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	colik.fgy1(tfgy 
	  , tree
	  , integrationMethod
	  , timeOfOriginBoundaryCondition
	  , maxHeight
	  , expmat 
	  , finiteSizeCorrection
	  , forgiveAgtY 
	  , AgtY_penalty
	  , returnTree
	  , solveODE 
	)
}

colik.fgy1 <- function(
	tfgy 
	  , tree
	  , integrationMethod='lsoda'
	  , timeOfOriginBoundaryCondition=T
	  , maxHeight=Inf
	  , expmat =F
	  , finiteSizeCorrection=F
	  , forgiveAgtY =1
	  , AgtY_penalty=1
	  , returnTree=F
	  , solveODE =0
){
	if (tfgy[[1]][1] < tfgy[[1]][2] ) stop('tfgy must be in order of decreasing time.')
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
		tree$lstates <- matrix(0, ncol = tree$n + tree$Nnode, nrow = tree$m)
		tree$lstates[1,] <- 1
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
	  worklist <- c(h0,h1)
	  numTries <- 0
	  L <- 0
	  while (length(worklist) > 1) {
	    hw0 <- worklist[1]
            hw1 <- worklist[2]
            worklist <- worklist[c(-1)]
            if (solveODE==0) {
			  out <- .solve.Q.A.L.deSolve(hw0, hw1, A0, L, tree, tfgy, integrationMethod=integrationMethod)
			} else if (solveODE==1){
				out <- .solve.Q.A.L.boost(hw0, hw1, A0, L, tree, tfgy)
			} else if (solveODE == 2){
				out <- .solve.Q.A.L.boost1(hw0, hw1, A0, L, tree, tfgy)
			} else{
				stop('Invalid *solveODE* value' )
			}
	    if (is(out,"warning")) {
              worklist <- c(hw0,(hw0+hw1)/2,worklist)
	      print("Try again")
	      next  # try again
        }	      	
	    Q <- out[[1]]
	    A <- out[[2]]
	    L <- L+out[[3]]

	     # clean output
	     if (is.nan(L)) {L <- Inf}
	     if (sum(is.nan(Q)) > 0) Q <- diag(length(A))
	     if (sum(is.nan(A)) > 0) A <- A0
		
	      #update mstates 
	      tree$mstates <- update_states1( tree$mstates , Q, extantLines)
              if (length(worklist)>1) {
	        if (length(extantLines) > 1 ){
	          A0 <- rowSums(tree$mstates[,extantLines]) #TODO faster compute this on fly
                } else if (length(extantLines)==1){
		  A0 <- tree$mstates[,extantLines]
		}
	      }
     } # end while worklist


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
					L <- (h1 - h0) * pa_corate[[2]]
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
				if (is.infinite( loglik)){
					if (returnTree){
						return(list( loglik = loglik
						 , tree = tree ))
					} else{
						return(loglik)
					}
				}
				
				# finite size corrections for lines not involved in coalescent
				if (finiteSizeCorrection){
					rco_finite_size_correction2(alpha, pa_corate[[1]], A, extantLines, tree$mstates )
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

