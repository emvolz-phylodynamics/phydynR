# This version breaks problem in to pieces to facilitate unit testing, but won't be as fast as pure cpp version
# used for direct comparison of old rcolgem methods and new phydynr methods

# NOTE this doesnt work, since src dir does not make it into package tarball
# instead , src these directly in testing code 
#~ sourceCpp(system.file('src/solveQALboost0.cpp', package='phydynR')) #TODO 
#~ sourceCpp(system.file('src/colikModular0.cpp', package='phydynR')) #TODO 

require(rcolgem)

##############################################################
.solve.Q.A.L.rcolgem <- function(h0, h1, A0, L0, tree, fgymat)
{ # uses C implementation
	Q0 <- diag(tree$m)
	parameters 	<- c(tree$m, tree$maxHeight, length(heights), sum(A0), as.vector(fgymat))
	y0 <- c( as.vector( Q0), A0,  L0 ) #
	o <- ode(y=y0, c(h0,h1), func = "dQAL", parms = parameters, dllname = system.file("src/rcolgem",package='rcolgem')
	  , initfunc = "initfunc", method='lsoda' )
	Q1 		<- t( matrix(  abs(o[nrow(o),2:(1 + tree$m^2)]) , nrow=tree$m) ) #NOTE the transpose
	A1 <- o[nrow(o), (1 + tree$m^2 + 1):(1 + tree$m^2 +  tree$m)]
	L1 <- o[nrow(o), ncol(o)]
	return ( list(  unname(Q1), unname(A1),  unname(L1) ) ) #
}

.solve.Q.A.L.boost <- function(h0, h1, A0, L0, tree, tfgy)
{ # uses C implementation
	solveQALboost0(tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
	 , h0
	 , h1
	 , L0
	 , A0
	 , tree$maxHeight)
}
################################################################
# helper function for updating ancestral node and likelihood terms
.update.alpha <- function(u,v, tree, fgy, A)
{ #CPP
	.F <- fgy$.F
	.G <- fgy$.G
	.Y <- fgy$.Y
	{
		.Y <- pmax(A, .Y)
		FklXpuk_Yk <- (.F * tree$ustates[u,]/.Y)
		FklXpvk_Yk <- (.F * tree$ustates[v,]/.Y)
		FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
		FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
		vk_Yk <- pmin(pmax(tree$ustates[v,]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
		uk_Yk <- pmin(pmax(tree$ustates[u,]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
		ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
	}
	
	coalescentRate <- max( sum(ratekl) , 0)
#~ 	if (is.nan(L))
#~ 	{
#~ 		L <- (h1 - h0) * sum(ratekl) 
#~ 	}
#~ 	logCoalescentSurvivalProbability <- -L
	
	if (sum(ratekl)==0) {ratekl <- rep(1/tree$m, tree$m) * 1e-6}
	# definitions of alpha state
	p_a <- as.vector( ratekl / sum(ratekl) )
	
	list( pa = p_a, corate =coalescentRate) 
}

.update.alpha.cpp <- function(u,v, tree, fgy, A)
{
#~ returns:  
#~ 	o["pa"] = pa;
#~ 	o["corate"] = corate;
	update_alpha0(tree$ustates[u,] 
	  , tree$ustates[v,]
	  , fgy$.F
	  , fgy$.Y
	)
}





################################################################################

#~ rcolgem::finite_size_correction

#~void finite_size_correction0(mat& mstates
#~   , int a
#~   , vec pa
#~   , uvec extantLines
#~   , vec A
#~ )

################################################################################
.extant.at.height <- function(h, tree)
{
	return( which( tree$heights <= h & tree$parentheight > h)  )
}

################################################################################
#~ update_states0...

update_states_rcolgem0 <- function( tree, Q, extantLines )
{
	if (length(extantLines) > 1)
	{
		nms <- t(rcolgem::update_mstates_arma( extantLines, t(Q), t(tree$mstates) ))
		tree$mstates <- nms
		# NOTE will have slightly different results if recomputing A here, influences fcs: 
		#A <- colSums(tree$mstates[extantLines,])
			# renormalise
		tree$mstates[extantLines,] <- pmax( tree$mstates[extantLines,], 0 )
		tree$mstates[extantLines,] <- tree$mstates[extantLines,] / rowSums( tree$mstates[extantLines,] )
	}
	else{
		tree$mstates[extantLines,] <- t( t(Q) %*% tree$mstates[extantLines,] )
		tree$mstates[extantLines,] <- pmax(tree$mstates[extantLines,],0) / sum(abs(tree$mstates[extantLines,]))
		#recalculate A
	}
	tree$mstates
}
################################################################################
colik.modular0 <- function(tree, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , AgtYboundaryCondition = TRUE # can also be numeric > 0
  , maxHeight = Inf 
  , expmat = FALSE
  , forgiveAgtY = .2 #can be NA
  , returnTree = FALSE
) 
{
	# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res, integrationMethod=integrationMethod) 
		
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		ih <- 1+floor( length(tfgy[[1]]) * h / (tree$maxSampleTime- tail(tfgy[[1]],1))  )
		list( .F = tfgy[[2]][[ih]]
		 , .G = tfgy[[3]][[ih]]
		 , .Y = tfgy[[4]][[ih]]
		 )
	}
	
	if (is.null( tree$n ) ) tree$n <- length( tree$sampleTimes)
	if (is.null(tree$m)) tree$m <- length( tfgy[[4]][[1]] )
	if (is.null( tree$lstates)) {
		tree$lstates <- matrix(0, ncol = tree$n + tree$Nnode, nrow = tree$m)
		tree$lstates[,1:tree$n ] <- t( tree$sampleStates )
	}
	tree$mstates <- tree$lstates
	if (is.null( tree$ustates)) tree$ustates <- matrix(0, ncol = tree$n + tree$Nnode, nrow = tree$m)
	
	# construct forcing timeseries for ode's, needed for rcolgem version
	heights <- sort( tree$maxSampleTime - tfgy[[1]] )
	heights <- heights[heights <= tree$maxHeight]
	heights <- heights[heights >= 0]
	fgymat <- t( sapply( heights, function(h) 
	  with(get.fgy(h), {
		c( as.vector(.F), as.vector(.G), as.vector(.Y) )
	  })
	) )
	fgymat <- pmax(fgymat, 0)
# TODO inadequate precision here causing confusion; 
#~ inhs <- tree$heights[(tree$n+1):(tree$n+tree$Nnode)]
#~ plot( ecdf( inhs))
	eventTimes <- unique( sort(tree$heights) )
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]
	}
	S <- 1
	L <- 0

# TODO some loss of precision in eventTimes & treeheights
	extantAtEvent_nodesAtHeight <- eventTimes2extant( eventTimes, tree$heights, tree$parentheight )
	extantAtEvent_list <- extantAtEvent_nodesAtHeight[[1]]
	nodesAtHeight <- extantAtEvent_nodesAtHeight[[2]]
	
	#debugging variables
	tree$A <- c() #h->A
	tree$lnS <- c()
	tree$lnr <- c()
	tree$ih <- c()
	
	for (ih in 1:(length(eventTimes)-1))
	{
		h0 <- eventTimes[ih]
		h1 <- eventTimes[ih+1]
		fgy <- get.fgy(h1)
		
		#get A0, process new samples, calculate state of new lines
		##extantLines <- .extant.at.height(h0, tree) #TODO would be faster to compute on fly
		extantLines <- extantAtEvent_list[[ih]]
		if (length(extantLines) > 1 ){
browser()
			A0 <- colSums(tree$mstates[extantLines,]) 
		} else if (length(extantLines)==1){ A0 <- tree$mstates[extantLines,] }
		out <- .solve.Q.A.L.rcolgem(h0, h1, A0, L, tree, fgymat)
		out2 <- .solve.Q.A.L.boost(h0, h1, A0, L, tree, tfgy) 
# TODO  
		Q <- out[[1]]
		A <- out[[2]]
		L <- out[[3]]

		# clean output
		if (is.nan(L)) {L <- Inf}
		if (sum(is.nan(Q)) > 0) Q <- diag(length(A))
		if (sum(is.nan(A)) > 0) A <- A0
		
		#update mstates 
		update_states0( tree$mstates , Q) #void
		mstates0 <- update_states_rcolgem0( tree, Q, extantLines )
#TODO 
		
		#if applicable: update ustate & calculate lstate of new line
		newNodes <- which( tree$heights == h1)
		newNodes <- newNodes[newNodes > length(tree$sampleTimes)] # <- internal nodes
		newNodes2 <- nodesAtHeight[[ih]]
		newNodes2 <- newNodes2[newNodes2 > tree$n] 
#TODO
		{
			if ( (sum(.Y) < length(extantLines)) & (!forgiveAgtY) ) { L <- Inf }
			else if ( (sum(.Y) < length(extantLines)) & (length(extantLines)/length(tree$tip.label)) > forgiveAgtY) { L <- Inf }
			
			#for (alpha in newNodes){
			if (length(newNodes)==1)
			{
				#? should rewrite to use helper functions
				alpha <- newNodes
				u <- tree$daughters[alpha,1]
				v <- tree$daughters[alpha,2]
				
				tree2 <- tree
				pa_corate2 <- .update.alpha.cpp(u,v, tree2, fgy, A)
				pa_corate <- .update.alpha(u,v, tree, fgy, A)
#TODO
				tree$coalescentRates[alpha] <- pa_corate[[2]] 
				if (is.nan(L))
				{
					warning('is.nan(L)')
					L <- (h1 - h0) * pa_corate[[2]]
				}
				tree$lstates[alpha,] <- pa_corate[[1]]
				tree$mstates[alpha,] <- pa_corate[[1]]
				
				# debug
				tree$A <- rbind( tree$A, A)
				tree$lnS <- c( tree$lnS, -L )
				tree$lnr <- c( tree$lnr, log(pa_corate[[2]]) )
				tree$ih <- c( tree$ih, h1)
				
				
				# update lik
				loglik <- loglik + log( pa_corate[[2]] ) - L ;
				
				# finite size corrections for lines not involved in coalescent
				mstates_rcolgem <- t( rcolgem::finite_size_correction( pa_corate[[1]], A, extantLines, t(tree$mstates) ) )
				mstates2 <- tree$mstates
				finite_size_correction0( mstates2, alpha, pa_corate[[1]], extantLines, A ) #void
				tree$mstates <- mstates_rcolgem
#TODO
			} else{
				stop("concurrent internal nodes")
			} 
		}
		
		# if coalescent occurred, reset cumulative hazard function
		if (length(newNodes) > 0) {
		  L <- 0 }
	}
	
	if (returnTree){
		return(list( loglik = loglik
		 , tree = tree ))
	}
	return( loglik)
} 

