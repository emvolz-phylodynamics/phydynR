# ancestral state estimation with a forward/backward algorithm
# inspired by http://www.biorxiv.org/content/biorxiv/early/2017/09/13/188516.full.pdf

ace <- function(tree, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , maxHeight = Inf 
  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
  , AgtY_penalty = 1 # penalises likelihood if A > Y
  , step_size_res = 10 # for adaptive ode solver, set to value < 1
)
{
	tfgy <- demographic.process.model( theta, x0, t0, tree$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	m <- nrow(tfgy[[2]][[1]] )
	
	bdt <- 	colik.pik.fgy(tree 
	  , tfgy
	  , timeOfOriginBoundaryCondition
	  , maxHeight
	  , forgiveAgtY 
	  , AgtY_penalty
	  , returnTree = TRUE
	  , step_size_res = step_size_res
	)[[2]]
	
	eventTimes <- unique( sort(tree$heights) )
	if (maxHeight) { 
		eventTimes <- eventTimes[eventTimes<=maxHeight]
	}
	# NOTE extant refers to node indices with parent lineage extant 
	extantAtEvent_nodesAtHeight <- eventTimes2extant( eventTimes, tree$heights, tree$parentheight ) #1-2 millisec
	extantAtEvent_list <- extantAtEvent_nodesAtHeight[[1]]
	nodesAtHeight <- extantAtEvent_nodesAtHeight[[2]]
	
	get.fgy <- function(h)
	{# NOTE tfgy needs to be in order of decreasing time, fist time point must correspond to most recent sample
		#ih <- 1+floor( length(tfgy[[1]]) * h / (tree$maxSampleTime- tail(tfgy[[1]],1))  )
		ih <- min( length(tfgy[[1]]), 1+floor( length(tfgy[[1]]) * h / (tree$maxSampleTime- tail(tfgy[[1]],1))  ) )
		list( .F = tfgy[[2]][[ih]]
		 , .G = tfgy[[3]][[ih]]
		 , .Y = tfgy[[4]][[ih]]
		 )
	}
	
	lstates <- bdt$lstates 
	acestates <- matrix( NA, nrow = nrow(lstates), ncol = ncol(lstates))
	# initialise ace of root + dgtrs of root 
	acestates[,bdt$root] <- lstates[,bdt$root] 
	fgy <- get.fgy( max( bdt$heights ))
	x <- t(fgy$.F) %*% acestates[, bdt$root ] 
	x <- x/ sum(x)
	x <- (x + acestates[,bdt$root])/2
	acestates[, bdt$daughters[bdt$root,] ] <- x
	
	for (ih in (length(eventTimes)):2)
	{
		h1 <- eventTimes[ih]
		h0 <- eventTimes[ih-1]
		
		#get A0, process new samples, calculate state of new lines
		extantLines <- extantAtEvent_list[[ih-1]]
		
		# acestates for nodes_h1 need to be initialised 
		nodes_h1 <- nodesAtHeight[[ih]]
		nodes_h0 <- nodesAtHeight[[ih-1]]
		
		# derive Q between h1 and h0 using forward equations
		if ( h1  > h0 ){
			QQ <- solveQfwd0(times=tfgy[[1]], Fs=tfgy[[2]], Gs=tfgy[[3]], Ys=tfgy[[4]], deaths=tfgy[[6]], m, h1, h0)
			acestates[, extantLines] <- QQ %*% acestates[, extantLines ]
		}
		
		## initialise ace states for next generation 
		fgy <- get.fgy(h0)
		for (u in nodes_h0 ){
			x <- acestates[, u] * lstates[, u] 
			acestates[, u] <- x / sum(x) 
			
			## catch NAs which may occur if maxHeight option is enabled
			if (any(is.na( acestates[,u] )) & !any(is.na(lstates[,u]))){
				# set acestate to lstate 
				acestates[,u] <- lstates[,u]
			}
			
			v0 <- bdt$daughters[u,1]
			v1 <- bdt$daughters[u,2]
			x <- t(fgy$.F) %*% acestates[, u ] 
			x <- x/ sum(x)
			x <- (x + acestates[,u])/2
			if( !is.na( v0 )){
				acestates[,v0] <- x
			}
			if( !is.na( v1 )){
				acestates[,v1] <- x
			}
		}
	}
	bdt$acestates <- acestates
	bdt
}
