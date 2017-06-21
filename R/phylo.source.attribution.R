## compute 1) matrix of infector probabilities and 2) probability that infector is in sample

require(deSolve)
require(Rcpp)

phylo.source.attribution <- function( trees, sampleTimes, f.t, Y.t, maxTMRCA = NULL, res = 1e3, treeErrorTol = 1e-3)
{
	if (is.null( maxTMRCA )) maxTMRCA <- Inf
	if ( class( trees) == 'phylo' ) trees <- list( trees )
	if (!is.function(f.t)) {
		ft <- f.t
		f.t <- function(t) ft
	}
	if (!is.function(Y.t)) {
		Yt <- Y.t
		Y.t <- function(t) Yt
	}
	
	mst <- max(sampleTimes )
	sampleHeights <- mst - sampleTimes
	n <- length(sampleTimes)
	
	A.t <- function( t, bdt )
	{
		h <- mst - t
		sum(sampleHeights <= h ) - sum( bdt$heights[(n+1):(n+bdt$Nnode)] <= h )
	}
	
	d.Lambda <- function( h, Lambda, parms , ... )
	{
		t <- mst - h
		f <- f.t( t )
		Y <- Y.t (t )
		A <- A.t( t, parms$bdt )
		list( c(Lambda = unname( f*(Y-A)/Y^2)) ) 
		#list( c(Lambda = unname( f/Y )) ) 
	}
	
	W <- matrix ( 0, nrow=n, ncol =n); 
	TIPLABELS <- trees[[1]]$tip.label
	rownames(W) = colnames(W) <- TIPLABELS
	ntrees <- length(trees)
	for (tree in trees )
	{
		bdt <- DatedTree( tree, sampleTimes , tol = treeErrorTol)
		heights <- seq(0, min( bdt$maxHeight, maxTMRCA), length.out = res )
		Lambda <- ode( c(Lambda = 0 ), times = heights , func = d.Lambda, parms = list( bdt = bdt ),method = 'lsoda' )[,2]
		dh <- heights[2] - heights[1]
		Lambda.h <- function(h) Lambda[min(length(Lambda), 1 + floor( h / dh ))] # fast interpolator
		
		# traverse nodes in postorder, updating psi & W along the way
		bdt <- reorder.phylo( bdt, order = 'postorder')
		ancs_postorder <- c(); for (a in bdt$edge[,1] ) if (!(a %in% ancs_postorder)) ancs_postorder <- c( ancs_postorder, a )
		# update node_daughters
		node_daughters <- t( sapply( 1:(bdt$Nnode+n), function(u) {
			i <- which( bdt$edge[,1] == u)
			if (length(i)==2)
			{
				return( bdt$edge[i,2] ) 
			} else{
				return( c(NA, NA))
			}
		}) )
		
		tips2Wcoords <- match( bdt$tip.label , TIPLABELS )
		#  check that reorder edges does not screw up heights vec
		PSI  <- matrix (0, nrow=n + bdt$Nnode, ncol = n )
		PSI[cbind(1:n, 1:n)] <- 1
		for (inode in 1:length(ancs_postorder) )
		{ # clade (a, (u,v))
			a <- ancs_postorder[inode]
			if ( bdt$heights[a] < maxTMRCA )
			{
				u <- node_daughters[a, 1]
				v <- node_daughters[a, 2]
				dpsi_au <- exp(-( Lambda.h( bdt$heights[a] ) - Lambda.h(bdt$heights[u]) )) 
				dpsi_av <- exp(-( Lambda.h( bdt$heights[a] ) - Lambda.h(bdt$heights[v]) )) 
				PSI[a,] <- dpsi_au * PSI[u,] + dpsi_av * PSI[v,] 
				utips <- which(PSI[u,] > 0)
				vtips <- which(PSI[v,] > 0)
				utipsW <- tips2Wcoords[utips]
				vtipsW <- tips2Wcoords[vtips]
				psi_a <- PSI[a,] 
				W <- updateWCpp( W, psi_a , utips, vtips, utipsW, vtipsW )
				if (F)
				{
					for (ut in utips) for (vt in vtips){
						uname <- bdt$tip.label[ut]
						vname <- bdt$tip.label[vt]
						W[uname, vname] = W[vname, uname] <- W[uname, vname] + PSI[a, ut] * PSI[a, vt] /2
					}
				}
				PSI[a,] <- PSI[a,] / 2
			}
		}
	}
	W <- W / ntrees
	missingInfector <- setNames( 1 - colSums( W), colnames(W))
	list( W = W, missingInfector = missingInfector )
}




phylo.source.attribution.multiDeme.model <- function( tree
  , sampleTimes
  , sampleStates
  , maxHeight
  , theta, demographic.process.model, x0, t0
  , res = 1e3
  , treeErrorTol = 1e-2
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = FALSE
) 
{
	bdt <- DatedTree( tree, sampleTimes , sampleStates = sampleStates, tol = treeErrorTol)
	
	tfgy <- demographic.process.model( theta, x0, t0, bdt$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	phylo.source.attribution.multiDeme.fgy( bdt
	  , maxHeight
	  , tfgy
	  , treeErrorTol = 1e-3
	  , timeOfOriginBoundaryCondition = FALSE
	) 
}



phylo.source.attribution.multiDeme.fgy <- function( dated_tree
  , maxHeight
  , tfgy
  , treeErrorTol = 1e-2
  , timeOfOriginBoundaryCondition = FALSE
  , mode = 2
) 
{
	bdt <- dated_tree
	mel <- min(bdt$edge.length)
	if (mel <= 0) stop('Tree has a zero edge length. Try setting minEdgeLength to a value greater than zero. \n')
	bdt <- reorder.phylo( bdt, 'postorder' )
	bdt$heights <- signif( bdt$heights, digits = floor( 1 / bdt$maxHeight /10 )  +  6 ) #TODO move to DatedTree
	
	times <- tfgy[[1]]
	Fs <- tfgy[[2]]
	Gs <- tfgy[[3]]
	Ys <- tfgy[[4]]
	
	m <- nrow(Fs[[1]])
	if (m < 2)  stop('Currently only models with at least two demes are supported')
	DEMES <- names(Ys[[1]] )
	# for unstructured models
	if (m == 2 & DEMES[2]=='V2' & ncol(bdt$sampleStates) == 1){
		stop('multiDeme version of source attribution should use sampleStates, but none provided.')
	}
	
	#demographic model should be in order of decreasing time:
	fgyi <- 1:length(Fs)
	if (times[2] > times[1])
		fgyi <- length(Fs):1
	#
	heights <- times[fgyi[1]] - times
	
	## note bdt postorder
	ndesc <- rep(0, bdt$n + bdt$Nnode )
	for (iedge in 1:nrow(bdt$edge)){
		a <- bdt$edge[iedge,1]
		u <- bdt$edge[iedge,2]
		ndesc[a] <- ndesc[a] + ndesc[u] + 1
	}
	
	## definition is complicated since multiple nodes can exist at a given height
	uniqhgts <- sort( unique(bdt$heights) )
	eventHeights <- rep(NA, (bdt$n+bdt$Nnode) )
	eventIndicatorNode <- rep(NA, bdt$n + bdt$Nnode )
	events <- rep(NA, bdt$n + bdt$Nnode )
	hgts2node <- lapply( uniqhgts, function(h) which( bdt$heights==h) )
	k <- 1
	l <- 1
	for (h in uniqhgts){
		us  <- hgts2node[[k]]
		if (k < length(hgts2node) | length(us) == 1)
		{ # do not count events for polytomous root (multiple instances of maxheight)
			if (length(us) > 1){
				i_us <- sort( index.return=TRUE, decreasing=FALSE , ndesc[us] )$ix
				for (u in us[i_us] ){
					eventHeights[l] <- h
					events[l] <- ifelse( u <= bdt$n, 0, 1 )
					eventIndicatorNode[l] <- u
					l <- l + 1
				}
			} else {
				eventHeights[l] <- h
				events[l] <- ifelse( us <= bdt$n, 0, 1 )
				eventIndicatorNode[l] <- us 
				l <- l + 1
			}
		}
		k <- k + 1
	}
	excl <- is.na(eventHeights) | is.na(events) | is.na( eventIndicatorNode )
	events <- events[!excl]
	eventIndicatorNode <- eventIndicatorNode[!excl]
	eventHeights <- eventHeights[!excl]
#~ browser()
print('start source attrib') 
print(date())
#~ browser()
if(mode == 1)
{
	W <- sourceAttribMultiDemeCpp(
	  heights
	  , Fs[fgyi]
	  , Gs[fgyi]
	  , Ys[fgyi]
	  , events
	  , eventIndicatorNode
	  , eventHeights
	  , t(bdt$sampleStates) # m X n 
	  , bdt$daughters
	  , bdt$n
	  , bdt$Nnode
	  , m
	  , TRUE
	  , maxHeight
	)
} else{
	W <- sourceAttribMultiDemeCpp2(
	  heights
	  , Fs[fgyi]
	  , Gs[fgyi]
	  , Ys[fgyi]
	  , events
	  , eventIndicatorNode
	  , eventHeights
	  , t(bdt$sampleStates) # m X n 
	  , bdt$daughters
	  , bdt$n
	  , bdt$Nnode
	  , m
	  , TRUE
	  , maxHeight
	  , 10
	)
}
	
	W$donor <- bdt$tip.label[ W$donor ]
	W$recip <- bdt$tip.label[ W$recip ]
print('source attrib complete')
print(date())
	W
}


# this variant of the SA function is tailored to clustering of HIV sequences (1 per host)
.cd42stage <- function(cd4)
{
	# from cori paper
	#~ k =1: CD4>=500
	#~ k = 2 : 350<=CD4<500. 
	#~ k = 3: 200<=CD4<350
	#~ k = 4 : CD4<200
	if (is.na(cd4) ) return (NA) 
	if (cd4  >= 500) return (2)
	if (cd4 >= 350) return(3)
	if (cd4 >= 200 ) return(4)
	return(5)
}
phylo.source.attribution.hiv.msm <- function( tree
  , sampleTimes # must use years
  , cd4s # named numeric vector, cd4 at time of sampling 
  , ehi # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  , numberPeopleLivingWithHIV # scalar
  , numberNewInfectionsPerYear # scalar 
  , maxHeight
  , res = 1e3
  , treeErrorTol = 1/12
  , minEdgeLength=1/52
  , mode = 2
) {
print("NOTE : sample times must be in units of years") 	
	# note time units in year
	stage_prog_yrs <- c( .5, 3.32, 2.7, 5.50, 5.06 ) #cori AIDS
	stageprog_rates <- setNames( 1 / (stage_prog_yrs) 
	 , c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')  )
	
	pstartstage <- c( 
	  pstartstage1 = 0.76
	 , pstartstage2 = 0.19
	 , pstartstage3 = 0.05
	 , pstartstage4 = 0
	)
	progparms <- c( stageprog_rates, pstartstage )
	
	DEMENAMES <- paste(sep='', 'stage', 0:4)
	
	# prior probability of each stage 
	# NOTE this would only be approx, since does not account for stage 0 going to stages > 1: 
	#pstage <- (1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')] ) / sum(1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')])
	prEverReach <- cumsum(  pstartstage )
	meanDurNRStages <- (1/progparms[c( 'gamma1', 'gamma2', 'gamma3', 'gamma4')] )
	pstage1 <- c( 1/progparms['gamma0'], prEverReach * meanDurNRStages)
	pstage <- pstage1 / sum(pstage1)
	pstage[1] <- (numberNewInfectionsPerYear/2) / (numberPeopleLivingWithHIV)
	pstage <- pstage / sum(pstage) 
	
	sampleStates <- matrix( 0, nrow = length(tree$tip.label), ncol = 5 )
	rownames(sampleStates) <- names(sampleTimes)
	colnames(sampleStates) <- DEMENAMES
	for ( tl in names(sampleTimes)){
		stage <- .cd42stage( cd4s[tl] )
		if (is.na( stage) ){
			sampleStates[tl, ] <- pstage 
		}
		else if (stage > 2){
			sampleStates[tl, stage] <- 1
		} else{
			sampleStates[tl, 1] <- pstage[1] / (pstage[1] + pstage[2] )
			sampleStates[tl, 2] <- pstage[2] / (pstage[1] + pstage[2] )
		}
		if (!is.na( ehi[tl] )){
			if (ehi[tl]) {
				sampleStates[tl,1] <- 1
				sampleStates[tl,2:5] <- 0
			}
		}
	}
	
	# make tfgy 
	times <- seq(0 , maxHeight+1, length.out = res)
	Ys <- lapply( 1:res, function(i) setNames( numberPeopleLivingWithHIV * pstage, DEMENAMES ) )
	Fs <- lapply( 1:res, function(i) {
		FF <- matrix ( 0, nrow = 5, ncol = 5)
		rownames(FF) = colnames(FF) <- DEMENAMES
		FF[, 1] <- numberNewInfectionsPerYear * pstage
		FF
	})
	Gs <- lapply( 1:res, function(i ){
		GG <- matrix(0, nrow = 5, ncol = 5 )
		rownames(GG) = colnames(GG) <- DEMENAMES
		GG[1,2:5] <- pstage[1] * numberPeopleLivingWithHIV * pstartstage
		GG[2,3]   <- pstage[2] * numberPeopleLivingWithHIV * progparms['gamma1'] 
		GG[3,4]   <- pstage[3] * numberPeopleLivingWithHIV * progparms['gamma2'] 
		GG[4,5]   <- pstage[4] * numberPeopleLivingWithHIV * progparms['gamma3'] 
		GG
	})
	tfgy <- list( times, Fs, Gs, Ys )
	
	bdt <- DatedTree( tree , sampleTimes , sampleStates = sampleStates , tol = treeErrorTol , minEdgeLength=minEdgeLength)
	mel <- min(bdt$edge.length)
	if (mel <= 0) stop('Tree has a zero edge length. Try setting minEdgeLength to a value greater than zero. \n')
	phylo.source.attribution.multiDeme.fgy( bdt
	  , maxHeight
	  , tfgy
	  , treeErrorTol = 1e-3
	  , timeOfOriginBoundaryCondition = FALSE
	  , mode = mode
	)
}

phylo.source.attribution.hiv.het <- function( tree
  , sampleTimes # must use years
  , sex # m or f for each sample, no NA allowed
  , cd4s # named numeric vector, cd4 at time of sampling 
  , ehi # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
  , numberPeopleLivingWithHIV_m # scalar
  , numberPeopleLivingWithHIV_f # scalar
  , numberNewInfectionsPerYear_m # scalar ; new infections *in* men, not *by* men
  , numberNewInfectionsPerYear_f # scalar 
  , maxHeight
  , res = 1e3
  , treeErrorTol = 1e-2
  , minEdgeLength=1/52
) {
cat("NOTE : sample times must be in units of years\n")
	if (length(sampleTimes)!=length(tree$tip.label)) stop('Sample time must be defined for all tips in tree')
	if(is.null(names(sampleTimes))) names(sampleTimes) <- tree$tip.label
	if (any(is.na(sex))) stop('sex variable must be defined for all samples')
	if (length(sex)!=length(sampleTimes)) stop('sex variable must be defined for all samples')
	if(is.null(names(sex))) names(sex) <- names(sampleTimes)
	
	if(is.null(names(cd4s))) names(cd4s) <- names(sampleTimes)
	
	if(is.null(names(ehi))) names(ehi) <- names(sampleTimes)
	
#~ 	if (length(setdiff( names(sampleTimes), names(sex) )
	# note time units in year
	stage_prog_yrs <- c( .5, 3.32, 2.7, 5.50, 5.06 ) #cori AIDS
	stageprog_rates <- setNames( 1 / (stage_prog_yrs) 
	 , c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')  )
	
	pstartstage <- c( 
	  pstartstage1 = 0.76
	 , pstartstage2 = 0.19
	 , pstartstage3 = 0.05
	 , pstartstage4 = 0
	)
	progparms <- c( stageprog_rates, pstartstage )
	
	DEMENAMES <- c( paste(sep='', 'm_stagem', 0:4)
	 , paste(sep='', 'f_stage', 0:4) )
	
	# prior probability of each stage 
	# NOTE this would only be approx, since does not account for stage 0 going to stages > 1: 
	#pstage <- (1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')] ) / sum(1/progparms[c('gamma0', 'gamma1', 'gamma2', 'gamma3', 'gamma4')])
	prEverReach <- cumsum(  pstartstage )
	meanDurNRStages <- (1/progparms[c( 'gamma1', 'gamma2', 'gamma3', 'gamma4')] )
	pstage1 <- c( 1/progparms['gamma0'], prEverReach * meanDurNRStages)
	pstage <- pstage1 / sum(pstage1)
	pstage[1] <- ((numberNewInfectionsPerYear_f+numberNewInfectionsPerYear_m)/2) / (numberPeopleLivingWithHIV_m + numberPeopleLivingWithHIV_f)
	pstage <- pstage / sum(pstage) 
	
	sampleStates <- matrix( 0, nrow = length(tree$tip.label), ncol = 5*2 )
	rownames(sampleStates) <- names(sampleTimes)
	colnames(sampleStates) <- DEMENAMES
	for ( tl in names(sampleTimes)){
		stage <- .cd42stage( cd4s[tl] )
		
		i <- ifelse( sex[tl]=='m' , 0 ,  5)
		if (is.na( stage) ){
			i <- 6:10
			if ( sex[tl]=='m' ){
				i <- 1:5
			}
			sampleStates[tl,i ] <- pstage 
			if (!is.na( ehi[tl] ) ){
				if (ehi[tl]) {
					i <- ifelse( sex[tl]=='m' , 0 ,  5)
					sampleStates[tl,1 + i] <- 1
					sampleStates[tl,2:5 + i] <- 0
				}
			}
		} else if (stage > 2){
			sampleStates[tl, stage+i] <- 1
		} else{
			sampleStates[tl, 1+i] <- pstage[1] / (pstage[1] + pstage[2] )
			sampleStates[tl, 2+i] <- pstage[2] / (pstage[1] + pstage[2] )
		}
		if (!is.na( ehi[tl] ) & !is.na(stage)){
			if (ehi[tl]) {
				sampleStates[tl,1 + i] <- 1
				sampleStates[tl,2:5 + i] <- 0
			}
		}
	}
	
	# make tfgy 
	times <- seq(0 , maxHeight+1, length.out = res)
	Ys <- lapply( 1:res, function(i) setNames( c(numberPeopleLivingWithHIV_m * pstage,numberPeopleLivingWithHIV_m * pstage) , DEMENAMES)  ) 
	Fs <- lapply( 1:res, function(i) {
		FF <- matrix ( 0, nrow = 5*2, ncol = 5*2)
		rownames(FF) = colnames(FF) <- DEMENAMES
		FF[1:5, 6] <- numberNewInfectionsPerYear_m * pstage
		FF[6:10, 1] <- numberNewInfectionsPerYear_f * pstage
		FF
	})
	Gs <- lapply( 1:res, function(i ){
		GG <- matrix(0, nrow = 5*2, ncol = 5*2 )
		rownames(GG) = colnames(GG) <- DEMENAMES
		GG[1,2:5] <- pstage[1] * numberPeopleLivingWithHIV_m * pstartstage
		GG[2,3]   <- pstage[2] * numberPeopleLivingWithHIV_m * progparms['gamma1'] 
		GG[3,4]   <- pstage[3] * numberPeopleLivingWithHIV_m * progparms['gamma2'] 
		GG[4,5]   <- pstage[4] * numberPeopleLivingWithHIV_m * progparms['gamma3'] 
		
		GG[5+1,5+2:5] <- pstage[1] * numberPeopleLivingWithHIV_f * pstartstage
		GG[5+2,5+3]   <- pstage[2] * numberPeopleLivingWithHIV_f * progparms['gamma1'] 
		GG[5+3,5+4]   <- pstage[3] * numberPeopleLivingWithHIV_f * progparms['gamma2'] 
		GG[5+4,5+5]   <- pstage[4] * numberPeopleLivingWithHIV_f * progparms['gamma3'] 
		
		GG
	})
	tfgy <- list( times, Fs, Gs, Ys )
	bdt <- DatedTree( tree , sampleTimes , sampleStates = sampleStates , tol = treeErrorTol, minEdgeLength=minEdgeLength )
	mel <- min(bdt$edge.length)
	if (mel <= 0) stop('Tree has a zero edge length. Try setting minEdgeLength to a value greater than zero. \n')
	phylo.source.attribution.multiDeme.fgy( bdt
	  , maxHeight
	  , tfgy
	  , treeErrorTol = 1e-3
	  , timeOfOriginBoundaryCondition = FALSE
	)
	
}

