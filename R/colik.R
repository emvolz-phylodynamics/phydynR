require(RcppArmadillo)
require(Rcpp)
#~ sourceCpp('colik.cpp') 

colik <- function(bdt, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='lsoda'
  , timeOfOriginBoundaryCondition = TRUE
  , AgtYboundaryCondition = TRUE # can also be numeric > 0
  , maxHeight = Inf 
  , expmat = FALSE
  , forgiveAgtY = .2 #can be NA
) 
{
	tfgy <- demographic.process.model( theta, x0, t0, bdt$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	# TODO need to c-ify this code; should probably go with likelihood and should check Y>A for each deme
	if (!is.na( forgiveAgtY )){
		if (forgiveAgtY < 0 | forgiveAgtY > 1 ){
			stop('forgiveAgtY must be in (0,1)')
		} 
		soinh <- sort( bdt$heights[(1+bdt$n):(bdt$n+bdt$Nnode)] ) 
		A.h <- function(h) {  sum(bdt$sortedSampleHeights <= h) - sum(soinh <=h)  }
		#for (i in 1:length(tfgy[[1]])){
		x <- sum( sapply( 1:length(tfgy[[1]]), function(i){
			t <- tfgy[[1]][i]
			h <- bdt$maxSampleTime - t
			if (h < 0) return(FALSE)
			if (sum(soinh < h)  > (1-forgiveAgtY)*bdt$Nnode) return(FALSE)
			y <- tfgy[[4]][[i]]
			if (sum(y) < A.h( h )) return(TRUE)
			return(FALSE)
		}))
#~ browser()
		if (sum(x) > 0 ) return(-Inf)
	}
	
	colik.fgy(bdt
	  , tfgy
	  , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
	  , AgtYboundaryCondition = AgtYboundaryCondition
      , maxHeight = maxHeight
      , expmat = expmat
    ) 
}

colik.fgy <- function (bdt, tfgy, timeOfOriginBoundaryCondition = TRUE, 
    AgtYboundaryCondition = TRUE, maxHeight = Inf, expmat = FALSE) 
{
    
    bdt <- reorder.phylo(bdt, "postorder")
    bdt$heights <- signif(bdt$heights, digits = floor(1/bdt$maxHeight/10) + 
        6)
    times <- tfgy[[1]]
    Fs <- tfgy[[2]]
    Gs <- tfgy[[3]]
    Ys <- tfgy[[4]]
    m <- nrow(Fs[[1]])
    if (m < 2) 
        stop("Currently only models with at least two demes are supported")
    DEMES <- names(Ys[[1]])
    if (m == 2 & DEMES[2] == "V2" & ncol(bdt$sampleStates) == 
        1) {
        bdt$sampleStates <- cbind(bdt$sampleStates, rep(0, bdt$n))
        bdt$sortedSampleStates <- cbind(bdt$sortedSampleStates, 
            rep(0, bdt$n))
    }
    fgyi <- 1:length(Fs)
    if (times[2] > times[1]) 
        fgyi <- length(Fs):1
    
    t0 <- times[fgyi[length(fgyi)]]
    if (timeOfOriginBoundaryCondition) {
        if ((bdt$maxSampleTime - t0) < bdt$maxHeight) 
            return(-Inf)
    }
    
    heights <- times[fgyi[1]] - times[fgyi]
    ndesc <- rep(0, bdt$n + bdt$Nnode)
    for (iedge in 1:nrow(bdt$edge)) {
        a <- bdt$edge[iedge, 1]
        u <- bdt$edge[iedge, 2]
        ndesc[a] <- ndesc[a] + ndesc[u] + 1
    }
    uniqhgts <- sort(unique(bdt$heights))
    eventHeights <- rep(NA, (bdt$n + bdt$Nnode))
    eventIndicatorNode <- rep(NA, bdt$n + bdt$Nnode)
    events <- rep(NA, bdt$n + bdt$Nnode)
    hgts2node <- lapply(uniqhgts, function(h) which(bdt$heights == 
        h))
    k <- 1
    l <- 1
    for (h in uniqhgts) {
        us <- hgts2node[[k]]
        if (k < length(hgts2node) | length(us) == 1) {
            if (length(us) > 1) {
                i_us <- sort(index.return = TRUE, decreasing = FALSE, 
                  ndesc[us])$ix
                for (u in us[i_us]) {
                  eventHeights[l] <- h
                  events[l] <- ifelse(u <= bdt$n, 0, 1)
                  eventIndicatorNode[l] <- u
                  l <- l + 1
                }
            }
            else {
                eventHeights[l] <- h
                events[l] <- ifelse(us <= bdt$n, 0, 1)
                eventIndicatorNode[l] <- us
                l <- l + 1
            }
        }
        k <- k + 1
    }
    excl <- is.na(eventHeights) | is.na(events) | is.na(eventIndicatorNode) | 
        (eventHeights > maxHeight)
    events <- events[!excl]
    eventIndicatorNode <- eventIndicatorNode[!excl]
    eventHeights <- eventHeights[!excl]
    if (!expmat) {
        ll <- colik2cpp(heights, Fs[fgyi], Gs[fgyi], Ys[fgyi], 
            events, eventIndicatorNode, eventHeights, t(bdt$sortedSampleStates), 
            bdt$daughters, bdt$n, bdt$Nnode, m, AgtYboundaryCondition)
    }
    else {
#~     stop('Error: *expmat* option is not yet implemented. ')
		warning('*expmat* option is not debugged ')
        ll <- colik3cpp(heights, Fs[fgyi], Gs[fgyi], Ys[fgyi], 
            events, eventIndicatorNode, eventHeights, t(bdt$sortedSampleStates), 
            bdt$daughters, bdt$n, bdt$Nnode, m, AgtYboundaryCondition)
    }
    ll
}


colik.rcolgem <- function(bdt, theta, demographic.process.model, x0, t0, res = 1e3
  , integrationMethod='rk4'
  , maxHeight = Inf 
  , forgiveAgtY = .2 #can be NA
){
	require(rcolgem)
	tfgy <- demographic.process.model( theta, x0, t0, bdt$maxSampleTime, res = res, integrationMethod=integrationMethod) 
	
	# TODO need to c-ify this code; should probably go with likelihood and should check Y>A for each deme
	if (!is.na( forgiveAgtY )){
		if (forgiveAgtY < 0 | forgiveAgtY > 1 ){
			stop('forgiveAgtY must be in (0,1)')
		} 
		soinh <- sort( bdt$heights[(1+bdt$n):(bdt$n+bdt$Nnode)] ) 
		A.h <- function(h) {  sum(bdt$sortedSampleHeights <= h) - sum(soinh <=h)  }
		#for (i in 1:length(tfgy[[1]])){
		x <- sum( sapply( 1:length(tfgy[[1]]), function(i){
			t <- tfgy[[1]][i]
			h <- bdt$maxSampleTime - t
			if (h < 0) return(FALSE)
			if (sum(soinh < h)  > (1-forgiveAgtY)*bdt$Nnode) return(FALSE)
			y <- tfgy[[4]][[i]]
			if (sum(y) < A.h( h )) return(TRUE)
			return(FALSE)
		}))
		if (sum(x) > 0 ) return(-Inf)
	}
	
	colik.rcolgem.fgy(bdt
	  , tfgy
      , maxHeight = maxHeight
      , forgiveAgtY = forgiveAgtY
	  , integrationMethod = integrationMethod
    ) 

}

colik.rcolgem.fgy <- function(bdt, tfgy, timeOfOriginBoundaryCondition = TRUE
    , maxHeight = Inf, expmat = TRUE
    , forgiveAgtY = .2
    , integrationMethod = 'rk4')
{
	require(rcolgem)
#~ browser()
	excl <- tfgy[[1]] > bdt$maxSampleTime
	tfgy[[1]] <- tfgy[[1]][!excl]
	tfgy[[2]] <- tfgy[[2]][!excl]
	tfgy[[3]] <- tfgy[[3]][!excl]
	tfgy[[4]] <- tfgy[[4]][!excl]
	if ( tfgy[[1]][2] > tfgy[[1]][1]){
		# should be in order of increasing height = decreasing time
		fgyi <- length(tfgy[[1]]):1
		tfgy[[1]] <- tfgy[[1]][fgyi]
		tfgy[[2]] <- tfgy[[2]][fgyi]
		tfgy[[3]] <- tfgy[[3]][fgyi]
		tfgy[[4]] <- tfgy[[4]][fgyi]
	}
	coalescent.log.likelihood.fgy2(bdt
	, tfgy
	, integrationMethod = integrationMethod
	,  censorAtHeight= ifelse( is.infinite( maxHeight ), FALSE, maxHeight )
	, forgiveAgtY=forgiveAgtY
	, returnTree=FALSE)
#~ browser()
#~ i <- length(tfgy[[1]]):1
#~ 	coalescent.log.likelihood.fgy2(bdt, tfgy[[1]], tfgy[[2]], tfgy[[3]], tfgy[[4]]
#~ 	, integrationMethod = integrationMethod
#~ 	,  censorAtHeight= FALSE
#~ 	, forgiveAgtY=forgiveAgtY
#~ 	, returnTree=FALSE)

}



det.solve.As <- function(times, Fs, Gs, Ys, bdt )
{
	step_size_multiplier <- 1
	m <- nrow(Fs[[1]] )
	deltah <- abs( times[2] - times[1] )
	
	# CO HEIGHTS 
	parms <- list(times = times, Fs = Fs, Gs = Gs, Ys = Ys , m = m, deltah = abs( times[2] - times[1] )
	  , times.size  = length(times ) )
	maxHeight <- bdt$maxSampleTime - times[length(times)]
	ih <- 1
	AL <- rep(0, m + 1)
	tAL <- c()
	h <- 0
	heights <- c()
	L <- c()
	coheights <- c()
	while (ih < length( bdt$sortedSampleHeights )){
		AL <- AL + c(bdt$sortedSampleStates[ih, ], 0)
		if (bdt$sortedSampleHeights[ih] > h){
			h1 <- bdt$sortedSampleHeights[ih+1]
			datimes <- seq(h, h1, length.out = max(2, ceiling( (h1 - h)/( deltah) ))  )
			o <- ode(y = AL, times = datimes, func = dAL, parms = parms,  method = 'lsoda')
			tAL <- rbind( tAL, o )
			AL <- o[nrow(o), 2:ncol(o)]
			heights <- c( heights , datimes )
			
			dL <- o[nrow(o), ncol(o)] - o[1, ncol(o)]
			
			L <- c( L, o[, ncol(o)] )
			h  <- bdt$sortedSampleHeights[ih + 1]
		} 
		ih <- ih + 1
	}
	# last interval 
	AL <- AL + c(bdt$sortedSampleStates[ih, ], 0)
	ntlA <- sum(AL[1:(length(AL)-1)]) # total A at beginning of last interval
	h1 <- maxHeight 
	datimes <- seq(h, h1, length.out = max(2, ceiling( (h1 - h)/(step_size_multiplier * deltah) ))  )
	o <- ode(y = AL, times = datimes, func = dAL, parms = parms,  method = 'lsoda')
	tAL <- rbind( tAL, o )
	AL <- o[nrow(o), 2:ncol(o)]
	heights <- c( heights , datimes )
	lastA <- sum(AL[1:(length(AL)-1)]) # A at end of last interval
	# /CO HEIGHTS
	
	Amat <- sapply( 1:m, function(k ){
		approx( tAL[,1], tAL[,1+k], xout = maxSampleTime - times  , rule = 2)$y
	})
	As <- lapply( 1:nrow(Amat), function(i){
		Amat[i, ]
	})
	
	return (As) # TODO only at times?

}



optim.colik <- function(tre
  , dm
  , start
  , est_pars
  , ic_pars
  , t0 = NULL
  , parm_lowerBounds = c()
  , parm_upperBounds = c()
  , timeOfOriginBoundaryCondition = FALSE
  , AgtYboundaryCondition = TRUE # can also be numeric
  , res = 1e3
  , integrationMethod = 'lsoda' 
  , control = list()
  ,  ... )
{
	theta <- start 
	est_ic_pars <- intersect( est_pars, ic_pars) 
	nonic_pars <- setdiff( est_pars, c( est_ic_pars, 't0') )
	.objfun <- function( of_theta)
	{
		for (pn in names(parm_lowerBounds)){
			if ( of_theta[ pn ] < parm_lowerBounds[pn]) return(Inf)
		}
		for (pn in names(parm_upperBounds)){
			if ( of_theta[ pn ] > parm_upperBounds[pn]) return(Inf)
		}
		x0 <- start[ic_pars]
		x0[ est_ic_pars] <- unname( of_theta[ est_ic_pars ]  )
		if (is.null(t0)){
			if (is.na(of_theta['t0'])){
				t0 <- 0
			} else{
				t0 <- of_theta['t0']
			}
		}
		print( of_theta )
		.theta <- theta
		.theta[ names(of_theta)] <- unname( of_theta )
		-colik(tre
		  , as.list( .theta )
		  , dm
		  , x0 = x0
		  , t0 = t0
		  , res = res
		  , integrationMethod = integrationMethod
		  , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
		  , AgtYboundaryCondition = AgtYboundaryCondition # important/necessary
		)
	}
	of_theta <- start[est_pars] 
	optim( par = of_theta, fn = .objfun, control = control )
}


##########################################
colik4.fgy <- function (bdt, tfgy, timeOfOriginBoundaryCondition = TRUE
    , AgtYboundaryCondition = FALSE #TODO
    , maxHeight = Inf
    , finiteSizeCorrection = TRUE
    , detailedOutput=FALSE) 
{
	if (AgtYboundaryCondition) warning("AgtYboundaryCondition not implemented")
    bdt <- reorder.phylo(bdt, "postorder")
    bdt$heights <- signif(bdt$heights, digits = floor(1/bdt$maxHeight/10) + 
        6)
    times <- tfgy[[1]]
#~     if (times[2] > times[1]) stop("Times must be decreasing")
    Fs <- tfgy[[2]]
    Gs <- tfgy[[3]]
    Ys <- tfgy[[4]]
    m <- nrow(Fs[[1]])
    if (m < 2) 
        stop("Currently only models with at least two demes are supported")
    DEMES <- names(Ys[[1]])
    if (m == 2 & DEMES[2] == "V2" & ncol(bdt$sampleStates) == 
        1) {
        bdt$sampleStates <- cbind(bdt$sampleStates, rep(0, bdt$n))
        bdt$sortedSampleStates <- cbind(bdt$sortedSampleStates, 
            rep(0, bdt$n))
    }
    fgyi <- 1:length(Fs)
    if (times[2] > times[1]) 
        fgyi <- length(Fs):1
    
    t0 <- times[fgyi[length(fgyi)]]
    if (timeOfOriginBoundaryCondition) {
        if ((bdt$maxSampleTime - t0) < bdt$maxHeight) 
            return(-Inf)
    }
    
#~     heights <- times[fgyi[1]] - times[fgyi] #TODO heights shoudl be rel to mst?
    heights <- bdt$maxSampleTime - times[fgyi] #TODO heights shoudl be rel to mst?
    ndesc <- rep(0, bdt$n + bdt$Nnode)
    for (iedge in 1:nrow(bdt$edge)) {
        a <- bdt$edge[iedge, 1]
        u <- bdt$edge[iedge, 2]
        ndesc[a] <- ndesc[a] + ndesc[u] + 1
    }
    uniqhgts <- sort(unique(bdt$heights))
    eventHeights <- rep(NA, (bdt$n + bdt$Nnode))
    eventIndicatorNode <- rep(NA, bdt$n + bdt$Nnode)
    events <- rep(NA, bdt$n + bdt$Nnode)
    hgts2node <- lapply(uniqhgts, function(h) which(bdt$heights == 
        h))
    k <- 1
    l <- 1
    for (h in uniqhgts) {
        us <- hgts2node[[k]]
        if (k < length(hgts2node) | length(us) == 1) {
            if (length(us) > 1) {
#~ browser()
                i_us <- sort(index.return = TRUE, decreasing = FALSE, 
                  ndesc[us])$ix
                for (u in us[i_us]) {
                  eventHeights[l] <- h
                  events[l] <- ifelse(u <= bdt$n, 0, 1)
                  eventIndicatorNode[l] <- u
                  l <- l + 1
                }
            }
            else {
                eventHeights[l] <- h
                events[l] <- ifelse(us <= bdt$n, 0, 1)
                eventIndicatorNode[l] <- us
                l <- l + 1
            }
        }
        k <- k + 1
    }
#~ browser()
#~ bdt$tip.label [ eventIndicatorNode[ events==0 ] ]
    excl <- is.na(eventHeights) | is.na(events) | is.na(eventIndicatorNode) | 
        (eventHeights > maxHeight)
    events <- events[!excl]
    eventIndicatorNode <- eventIndicatorNode[!excl]
    eventHeights <- eventHeights[!excl]

#TODO bug? sum(events)!=Nnode? can appear if time res is small
#~ 	sevents <- events[ events==0 ]
#~ browser()
    o <- colik4cpp( times[fgyi], Fs[fgyi], Gs[fgyi], Ys[fgyi]
	  , bdt$sortedSampleHeights,  t(bdt$sortedSampleStates)
     , events, eventIndicatorNode
     , eventHeights, bdt$daughters
     , bdt$maxSampleTime
     ,  m
     , finiteSizeCorrection
     , maxHeight
    )
#~     print('Complete')
#~     print(o$loglik)
#~ browser()
	if (detailedOutput){return(o)}
    o$loglik
}
