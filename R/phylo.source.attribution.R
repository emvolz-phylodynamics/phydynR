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
		Lambda <- ode( c(Lambda = 0 ), times = heights , func = d.Lambda, parms = list( bdt = bdt ),method = 'adams' )[,2]
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
