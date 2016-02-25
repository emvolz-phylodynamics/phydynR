require(phydynR)

births <- c(I = 'parms$beta * I' )
deaths <- c(I = 'parms$gamma * I' )

dm <- build.demographic.process(births, deaths = deaths
  , parameterNames=c('beta', 'gamma') , rcpp=FALSE, sde = TRUE)

show.demographic.process( dm
 , theta = list( beta = 1.5, gamma = 1 )
 , x0  = c( I = 1 )
 , t0 = 0
 , t1 = 10 
) #TODO legend labels shows V2


dm.det <- build.demographic.process(births = c(I = 'beta*I')
  , deaths = c(I = 'gamma * I')
  #, nonDemeDynamics = NULL 
  , parameterNames=c('beta', 'gamma') 
  , rcpp=TRUE, sde = FALSE)

show.demographic.process( dm.det
 , theta = list( beta = 1.5, gamma = 1 )
 , x0  = c( I = 1 )
 , t0 = 0
 , t1 = 10 
 , log = 'y'
) #TODO legend labels shows V2


tre <- sim.co.tree(   list( beta = 1.5, gamma = 1 )
  , dm.det
  , x0  = c(I = 1 )
  , t0 = 0
  , sampleTimes = seq(5, 15, length.out=100)
  , res = 100 
) #TODO warnings




objfun <- function( of_theta, .tre)
{
	lnbeta <- of_theta['lnbeta']
	lnI0 <- of_theta['lnI0']
	b <- exp(lnbeta)
	I0 <- exp(lnI0)
	print( c( b, I0))
	-colik(.tre
	  , list( beta = unname(b), gamma = 1)
	  , dm.det
	  , x0 = c( I = unname(I0) )
	  , t0 = -1
	  , res = 1e3
	  , timeOfOriginBoundaryCondition = FALSE
	  , AgtYboundaryCondition = TRUE # important/necessary
	)
}



#~ fit <- optim( par = c( lnbeta = log(2), lnI0 = log(1))
#~  , fn = objfun
#~  , control = list( trace = 6 )
#~ )

if (F)
{
	nfits <- 10
	fits <- lapply( 1:nfits, function(ifit){
		.tre <- sim.co.tree(   list( beta = 1.5, gamma = 1 )
		  , dm.det
		  , x0  = c(I = 1 )
		  , t0 = 0
		  , sampleTimes = seq(5, 15, length.out=100)
		  , res = 100 
		) #TODO warnings
		
		optim( par = c( lnbeta = log(2), lnI0 = log(1))
		  , fn = objfun
		  , control = list( trace = 6 )
		  , .tre = .tre
		)
	})
	betas <- unname( sapply(fits, function(fit) exp(fit$par['lnbeta']) ) )
}

if (F)
{
	.objfun <- function( of_theta, .tre)
	{
		b <- (of_theta['beta'])
		I0 <- (of_theta['I0'])
		if (b < 0 ) return(Inf)
		if (I0 < 0) return (Inf)
		print( c( b, I0))
		-colik(.tre
		  , list( beta = unname(b), gamma = 1)
		  , dm.det
		  , x0 = c( I = unname(I0) )
		  , t0 = -1
		  , res = 1e3
		  , timeOfOriginBoundaryCondition = FALSE
		  , AgtYboundaryCondition = TRUE # important/necessary
		)
	}
		.tre <- sim.co.tree(   list( beta = 1.5, gamma = 1 )
		  , dm.det
		  , x0  = c(I = 1 )
		  , t0 = 0
		  , sampleTimes = seq(5, 15, length.out=100)
		  , res = 100 
		) 
		
		.o <- optim( par = c( beta = (2), I0 = (.1))
		  , fn = .objfun
		  , control = list( trace = 6 )
		  , .tre = .tre
		)

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
  , AgtYboundaryCondition = TRUE
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
		t0 <- ifelse(is.null(t0), 0, 
		    ifelse( is.na( of_theta['t0'] ), 0, of_theta['t0'] )
		  )
		print( of_theta )
		.theta <- theta
		.theta[ names(of_theta)] <- unname( of_theta )
		-colik(tre
		  , as.list( .theta )
		  , dm
		  , x0 = x0
		  , t0 = t0
		  , res = 1e3
		  , timeOfOriginBoundaryCondition = timeOfOriginBoundaryCondition
		  , AgtYboundaryCondition = AgtYboundaryCondition # important/necessary
		)
	}
	of_theta <- start[est_pars] 
	optim( par = of_theta, fn = .objfun, control = control )
}


if (T)
{
	fit <- optim.colik(tre
	  , dm.det
	  , start = c( beta = 2, gamma = 1, I = 1, t0 = -1)
	  , est_pars = c( 'beta', 'I')
	  , ic_pars = c('I')
	  , t0 = -1
	  , parm_lowerBounds = c( I = 0, beta = 0)
	  , timeOfOriginBoundaryCondition = FALSE
	  , AgtYboundaryCondition = TRUE
	  , control = list()
	)
}

if (F){
source('../R/colik.R')
sourceCpp( '../src/colik.cpp' )
b <- 1.5
I0 <- 1
print(system.time( { print(
colik(tre
	  , c( beta = unname(b), gamma = 1)
	  , dm.det
	  , x0 = c( I = unname(I0) )
	  , t0 = -1
	  , res = 1e3
	  , timeOfOriginBoundaryCondition = FALSE
	  , AgtYboundaryCondition = FALSE
	)
)}))
}
