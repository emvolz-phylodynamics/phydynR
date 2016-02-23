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




objfun <- function( of_theta)
{
	lnbeta <- of_theta['lnbeta']
	lnI0 <- of_theta['lnI0']
	b <- exp(lnbeta)
	I0 <- exp(lnI0)
	print( c( b, I0))
	-colik(tre
	  , list( beta = unname(b), gamma = 1)
	  , dm.det
	  , x0 = c( I = unname(I0) )
	  , t0 = -1
	  , res = 1e3
	  , timeOfOriginBoundaryCondition = FALSE
	  , AgtYboundaryCondition = TRUE # important/necessary
	)
}



fit <- optim( par = c( lnbeta = log(2), lnI0 = log(1))
 , fn = objfun
 , control = list( trace = 6 )
)

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
