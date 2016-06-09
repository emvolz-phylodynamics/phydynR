\name{colik}
\alias{colik}
\title{ Compute a coalescent log likelihood. }
\description{Computes the log likelihood using a coalescent (or structured coalescent) genealogical model based on a user-supplied demographic process. }
\usage{
\method{colik}{default}(bdt, theta, demographic.process.model, x0, t0, res = 1000,  integrationMethod = "adams", timeOfOriginBoundaryCondition = TRUE,  AgtYboundaryCondition = TRUE) 
}
\arguments{
\item{bdt}{A DatedTree object. }
\item{theta}{ A named numeric vector or named list of parameter values used by the demographic model.  }
\item{demographic.process.model}{ See ?build.demographic.process  }
\item{x0}{ A named vector of initial conditions required by the model. This includes demes and any other dynamic variables. }
\item{t0}{ The time of origin of the process. Should predate the root of the tree. }
\item{res}{ Integer number of time steps to use when simulating model.  }
\item{integrationMethod}{If simulating an ODE model, this provides the integration routine corresponding to options in deSolve. }
\item{timeOfOriginBoundaryCondition }{ If TRUE, will return -Inf if the root of the tree precedes the time of origin. }
\item{maxHeight}{ Will only count internode intervals in the likelihood that occur after maxHeight years before present. Useful for large trees and when you don't want to model the entire demographic history.   }
\item{expmat}{ if TRUE, the likelihood will be computed using matrix exponentiation rather than ODE solvers. This will improve speed but decrease accuracy. NOT YET IMPLEMENTED.   }
\item{finiteSizeCorrection}{ If TRUE, update states of lineages not involved in coalescent event. }
\item{forgiveAgtY}{If number of extant lineages exceeds simulated population size, return -Inf if this value is zero, or forgive the discrepancy if zero. If between zero and one, only forgive the discrepancy if this proportion of lineages is less than the given value. }
\item{AgtY_penalty}{If number of extant lineages exceeds simulated population size, penalise likelihood with value L*AgtY_penalty where L is the cumulative coalescent rate within the given internode interval. 0<= AgtY_penalty <= Inf.}
\item{returnTree}{If TRUE, a copy of the tree is also returned, which includes the inferred states of lineages and likelihood terms at each internal node.}
}
\value{
The coalescent (or structured coalescent) log likelihood (numeric). 
}
\author{Erik Volz}
\examples{
# A simple exponential growth model with birth rates beta and death rates gamma: 
dm <- build.demographic.process(births=c( I = 'parms$beta * I' )
  , deaths = c( I = 'parms$gamma * I' )
  , parameterNames=c('beta', 'gamma') 
  , rcpp=FALSE
  , sde=FALSE)
# simulate a tree based on the model: 
tre <- sim.co.tree(   list( beta = 1.5, gamma = 1 )
  , dm
  , x0  = c(I = 1 )
  , t0 = 0
  , sampleTimes = seq(10, 15, length.out=50)
  , res = 1000
)
# Compute a likelihood
colik(tre
  , list( beta = 1.5, gamma = 1)
  , dm
  , x0 = c( I = 1 )
  , t0 = -1
  , res = 1e3
)
}
\keyword{}
