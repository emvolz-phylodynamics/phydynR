\name{build.demographic.process}
\alias{build.demographic.process}
\title{Build a demographic process model (ODE or SDE)}
\description{A model is constructed by supplying birth rates, migration rates, and death rates. The model can be used to simulate demographic histories or as an input to coalescent simulation or likelihood calculation. }
\usage{
\method{build.demographic.process}{default}(births, nonDemeDynamics = NA, migrations = NA, deaths = NA, parameterNames = c(), rcpp = TRUE, sde = FALSE) 
}
\arguments{
\item{births}{ A named character vector or matrix of mathematical expressions that describe model birth rates. Names correspond to each deme. If there is more than one deme, a matrix must be supplied describing birth rates between each pair of demes. Either R or C code can be supplied (see rcpp option). }
\item{nonDemeDynamics}{ A named character vector of mathematical expressions that describe rate of change for dynamic variables which are not demes; for example: an equation for the number susceptible in an SIR model.  }
\item{migrations}{ A named character vector or matrix of mathematical expressions that describe model migration rates. Names correspond to each deme. If there is more than one deme, a matrix must be supplied describing migration rates between each pair of demes. Either R or C code can be supplied (see rcpp option). }
\item{deaths}{ A named character vector of mathematical expressions that describe death rates in each deme.}
\item{parameterNames}{ A character vector providing the names of static parameters used in the model.  }
\item{rcpp}{ If TRUE, the expressions are interpreted as C code using the Rcpp package. }
\item{sde}{ If TRUE, a stochastic differential equation model is constructed; if FALSE, an ordinary differential equation model is constructed. }
}
\value{
An object of class \code{demographic.process}, which is a function that can be used to simulate the model. The model can also be used as an input to tree simulation or likelihood calculation. 
}
\author{Erik Volz}
\examples{
# A simple exponential growth model with birth rates beta and death rates gamma: 
dm <- build.demographic.process(births=c( I = 'parms$beta * I' )
  , deaths = c( I = 'parms$gamma * I' )
  , parameterNames=c('beta', 'gamma') 
  , rcpp=FALSE
  , sde = TRUE
)
# do a simulation and plot the trajectory: 
show.demographic.process( dm
 , theta = list( beta = 1.5, gamma = 1 )
 , x0  = c( I = 1 )
 , t0 = 0
 , t1 = 10 
) 
}
\keyword{}
