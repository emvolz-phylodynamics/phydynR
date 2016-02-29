\name{sim.co.tree}
\alias{sim.co.tree}
\title{Simulate a coalescent genealogy based on a demographic process. }
\description{A tree is simulated based on a demographic process and user-supplied parameters and initial conditions. The times and types of samples must also be supplied. }
\usage{
\method{sim.co.tree}{default}(theta, demographic.process.model, x0, t0, sampleTimes, sampleStates = NULL, res = 1000, integrationMethod = "adams")
}
\arguments{
\item{theta}{ A named vector of parameter values required by the model. }
\item{demographic.process.model}{ An object of class demographic.process. See build.demographic.process. }
\item{x0}{ A named vector of initial conditions required by the model. This includes demes and any other dynamic variables. }
\item{t0}{ The time of origin of the process }
\item{sampleTimes }{ A numeric vector providing the times of samples. If named, taxon labels will be based on names of the corresponding sample times. }
\item{sampleStates}{  For models with more than one deme, a matrix must be supplied describing the probability of sampling from each deme. Each row corresponds to a sample in the same order as sampleTimes. Each column corresponds to probability of sampling each deme. Column names should be defined which correspond to deme names in the model. Rows should sum to one. }
\item{res}{ Integer number of time steps to use when simulating model.  }
\item{integrationMethod}{If simulating an ODE model, this provides the integration routine corresponding to options in deSolve. }
}
\value{
An object of class DatedTree, which subclasses ape::phylo. The $heights attribute provides the time of each node in the tree. The $maxHeight attribute provides the time of the root. 
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
# plot the tree 
plot( tre )
}
\keyword{}
