% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeSimulatorCpp2.R
\name{show.demographic.process}
\alias{show.demographic.process}
\title{Show demographic process}
\usage{
show.demographic.process(
  demo.model,
  theta,
  x0,
  t0,
  t1,
  res = 1000,
  integrationMethod = "lsoda",
  legend_position = "bottomright",
  ...
)
}
\arguments{
\item{demo.model}{An object of class \code{demographic.process}, which is a function
that can be used to simulate the model.}

\item{theta}{A named numeric vector or named list of parameter values used by
the demographic model}

\item{x0}{A named vector of initial conditions required by the model. This
includes demes and any other dynamic variables.}

\item{t0}{Initial time for the simulation (to show the demographic process)}

\item{t1}{Final time for tge simulation (to show the demographic process)}

\item{res}{Integer number of time steps to use when simulating model.}

\item{integrationMethod}{If simulating an ODE (ordinary differential equation)
model, this provides the integration routine corresponding to options in
deSolve.}

\item{legend_position}{String for position of legend in final plot. Default
is set to "bottomright"}

\item{...}{Additional arguments that can be passed to the function, 
such as graphical parameters.}
}
\description{
Shows a graphical representation of the demographic process given values 
for each parameter of interest, and initial population sizes. The user should
also provide the initial size for the population of interest in the their model.
}
\examples{
# A simple exponential growth model with birth rates beta and death rates gamma:
# I is the number of infected individuals.
dm <- build.demographic.process(births=c(I = 'parms$beta * I'),
                                deaths = c(I = 'parms$gamma * I'),
                                parameterNames=c('beta', 'gamma'),
                                rcpp=FALSE,
                                sde = TRUE)
# Do a simulation and plot the trajectory:
show.demographic.process(dm,
                         theta = list(beta = 1.5, gamma = 1),
                         x0  = c(I = 1),
                         t0 = 0,
                         t1 = 10)
}
\seealso{
\code{\link[phydynR]{build.demographic.process}}
}
