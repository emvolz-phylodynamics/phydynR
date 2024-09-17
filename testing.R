library(phydynR)

tree <- read.tree(system.file("extdata/sirModel0.nwk", package = "phydynR"))
tree$tip.label <- unlist(lapply(tree$tip.label, function(x) paste("seq", x, sep = "_")))


parms_truth <- list( beta = 0.00020002, 
                     gamma = 1, 
                     S0 = 9999, 
                     t0 = 0 )

sampleTimes <- rep(12, 75)
sampleTimes <- setNames(sampleTimes, tree$tip.label)

bdt <- DatedTree( phylo = tree, sampleTimes = sampleTimes)
bdt

names(bdt)
bdt$sampleTimes
bdt$sampleStates


births <- c( I = "parms$beta * S * I" )
deaths <- c( I = "parms$gamma * I" )
nonDemeDynamics <- c("S = -parms$beta * S * I")

# initial number of I and S
x0 <- c(I = 1, S = unname(parms_truth$S) )

# initial t0 (time of origin of the process)
t0 <- bdt$maxSampleTime - max(bdt$heights) - 1
t0 <- 13


dm <- build.demographic.process(births = births,
                                deaths = deaths,
                                parameterNames = names(parms_truth),
                                rcpp = FALSE,
                                sde = FALSE)

colik(tree = bdt,
      theta = parms_truth,
      demographic.process.model = dm,
      x0 = x0,
      t0 = t0,
      res = 1000,
      integrationMethod = "rk4")
