# phydynR <img src="man/figures/phydynr.png" align="right" width="200"/>

phydynR is an [R](https://www.r-project.org/) package based on the coalescent
theory in which you can

* define infectious disease models or ecological process models in terms of 
ordinary differential equations (ODEs) or stochastic differential equations (SDEs);

* simulate genealogies conditional on a process model;

* compute likelihoods of phylogenetic trees also conditional on a process model.



## Installation

```r
# You will need to install the R package devtools 
# (https://github.com/r-lib/devtools)

install.packages("devtools")
devtools::install_github("emvolz-phylodynamics/phydynR")
```



## Tutorials

* We recommend that you read the [Get started](http://emvolz-phylodynamics.github.io/phydynR/articles/phydynR.html) to 
understand the basic functions in phydynR.

* You can later explore other tutorials:
  
  - [Estimating transmission rates using the SIR model.](http://emvolz-phylodynamics.github.io/phydynR/articles/sir_model.html)
  - [Estimating transmission rates using a slightly more complex model for HIV.](http://emvolz-phylodynamics.github.io/phydynR/articles/HIV_epidemics.html)
  - [Estimating genealogies with an epidemiological coalescent model.](http://emvolz-phylodynamics.github.io/phydynR/articles/simulate_genealogies.html)
  - [You can try a more complex example of estimating number of infected individuals and transmission rates for HIV in Senegal.](http://emvolz-phylodynamics.github.io/phydynR/articles/SenegalHIVmodel.html)



## Author

phydynR has been developed by [Erik Volz](https://profiles.imperial.ac.uk/e.volz)



## Related softwares

phydynR works on a fixed phylogenetic tree and therefore is *not* a program that 
will estimate the tree for you.

[PhyDyn](https://github.com/mrc-ide/PhyDyn): If you would like to take into 
consideration the uncertainty on the tree estimates, check out our other software 
PhyDyn implemented in [BEAST 2](https://www.beast2.org/). For details on how to use
it start [here](https://github.com/mrc-ide/PhyDyn/wiki). PhyDyn is substantially 
slow to run.

[Coalescent.jl](https://emvolz.github.io/Coalescent.jl/dev/intro/): implements
similar specification of demographic process but it is substantially faster than
phydynR. Coalescent.jl is based on the [Julia](https://julialang.org/) programming language.



## References

If you would like to understand more about the models implemented in phydynR,
check out [Volz, 2012](http://www.genetics.org/content/190/1/187)
