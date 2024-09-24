# phydynR

phydynR is an R package that provides functions for defining infectious disease models or ecological process models in terms of ordinary differential equations (ODEs) or stochastic differential equations (SDEs).

phydynR also provides functions for simulating genealogies conditional on process models and functions for computing likelihoods of trees conditional on a process model.

## Installation

```r
# You will need to install the R package devtools 
# (https://github.com/r-lib/devtools)

install.packages("devtools")
devtools::install_github("emvolz-phylodynamics/phydynR")
```

## How to use it?

* We recommend that you read the [Get started](articles/phydynR.html) link to 
understand a simple example using phydynR.

* You can later explore other tutorials:
  
  - [Estimating transmission rates using the SIR model](articles/sir_model.html)
  - [Estimating transmission rates using a slightly more complex model for HIV](articles/HIV_epidemics.html)
  - [Estimating genealogies with an epidemiological coalescent model](articles/simulate_genealogies.html)
  - You can try a more complex example of using phydynR using this 
  [tutorial](articles/SenegalHIVmodel.html).

## Authors

phydynR has been developed by [Erik M. Volz](https://profiles.imperial.ac.uk/e.volz)

## References

If you would like to understand more about the models implemented in phydynR,
check out this paper [Volz, 2012](http://www.genetics.org/content/190/1/187)
