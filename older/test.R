#' @title phydynR
#' @name phydynR
#' @docType package
#' @useDynLib phydynR
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @description \packageDescription{phydynR} Phylodynamic inference and 
#'    simulation in R.
#' @details This package provides an interface to building demographic process 
#'    models (ODEs or SDEs) and simulating genealogies based on these models. 
#'    Low level functions can be used to construct more complex models. 
#'    A structured coalescent likelihood allows for inference of model 
#'    parameters provided a time-scaled genealogy. 
#' @author Erik Volz
#' @section Maintainer: \packageMaintainer{phydynR}
#' @references 
#' 1. Volz, Erik M. Complex population dynamics and the coalescent under 
#'   neutrality. Genetics 190.1 (2012): 187-201.
#'   
#' 2. Volz, Erik M., Katia Koelle, and Trevor Bedford. Viral phylodynamics. 
#'   PLoS Comput Biol 9.3 (2013): e1002947.
#' @examples 
#' See the README for common usage and examples.
"_PACKAGE"