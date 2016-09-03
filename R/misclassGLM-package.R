#' misclassGLM
#'
#' Estimates models that extend the standard GLM to take
#' misclassification into account. The models require side information from a secondary data set
#' on the misclassification process, i.e. some sort of misclassification
#' probabilities conditional on some common covariates.
#' A detailed description of the algorithm can be found in
#' Dlugosz, Mammen and Wilke (2015) \url{http://www.zew.de/PU70410}.
#'
#' The two main functions are \code{\link{misclassGLM}} and \code{\link{misclassMlogit}}.
#'
#' @docType package
#' @name misclassGLM
#' @import stats
#' @import Matrix
#' @import MASS
#' @import ucminf
#' @import numDeriv
#' @import bigmemory
#' @import foreach
#' @import mlogit
#' @useDynLib misclassGLM
NULL
