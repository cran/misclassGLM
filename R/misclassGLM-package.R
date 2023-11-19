#' @keywords internal
#' @aliases misclassGLM-package
"_PACKAGE"

#' Estimates models that extend the standard GLM to take
#' misclassification into account. The models require side information from a secondary data set
#' on the misclassification process, i.e. some sort of misclassification
#' probabilities conditional on some common covariates.
#' A detailed description of the algorithm can be found in
#' Dlugosz, Mammen and Wilke (2015) \url{https://www.zew.de/publikationen/generalised-partially-linear-regression-with-misclassified-data-and-an-application-to-labour-market-transitions}.
#'
#' The two main functions are \code{\link{misclassGLM}} and \code{\link{misclassMlogit}}.
#'
#' @import stats
#' @import Matrix
#' @import MASS
#' @import ucminf
#' @import numDeriv
#' @import foreach
#' @import mlogit
#' @useDynLib misclassGLM
NULL
