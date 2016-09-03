cfmlogitValidation <- function(par, Y, X, W, SetM) .Call("cfmlogitValidation", par, Y, X, W, SetM)
cgmlogitValidation <- function(par, Y, X, W, SetM) .Call("cgmlogitValidation", par, Y, X, W, SetM)


#predict.mlogit <- function(object, newdata, ...) {
#  if (!inherits(object, "mlogit")) stop("estimate is no mlogit result")
#
#  type <- match.arg(type)
#
#  n <- nrow(newdata)
#  k <- length(object$coefficients) / (ncol(newdata) + 1)
#
#  if (k != floor(k)) stop("k does not fit")
#
#  get.eta <- function(X, beta, k) {
#    d <- ncol(X)
#    ret <- do.call(cbind, lapply(1:k, function(l) {
#      X %*% beta[((0:d) * k) + l][-1] + beta[l]
#    }))
#    return(ret)
#  }
#
#  eta <- get.eta(as.matrix(newdata[,, drop = FALSE]), object$coefficients, k)
#
#  summen <- rowSums(exp(eta)) + 1
#
#  ret <- cbind(exp(eta), rep.int(1, n)) / summen
#
#  tmp <- c()
#  for (j in 1:k) {
#    tmp <- c(tmp, paste("M", j, concat = ""))
#  }
#  dimnames(ret)[[2]] <- c(tmp, "M0")
#
#  return(ret)
#}
