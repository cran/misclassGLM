fgaussValidation <- function(par, Y, X, W, SetM) .Call(cfgaussValidation, par, Y, X, W, SetM)
ggaussValidation <- function(par, Y, X, W, SetM) .Call(cggaussValidation, par, Y, X, W, SetM)
flogitValidation <- function(par, Y, X, W, SetM) .Call(cflogitValidation, par, Y, X, W, SetM)
glogitValidation <- function(par, Y, X, W, SetM) .Call(cglogitValidation, par, Y, X, W, SetM)

fprobitValidation <- function(par, Y, X, W, SetM) {
  n <- nrow(X)
  d <- ncol(X)
  p <- d + ncol(SetM) + 1

  lik <- rep_len(0.0, n)

  Xbeta <- rep.int(par[1], n) + X %*% par[2:(d + 1)]

  for (k in 1:nrow(SetM)) {
    tmp <- qnorm(Xbeta + rep.int(SetM[k,] %*% par[(d + 2):p], n), sd = 1)
    idx <- (Y == TRUE)
    tmp[!idx] <- 1 - tmp[!idx]
    lik <- lik + tmp * W[, k]
  }

  return(-1 * sum(log(lik)))
}

gprobitValidation <- function(par, Y, X, W, SetM) {
  # @TODO!!
  n <- nrow(X)
  ms <- nrow(SetM)
  pm <- ncol(SetM)
  px <- ncol(X)

  p <- px + pm + 1

  predicts <- matrix(rep.int(0, n * ms), nrow = n)
  exppredicts <- matrix(rep.int(0, n * ms), nrow = n)

  ret <- rep.int(0.0, p + 1)

  tmp <- rep.int(par[1], n) + X %*% par[2:(px + 1)]

  for (j in 1:ms) {
    predicts[, j] <- tmp + rep.int(SetM[j, ] %*% par[(px + 2):p], n)
    exppredicts[, j] <- exp(predicts[, j] * predicts[, j]) * W[, j]
  }

  tempsumm <- rowSums(exppredicts)

  tmp <- rowSums(predicts * exppredicts) / tempsumm

  ret[1] <- sum(tmp)

  for (k in 1:px) {
    ret[k + 1] <- sum(tmp * X[, k])
  }

  for (k in 1:pm) {
    ret[k + px + 1] <- sum(((predicts * exppredicts) %*% SetM[, k]) / tempsumm)
  }

  ret[p + 1] <- sum(rowSums(((predicts * predicts) - 1) * exppredicts) / tempsumm)

  return(-1 * ret)
}
