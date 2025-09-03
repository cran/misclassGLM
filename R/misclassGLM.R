#' GLM estimation under misclassified covariate
#'
#' \code{misclassGLM} computes estimator for a GLM with a misclassified covariate
#' using additional side information on the misclassification process
#'
#' @param Y a vector of integers or numerics. This is the dependent variable.
#' @param X a matrix containing the independent variables.
#' @param P probabilities corresponding to each of the potential pattern conditional on the other
#'        covariates denoted in x.
#' @param na.action how to treat NAs
#' @param family a description of the error distribution and link function to be used in the model.
#'        This can be a character string naming a family function, a family function or the result
#'        of a call to a family function. (See \code{\link{family}} for details of family functions.)
#' @param control options for the optimization procedure (see \code{\link{optim}},
#'        \code{\link[ucminf]{ucminf}} for options and details).
#' @param par (optional) starting parameter vector
#' @param setM (optional) matrix, rows containing potential patterns for a misclassified (latent) covariate M
#'        in any coding for a categorical independent variable, e.g. dummy coding (default: Identity).
#' @param x logical, add covariates matrix to result?
#' @param robust logical, if true the computed asymptotic standard errors are replaced by their
#'        robust counterparts.
#' @examples
#' ## simulate data
#' \donttest{
#' data <- simulate_GLM_dataset()
#' }
#'
#' ## estimate model without misclassification error
#' \donttest{
#' summary(lm(Y ~ X + M2, data))
#' }
#'
#' ## estimate model with misclassification error
#' \donttest{
#' summary(lm(Y ~ X + M, data))
#' }
#'
#' ## estimate misclassification probabilities
#' \donttest{
#' Pmodel <- glm(M2 ~ M + X, data = data, family = binomial("logit"))
#' summary(Pmodel)
#' }
#'
#' ## construct a-posteriori probabilities from Pmodel
#' \donttest{
#' P <- predict(Pmodel, newdata = data, type = "response")
#' P <- cbind(1 - P, P)
#' dimnames(P)[[2]] <- c("M0", "M1") ## speaking names
#' }
#'
#' ## estimate misclassGLM
#' \donttest{
#' est <- misclassGLM(Y = data$Y,
#'                    X = as.matrix(data[, 2, drop = FALSE]),
#'                    setM = matrix(c(0, 1), nrow = 2),
#'                    P = P)
#' summary(est)
#' }
#'
#' ## and bootstrapping the results from dataset
#' \dontrun{
#'   summary(boot.misclassGLM(est,
#'                            Y = data$Y,
#'                            X = data.matrix(data[, 2, drop = FALSE]),
#'                            Pmodel = Pmodel,
#'                            PX = data,
#'                            repetitions = 100))
#' }
#'
#' @export
misclassGLM <- function(Y, X, setM, P,
                        na.action = na.omit,
                        family = gaussian(link = "identity"),
                        control = list(),
                        par = NULL,
                        x = FALSE,
                        robust = FALSE) {
  if (!is.matrix(X)) stop("X is no matrix")
  if (!is.matrix(setM)) stop("setM is no matrix")
  if (!is.vector(Y)) stop("Y must be a vector")
  if (length(Y) != nrow(X)) stop("dimensions of Y and X do not match")
  if (!is.matrix(P)) stop("P is no matrix")
  if (nrow(setM) != ncol(P)) stop("dimensions of P and setM do not match")
  if (nrow(P) != nrow(X)) stop("dimensions of P and X do not match")

  g <- NULL

  control <- do.call("make.control", control)

  # log-likelihood definieren

    # Validation
    if (family$family == "gaussian") {
      f <- fgaussValidation
      if (control$method == "BFGS" || control$method == "BFGS2" || control$method == "CG" ||
          control$method == "nlm") {
        g <- ggaussValidation
      }
    } else if (family$family == "binomial") {
      if (family$link == "logit") {
        f <- flogitValidation
        if (control$method == "BFGS" || control$method == "BFGS2" || control$method == "CG" ||
            control$method == "nlm") {
          g <- glogitValidation
        }
      } else if (family$link == "probit") {
        f <- fprobitValidation
        if (control$method == "BFGS" || control$method == "BFGS2" || control$method == "CG" ||
            control$method == "nlm") {
          g <- gprobitValidation
        }
      } else {
        stop("unknown link")
      }
    } else {
      stop("unknown family")
    }

    # starting parameters
    length.par <- ncol(X) + ncol(setM) + 1
    if (is.null(par)) {
      cat("estimate starting parameters...")
      par <- rep.int(0, length.par)
      if (family$family == "gaussian") {
        par <- c(par, 0.5)
        par[1:(ncol(X) + 1)] <- lm.fit(x = cbind(rep.int(1, nrow(X)), X), y = Y)$coefficients
      } else {
        par[1:(ncol(X) + 1)] <- glm.fit(x = cbind(rep.int(1, nrow(X)), X), y = Y,
                                        family = family)$coefficients
      }
    } else {
      if (family$family == "gaussian") length.par <- length.par + 1
      if (length(par) != length.par) stop("starting parameter par is not valid.")
    }

  # optimise
  cat("optimize...")
  erg <- fit.algorithm(par, f, g, Y, X, P, setM, control)
  cat("ok.\n")

  # fertig
  ret <- list()
  ret$family <- family
  ret$setM <- setM
  p <- ncol(X) + ncol(setM) + 1
  ret$coefficients <- erg$par[1:p]
  names(ret$coefficients) <- c("(Intercept)", dimnames(X)[[2]], dimnames(P)[[2]][-1])

  ret$y <- Y
  ret$na.action <- na.action
  ret$df.residual <- nrow(X) - length(par)
  ret$optim <- erg
  if (!is.null(g)) ret$optim$gradient <- g(erg$par, Y, X, P, setM)
  if (control$method == "BFGS2") {
    ret$iter <- erg$info[4]
    ret$optim$call <- "ucminf"
  } else {
    ret$optim$call <- "optim"
    #ret$optim$hessian <- hessian(f, erg$par, Y=Y, X=X, W=P, SetM=setM)
    ret$iter <- erg$counts
  }
  if (x) ret$x <- cbind(rep.int(1, nrow(X)), X, P[, -1, drop = FALSE])

  ret$SEtype <- "standard"
  ret$setM <- setM
  ret$control <- control

  class(ret) <- c("misclassGLM", "glm", "lm")

  ret$fitted.values <- predict(ret, X, P, type = "response")
  ret$prior.weights <- rep_len(1, nrow(X))

  if (robust) ret <- estimateRobustSE(ret, Y, X, P)

  return(ret)
}


#' Predict Method for \code{misclassGLM} Fits
#'
#' Obtains predictions
#'
#' @usage ## S3 method for class 'misclassGLM'
#'        \method{predict}{misclassGLM}(object, X, P = NULL, type = c("link", "response"),
#'                                      na.action = na.pass, ...)
#' @param object a fitted object of class inheriting from 'misclassGLM'.
#' @param X matrix of fixed covariates
#' @param P a-posteriori probabilities for the true values of the misclassified variable.
#'        If provided, the conditional expectation on X,P is computed, otherwise a set of marginal
#'        predictions is provided, one for each alternative.
#' @param type the type of prediction required. The default is on the scale of the linear predictors;
#'        the alternative "response" is on the scale of the response variable.
#'        Thus for a default binomial model the default predictions are of log-odds
#'        (probabilities on logit scale) and type = "response" gives the predicted probabilities.
#'
#'        The value of this argument can be abbreviated.
#' @param na.action function determining what should be done with missing values in \code{newdata}.
#'        The default is to predict NA.
#' @param \dots additional arguments (not used at the moment)
#' @seealso \code{\link{misclassGLM}}
#' @export
predict.misclassGLM <- function(object, X, P = NULL, type = c("link", "response"),
                                na.action = na.pass, ...) {

  type <- match.arg(type)

    ret <- matrix(rep_len(0.0, nrow(object$setM) * nrow(X)), nrow = nrow(X))
    for (m in nrow(object$setM)) {
      ret[, m] <- cbind(rep_len(1, nrow(X)), X,
                        matrix(rep.int(object$setM[m,], nrow(X)), nrow = nrow(X), byrow = TRUE)) %*%
        object$coefficients
      if (type == "response") family(object)$linkinv(ret[, m])
    }
    if (type == "response" & !is.null(P)) ret <- rowSums(ret * P)

  return(ret)
}


#' Compute Bootstrapped Standard Errors for \code{misclassGLM} Fits
#'
#' Obtain bootstrapped standard errors.
#'
#' @param ret a fitted object of class inheriting from 'misclassGLM'.
#' @param Y a vector of integers or numerics. This is the dependent variable.
#' @param X a matrix containing the independent variables.
#' @param Pmodel a fitted model (e.g. of class 'GLM' or 'mlogit') to implicitly produce variations
#'        of the predicted true values probabilities. (Usually conditional on the observed
#'        misclassified values and additional covariates.)
#' @param PX covariates matrix suitable for predicting probabilities from \code{Pmodel},
#'        usually including the mismeasured covariate.
#' @param boot.fraction fraction of sample to be used for estimating the bootstrapped standard
#'        errors, for speedup.
#' @param repetitions number of bootstrap samples to be drown.
#' @seealso \code{\link{misclassGLM}}
#' @export
boot.misclassGLM <- function(ret, Y, X, Pmodel, PX,
                             boot.fraction = 1, repetitions = 1000) {
  if (!("misclassGLM" %in% class(ret))) stop("object is no misclassGLM result")
  family <- ret$family
  if (boot.fraction > 1) stop("boot.fraction must be <1")
  if (boot.fraction <= 0) stop("boot.fraction must be >0")

  if (!is.matrix(X)) stop("X is no matrix")
  if (!is.vector(Y)) stop("Y must be a vector")
  if (length(Y) != nrow(X)) stop("dimensions of Y and X do not match")

  g <- NULL

  # log-likelihood definieren
  if (family$family == "gaussian") {
    f <- cfgaussValidation
    if (ret$control$method == "BFGS" || ret$control$method == "BFGS2" ||
        ret$control$method == "CG" || ret$control$method == "nlm") {
      g <- cggaussValidation
    }
  } else if (family$family == "binomial") {
    if (family$link == "logit") {
      f <- cflogitValidation
      if (ret$control$method == "BFGS" || ret$control$method == "BFGS2" ||
          ret$control$method == "CG" || ret$control$method == "nlm") {
        g <- cglogitValidation
      }
    } else if (family$link == "probit") {
      f <- fprobitValidation
      if (ret$control$method == "BFGS" || ret$control$method == "BFGS2" ||
          ret$control$method == "CG" || ret$control$method == "nlm") {
        g <- gprobitValidation
      }
    } else {
      stop("unknown link")
    }
  } else {
    stop("unknown family")
  }

  return(bootstrapping(ret, Y, X, Pmodel, PX, f, g, boot.fraction, repetitions))
}


estimateRobustSE <- function(object, Y, X, P) {
  if (!inherits(object, "misclassGLM")) stop("object is not of class misclassGLM.")

  if (object$SEtype != "standard") stop("robust errors require usual computation of SE.")

  ok <- try(covmat <- chol2inv(chol(object$optim$hessian)), silent = TRUE)
  if (inherits(ok, "try-error")) {
    print(ok[1])
    return(NULL)
  }

  if (object$family$family == "gaussian") {
    g <- cggaussValidation
  } else if (object$family$family == "binomial") {
    if (object$family$link == "logit") {
      g <- cglogitValidation
    } else if (object$family$link == "probit") {
      g <- gprobitValidation
    } else {
      stop("unknown link")
    }
  } else {
    stop("unknown family")
  }
  Lstrich <- matrix(rep.int(0, length(object$optim$par)^2), ncol = length(object$optim$par))
  for (i in 1:nrow(X)) {
    if (inherits(Y, "matrix")) {
      Lstrich <- Lstrich + tcrossprod(g(object$optim$par, Y[i,, drop = FALSE], X[i,, drop = FALSE],
                                        P[i,, drop = FALSE], object$setM))
    } else {
      Lstrich <- Lstrich + tcrossprod(g(object$optim$par, Y[i], X[i,, drop = FALSE],
                                        P[i,, drop = FALSE], object$setM))
    }
  }
  object$CoVar <- covmat %*% Lstrich %*% covmat

  return(object)
}

#' @export
summary.misclassGLM <- function(object, ...) {
  coef.p <- object$coefficients
  p <- length(coef.p)

  if (is.null(object$CoVar)) {
    ok <- try(covmat <- chol2inv(chol(object$optim$hessian)), silent = TRUE)
    if (inherits(ok, "try-error")) {
      print(ok[1])
      return(NULL)
    }
  } else {
    covmat <- object$CoVar
  }
  var.cf <- diag(covmat[1:p, 1:p])
  s.err <- sqrt(var.cf)
  tvalue <- coef.p / s.err
  if (object$SEtype == "standard") {
    df.r <- object$df.residual
    pvalue <- 2 * pt(-abs(tvalue), df.r)
  } else if (object$SEtype == "robust") {
    df.r <- object$df.residual
    pvalue <- 2 * pt(-abs(tvalue), df.r)
  } else if (object$SEtype == "bootstrap") {
    pvalue <- rep_len(0.0, p)
    for (j in 1:p) {
      tmp <- object$boot.beta[, j]
      pvalue[j] <- 2 * min(length(which(tmp <= 0)),
                           length(which(tmp >= 0))) / (nrow(object$boot.beta))
    }
  }
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

  keep <- match(c(#"call",
          "family", #"deviance", "aic",
          #"contrasts",
          "df.residual", #"null.deviance", "df.null",
          "iter", "na.action", "SEtype"), names(object), 0L)
  ret <- c(object[keep],
      list(deviance.resid = residuals(object, type = "deviance"),
          coefficients = coef.table, #aliased = aliased,
          #dispersion = dispersion, df = c(object$rank, df.r, df.f),
          #cov.unscaled = covmat.unscaled,
          cov.scaled = covmat
          ))

  class(ret) <- c("summary.misclassGLM", "summary.glm")
  return(ret)
}

#' @export
print.summary.misclassGLM <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  class(x) <- setdiff(class(x), "summary.misclassGLM")
  print(x = x, digits = digits, ...)

  if (x$SEtype == "robust") {
    cat("Results show robust standard errors.\n")
  } else if (x$SEtype == "bootstrap") {
    cat("Results based on ")
    cat(nrow(x$boot.beta))
    cat(" bootstrap samples.\n\n")
  }
}

#' @export
vcov.misclassGLM <- function(object, ...) {
  summary.misclassGLM(object)$cov.scaled
}

#' @export
df.residual.misclassGLM <- function(object, ...) {
  object$df.residuals
}

#' Compute Marginal Effects for \code{misclassGLM} Fits
#'
#' Obtain marginal Effects.
#'
#' @param w a fitted object of class inheriting from 'misclassGLM'.
#' @param x.mean logical, if true computes marginal effects at mean, otherwise average marginal
#'               effects.
#' @param rev.dum logical, if true, computes differential effects for switch from 0 to 1.
#' @param digits number of digits to be presented in output.
#' @param \dots further arguments passed to or from other functions.
#' @seealso \code{\link{misclassGLM}}
#' @export
mfx.misclassGLM <- function(w, x.mean = TRUE, rev.dum = TRUE, digits = 3, ...) {
    if (!("glm" %in% class(w))) {
      stop("Please provide an object from 'glm()' or 'mlogit()'.\n")
    }
    if ("glm" %in% class(w)) link <- w$family$link
    if (link != "probit" & link != "logit") {
      stop("Need a binary probit, a logit or mlogit model'.")}
    if (is.null(dim(w$x))) {
      stop("Please specify 'x = TRUE' in glm().\n")}

    x.bar <- as.matrix(colMeans(w$x))
    b.est <- as.matrix(coef(w))
    K <- nrow(x.bar)

    xb <- t(x.bar) %*% b.est

    if (!x.mean) {
      xb2 <- as.matrix(w$x) %*% b.est
      if (link == "probit") f.xb <- mean(dnorm(xb2))
      if (link == "logit")  f.xb <- mean(dlogis(xb2))
    }
    if (link == "probit") f.xb <- dnorm(xb)
    if (link == "logit" ) f.xb <- dlogis(xb)
    me <- f.xb * coef(w)

    bx <- b.est %*% t(x.bar)

    if (link == "probit") {
      dr <- diag(1, K, K) - as.numeric(xb) * bx
      va <- as.numeric(f.xb)^2 * dr %*% vcov(w) %*% t(dr)
    } else if (link == "logit") {
      pg <- as.numeric(plogis(xb))
      dr <- diag(1, K, K) + (1 - 2 * pg) * bx
      va <- (pg * (1 - pg))^2 * dr %*% vcov(w) %*% t(dr)
    }
    se <- sqrt(diag(va))

    if (rev.dum) {
      for (i in 1:ncol(w$x)) {
        if (identical(sort(unique(w$x[,i])), c(0, 1)) ||
            ("misclassGLM" %in% class(w) && i >= ncol(w$x) - dim(w$setM)[2]) ) {
          x.d1 <- x.bar; x.d1[i, 1] <- 1
          x.d0 <- x.bar; x.d0[i, 1] <- 0
          if (link == "probit") {
            me[i] <- pnorm(t(x.d1) %*% b.est) -
              pnorm(t(x.d0) %*% b.est)
            dr2 <- dnorm(t(x.d1) %*% b.est) %*%  t(x.d1) -
              dnorm(t(x.d0) %*% b.est) %*%  t(x.d0)
            va2 <- dr2 %*% vcov(w) %*% t(dr2)
          } else if (link == "logit") {
            me[i] <- plogis(t(x.d1) %*% b.est) -
              plogis(t(x.d0) %*% b.est)
            dr2 <- dlogis(t(x.d1) %*% b.est) %*%  t(x.d1) -
              dlogis(t(x.d0) %*% b.est) %*%  t(x.d0)
            va2 <- dr2 %*% vcov(w) %*% t(dr2)
          }
          se[i] <- sqrt(as.numeric(va2))
        }
      }
    }
    out <- data.frame(effect = me, error = se)
    out$t.value <- out$effect / out$error
    out$p.value <- 2 * (1 - pt(abs(out[, 3]), w$df.residual))
    out <- round(out, digits = digits)
    result <- list(link = link, f.xb = f.xb, w = w, out = out)
    class(result) <- "mfx"
    return(result)

}


#' Simulate a Data Set to Use With \code{misclassGLM}
#'
#' simulates a data set with
#' - one continuous variable X drawn from a Gaussian distribution,
#' - a binary or trinary variable M with misclassification (M2)
#' - a dependent variable either with added Gaussian noise or drawn from a logit distribution
#'
#' This can be used to demonstrate the abilities of \code{\link{misclassGLM}}. For an example
#' see \code{\link{misclassGLM}}.
#'
#' @param n number observations
#' @param const constant
#' @param alpha parameter for X
#' @param beta parameter for M(1)
#' @param beta2 parameter for M2, if NULL, M is a binary covariate, otherwise a three-valued
#'              categorical
#' @param logit logical, if true logit regression, otherwise Gaussian regression
#' @seealso \code{\link{misclassGLM}}
#' @export
simulate_GLM_dataset <- function(n = 50000, const = 0, alpha = 1, beta = -2, beta2 = NULL,
                                 logit = FALSE) {
  set.seed(30)

  X <- rnorm(n, 0, 2)

  if (is.null(beta2)) {
    Mmisclass <- rbinom(n, 1, 0.5)

    M <- rep.int(0, n)
    for (i in 1:n) {
      temp <- exp(X[i] + Mmisclass[i])
      M[i] <- rbinom(1, 1, temp / (1 + temp))
    }

    Y <- (alpha * X) + beta * M + const
  } else {
    M <- matrix(rep.int(0, 2 * n), ncol = 2)
    for (i in 1:n) {
      tmp1 <- exp(Mmisclass[i, 1] - Mmisclass[i, 2] - (1 - Mmisclass[i, 1]) * (1 - Mmisclass[i, 2]))
      tmp2 <- exp(Mmisclass[i, 2] - Mmisclass[i, 1] - (1 - Mmisclass[i, 1]) * (1 - Mmisclass[i, 2]))
      s <- tmp1 + tmp2 + 1
      M[i,] <- rmultinom(1, 1, c(tmp1, tmp2, 1) / s)[-3]
    }

    Y <- alpha * X + beta * M[, 1] + beta2 * M[, 2] + const
  }

  if (logit) {
    Y <- exp(Y)
    for (i in 1:n) {
      Y[i] <- rbinom(1, 1, Y[i] / (1 + Y[i]))
    }
  } else {
    E <- rnorm(n, 0, 0.5)
    Y <- Y + E
  }

  return(data.frame(Y = Y, X = X, M = Mmisclass, M2 = M))
}

