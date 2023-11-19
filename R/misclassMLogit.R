#' Mlogit estimation under misclassified covariate
#'
#' \code{misclassMLogit} computes estimator for a GLM with a misclassified covariate
#' using additional side information on the misclassification process
#'
#' @param Y a matrix of 0s and 1s, indicating the target class. This is the dependent variable.
#' @param X a matrix containing the independent variables
#' @param setM matrix, rows containing potential patterns for a misclassed (latent) covariate M in
#'        any coding for a categorical independent variable, e.g. dummy coding.
#' @param P probabilities corresponding to each of the potential pattern conditional on the other
#'        covariates denoted in x.
#' @param na.action how to treat NAs
#' @param control options for the optimization procedure (see \code{\link{optim}},
#'        \code{\link{ucminf}} for options and details).
#' @param par (optional) starting parameter vector
#' @param baseoutcome reference outcome class
#' @param x logical, add covariates matrix to result?
#' @examples
#' ## simulate data
#' \donttest{
#' data <- simulate_mlogit_dataset()
#' }
#'
#' ## estimate model without misclassification error
#' \donttest{
#' library(mlogit)
#' data2 <- mlogit.data(data, varying = NULL, choice = "Y", shape = "wide")
#' summary(mlogit(Y ~ 1 | X + M2, data2, reflevel = "3"))
#' }
#'
#' ## estimate model with misclassification error
#' \donttest{
#' summary(mlogit(Y ~ 1 | X + M, data2, reflevel = "3"))
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
#' Yneu <- matrix(rep.int(0, nrow(data) * 3), ncol = 3)
#' for (i in 1:nrow(data)) Yneu[i, data$Y[i]] <- 1
#' est <- misclassMlogit(Y = Yneu,
#'                       X = as.matrix(data[, 2, drop = FALSE]),
#'                       setM = matrix(c(0, 1), nrow = 2),
#'                       P = P)
#' summary(est)
#' }
#'
#' ## and bootstrapping the results from dataset
#' \dontrun{
#' summary(boot.misclassMlogit(est,
#'                          Y = Yneu,
#'                          X = data.matrix(data[, 2, drop = FALSE]),
#'                          Pmodel = Pmodel,
#'                          PX = data,
#'                          repetitions = 100))
#' }
#'
#' @export
misclassMlogit <- function(Y, X, setM, P,
                        na.action = na.omit,
                        control = list(),
                        par = NULL,
                        baseoutcome = NULL,
                        x = FALSE) {
  if (!is.matrix(X)) stop("X is no matrix")
  if (!is.matrix(setM)) stop("setM is no matrix")
  if (!is.matrix(Y)) stop("Y must be a matrix")
  if (nrow(Y) != nrow(X)) stop("dimensions of Y and X do not match")
  if (is.null(baseoutcome)) baseoutcome <- which.max(colSums(Y))
  if (baseoutcome < ncol(Y)) Y <- cbind(Y[, -baseoutcome], Y[, baseoutcome])

  if (!is.matrix(P)) stop("P is no matrix")
  if (nrow(setM) != ncol(P)) stop("dimensions of P and setM do not match")
  if (nrow(P) != nrow(X)) stop("dimensions of P and X do not match")

  g <- NULL

  control <- do.call("make.control", control)

  # log-likelihood definieren
    f <- fmlogitValidation
    if (control$method == "BFGS" | control$method == "BFGS2" | control$method == "CG" | control$method == "nlm") {
      g <- gmlogitValidation
    }

    # starting parameters
    length.par <- (ncol(X) + ncol(setM) + 1) * (ncol(Y) - 1)
    if (is.null(par)) {
      cat("estimate starting parameters...")
      par <- rep.int(0.0, length.par)
    } else {
      if (length(par) != length.par) stop("starting parameter par is not valid.")
    }

  # optimise
  cat("optimize...")
  erg <- fit.algorithm(par, f, g, Y, X, P, setM, control)
  cat("ok.\n")

  # fertig
  ret <- list()
  ret$baseoutcome <- baseoutcome
  ret$family <- family
  ret$setM <- setM
  p <- ncol(X) + ncol(setM) + 1

  ret$coefficients <- erg$par[1:(p * (ncol(Y) - 1))]
  for (l in 1:(ncol(Y) - 1)) {
    names(ret$coefficients)[(l - 1) * p + 1:p] <- c(paste0("alt", l, ":(Intercept)"),
    unlist(lapply(dimnames(X)[[2]], function(n) paste0("alt", l, ":", n))),
    unlist(lapply(dimnames(P)[[2]][-1], function(p) paste0("alt", l, ":", p))))
  }

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

  class(ret) <- c("misclassMlogit", "mlogit")

  ret$prior.weights <- rep_len(1, nrow(X))
  ret$fitted.values <- predict(ret, X, P, type = "response")
  ret$control <- control

  ret$logLik <- erg$value

  ret$formula <- as.formula(paste0("y ~ ", paste(dimnames(X)[[2]], collapse = " + "), " + M"))

  return(ret)
}


#' Predict Method for \code{misclassMlogit} Fits
#'
#' Obtains predictions
#'
#' @usage ## S3 method for class 'misclassMlogit'
#'        \method{predict}{misclassMlogit}(object, X, P = NULL, type = c("link", "response"),
#'        na.action = na.pass, ...)
#' @param object a fitted object of class inheriting from 'misclassMlogit'.
#' @param X matrix of fixed covariates.
#' @param P a-posteriori probabilities for the true values of the misclassified variable. If provided,
#'        the conditional expectation on X,P is computed, otherwise a set of marginal predictions is
#'        provided, one for each alternative.
#' @param type the type of prediction required. The default is on the scale of the linear predictors;
#'        the alternative "response" is on the scale of the response variable.
#'        Thus for a default binomial model the default predictions are of log-odds (probabilities
#'        on logit scale) and type = "response" gives the predicted probabilities.
#'
#'        The value of this argument can be abbreviated.
#' @param na.action function determining what should be done with missing values in \code{newdata}.
#'        The default is to predict NA.
#' @param \dots additional arguments (not used at the moment)
#' @seealso \code{\link{misclassMlogit}}
#' @export
predict.misclassMlogit <- function(object, X, P = NULL, type = c("link", "response"),
                                   na.action = na.pass, ...) {

  type <- match.arg(type)

    k <- length(object$coefficients) / (ncol(X) + 1 + ncol(object$setM))
    ret <- list(nrow(object$setM))
    for (m in 1:nrow(object$setM)) {
      eta <- get.eta(cbind(X, object$setM[m,]), object$coefficients, k)
      if (type == "response") {
        summen <- rowSums(exp(eta)) + 1
        ret[[m]] <- cbind(exp(eta), rep.int(1, nrow(eta))) / summen
        dimnames(ret[[m]])[[2]] <- c(paste0("M", 1:k), "M0")
      } else {
        ret[[m]] <- eta
        dimnames(ret[[m]])[[2]] <- paste0("M", 1:k)
      }
    }
    if (type == "response" & !is.null(P)) {
      tmp <- ret
      ret <- tmp[[1]] * P[, 1]
      for (m in 2:nrow(object$setM)) {
        ret <- ret + tmp[[m]] * P[, m]
      }
    }

  return(ret)
}


get.eta <- function(X, beta, k) {
  d <- ncol(X)
  ret <- do.call(cbind, lapply(1:k, function(l) {
    X %*% beta[(1:d) + (d + 1)*(l - 1) + 1] + beta[1 + (d + 1)*(l - 1)]
  }))
  return(ret)
}


#' Compute Bootstrapped Standard Errors for \code{misclassMlogit} Fits
#'
#' Obtain bootstrapped standard errors.
#'
#' @param ret a fitted object of class inheriting from 'misclassMlogit'.
#' @param Y a matrix of 0s and 1s, indicating the target class. This is the dependent variable.
#' @param X a matrix containing the independent variables.
#' @param Pmodel a fitted model (e.g. of class 'GLM' or 'mlogit') to implicitly produce variations
#'        of the predicted true values probabilities. (Usually conditional on the observed
#'        misclassified values and additional covariates.)
#' @param PX covariates matrix suitable for predicting probabilities from \code{Pmodel},
#'        usually including the mismeasured covariate.
#' @param boot.fraction fraction of sample to be used for estimating the bootstrapped standard
#'        errors, for speedup.
#' @param repetitions number of bootstrap samples to be drown.
#' @seealso \code{\link{misclassMlogit}}
#' @export
boot.misclassMlogit <- function(ret, Y, X, Pmodel, PX,
                             boot.fraction = 1, repetitions = 1000) {
  if (!("misclassMlogit" %in% class(ret))) stop("object is no misclassMlogit result")
  baseoutcome <- ret$baseoutcome
  if (boot.fraction > 1) stop("boot.fraction must be <1")
  if (boot.fraction <= 0) stop("boot.fraction must be >0")

  if (is.matrix(Y)) {
    k <- ncol(Y) - 1
    if (k == 0) k <- 1
  } else {
    k <- 1
  }

  if (!is.matrix(X)) stop("X is no matrix")

  # letzte Kategorie als Referenz
  # daher: umsortieren!
  if (!is.matrix(Y)) stop("Y must be a matrix")
  if (nrow(Y) != nrow(X)) stop("dimensions of Y and X do not match")
  if (is.null(baseoutcome)) baseoutcome <- which.max(colSums(Y))
  if (baseoutcome < k) Y <- cbind(Y[, -baseoutcome], Y[, baseoutcome])

  g <- NULL

  # log-likelihood definieren
  f <- cfmlogitValidation
  if (ret$control$method == "BFGS" || ret$control$method == "BFGS2" ||
    ret$control$method == "CG" || ret$control$method == "nlm") {
    g <- cgmlogitValidation
  }

  return(bootstrapping(ret, Y, X, Pmodel, PX, f, g, boot.fraction, repetitions))
}

#' @export
summary.misclassMlogit <- function(object, ...) {
  b <- coef(object, fixed = TRUE)
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  print(CoefTable)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)")
  object$CoefTable <- CoefTable

  #object$lratio <- lratio(object)
  #object$mfR2 <- mfR2(object)

  #if (!is.null(object$rpar)) {
  #  rpar <- object$rpar
  #  object$summary.rpar <- t(sapply(rpar, summary))
  #}
  class(object) <- c("summary.misclassMlogit", "summary.mlogit", "mlogit")
  return(object)
}


#' @export
print.summary.misclassMlogit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  class(x) <- setdiff(class(x), "summary.misclassMlogit")
  #print(x = x, digits = digits, ...)

    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat("Frequencies of alternatives:")
    print(prop.table(x$freq), digits = digits)
    cat("\n")
    print(x$est.stat)
    cat("\nCoefficients :\n")
    printCoefmat(x$CoefTable, digits = digits)
    cat("\n")
    cat(paste("Log-Likelihood: ", signif(x$logLik, digits), "\n",
              sep = ""))
    #if (has.intercept(x$formula)) {
    #  cat("McFadden R^2: ", signif(x$mfR2, digits), "\n")
    #  cat("Likelihood ratio test : ", names(x$lratio$statistic),
    #      " = ", signif(x$lratio$statistic, digits), " (p.value = ",
    #      format.pval(x$lratio$p.value, digits = digits), ")\n",
    #      sep = "")
    #}
    #if (!is.null(x$summary.rpar)) {
    #  cat("\nrandom coefficients\n")
    #  print(x$summary.rpar)
    #}
    #invisible(x)

  if (x$SEtype == "robust") {
    cat("Results show robust standard errors.\n")
  } else if (x$SEtype == "bootstrap") {
    cat("Results based on ")
    cat(nrow(x$boot.beta))
    cat(" bootstrap samples.\n\n")
  }
}


#' @export
vcov.misclassMlogit <- function(object, ...) {
  if (is.null(object$CoVar)) {
    ok <- try(covmat <- chol2inv(chol(object$optim$hessian)), silent = TRUE)
    if (inherits(ok, "try-error")) {
      print(ok[1])
      return(NULL)
    }
  } else {
    covmat <- object$CoVar
  }
  return(covmat)
}


#' Compute Marginal Effects for 'misclassMlogit' Fits
#'
#' Obtain marginal effects.
#'
#' @param w a fitted object of class inheriting from 'misclassMlogit'.
#' @param x.mean logical, if true computes marginal effects at mean, otherwise average marginal
#'        effects.
#' @param rev.dum logical, if true, computes differential effects for switch from 0 to 1.
#' @param outcome for which the ME should be computed.
#' @param baseoutcome base outcome, e.g. reference class of the model.
#' @param digits number of digits to be presented in output.
#' @param \dots further arguments passed to or from other functions.
#' @seealso \code{\link{misclassMlogit}}
#' @export
mfx.misclassMlogit <- function(w, x.mean = TRUE, rev.dum = TRUE, outcome = 2, baseoutcome = 1,
    digits = 3, ...) {
  if (outcome == baseoutcome) return(NULL)
  if (!("mlogit" %in% class(w))) {
    stop("Please provide an object from 'mlogit()'.\n")
  }
  if (is.null(dim(w$x))) {
    stop("Please specify 'x = TRUE' in misclassMlogit().\n")}

  x.bar <- as.matrix(colMeans(w$x))
  b.est <- as.matrix(coef(w))
  K <- nrow(x.bar)

  alt <- length(coef(w)) / K
  xb <- rep.int(0, alt)
  for (i in 1:alt) {
    xb[i] <- t(x.bar) %*% b.est[(0:(K - 1)) * alt + i,]
  }
  xb <- exp(xb)

  if (!x.mean) stop('not implemented')

  if (outcome < baseoutcome) {
    f.xb <- xb[outcome] / ((1 + sum(xb)))^2
    me <- f.xb * coef(w)[(0:(K - 1)) * alt + outcome]
  } else if (outcome > baseoutcome) {
    f.xb <- xb[outcome - 1] / ((1 + sum(xb)))^2
    me <- f.xb * coef(w)[(0:(K - 1)) * alt + outcome - 1]
  }

  bx <- list()
  for (i in 1:alt) {
    bx[[i]] <- b.est[(0:(K - 1)) * alt + i,] %*% t(x.bar)
  }

  if (baseoutcome > outcome) {
    pg <- xb[outcome] / (1 + sum(xb))
    dr <- diag(1, K, K) + (1 - 2 * pg) * bx[[outcome]]
    va <- (pg - pg^2)^2 * dr %*% vcov(w)[(0:(K - 1)) * alt + outcome,
                                         (0:(K - 1)) * alt + outcome] %*% t(dr)
  } else if (baseoutcome < outcome) {
    pg <- xb[outcome - 1] / (1 + sum(xb))
    dr <- diag(1, K, K) + (1 - 2 * pg) * bx[[outcome - 1]]
    va <- (pg - pg^2)^2 * dr %*% vcov(w)[(0:(K - 1)) * alt + outcome - 1,
                                         (0:(K - 1)) * alt + outcome - 1] %*% t(dr)
  }

  se <- sqrt(diag(va))

  if (rev.dum) {
    for (i in 1:ncol(w$x)) {
      if (identical(sort(unique(w$x[,i])), c(0, 1)) || (i >= ncol(w$x) - dim(w$setM)[2]) ) {
        x.d1 <- x.bar; x.d1[i, 1] <- 1
        x.d0 <- x.bar; x.d0[i, 1] <- 0
        xb0 <- rep.int(0, alt)
        xb1 <- rep.int(0, alt)
        for (j in 1:alt) {
          xb0[j] <- t(x.d0) %*% b.est[(0:(K - 1)) * alt + j,]
          xb1[j] <- t(x.d1) %*% b.est[(0:(K - 1)) * alt + j,]
        }
        xb1 <- exp(xb1)
        xb0 <- exp(xb0)
        if (baseoutcome > outcome) {
          pr1 <- xb1[outcome] / (1 + sum(xb1))
          pr0 <- xb0[outcome] / (1 + sum(xb0))
          dr2 <- (pr1 - pr1^2) %*% t(x.d1) - (pr0 - pr0^2) %*% t(x.d0)
          va2 <- dr2 %*% vcov(w)[(0:(K - 1)) * alt + outcome,
                                 (0:(K - 1)) * alt + outcome] %*% t(dr2)
        } else if (baseoutcome < outcome) {
          pr1 <- xb1[outcome - 1] / (1 + sum(xb1))
          pr0 <- xb0[outcome - 1] / (1 + sum(xb0))
          dr2 <- (pr1 - pr1^2) %*% t(x.d1) - (pr0 - pr0^2) %*% t(x.d0)

          va2 <- dr2 %*% vcov(w)[(0:(K - 1)) * alt + outcome - 1,
                                 (0:(K - 1)) * alt + outcome - 1] %*% t(dr2)
        }
        me[i] <- pr1 - pr0

        se[i] <- sqrt(as.numeric(va2))
      }
    }
  }
  out <- data.frame(effect = me, error = se)
  out$t.value <- out$effect / out$error
  out$p.value <- 2 * (1 - pt(abs(out[, 3]), w$df.residual))
  out <- round(out, digits = digits)
  result <- list(f.xb = f.xb, w = w, out = out)
  class(result) <- c("mfx.misclassMlogit", "mfx")
  return(result)
}


#' Simulate a Data Set to Use With \code{misclassMlogit}
#'
#' simulates a data set with
#' - one continuous variable X drawn from a Gaussian distribution,
#' - a binary or trinary variable M with misclassification (M2)
#' - a dependent variable drawn from a multionomial distribution dependent on X and M.
#'
#' This can be used to demonstrate the abilities of misclassMlogit. For an example
#' see \code{misclassMlogit}.
#'
#' @param n number observations
#' @param const constants
#' @param alpha parameters for X
#' @param beta parameters for M(1)
#' @param beta2 parameters for M2, if NULL, M is a binary covariate, otherwise a three-valued
#'        categorical.
#' @seealso \code{\link{misclassMlogit}}
#' @export
simulate_mlogit_dataset <- function(n = 1000, const = c(0, 0), alpha = c(1, 2),
                                    beta = -2 * c(1, 2), beta2 = NULL) {
  set.seed(30)

  X <- rnorm(n, 0, 2)

  if (is.null(beta2)) {
    Mmisclass <- rbinom(n, 1, 0.5)

    M <- rep.int(0, n)
    for (i in 1:n) {
      temp <- exp(X[i] + Mmisclass[i])
      M[i] <- rbinom(1, 1, temp / (1 + temp))
    }

    Y1 <- (alpha[1] * X) + beta[1] * M + const[1]
    Y2 <- (alpha[2] * X) + beta[2] * M + const[2]
  } else {
    M <- matrix(rep.int(0, 2 * n), ncol = 2)
    for (i in 1:n) {
      if (sum(Mmisclass[i,]) == 2) Mmisclass[i,] <- c(0, 0)
      temp <- exp(X[i])
      misclass <- rbinom(1, 1, temp / (1 + temp))
      if (rbinom(1, 1, 0.5)) {
        M[i, 1] <- ifelse(misclass, 1 - Mmisclass[i, 1], Mmisclass[i, 1])
      } else {
        M[i, 2] <- ifelse(misclass, 1 - Mmisclass[i, 2], Mmisclass[i, 2])
      }
      if (sum(M[i,]) == 2) M[i,] <- c(0, 0)
    }

    Y1 <- alpha[1] * X + beta[1] * M[, 1] + beta2[1] * M[, 2] + const[1]
    Y2 <- alpha[2] * X + beta[2] * M[, 1] + beta2[2] * M[, 2] + const[2]
  }
  Y1 <- exp(Y1)
  Y2 <- exp(Y2)
  tmp <- rep.int(1, n) + Y1 + Y2
  Y <- matrix(rep.int(0, 3 * n), ncol = 3)
  for (i in 1:n) {
    prob <- c(Y1[i] / tmp[i], Y2[i] / tmp[i], 1 / tmp[i])
    Y[i,] <- t(rmultinom(1, 1, prob))
  }

  return(data.frame(Y = as.factor(Y[, 1] + 2 * Y[, 2] + 3 * (1 - Y[, 1] - Y[, 2])),
                    X = X,
                    M = Mmisclass,
                    M2 = M))
}
