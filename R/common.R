
make.control <- function(max.retries = 10, method = "BFGS", ...) {
  if (max.retries <= 0) stop("maximum number of retries must be > 0")
  list(max.retries = max.retries, method = method, ...)
}

fit.algorithm <- function(par, f, g, Y, X, P, setM, 
                          control, silent = FALSE) {
  
  method <- control$method
  max.retries <- control$max.retries
  control$method <- NULL
  control$max.retries <- NULL
  retry <- TRUE
  retry_count <- 0
  
  while (retry) {
    if (method == "BFGS2") {
      erg <- ucminf(par = par, fn = f, gr = g, Y = Y, X = X, W = P, SetM = setM, 
                    control = control, hessian = TRUE)
    } else {
      erg <- optim(par = par, fn = f, gr = g, Y = Y, X = X, W = P, SetM = setM, 
                   method = method, control = control, hessian = TRUE)
    }
    # perform check
    erg$hessian <- (erg$hessian + t(erg$hessian)) / 2
    retry <- class(try(chol(erg$hessian))) == "try-error"
    if (retry) {
      if (!silent) cat("wrong hessian...resetting parameter ")
      #PD <- Matrix::nearPD(erg$hessian)
      #print(PD$normF)
      #stop()
      retry_count <- retry_count + 1
      if (retry_count > max.retries) stop("too many retries... hessian is not positive definite!")
      #print(erg$par)
      i <- which.max(abs(erg$par))
      if (!silent) cat(i)
      if (!silent) cat(" (")
      if (!silent) cat(erg$par[i])
      if (!silent) cat(")...")
      par <- erg$par
      par[i] <- 0
    } 
  }
  
  return(erg)
}


bootstrapping <- function(ret, Y, X, Pmodel, PX, f, g,
                                boot.fraction = 1, repetitions = 1000) {

  if (is.matrix(Y)) {
    k <- ncol(Y) - 1
    if (k == 0) k <- 1
  } else {
    k <- 1
  }
  
  simuliere.predict <- function(Pmodel, newdata, u, prep, s) {
    Pmodel$coefficients <- prep
    P <- predict(Pmodel, newdata = newdata[s,], type = "response")
    if (u == 2) {
      P0 <- 1 - P
      P <- cbind(P0, P)
      dimnames(P)[[2]] <- c("M0", "M1")
    }
    
    return(P)
  }
  
  u <- nrow(ret$setM)
  n <- nrow(X)
  m <- floor(boot.fraction * n)
  par <- ret$optim$par
  boot.beta <- matrix(0.0, nrow = repetitions, ncol = length(par))
    
  cat("bootstrapping")
  
  s <- as.big.matrix(matrix(as.integer(0), nrow = m, ncol = repetitions), type = "integer")#, backingfile="sbacking.bigmemory")
  prep <- as.big.matrix(matrix(0.0, nrow = repetitions,
                               ncol = length(Pmodel$coefficients)))#, backingfile="prepbacking.bigmemory")
  sdescription <- describe(s)
  prepdescription <- describe(prep)
  Xdescription <- NULL
  PXdescription <- NULL
  if (is.big.matrix(X)) Xdescription <- describe(X)
  if (is.big.matrix(PX)) PXdescription <- describe(PX)
  
  for (i in 1:repetitions) {
    # ziehen
    s[, i] <- sample.int(n, m, replace = TRUE)
    prep[i,] <- mvrnorm(n = 1, mu = coef(Pmodel), Sigma = vcov(Pmodel))
  }
  
  `%op%` <- ifelse(getDoParWorkers() > 1, `%dopar%`, `%do%`)
  
  tmp <- foreach(i = 1:repetitions, .inorder = FALSE) %op% {
    s <- attach.big.matrix(sdescription)
    prep <- attach.big.matrix(prepdescription)
    if (!is.null(Xdescription)) X <- attach.big.matrix(Xdescription)
    if (!is.null(PXdescription)) PX <- attach.big.matrix(PXdescription)
    si <- as.vector(sub.big.matrix(s, firstCol = i, lastCol = i)[,])
    
    # P simulieren
    P <- simuliere.predict(Pmodel, newdata = PX, u, prep[i,, drop = TRUE], si)
    
    # single step
    if (k > 1) {
      tmp2 <- fit.algorithm(par, f, g, Y[si,, drop = FALSE], X[si,, drop = FALSE], P, ret$setM,
                                 ret$control, silent = TRUE)$par
    } else {
      tmp2 <- fit.algorithm(par, f, g, Y[si], X[si,, drop = FALSE], P, ret$setM,
                                 ret$control, silent = TRUE)$par
    }
    cat(".")
    
    rm(si)
    rm(P)
    gc()
    
    tmp2
  }
  
  for (i in 1:repetitions) {
    boot.beta[i,] <- tmp[[i]]
  }
  rm(s)
  rm(prep)
  
  cat("done.\n")
  
  ret$boot.beta <- boot.beta
  ret$CoVar <- m / n * cov(boot.beta)
  ret$boot.fraction <- boot.fraction
  ret$SEtype <- "bootstrap"
  
  return(ret)
}
