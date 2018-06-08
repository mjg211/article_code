###############################################################################
# Name:        gs.sw                                                          #
# Description: Determines group sequential SW-CRT designs analysed with the   #
#              Hussey and Hughes model.                                       #
###############################################################################
# Inputs                                                                      #
# ######                                                                      #
#                                                                             #
#        X - Matrix of the binary X_{ij} that indicate whether a cluster      #
#            receives the intervention in a particular time period.           #
#  sigma.e - Residual standard deviaiton.                                     #
#  sigma.c - Between cluster standard deviation.                              #
#    alpha - Desired type-I error rate.                                       #
#     beta - Desired type-II error rate.                                      #
#    delta - Clinically relevant difference to power for.                     #
#    set.T - Vector of time periods after which interim analyses are to be    #
#            performed.                                                       #
# stopping - Character string indicating if early stopping is to be included  #
#            for futility, efficacy, or futility and efficacy.                #
#  gamma.f - Power in the spending function f.                                #
#  gamma.e - Power in the spending function e.                                #
#  summary - Logical; indicates whether a summary of the programmes progress  #
#            should be printed.                                               #
###############################################################################
# Outputs                                                                     #
# #######                                                                     #
#                                                                             #
# A list is returned containing the inputs and the following factors:         #
#                                                                             #
#       e - Vector of efficacy stopping boundaries.                           #
#       f - Vector of futility stopping boundaries.                           #
#       I - Vector of information levels (across the interim analyses).       #
#  Lambda - Covariance matrix of the test statistics (across the interim      #
#           analyses).                                                        #
#       n - Required per cluster per period sample size.                      #
#    Nvec - Vector of sample sizes required (across the interim analyses).    #
#    n.sw - Per cluster per period sample size required by a corresponding    #
#           fixed sample SW-CRT design.                                       #
# perf.H0 - Performance of the GS SW-CRT design under H0 (type-I error rate   #
#           and ESS).                                                         #
# perf.H1 - Performance of the GS SW-CRT design under H1 (power and ESS).     #
###############################################################################
library(mvtnorm)

gs.sw <- function(X, sigma.e = sqrt(0.51), sigma.c = sqrt(0.02), alpha = 0.05,
                  beta = 0.1, delta = 0.2, set.T = 2:5, stopping = "both",
                  gamma.f = 0.5, gamma.e = 0.5, summary = TRUE){

  ##### ERROR CHECKING ########################################################

  if (!is.matrix(X) | ncol(X) == 1 | nrow(X) == 1 | !all(X %in% c(0, 1))){
    stop("X must be an indicator matrix with at least 2 rows and 2 columns.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if (sigma.c <= 0){
    stop("Between cluster standard deviation sigma.e must be strictly positive.")
  }
  if ((alpha <= 0) | (alpha >= 1)){
    stop("Type-I error rate alpha must be strictly between 0 and 1.")
  }
  if ((beta <= 0) | (beta >= 1)){
    stop("Type-II error rate beta must be strictly between 0 and 1.")
  }
  if (delta <= 0){
    stop("Clinically relevant difference delta to power for must be strictly positive.")
  }
  if (!is.vector(set.T) | !all(set.T %in% 1:ncol(X)) | sum(X[, 1:set.T[1]]) == 0){
    stop("set.T must be a vector with elements in (min_j sum(X[, 1:j]) > 0):ncol(X).")
  }
  if (!(stopping %in% c("efficacy", "futility", "both"))){
    stop("stopping must be set to one of \"efficacy\", \"futility\", or \"both\".")
  }
  if (gamma.f <= 0){
    stop("gamma.f must be greater than 0.")
  }
  if (gamma.e <= 0){
    stop("gamma.e must be greater than 0.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  informationCRCT <- function(n, X, set.T, sigma.e, sigma.c){
    if (length(set.T) == 1 && set.T == 1){
      C      <- length(X)
      sigma2 <- sigma.e^2/n
      U      <- sum(X)
      V      <- sum(X^2)
      W      <- sum(X)^2
      I      <- ((C*U - W)*sigma2 + (U^2 + C*U - W - C*V)*sigma.c^2)/
        (C*sigma2*(sigma2 + sigma.c^2))
    } else {
      C      <- nrow(X)
      sigma2 <- sigma.e^2/n
      I      <- numeric(length(set.T))
      for (i in 1:length(set.T)){
        U    <- sum(X[, 1:set.T[i]])
        if (set.T[i] > 1){
          V  <- sum(rowSums(X[, 1:set.T[i]])^2)
          W  <- sum(colSums(X[, 1:set.T[i]])^2)
        } else {
          V  <- sum(X[, set.T[i]]^2)
          W  <- sum(X[, set.T[i]])^2
        }
        I[i] <- ((C*U - W)*sigma2 + (U^2 + C*set.T[i]*U - set.T[i]*W - C*V)*sigma.c^2)/
          (C*sigma2*(sigma2 + set.T[i]*sigma.c^2))
      }
    }
    return(I)
  }

  boundaryF <- function(crit, pi2l, prefbounds, preebounds, currSigma, currI,
                        delta){
    integral <- pmvnorm(lower = c(prefbounds, -Inf),
                        upper = c(preebounds, crit),
                        mean = rep(delta, length(currI))*sqrt(currI),
                        sigma = currSigma)[1]
    return((pi2l - integral)^2)
  }

  boundaryE <- function(crit, pi1l, prefbounds, preebounds, currSigma){
    integral <- pmvnorm(lower = c(prefbounds, crit),
                        upper = c(preebounds, Inf),
                        mean = numeric(length(prefbounds) + 1),
                        sigma = currSigma)[1]
    return((pi1l - integral)^2)
  }

  perfF <- function(theta, L, f, I, nvec, Sigma){
    PF <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PF[l] <- pmvnorm(lower = -Inf, upper = f[l],
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PF[l] <- pmvnorm(lower = c(f[1:(l - 1)], -Inf),
                         upper = c(rep(Inf, l - 1), f[l]),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    PR <- 1 - sum(PF)
    EN <- sum(nvec[1:(L - 1)]*PF[1:(L - 1)]) +
      nvec[L]*(1 - sum(PF[1:(L - 1)]))
    return(c(PR, EN))
  }

  perfE <- function(theta, L, e, I, nvec, Sigma){
    PE <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PE[l] <- pmvnorm(lower = e[l], upper = Inf,
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PE[l] <- pmvnorm(lower = c(rep(-Inf, l - 1), e[l]),
                         upper = c(e[1:(l - 1)], Inf),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    PR <- sum(PE)
    EN <- sum(nvec[1:(L - 1)]*PE[1:(L - 1)]) +
      nvec[L]*(1 - sum(PE[1:(L - 1)]))
    return(c(PR, EN))
  }

  perfEF <- function(theta, L, e, f, I, nvec, Sigma){
    PE <- numeric(L)
    PF <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PE[l] <- pmvnorm(lower = e[l], upper = Inf,
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
        PF[l] <- pmvnorm(lower = -Inf, upper = f[l],
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PE[l] <- pmvnorm(lower = c(f[1:(l - 1)], e[l]),
                         upper = c(e[1:(l - 1)], Inf),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
        PF[l] <- pmvnorm(lower = c(f[1:(l - 1)], -Inf),
                         upper = c(e[1:(l - 1)], f[l]),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    PR <- sum(PE)
    EN <- sum(nvec*(PE + PF))
    return(c(PR, EN))
  }

  covariance <- function(L, I){
    Sigma <- diag(1, L, L)
    for (l1 in 2:L){
      for (l2 in 1:(l1 - 1)){
        Sigma[l1, l2] <- sqrt(I[l2]/I[l1])
        Sigma[l2, l1] <- Sigma[l1, l2]
      }
    }
    return(Sigma)
  }

  family <- function(L, I, spend, gamma){
    pi      <- numeric(L)
    pi[1]   <- min(spend*(I[1]/I[L])^gamma, spend)
    for (l in 2:L){
      pi[l] <- min(spend*(I[l]/I[L])^gamma, spend) -
        min(spend*(I[l - 1]/I[L])^gamma, spend)
    }
    return(pi)
  }

  nSeqCRCTF <- function(n, X, delta, alpha, beta, sigma.e, sigma.c, set.T, gamma.f){
    C     <- nrow(X)
    L     <- length(set.T)
    I     <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
    Sigma <- covariance(L, I)
    pi2   <- family(L, I, beta, gamma.f)
    f     <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
      } else if ((l > 1) & (l < L)){
        f[l] <- suppressWarnings(optim(par = qnorm(pi2[l],
                                                   mean = delta*sqrt(I[l])),
                                       fn = boundaryF, pi2l = pi2[l],
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = rep(Inf, l - 1),
                                       currSigma = Sigma[1:l, 1:l],
                                       currI = I[1:l], delta = delta)$par)
      } else {
        f[l] <- suppressWarnings(optim(par = qnorm(1 - alpha/L),
                                       fn = boundaryE, pi1l = alpha,
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = rep(Inf, L - 1),
                                       currSigma = Sigma)$par)
      }
    }
    nvec <- n*C*set.T
    PR   <- perfF(delta, L, f, I, nvec, Sigma)[1]
    return((beta - (1 - PR))^2)
  }

  nSeqCRCTE <- function(n, X, delta, alpha, beta, sigma.e, sigma.c, set.T, gamma.e){
    C     <- nrow(X)
    L     <- length(set.T)
    I     <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
    Sigma <- covariance(L, I)
    pi1   <- family(L, I, alpha, gamma.e)
    e     <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        e[l] <- qnorm(1 - pi1[l])
      } else {
        e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                       fn = boundaryE, pi1l = pi1[l],
                                       prefbounds = rep(-Inf, l - 1),
                                       preebounds = e[1:(l - 1)],
                                       currSigma = Sigma[1:l, 1:l])$par)
      }
    }
    nvec <- n*C*set.T
    PR   <- perfE(delta, L, e, I, nvec, Sigma)[1]
    return((beta - (1 - PR))^2)
  }

  nSeqCRCTEF <- function(n, X, delta, alpha, beta, sigma.e, sigma.c, set.T,
                         gamma.f, gamma.e){
    C     <- nrow(X)
    L     <- length(set.T)
    I     <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
    Sigma <- covariance(L, I)
    pi1   <- family(L, I, alpha, gamma.e)
    pi2   <- family(L, I, beta, gamma.f)
    e     <- numeric(L)
    f     <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
        e[l] <- qnorm(1 - pi1[l])
      } else if ((l > 1) & (l < L)){
        f[l] <- suppressWarnings(optim(par = qnorm(pi2[l],
                                                   mean = delta*sqrt(I[l])),
                                       fn = boundaryF, pi2l = pi2[l],
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = e[1:(l - 1)],
                                       currSigma = Sigma[1:l, 1:l],
                                       currI = I[1:l], delta = delta)$par)
        e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                       fn = boundaryE, pi1l = pi1[l],
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = e[1:(l - 1)],
                                       currSigma = Sigma[1:l, 1:l])$par)
      } else {
        e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                       fn = boundaryE, pi1l = pi1[l],
                                       prefbounds = f[1:(l - 1)],
                                       preebounds = e[1:(l - 1)],
                                       currSigma = Sigma[1:l, 1:l])$par)
        f[l] <- e[l]
      }
    }
    nvec <- n*C*set.T
    PR   <- perfEF(delta, L, e, f, I, nvec, Sigma)[1]
    return((beta - (1 - PR))^2)
  }

  nCRCT <- function(n, X, delta, beta, sigma.e, sigma.c, e){
    I     <- informationCRCT(n, X, ncol(X), sigma.e, sigma.c)
    PnotR <- pnorm(q = e, mean = delta*sqrt(I))
    Score <- (beta - PnotR)^2
    return(Score)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  n.sw     <- suppressWarnings(optim(par = 10^-2, fn = nCRCT, X = X,
                                     delta = delta, beta = beta,
                                     sigma.e = sigma.e, sigma.c = sigma.c,
                                     e = qnorm(1 - alpha))$par)

  if (summary == TRUE){
    print("Determining exact GS design...")
  }

  if (length(set.T) > 1){
    if (stopping == "futility"){
      n.seq.sw <- suppressWarnings(optim(par = n.sw, fn = nSeqCRCTF, X = X,
                                         delta = delta, alpha = alpha,
                                         beta = beta, sigma.e = sigma.e,
                                         sigma.c = sigma.c, set.T = set.T,
                                         gamma.f = gamma.f)$par)
    } else if (stopping == "efficacy"){
      n.seq.sw <- suppressWarnings(optim(par = n.sw, fn = nSeqCRCTE, X = X,
                                         delta = delta, alpha = alpha,
                                         beta = beta, sigma.e = sigma.e,
                                         sigma.c = sigma.c, set.T = set.T,
                                         gamma.e = gamma.e)$par)
    } else if (stopping == "both"){
      n.seq.sw <- suppressWarnings(optim(par = n.sw, fn = nSeqCRCTEF, X = X,
                                         delta = delta, alpha = alpha,
                                         beta = beta, sigma.e = sigma.e,
                                         sigma.c = sigma.c, set.T = set.T,
                                         gamma.f = gamma.f,
                                         gamma.e = gamma.e)$par)
    }
  } else {
    n.seq.sw <- ceiling(n.sw)
  }

  C        <- nrow(X)
  L        <- length(set.T)
  n        <- ceiling(n.seq.sw)
  nvec     <- n*C*set.T
  I        <- informationCRCT(n, X, set.T, sigma.e, sigma.c)
  if (L > 1){
    Lambda <- covariance(L, I)
  } else {
    Lambda <- 1
  }
  if (L > 1){
    if (stopping == "futility"){
      pi2      <- family(L, I, beta, gamma.f)
      f        <- numeric(L)
      for (l in 1:L){
        if (l == 1){
          f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
        } else if ((l > 1) & (l < L)){
          f[l] <- suppressWarnings(optim(par = qnorm(pi2[l],
                                                     mean = delta*sqrt(I[l])),
                                         fn = boundaryF, pi2l = pi2[l],
                                         prefbounds = f[1:(l - 1)],
                                         preebounds = rep(Inf, l - 1),
                                         currSigma = Lambda[1:l, 1:l],
                                         currI = I[1:l], delta = delta)$par)
        } else {
          f[l] <- suppressWarnings(optim(par = qnorm(1 - alpha/L), fn = boundaryE,
                                         pi1l = alpha, prefbounds = f[1:(l - 1)],
                                         preebounds = rep(Inf, L - 1),
                                         currSigma = Lambda)$par)
        }
      }
      e <- c(rep(Inf, L - 1), f[L])
      perf.H0        <- perfF(0, L, f, I, nvec, Lambda)
      perf.H1        <- perfF(delta, L, f, I, nvec, Lambda)
    } else if (stopping == "efficacy"){
      pi1   <- family(L, I, alpha, gamma.e)
      e     <- numeric(L)
      for (l in 1:L){
        if (l == 1){
          e[l] <- qnorm(1 - pi1[l])
        } else {
          e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                         fn = boundaryE, pi1l = pi1[l],
                                         prefbounds = rep(-Inf, l - 1),
                                         preebounds = e[1:(l - 1)],
                                         currSigma = Lambda[1:l, 1:l])$par)
        }
      }
      f <- c(rep(-Inf, L - 1), e[L])
      perf.H0        <- perfE(0, L, e, I, nvec, Lambda)
      perf.H1        <- perfE(delta, L, e, I, nvec, Lambda)
    } else if (stopping == "both"){
      pi1   <- family(L, I, alpha, gamma.e)
      pi2   <- family(L, I, beta, gamma.f)
      e     <- numeric(L)
      f     <- numeric(L)
      for (l in 1:L){
        if (l == 1){
          f[l] <- qnorm(pi2[l], mean = delta*sqrt(I[l]))
          e[l] <- qnorm(1 - pi1[l])
        } else if ((l > 1) & (l < L)){
          f[l] <- suppressWarnings(optim(par = qnorm(pi2[l],
                                                     mean = delta*sqrt(I[l])),
                                         fn = boundaryF, pi2l = pi2[l],
                                         prefbounds = f[1:(l - 1)],
                                         preebounds = e[1:(l - 1)],
                                         currSigma = Sigma[1:l, 1:l],
                                         currI = I[1:l], delta = delta)$par)
          e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                         fn = boundaryE, pi1l = pi1[l],
                                         prefbounds = f[1:(l - 1)],
                                         preebounds = e[1:(l - 1)],
                                         currSigma = Sigma[1:l, 1:l])$par)
        } else {
          e[l] <- suppressWarnings(optim(par = qnorm(1 - pi1[l]),
                                         fn = boundaryE, pi1l = pi1[l],
                                         prefbounds = f[1:(l - 1)],
                                         preebounds = e[1:(l - 1)],
                                         currSigma = Sigma[1:l, 1:l])$par)
          f[l] <- e[l]
        }
      }
      perf.H0        <- perfEF(0, L, e, f, I, nvec, Lambda)
      perf.H1        <- perfEF(delta, L, e, f, I, nvec, Lambda)
    }
  } else {
    e <- qnorm(1 - alpha)
    f <- qnorm(1 - alpha)
    perf.H0 <- c(pnorm(q = e, mean = 0, lower.tail = F), nvec)
    perf.H1 <- c(pnorm(q = e, mean = delta*sqrt(I), lower.tail = F), nvec)
  }


  if (summary == TRUE){
    print("Outputting...")
  }

  output <- list(alpha = alpha, beta = beta, delta = delta, e = e, f = f,
                 gamma.f = gamma.f, gamma.e = gamma.e, I = I, Lambda = Lambda,
                 n = n, Nvec = nvec, n.sw = n.sw, perf.H0 = perf.H0, perf.H1 = perf.H1,
                 set.T = set.T, sigma.c = sigma.c, sigma.e = sigma.e, X = X)
  class(output) <- "gsSW"
  return(output)
}
