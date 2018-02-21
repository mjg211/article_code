################################################################################
# Function: optimal.gs.sw.R                                                    #
# Author:   Michael J. Grayling (mjg211@cam.ac.uk)                             #
# Last Edited: 21/02/2018                                                      #
################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
#       C - The number of clusters.                                            #
#      Ti - The number of time periods.                                        #
# sigma.c - The between cluster standard deviation.                            #
# sigma.e - The residual standard deviation.                                   #
#   alpha - The desired type-I error-rate.                                     #
#    beta - The desired type-II error-rate.                                    #
#   delta - The treatment effect to power for.                                 #
#   set.T - The set of time periods after which an interim analysis should be  #
#           conducted.                                                         #
#       w - A vector of weights to use in the objective function.              #
#    N.SW - The penalty to use in the objective function.                      #
#   n.max - The maximum allowed value for n.                                   #
#  rho.CE - The rarity parameter to use in CEoptim.                            #
#    N.CE - The number of samples to use in CEoptim.                           #
#    seed - The random number seed to use.                                     #
# summary - Allows printing of progress to the console.                        #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# A list containing all input parameters and the following objects:            #
#                                                                              #
#        e - The optimal vector of efficacy boundaries.                        #
#        f - The optimal vector of futility boundaries.                        #
#        I - The vector of information levels for the optimal design.          #
#   Lambda - The covariance matrix of the standard test statistics for the     #
#            optimal design.                                                   #
#        O - The value of the objective function for the optimal design.       #
#     perf - A vector of length 4 containing P(0), P(delta), ESS(0), and       #
#            ESS(delta) for the optimal design.                                #
# run.time - The run time of the executed command (time for optimal design     #
#            determination).                                                   #
#        X - The matrix of binary treatment indicators for the optimal         #
#            design.                                                           #
################################################################################

library(mvtnorm)
library(CEoptim)

optimal.gs.sw <- function(C = 4, Ti = 5, sigma.c = sqrt(0.02),
                          sigma.e = sqrt(0.51), alpha = 0.05, beta = 0.1,
                          delta = 0.2, set.T = 2:5, w = rep(1/3, 3),
                          N.SW = 1400, n.max = 400, rho.CE = 0.01,
                          N.CE = 10000*(C + 2*length(set.T)),
                          seed = Sys.time(), summary = TRUE){
  
  # Set seed for repeatability
  set.seed(seed)
  
  # Facilitate storing run time for reference
  start.time <- Sys.time()
  
  ##### ERROR CHECKING ########################################################

  if ((C%%1 != 0) | (C < 2)){
    stop("C must be a whole number greater than or equal to 2.")
  }
  if ((Ti%%1 != 0) | (Ti < 2)){
    stop("Ti must be a whole number greater than or equal to 2.")
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
  if ((length(set.T) == 1) | (any(set.T%%1 != 0)) | (any(set.T <= 0)) |
        (sum(set.T[2:length(set.T)] >=
               set.T[1:(length(set.T) - 1)]) != length(set.T) - 1) |
        (set.T[length(set.T)] != Ti)){
    stop("set.T must be a vector of positive integers of length at least 2, where each element is strictly larger than the last, and the final element is Ti.")
  }
  if ((length(w) != 3) | (any(w < 0)) | (sum(w[1:2]) == 0)){
    stop("w must be a vector of length 3, consisting of elements greater than or equal to 0, with at least one of the first 2 elements strictly positive.")
  }
  if (N.SW <= 0){
    stop("N.SW must be strictly positive.")
  }
  if ((n.max%%1 != 0) | (n.max < 2)){
    stop("n.max must be a whole number greater than or equal to 2.")
  }
  if (rho.CE <= 0){
    stop("rho.CE must be strictly positive.")
  }
  if ((N.CE%%1 != 0) | (N.CE <= 0)){
    stop("N.CE must be a positive whole number.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }

  ##### FUNCTION INITIALISATION ###############################################

  # Function to determine performance of chosen design
  perfEF <- function(theta, L, e, f, I, Nvec, Sigma){
    PE        <- numeric(L)
    PF        <- numeric(L)
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
    PR        <- sum(PE)
    EN        <- sum(Nvec*(PE + PF))
    return(c(PR, EN))
  }
  
  # Function to return information at each analysis
  informationCRT <- function(n, X, set.T, sigma.e, sigma.c){
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
        I[i] <- ((C*U - W)*sigma2 +
                   (U^2 + C*set.T[i]*U - set.T[i]*W - C*V)*sigma.c^2)/
                  (C*sigma2*(sigma2 + set.T[i]*sigma.c^2))
      }
    }
    return(I)
  }
  
  # Internal function for evaluating objective function
  internalObjFn <- function(f.r, C, Ti, n.S, delta, alpha, beta,
                            sigma.c, sigma.e, set.T, w, N.SW, L){
    n <- n.S[1]
    S <- n.S[2:(C + 1)]
    if (length(unique(S)) == 1){
      Score             <- Inf
    } else {
      X                 <- matrix(0, nrow = C, ncol = Ti)
      for (c in 1:C){
        if (S[c] != Ti + 1){
          X[c, S[c]:Ti] <- 1
        }
      }
      f                 <- f.r[1:L]
      e                 <- c(f[1:(L - 1)] + f.r[(L + 1):(2*L - 1)],
                             f[L])
      Nvec              <- n*C*set.T
      I                 <- informationCRT(n, X, set.T, sigma.e, sigma.c)
      Lambda <- diag(1, L, L)
      for (i in 2:L){
        for (j in 1:(i - 1)){
          Lambda[i, j]  <- sqrt(I[j]/I[i])
          Lambda[j, i]  <- Lambda[i, j]
        }
      }
      if (any(as.vector(is.nan(Lambda)))){
        Score   <- Inf
      } else {
        perf.H0 <- perfEF(0, L, e, f, I, Nvec, Lambda)
        perf.H1 <- perfEF(delta, L, e, f, I, Nvec, Lambda)
        Score   <- sum(w*c(perf.H0[2], perf.H1[2], Nvec[L])) + 
                     N.SW*(as.numeric(perf.H0[1] > alpha)*
                             (perf.H0[1] - alpha)/alpha + 
                           as.numeric(1 - perf.H1[1] > beta)*
                             (1 - perf.H1[1] - beta)/beta)
      }
    }
    return(Score)
  }
  
  # Objective function evaluator
  objFn <- function(f.r, n.S, C, Ti, delta, alpha, beta, sigma.c, sigma.e,
                    set.T, w, N.SW, L){
    n.S.int <- n.S + 1L
    Score   <- internalObjFn(f.r, C, Ti, n.S.int, delta, alpha, beta,
                             sigma.c, sigma.e, set.T, w, N.SW, L)
    return(Score)
  }

  ##### MAIN COMPUTATIONS #####################################################

  if (summary == TRUE){
    print("Initialising all required variables...")
  }

  # Initialise all required parameters
  L            <- length(set.T)
  contMean     <- rep(0, 2*L - 1)
  contSD       <- rep(10, 2*L - 1)
  contConstMat <- matrix(0, 2*L - 1, 2*L - 1)
  for (i in (L + 1):(2*L - 1)){
    contConstMat[i, i] <- -1
  }
  contConstVec <- numeric(2*L - 1)
  discCat      <- as.integer(c(n.max, set.T[1], rep(Ti + 1, C - 1)))
  discSmooth   <- 0.5

  if (summary == TRUE){
    print("Searching for optimal GS design...")
  }

  # Search for the optimal design
  optimal.design <- CEoptim(f = objFn,
                            f.arg = list(C = C, Ti = Ti, delta = delta,
                                         alpha = alpha, beta = beta,
                                         sigma.c = sigma.c, sigma.e = sigma.e,
                                         set.T = set.T, w = w, N.SW = N.SW,
                                         L = L),
                            continuous = list(mean = contMean, sd = contSD,
                                              conMat = contConstMat,
                                              conVec = contConstVec,
                                              smoothMean = 0.5,
                                              smoothSd = 0.5),
                            discrete = list(categories = discCat,
                                            smoothProb = discSmooth),
                            N = as.integer(N.CE), rho = rho.CE,
                            verbose = summary, noImproveThr = 20)
  
  # Return the parameters of the identified design and its performance
  O <- optimal.design$optimum
  n <- optimal.design$optimizer$discrete[1] + 1
  S <- optimal.design$optimizer$discrete[2:(C + 1)] + 1
  f <- optimal.design$optimizer$continuous[1:L]
  e <- c(f[1:(L - 1)] + optimal.design$optimizer$continuous[(L + 1):(2*L - 1)],
         f[L])
  
  X                 <- matrix(0, nrow = C, ncol = Ti)
  for (c in 1:C){
    if (S[c] != Ti + 1){
      X[c, S[c]:Ti] <- 1
    }
  }
  Nvec              <- n*C*set.T
  I                 <- informationCRT(n, X, set.T, sigma.e, sigma.c)
  Lambda            <- diag(1, L, L)
  for (i in 2:L){
    for (j in 1:(i - 1)){
      Lambda[i, j]  <- sqrt(I[j]/I[i])
      Lambda[j, i]  <- Lambda[i, j]
    }
  }
  perf.H0           <- perfEF(0, L, e, f, I, Nvec, Lambda)
  perf.H1           <- perfEF(delta, L, e, f, I, Nvec, Lambda)
  perf              <- c(perf.H0[1], perf.H1[1], perf.H0[2], perf.H1[2])
  names(perf)       <- c("P(0)", "P(delta)", "ESS(0)", "ESS(delta)")
  
  end.time          <- Sys.time()
  run.time          <- as.numeric(difftime(end.time, start.time, units = "min"))

  if (summary == TRUE){
    print("Outputting...")
  }

  output           <- list(alpha = alpha, beta = beta, C = C, delta = delta,
                           e = e, f = f, I = I, Lambda = Lambda, N.CE = N.CE,
                           N.SW = N.SW, n = n, O = O, perf = perf,
                           rho.CE = rho.CE, run.time = run.time, seed = seed,
                           set.T = set.T, sigma.c = sigma.c, sigma.e = sigma.e,
                           summary = summary, Ti = Ti, X = X)
  return(output)
}
