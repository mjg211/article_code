################################################################################
# Function:    optimal.gs.sw.R                                                 #
# Author:      Michael J. Grayling (mjg211@cam.ac.uk)                          #
# Last Edited: 21/02/2018                                                      #
################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
#     design - The design of the group sequential SW-CRT under consideration.  #
#              Should be a list containing at least the following components:  #
#              X, Ti, C, n, set.T, e, f, sigma.c, sigma.e, alphanumber of      #
#              clusters.                                                       #
#        tau - The treatment effects to assess performance at.                 #
#         mu - The mean response from an individual in period 1 when in the    #
#              control condition. Used in simulating data.                     #
#         pi - The period effects to use in simulating data.                   #
# replicates - The number of replicate simulations to use in assessing average #
#              peformance.                                                     #
#       seed - The random number seed to use.                                  #
#    summary - Allows printing of progress to the console.                     #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# A list containing the following objects:                                     #
#                                                                              #
#            power - Rejection rate for each value in tau.                     #
#       bias.naive - Average bias of the naive estimator for each value in     #
#                    tau.                                                      #
#         bias.adj - Average bias of the adjusted estimator for each value in  #
#                    tau.                                                      #
#       rmse.naive - RMSE of the naive estimator for each value in tau.        #
#         rmse.adj - RMSE of the adjusted estimator for each value in tau.     #
# p.value.naive.25 - 25th percentile of the p-values derived from the naive    #
#                    estimator for each value in tau.                          #
# p.value.naive.50 - 50th percentile of the p-values derived from the naive    #
#                    estimator for each value in tau.                          #
# p.value.naive.75 - 75th percentile of the p-values derived from the naive    #
#                    estimator for each value in tau.                          #
#   p.value.adj.25 - 25th percentile of the p-values derived from the          #
#                    adjusted estimator for each value in tau.                 #
#   p.value.adj.50 - 50th percentile of the p-values derived from the          #
#                    adjusted estimator for each value in tau.                 #
#   p.value.adj.75 - 75th percentile of the p-values derived from the          #
#                    adjusted estimator for each value in tau.                 #
#   coverage.naive - Coverage of the naive confidence interval for each value  #
#                    in tau.                                                   #
#     coverage.adj - Coverage of the adjusted confidence interval for each     #
#                    value in tau.                                             #
################################################################################

library(mvtnorm)
library(Matrix)

analyse.gs.sw <- function(design, tau = c(0, design$delta), mu = 0,
                          pi = rep(0, ncol(design$X)), replicates = 100000,
                          summary = TRUE, seed = Sys.time()){
  
  set.seed(seed)
  
  pValue <- function(theta, L, e, f, I, Sigma){
    PE        <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PE[l] <- pmvnorm(lower = e[l], upper = Inf,
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PE[l] <- pmvnorm(lower = c(f[1:(l - 1)], e[l]),
                         upper = c(e[1:(l - 1)], Inf),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    return(sum(PE))
  }
  
  pValueRoot <- function(theta, L, u.b, l.b, I, Sigma, h){
    PE        <- numeric(L)
    for (l in 1:L){
      if (l == 1){
        PE[l] <- pmvnorm(lower = u.b[l], upper = Inf,
                         mean = theta*sqrt(I[l]), sigma = 1)[1]
      } else {
        PE[l] <- pmvnorm(lower = c(l.b[1:(l - 1)], u.b[l]),
                         upper = c(u.b[1:(l - 1)], Inf),
                         mean = theta*sqrt(I[1:l]),
                         sigma = Sigma[1:l, 1:l])[1]
      }
    }
    return(sum(PE) - h)
  }
  
  singleTrial <- function(rep, Ti, C, n, means, periods, clusters,
                          treatments, chol.Sigma.22, chol.Sigma.c.bars,
                          Sigma.fact.mats, set.T, e, f, I,
                          Lambda, alpha, gls.factors, summary){
    if (rep%%1000 == 0 & summary == TRUE){
      print(paste("...replicate ", rep, "completed..."))
    }
    responses     <- NULL
    response.c    <- list()
    mean.c        <- list()
    for (t in 1:Ti){
      response.t  <- numeric(C*n)
      mean.t      <- means[(1 + (t - 1)*n*C):(t*n*C)]
      if (t == 1){
        for (c in 1:C){
          mean.c[[c]]     <- mean.t[(1 + (c - 1)*n):(c*n)]
          resp.c          <- as.numeric(rnorm(n)%*%chol.Sigma.22) +
                               mean.t[(1 + (c - 1)*n):(c*n)]
          response.c[[c]] <- resp.c
          response.t[(1 + (c - 1)*n):(c*n)] <- resp.c
        }
      } else {
        for (c in 1:C){
          resp.c          <- as.numeric(rnorm(n)%*%chol.Sigma.c.bars[[t - 1]] +
                                          mean.t[(1 + (c - 1)*n):(c*n)] +
                                          t(Sigma.fact.mats[[t - 1]]%*%
                                              (response.c[[c]] - mean.c[[c]])))
          mean.c[[c]]     <- c(mean.c[[c]], mean.t[(1 + (c - 1)*n):(c*n)])
          response.c[[c]] <- c(response.c[[c]], resp.c)
          response.t[(1 + (c - 1)*n):(c*n)]  <- resp.c
        }
      }
      responses     <- c(responses, response.t)
      if (t %in% set.T){
        df.analysis <- data.frame(Period = periods[1:(t*C*n)],
                                  Cluster = clusters[1:(t*C*n)],
                                  Treatment = treatments[1:(t*C*n)],
                                  Response = responses)
        rearrange   <- NULL
        for (c in 1:C){
          rearrange <- rbind(rearrange,
                             df.analysis[which(df.analysis$Cluster == c), ])
        }
        beta.gls <- gls.factors[[t]]%*%matrix(rearrange$Response, n*C*t, 1)
        mle.tau  <- beta.gls[t + 1]
        Z.t      <- mle.tau*sqrt(I[which(set.T == t)])
        if (Z.t > e[which(set.T == t)] | Z.t <= f[which(set.T == t)]){
          L             <- which(set.T == t)
          p.value.naive <- pnorm(Z.t, 0, lower.tail = F)
          lci.naive     <- mle.tau - qnorm(1 - alpha)*sqrt(1/I[L])
          if (L > 1){
            p.value.adj <- pValue(0, L, c(e[1:(L - 1)], Z.t), f[1:L], I[1:L],
                                  Lambda[1:L, 1:L])
            med.tau     <- uniroot(pValueRoot, c(-10000, 10000), L = L,
                                   u.b = c(e[1:(L - 1)], Z.t),
                                   l.b = f[1:L], I = I[1:L],
                                   Sigma = Lambda[1:L, 1:L], h = 0.5)$root
            lci.adj     <- uniroot(pValueRoot, c(-10000, 10000), L = L,
                                   u.b = c(e[1:(L - 1)], Z.t),
                                   l.b = f[1:L], I = I[1:L],
                                   Sigma = Lambda[1:L, 1:L], h = alpha)$root
          } else {
            p.value.adj <- pnorm(Z.t, 0, lower.tail = F)
            med.tau     <- uniroot(pValueRoot, c(-10000, 10000), L = L,
                                   u.b = Z.t, l.b = f[1:L], I = I[1:L],
                                   Sigma = Lambda[1:L, 1:L], h = 0.5)$root
            lci.adj     <- uniroot(pValueRoot, c(-10000, 10000), L = L,
                                   u.b = Z.t, l.b = f[1:L], I = I[1:L],
                                   Sigma = Lambda[1:L, 1:L], h = alpha)$root
          }
          return(c(Z.t > e[which(set.T == t)], n*C*t, mle.tau, med.tau,
                   p.value.naive, p.value.adj, lci.naive, lci.adj))
        }
      }
    }
  }
  
  wrapper <- function(rep, Ti, C, n, means, periods, clusters,
                      treatments, chol.Sigma.22, chol.Sigma.c.bars,
                      Sigma.fact.mats, set.T, e, f, I,
                      Lambda, alpha, gls.factors, summary){
    result <- singleTrial(rep, Ti, C, n, means, periods, clusters,
                          treatments, chol.Sigma.22, chol.Sigma.c.bars,
                          Sigma.fact.mats, set.T, e, f, I,
                          Lambda, alpha, gls.factors, summary)
    return(result)
  }
  
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
  
  ##### MAIN COMPUTATIONS ######################################################
  
  X                <- design$X
  Ti               <- ncol(X)
  C                <- nrow(X)
  n                <- design$n
  set.T            <- design$set.T
  L                <- length(set.T)
  e                <- design$e
  f                <- design$f
  sigma.c          <- design$sigma.c
  sigma.e          <- design$sigma.e
  alpha            <- design$alpha
  I                <- informationCRT(n, X, set.T, sigma.e, sigma.c)
  Lambda           <- diag(1, L, L)
  for (i in 2:L){
    for (j in 1:(i - 1)){
      Lambda[i, j] <- sqrt(I[j]/I[i])
      Lambda[j, i] <- Lambda[i, j]
    }
  }
  
  periods       <- numeric(C*Ti*n)
  treatments    <- numeric(C*Ti*n)
  clusters      <- numeric(C*Ti*n)
  for (j in 1:Ti){
    periods[(1 + (j - 1)*n*C):(j*n*C)]    <- rep(j, n*C)
    treatment.j                           <- numeric(n*C)
    cluster.j                             <- numeric(n*C)
    for (c in 1:C){
      treatment.j[(1 + (c - 1)*n):(c*n)]  <- rep(X[c, j], n)
      cluster.j[(1 + (c - 1)*n):(c*n)]    <- rep(c, n)
    }
    treatments[(1 + (j - 1)*C*n):(j*C*n)] <- treatment.j
    clusters[(1 + (j - 1)*C*n):(j*C*n)]   <- cluster.j
  }
  
  if (summary == TRUE){
    print("Computing all factors required for trial simulation...")
  }
  
  chol.Sigma.c.bars <- list()
  Sigma.fact.mats   <- list()
  for (t in 1:(Ti - 1)){
    Sigma.11               <- as(matrix(sigma.c^2, n, n), "dspMatrix") +
                                Diagonal(x = rep(sigma.e^2, n))
    Sigma.12               <- as(matrix(sigma.c^2, ncol = n*t, nrow = n),
                                 "dgeMatrix")
    Sigma.22.inv           <- Diagonal(n*t)/(sigma.e^2) -
                                ((sigma.c^2)/
                                   (sigma.e^4 + (sigma.c*sigma.e)^2*(1 + n*t)))*
                                      as(matrix(1, nrow = n*t, ncol = n*t),
                                         "dspMatrix")
    Sigma.fact             <- Sigma.12%*%Sigma.22.inv
    Sigma.c.bar            <- as(Sigma.11 - Sigma.fact%*%t(Sigma.12),
                                 "dspMatrix")
    chol.Sigma.c.bar       <- matrix(0, nrow = n, ncol = n)
    chol.Sigma.c.bar[upper.tri(chol.Sigma.c.bar,
                               diag = TRUE)] <- chol(Sigma.c.bar)@x
    chol.Sigma.c.bars[[t]] <- as(chol.Sigma.c.bar, "dtpMatrix")
    Sigma.fact.mats[[t]]   <- Sigma.fact
  }
  remove(Sigma.11)
  remove(Sigma.12)
  remove(Sigma.22.inv)
  remove(Sigma.c.bar)
  
  Sigma.22      <- as(matrix(sigma.c^2, n, n),
                 "dspMatrix") + Diagonal(x = rep(sigma.e^2, n))
  chol.Sigma.22 <- matrix(0, nrow = n, ncol = n)
  chol.Sigma.22[upper.tri(chol.Sigma.22, diag = TRUE)] <- chol(Sigma.22)@x
  chol.Sigma.22 <- as(chol.Sigma.22, "dtpMatrix")
  
  gls.factors <- list()
  for (t in 1:Ti){
    D.t       <- matrix(0, nrow = n*C*t, ncol = 1 + t)
    D.t[, 1]  <- 1
    for (i in 1:C){
      for (j in 1:t){
        if (j > 1){
          D.t[(1 + (i - 1)*n*t + n*(j - 1)):((i - 1)*n*t + n*j), j]     <- 1
        }
        if (X[i, j] == 1){
          D.t[(1 + (i - 1)*n*t + n*(j - 1)):((i - 1)*n*t + n*j), t + 1] <- 1
        }
      }
    }
    Sigma                              <- matrix(0, n*C*t, n*C*t)
    for (i in 1:C){
      Sigma[(1 + (i - 1)*n*t):(i*n*t),
            (1 + (i - 1)*n*t):(i*n*t)] <- matrix(sigma.c^2, n*t, n*t) +
                                            diag(sigma.e^2, n*t, n*t)
    }
    gls.factors[[t]] <- ginv(t(D.t)%*%ginv(Sigma)%*%D.t)%*%t(D.t)%*%ginv(Sigma)
  }
  
  power            <- numeric(length(tau))
  bias.naive       <- numeric(length(tau))
  bias.adj         <- numeric(length(tau))
  rmse.naive       <- numeric(length(tau))
  rmse.adj         <- numeric(length(tau))
  p.value.naive.25 <- numeric(length(tau))
  p.value.naive.50 <- numeric(length(tau))
  p.value.naive.75 <- numeric(length(tau))
  p.value.adj.25   <- numeric(length(tau))
  p.value.adj.50   <- numeric(length(tau))
  p.value.adj.75   <- numeric(length(tau))
  coverage.naive   <- numeric(length(tau))
  coverage.adj     <- numeric(length(tau))
  for (i in 1:length(tau)){
    if (summary == TRUE){
      print(paste("...beginning simulations for tau =", tau[i], "..."))
      print(paste("...current system time is", Sys.time(), "..."))
    }
    means                              <- numeric(C*Ti*n)
    for (j in 1:Ti){
      mean.j                           <- numeric(n*C)
      for (c in 1:C){
        mean.j[(1 + (c - 1)*n):(c*n)]  <- rep(tau[i]*X[c, j] + pi[j], n)
      }
      means[(1 + (j - 1)*C*n):(j*C*n)] <- mean.j
    }
    results <- lapply(1:replicates, wrapper, Ti = Ti, C = C, n = n,
                      means = means, periods = periods, clusters = clusters,
                      treatments = treatments, chol.Sigma.22 = chol.Sigma.22,
                      chol.Sigma.c.bars = chol.Sigma.c.bars,
                      Sigma.fact.mats = Sigma.fact.mats, set.T = set.T, e = e,
                      f = f, I = I, Lambda = Lambda, alpha = alpha,
                      gls.factors = gls.factors, summary = summary)
    results             <- matrix(unlist(results), ncol = 8, byrow = TRUE)
    power[i]            <- mean(results[, 1])
    bias.naive[i]       <- mean(results[, 3]) - tau[i]
    bias.adj[i]         <- mean(results[, 4]) - tau[i]
    rmse.naive[i]       <- sqrt(sum((results[, 3] - tau[i])^2)/replicates)
    rmse.adj[i]         <- sqrt(sum((results[, 4] - tau[i])^2)/replicates)
    p.value.naive.25[i] <- quantile(results[, 5], 0.25)
    p.value.naive.50[i] <- quantile(results[, 5], 0.50)
    p.value.naive.75[i] <- quantile(results[, 5], 0.75)
    p.value.adj.25[i]   <- quantile(results[, 6], 0.25)
    p.value.adj.50[i]   <- quantile(results[, 6], 0.50)
    p.value.adj.75[i]   <- quantile(results[, 6], 0.75)
    coverage.naive[i]   <- length(which(results[, 7] <= tau[i]))/replicates
    coverage.adj[i]     <- length(which(results[, 8] <= tau[i]))/replicates
  }
  if (summary == TRUE){
    print("...returning results...")
  }
  return(list(power = power, bias.naive = bias.naive, bias.adj = bias.adj,
              rmse.naive = rmse.naive, rmse.adj = rmse.adj,
              p.value.naive.25 = p.value.naive.25,
              p.value.naive.50 = p.value.naive.50,
              p.value.naive.75 = p.value.naive.75,
              p.value.adj.25 = p.value.adj.25, p.value.adj.50 = p.value.adj.50,
              p.value.adj.75 = p.value.adj.75, coverage.naive = coverage.naive,
              coverage.adj = coverage.adj))
}
