##### Load required packages ###################################################

library(iterpc)
library(MASS)
library(snowfall)
library(CEoptim)
library(mvtnorm)
library(ggplot2)
library(plyr)
library(tidyr)
library(Matrix)

##### admissible.ma.sw.crt(): The main function for computing the designs ######

################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
# D            - The number of treatment arms.                                 #
# poss.mCT     - A matrix with three columns, with each row corresponding to   #
#                an allowed combination of m, C, and T in that order. If left  #
#                NULL set.T, set.C, and set.m are used instead.                #
# set.T        - The set of allowed values of T.                               #
# set.C        - The set of allowed values of C.                               #
# set.m        - The set of allowed values of m.                               #
# alpha        - The desired type-I error rate.                                #
# beta         - The desired type-II error rate.                               #
# sigma.c      - The s.d. of the cluster level random effects.                 #
# sigma.pi     - The s.d. of the cluster-period random effects.                #
# sigma.s      - The s.d. of the individual level random effects.              #
# sigma.eps    - The residual s.d.                                             #
# design       - Either "cross-sectional" or "cohort".                         #
# delta        - The vector of treatment effects at which we power the trial.  #
# w            - The weight for the optimality criterion.                      #
# optimality   - The choice between D, E, or A optimality.                     #
# exhaustive   - Logical. Should an exhaustive or stochastic search be         #
#                utilised?                                                     #
# power.for    - Is "individual" or "combined" power desired?                  #
# mcc          - Should the bonferonni ("bonf") or no ("none") multiple        #
#                comparison correction be used?                                #
# parallel     - Should parallel computation be used for the exhaustive        #
#                search?                                                       #
# cpus         - The umber of cpus to use in parallel computation.             #
# N            - The integer value of N passed to CEoptim. See ?CEoptim.       #
# rho          - Value of rho passed to CEoptim. See ?CEoptim.                 #
# smoothProb   - Smoothing parameters passed to CEoptim. See ?CEoptim.         #
# noImproveThr - No improvement threshold passed to CEoptim. See ?CEoptim.     #
# summary      - Logical. Should a summary of the functions progress be        #
#                printed?                                                      #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# Returns the optimal X, and the optimal m if relevant, in a list.             #
################################################################################

admissible.ma.sw.crt <- function(D = 3, poss.mCT = NULL, set.T = 2:6,
                                 set.C = 2:6, set.m = 2:100, set.X = "all",
                                 alpha = 0.05, beta = 0.2, sigma.c = 1,
                                 sigma.pi = 0, sigma.s = 0, sigma.eps = 1,
                                 design = "cross_sectional",
                                 delta = rep(0.41, D - 1), w = 0,
                                 optimality = "D", exhaustive = TRUE,
                                 power.for = "individual", mcc = "bonf",
                                 parallel = TRUE, cpus = 6, N = 10000L,
                                 rho = 0.001, smoothProb = 0.5,
                                 noImproveThr = 5, summary = TRUE) {
  
  if (all(!exhaustive, design == "cohort")) {
    stop("Stochastic searchs for cohort designs not supported.")
  }
  
  ##### Function Initialisation ################################################
  
  check.rank <- function(index) {
    if (length(unique(all.combs[index, ])) > 1) {
      X                <- all.rows[all.combs[index, ], ]
      if (length(unique(as.vector(X))) == D){
        A              <- matrix(0, nrow = C*Ti,
                                 ncol = 1 + (Ti - 1) + (D - 1))
        A[, 1]         <- 1
        for (k in 1:C) {
          for (l in 1:Ti) {
            if (X[k, l] != 0) {
              A[(1 + (k - 1)*Ti + (l - 1)):((k - 1)*Ti + l),
                (Ti + 1):(Ti + X[k, l])]                        <- 1
            }
            if (l > 1){
              A[(1 + (k - 1)*Ti + (l - 1)):((k - 1)*Ti + l), l] <- 1
            }
          }
        }
        check          <- 0
        if (rankMatrix(A)[1] < ncol(A)) {
          return(F)
        } else {
          return(T)
        }
      } else {
        return(F)
      }
    } else {
      return(F)
    }
  }
  
  ind.design <- function(index) {
    Ti  <- poss.designs[index, 1]
    C   <- poss.designs[index, 2]
    m   <- poss.designs[index, 3]
    X   <- poss.switching.matrices[[Ti]][[C]][[poss.designs[index, 4]]]
    if (D == 2) {
      U <- sum(X)
      V <- sum(rowSums(X)^2)
      W <- sum(colSums(X)^2)
      if (design == "cohort") {
        sigma2  <- sigma.c^2 + sigma.pi^2 + sigma.s^2 + sigma.eps^2
        rho0    <- (sigma.c^2 + sigma.pi^2)/sigma2
        rho1    <- sigma.c^2/sigma2
        rho2    <- (sigma.c^2 + sigma.s^2)/sigma2
        psi     <- 1 + (m - 1)*rho0 - (m - 1)*rho1 - rho2
        xi      <- (m - 1)*rho1 + rho2
        gamma   <- psi + Ti*xi
        cov.tau <- (C*sigma2*gamma*psi/m)/((C*U - W)*gamma + (U^2 - C*V)*xi)
      } else {
        sigma2  <- sigma.eps^2/m
        cov.tau <- (C*sigma2*(sigma2 + Ti*sigma.c^2))/
                     ((C*U - W)*sigma2 + (U^2 + C*Ti*U - Ti*W - C*V)*sigma.c^2)
      }
    } else {
      Sigma.inv   <- matrix(0, nrow = m*C*Ti, ncol = m*C*Ti)
      Sigma.inv.c <- matrix(Sigma.inv.c.mats[[Ti]][[m]][2], m*Ti, m*Ti) +
                       diag(Sigma.inv.c.mats[[Ti]][[m]][1] -
                              Sigma.inv.c.mats[[Ti]][[m]][2], m*Ti, m*Ti)
      for (c in 1:C) {
        Sigma.inv[(1 + (c - 1)*m*Ti):(c*m*Ti),
                  (1 + (c - 1)*m*Ti):(c*m*Ti)] <- Sigma.inv.c
      }
      A      <- matrix(0, nrow = m*C*Ti, ncol = 1 + (Ti - 1) + (D - 1))
      A[, 1] <- 1
      for (i in 1:C) {
        for (j in 1:Ti) {
          if (X[i, j] != 0) {
            A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m),
              (Ti + 1):(Ti + X[i, j])]                                <- 1
          }
          if (j > 1) {
            A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m), j] <- 1
          }
        }
      }
      cov.beta <- ginv(t(A)%*%Sigma.inv%*%A)
      cov.tau  <- cov.beta[(Ti + 1):(Ti + D - 1), (Ti + 1):(Ti + D - 1)]
    }
    if (mcc == "none") {
      e <- qnorm(1 - alpha)
    } else if (mcc == "bonf") {
      e <- qnorm(1 - alpha/(D - 1))
    } else if (mcc == "dunnett") {
      e <- qmvnorm(1 - alpha, sigma = cov.tau)$quantile
    }
    obs     <- m*C*Ti
    if (D > 2) {
      opt.D <- det(cov.tau)
      opt.A <- mean(diag(cov.tau))
      opt.E <- max(diag(cov.tau))
    } else {
      opt.D <- opt.A <- opt.E <- cov.tau
    }
    if (D > 2) {
      cov.z    <- diag(1/sqrt(diag(cov.tau)))%*%cov.tau%*%
                    diag(1/sqrt(diag(cov.tau)))
    } else {
      cov.z    <- 1
    }
    type.I     <- 1 - pmvnorm(lower = rep(-Inf, D - 1), upper = rep(e, D - 1),
                              sigma = cov.z, mean = rep(0, D - 1))[1]
    power      <- numeric(1 + (D - 1))
    if (D > 2) {
      power[1] <- 1 - pmvnorm(lower = rep(-Inf, D - 1), upper = rep(e, D - 1),
                              sigma = cov.z,
                              mean = delta/sqrt(diag(cov.tau)))[1]
      for (d in 1:(D - 1)) {
        power[d + 1] <- pnorm(e, mean = delta[d]/sqrt(diag(cov.tau)[d]),
                              lower.tail = FALSE)
      }
    } else {
      power    <- rep(pnorm(e, mean = delta/sqrt(cov.tau),
                            lower.tail = FALSE), 2)
    }
    return(c(e, type.I, power, obs, opt.D, opt.A, opt.E))
  }
  
  ind.design.stoch.with.n <- function(par, D, Ti, C, set.m, all.switches,
                                      Sigma.inv.c.mats, optimality) {
    m        <- set.m[par[1] + 1]
    switches <- par[2:(C + 1)] + 1
    if (length(unique(switches)) == 1) {
      return(Inf)
    }
    X        <- matrix(0, nrow = C, ncol = Ti)
    for (c in 1:C) {
      X[c, ] <- all.switches[switches[c], ]
    }
    if (length(unique(as.vector(X))) < D) {
      return(Inf)
    }
    A      <- matrix(0, nrow = m*C*Ti, ncol = 1 + (Ti - 1) + (D - 1))
    A[, 1] <- 1
    for (i in 1:C) {
      for (j in 1:Ti) {
        if (X[i, j] != 0) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m),
            (Ti + 1):(Ti + X[i, j])]                                <- 1
        }
        if (j > 1) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m), j] <- 1
        }
      }
    }
    if (rankMatrix(A)[1] < ncol(A)) {
      return(Inf)
    }
    Sigma.inv   <- matrix(0, nrow = m*C*Ti, ncol = m*C*Ti)
    Sigma.inv.c <- matrix(Sigma.inv.c.mats[[m]][2], m*Ti, m*Ti) +
                     diag(Sigma.inv.c.mats[[m]][1] -
                            Sigma.inv.c.mats[[m]][2], m*Ti, m*Ti)
    for (c in 1:C) {
      Sigma.inv[(1 + (c - 1)*m*Ti):(c*m*Ti),
                (1 + (c - 1)*m*Ti):(c*m*Ti)] <- Sigma.inv.c
    }
    cov.beta   <- ginv(t(A)%*%Sigma.inv%*%A)
    cov.tau    <- cov.beta[(Ti + 1):(Ti + D - 1), (Ti + 1):(Ti + D - 1)]
    if (optimality == "D") {
      if (D > 2) {
        obj.fn <- det(cov.tau)
      } else {
        obj.fn <- cov.tau
      }
    } else if (optimality == "A") {
      if (D > 2) {
        obj.fn <- mean(diag(cov.tau))
      } else {
        obj.fn <- cov.tau
      }
    } else if (optimality == "E") {
      if (D > 2) {
        obj.fn <- max(diag(cov.tau))
      } else {
        obj.fn <- cov.tau
      }
    }
    return(obj.fn)
  }
  
  ind.design.stoch.without.n <- function(par, D, Ti, m, C, all.switches,
                                         Sigma.inv.c.mats, optimality) {
    switches <- par + 1
    switches <- as.numeric(switches)
    if (length(unique(switches)) == 1) {
      return(Inf)
    }
    X        <- matrix(0, nrow = C, ncol = Ti)
    for (c in 1:C) {
      X[c, ] <- all.switches[switches[c], ]
    }
    if (length(unique(as.vector(X))) < D) {
      return(Inf)
    }
    A      <- matrix(0, nrow = m*C*Ti, ncol = 1 + (Ti - 1) + (D - 1))
    A[, 1] <- 1
    for (i in 1:C) {
      for (j in 1:Ti) {
        if (X[i, j] != 0) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m),
            (Ti + 1):(Ti + X[i, j])]                                <- 1
        }
        if (j > 1) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m), j] <- 1
        }
      }
    }
    if (rankMatrix(A)[1] < ncol(A)) {
      return(Inf)
    }
    Sigma.inv   <- matrix(0, nrow = m*C*Ti, ncol = m*C*Ti)
    Sigma.inv.c <- matrix(Sigma.inv.c.mats[[m]][2], m*Ti, m*Ti) +
                     diag(Sigma.inv.c.mats[[m]][1] -
                            Sigma.inv.c.mats[[m]][2], m*Ti, m*Ti)
    for (c in 1:C) {
      Sigma.inv[(1 + (c - 1)*m*Ti):(c*m*Ti),
                (1 + (c - 1)*m*Ti):(c*m*Ti)] <- Sigma.inv.c
    }
    cov.beta   <- ginv(t(A)%*%Sigma.inv%*%A)
    cov.tau    <- cov.beta[(Ti + 1):(Ti + D - 1), (Ti + 1):(Ti + D - 1)]
    if (optimality == "D") {
      if (D > 2) {
        obj.fn <- det(cov.tau)
      } else {
        obj.fn <- cov.tau
      }
    } else if (optimality == "A") {
      if (D > 2) {
        obj.fn <- mean(diag(cov.tau))
      } else {
        obj.fn <- cov.tau
      }
    } else if (optimality == "E") {
      if (D > 2) {
        obj.fn <- max(diag(cov.tau))
      } else {
        obj.fn <- cov.tau
      }
    }
    return(obj.fn)
  }
  
  ##### Main Computations ######################################################
  
  if (exhaustive) {
    if (summary) {
      message("Determining X_{C,T} matrices for each C,T combination...")
    }
    if (is.null(poss.mCT)) {
      poss.mCT                              <- expand.grid(set.T, set.C, set.m)
      poss.switching.matrices               <- list()
      for (i in 1:length(set.T)){
        poss.switching.matrices[[set.T[i]]] <- list()
      }
      poss.CT                               <- expand.grid(set.T, set.C)
    } else {
      poss.mCT            <- as.matrix(poss.mCT)
      colnames(poss.mCT)  <- NULL
      set.T               <- unlist(unique(poss.mCT[, 1]))
      poss.CT             <- matrix(poss.mCT[1, 1:2], nrow = 1, ncol = 2)
      if (nrow(poss.mCT) > 1) {
        for (i in 2:nrow(poss.mCT)) {
          local.CT        <- unlist(poss.mCT[i, 1:2])
          check           <- 0
          for (j in 1:nrow(poss.CT)) {
            if (sum(poss.CT[j, ] == local.CT) == 2) {
              check       <- 1
              break
            }
          }
          if (check == 0) {
            poss.CT       <- rbind(poss.CT, local.CT)
          }
        }
        rownames(poss.CT) <- NULL
      }
      poss.switching.matrices               <- list()
      for (i in 1:length(set.T)){
        poss.switching.matrices[[set.T[i]]] <- list()
      }
    }
    num.poss.switching.matrices <- numeric(nrow(poss.mCT))
    for (i in 1:nrow(poss.CT)) {
      Ti                   <- as.numeric(poss.CT[i, 1])
      C                    <- as.numeric(poss.CT[i, 2])
      all.rows             <- getall(iterpc(D, Ti, 0:(D - 1), replace = T))
      if (set.X == "sw") {
        keep <- rep(T, nrow(all.rows))
        for (j in 1:nrow(all.rows)) {
          if (length(unique(all.rows[j, ])) < D) {
            keep[j] <- F
          }
        }
        all.rows <- all.rows[keep, ]
      }
      
      all.combs            <- getall(iterpc(nrow(all.rows), C, replace = T))
      if (set.X == "equal") {
        keep               <- rep(T, nrow(all.combs))
        for (j in 1:nrow(all.combs)) {
          if (length(unique(table(all.combs[j, ]))) != 1) {
            keep[j]        <- F
          }
        }
        all.combs          <- all.combs[keep, ]
      }
      suppressMessages(sfInit(parallel = parallel, cpus = cpus))
      suppressMessages(sfLibrary(Matrix))
      sfExport("all.rows", "all.combs", "C", "Ti", "D")
      results <- sfLapply(1:nrow(all.combs), check.rank)
      suppressMessages(sfStop())
      results <- unlist(results)
      all.combs <- all.combs[results, , drop = F]
      switching.matrices.i <- list()
      if (nrow(all.combs) > 0) {
        for (j in 1:nrow(all.combs)) {
          switching.matrices.i[[j]] <- all.rows[all.combs[j, ], ]
        }
        poss.switching.matrices[[Ti]][[C]] <- switching.matrices.i
        Ti.C.pairs                         <- which(poss.mCT[, 1] == Ti &
                                                      poss.mCT[, 2] == C)
        for (j in Ti.C.pairs) {
          num.poss.switching.matrices[j]   <- nrow(all.combs)
        }
      }
      if (summary) {
        message("Identified designs for C = ", poss.CT[i, 2], ", T = ", poss.CT[i, 1])
      }
    }
    if (summary) {
      message("Determining complete set of possible designs...")
    }
    names(poss.mCT)              <- NULL
    poss.mCT                     <- as.matrix(poss.mCT)
    poss.designs                 <- matrix(0,
                                           nrow =
                                             sum(num.poss.switching.matrices),
                                           ncol = 11 + D - 1 + 4)
    count                        <- 1
    for (i in 1:nrow(poss.mCT)) {
      if (num.poss.switching.matrices[i] > 0) {
        for (j in 1:num.poss.switching.matrices[i]) {
          poss.designs[count, 1:4] <- unlist(c(poss.mCT[i, ], j))
          count                    <- count + 1
        }
      }
    }
    if (summary) {
      message("Determined number of possible designs is ", count - 1, "...")
    }
    if (summary) {
      message("Storing parameters for use later...")
    }
    Sigma.inv.c.mats                  <- list()
    if (D > 2) {
      for (i in 1:length(set.T)) {
        Ti                            <- set.T[i]
        Sigma.inv.c.mats[[Ti]]        <- list()
        poss.m.for.T                  <- unique(poss.mCT[which(poss.mCT[, 1] ==
                                                                 Ti), 3])
        for (j in 1:length(poss.m.for.T)) {
          m                           <- as.numeric(poss.m.for.T[j])
          Sigma.c                     <- matrix(sigma.c^2, m*Ti, m*Ti) +
                                           diag(sigma.eps^2, m*Ti, m*Ti)
          Sigma.c.inv                 <- ginv(Sigma.c)
          Sigma.inv.c.mats[[Ti]][[m]] <- c(Sigma.c.inv[1, 1], Sigma.c.inv[1, 2])
        }
      }
    }
    if (summary) {
      message("Determining performance of each possible design...")
    }
    sink("NULL")
    suppressMessages(sfInit(parallel = parallel, cpus = cpus))
    suppressMessages(sfLibrary(mvtnorm))
    suppressMessages(sfLibrary(MASS))
    suppressMessages(sfLibrary(stats))
    sfExport("D", "poss.mCT", "poss.designs", "poss.switching.matrices",
             "Sigma.inv.c.mats", "alpha", "beta", "delta", "mcc", "sigma.c",
             "sigma.pi", "sigma.s", "sigma.eps", "design")
    results <- sfLapply(1:nrow(poss.designs), ind.design)
    suppressMessages(sfStop())
    sink()
    results <- matrix(unlist(results), ncol = 6 + D, byrow = T)
    if (summary) {
      message("Determining feasible and optimal designs...")
    }
    poss.designs[, 5:(12 + D - 2)]     <- results
    for (i in 1:4) {
      if (max(poss.designs[, 6 + D + i]) > min(poss.designs[, 6 + D + i])) {
        poss.designs[, 12 + D - 2 + i] <- (poss.designs[, 6 + D + i] -
                                             min(poss.designs[, 6 + D + i]))/
                                             (max(poss.designs[, 6 + D + i]) -
                                                min(poss.designs[, 6 + D + i]))
      } else {
        poss.designs[, 12 + D - 2 + i] <- 0
      }
    }
    colnames(poss.designs) <- c("T", "C", "m", "{C,T} Design Index", "e",
                                "type.I", "power.any",
                                paste("power", 1:(D - 1), sep = ""),
                                "obs", "opt.D", "opt.A", "opt.E", "scaled.obs",
                                "scaled.opt.D", "scaled.opt.A", "scaled.opt.E")
    if (power.for == "individual") {
      feasible.designs            <- rep(FALSE, nrow(poss.designs))
      for (i in 1:nrow(poss.designs)) {
        if (all(poss.designs[i, 8:(8 + D - 2)] >= 1 - beta)) {
          feasible.designs[i]        <- TRUE
        }
      }
      feasible.designs            <- poss.designs[feasible.designs, ]
    } else if (power.for == "combined") {
      feasible.designs            <- poss.designs[which(poss.designs[, 7] >=
                                                          1 - beta), ]
    }
    if (nrow(feasible.designs) > 0) {
      optimal.designs             <- list()
      for (i in 1:3) {
        optimal.designs[[i]]      <- list()
      }
      for (i in 1:length(w)) {
        obj.fn.i.D      <- w[i]*feasible.designs[, 12 + D - 2 + 1] +
                             (1 - w[i])*feasible.designs[, 12 + D - 2 + 1 + 1]
        obj.fn.i.A      <- w[i]*feasible.designs[, 12 + D - 2 + 1] +
                             (1 - w[i])*feasible.designs[, 12 + D - 2 + 1 + 2]
        obj.fn.i.E      <- w[i]*feasible.designs[, 12 + D - 2 + 1] +
                             (1 - w[i])*feasible.designs[, 12 + D - 2 + 1 + 3]
        opt.designs.i.D <- which(obj.fn.i.D == min(obj.fn.i.D))
        optimal.designs[[1]][[i]] <- feasible.designs[opt.designs.i.D, ]
        opt.designs.i.A <- which(obj.fn.i.A == min(obj.fn.i.A))
        optimal.designs[[2]][[i]] <- feasible.designs[opt.designs.i.A, ]
        opt.designs.i.E <- which(obj.fn.i.E == min(obj.fn.i.E))
        optimal.designs[[3]][[i]] <- feasible.designs[opt.designs.i.E, ]
      }
    } else {
      feasible.designs <- NULL
      optimal.designs  <- NULL
    }
    return(list(poss.designs = poss.designs,
                feasible.designs = feasible.designs,
                optimal.designs = optimal.designs,
                switching.matrices = poss.switching.matrices))
  } else {
    if (summary) {
      message("Determining possible treatment switching matrices...")
    }
    Ti           <- set.T
    C            <- set.C
    all.switches <- getall(iterpc(D, Ti, 0:(D - 1), replace = TRUE))
    if (summary) {
      message("Storing parameters for use later...")
    }
    Sigma.inv.c.mats        <- list()
    for (j in 1:length(set.m)){
      m                     <- set.m[j]
      Sigma.c               <- matrix(sigma.c^2, m*Ti, m*Ti) +
                                 diag(sigma.eps^2, m*Ti, m*Ti)
      Sigma.c.inv           <- ginv(Sigma.c)
      Sigma.inv.c.mats[[m]] <- c(Sigma.c.inv[1, 1], Sigma.c.inv[1, 2])
    }
    if (summary) {
      message("Performing stochastic search for optimal design...")
    }
    if (length(set.m) > 1) {
      discCat        <- as.integer(c(length(set.m) - 1,
                                     rep(nrow(all.switches) - 1, C)))
      optimal.design <- CEoptim(f = ind.design.stoch.with.n,
                                f.arg = list(D = D, Ti = Ti, C = C,
                                             set.m = set.m,
                                             all.switches = all.switches,
                                             Sigma.inv.c.mats =
                                               Sigma.inv.c.mats,
                                             optimality = optimality),
                                discrete = list(categories = discCat,
                                                smoothProb = 0.5), N = 1000L,
                                rho = 0.01, verbose = T)
      score          <- optimal.design$optimum
      m              <- optimal.design$optimizer$discrete[1]
      X              <-
        all.switches[optimal.design$optimizer$discrete[2:(C + 1)], ]
      return(list(m = m, X = X, score = score))
    } else {
      m              <- set.m
      discCat        <- as.integer(rep(nrow(all.switches) - 1, C))
      optimal.design <- CEoptim(f = ind.design.stoch.without.n,
                                f.arg = list(D = D, Ti = Ti, C = C, m = m,
                                             all.switches = all.switches,
                                             Sigma.inv.c.mats =
                                               Sigma.inv.c.mats,
                                             optimality = optimality),
                                discrete = list(categories = discCat,
                                                smoothProb = smoothProb), N = N,
                                rho = rho, verbose = TRUE,
                                noImproveThr = noImproveThr)
      score           <- optimal.design$optimum
      X               <- all.switches[optimal.design$optimizer$discrete + 1, ]
      return(list(X = X, score = score))
    }
  }
}

##### ind.design.eval(): Short function used to evaluate performance of a ######
##### particular design ########################################################

ind.design.eval <- function(X, m, C, Ti, sigma.c, sigma.eps, sigma.s, sigma.pi,
                            design, Sigma.inv.c.mats, alpha, mcc, delta) {
  D <- length(unique(as.vector(X)))
  if (D == 2) {
    U <- sum(X)
    V <- sum(rowSums(X)^2)
    W <- sum(colSums(X)^2)
    if (design == "cohort") {
      sigma2  <- sigma.c^2 + sigma.pi^2 + sigma.s^2 + sigma.eps^2
      rho0    <- (sigma.c^2 + sigma.pi^2)/sigma2
      rho1    <- sigma.c^2/sigma2
      rho2    <- (sigma.c^2 + sigma.s^2)/sigma2
      psi     <- 1 + (m - 1)*rho0 - (m - 1)*rho1 - rho2
      xi      <- (m - 1)*rho1 + rho2
      gamma   <- psi + Ti*xi
      cov.tau <- (C*sigma2*gamma*psi/m)/((C*U - W)*gamma + (U^2 - C*V)*xi)
    } else {
      sigma2  <- sigma.eps^2/m
      cov.tau <- (C*sigma2*(sigma2 + Ti*sigma.c^2))/
        ((C*U - W)*sigma2 + (U^2 + C*Ti*U - Ti*W - C*V)*sigma.c^2)
    }
  } else {
    Sigma.inv   <- matrix(0, nrow = m*C*Ti, ncol = m*C*Ti)
    Sigma.inv.c <- matrix(Sigma.inv.c.mats[[Ti]][[m]][2], m*Ti, m*Ti) +
      diag(Sigma.inv.c.mats[[Ti]][[m]][1] -
             Sigma.inv.c.mats[[Ti]][[m]][2], m*Ti, m*Ti)
    for (c in 1:C) {
      Sigma.inv[(1 + (c - 1)*m*Ti):(c*m*Ti),
                (1 + (c - 1)*m*Ti):(c*m*Ti)] <- Sigma.inv.c
    }
    A      <- matrix(0, nrow = m*C*Ti, ncol = 1 + (Ti - 1) + (D - 1))
    A[, 1] <- 1
    for (i in 1:C) {
      for (j in 1:Ti) {
        if (X[i, j] != 0) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m),
            (Ti + 1):(Ti + X[i, j])]                                <- 1
        }
        if (j > 1) {
          A[(1 + (i - 1)*m*Ti + (j - 1)*m):((i - 1)*m*Ti + j*m), j] <- 1
        }
      }
    }
    cov.beta <- ginv(t(A)%*%Sigma.inv%*%A)
    cov.tau  <- cov.beta[(Ti + 1):(Ti + D - 1), (Ti + 1):(Ti + D - 1)]
  }
  if (mcc == "none") {
    e <- qnorm(1 - alpha)
  } else if (mcc == "bonf") {
    e <- qnorm(1 - alpha/(D - 1))
  } else if (mcc == "dunnett") {
    e <- qmvnorm(1 - alpha, sigma = cov.tau)$quantile
  }
  obs     <- m*C*Ti
  if (D > 2) {
    opt.D <- det(cov.tau)
    opt.A <- mean(diag(cov.tau))
    opt.E <- max(diag(cov.tau))
  } else {
    opt.D <- opt.A <- opt.E <- cov.tau
  }
  if (D > 2) {
    cov.z    <- diag(1/sqrt(diag(cov.tau)))%*%cov.tau%*%
      diag(1/sqrt(diag(cov.tau)))
  } else {
    cov.z    <- 1
  }
  type.I     <- 1 - pmvnorm(lower = rep(-Inf, D - 1), upper = rep(e, D - 1),
                            sigma = cov.z, mean = rep(0, D - 1))[1]
  power      <- numeric(1 + (D - 1))
  if (D > 2) {
    power[1] <- 1 - pmvnorm(lower = rep(-Inf, D - 1), upper = rep(e, D - 1),
                            sigma = cov.z,
                            mean = delta/sqrt(diag(cov.tau)))[1]
    for (d in 1:(D - 1)) {
      power[d + 1] <- pnorm(e, mean = delta[d]/sqrt(diag(cov.tau)[d]),
                            lower.tail = FALSE)
    }
  } else {
    power    <- rep(pnorm(e, mean = delta/sqrt(cov.tau),
                          lower.tail = FALSE), 2)
  }
  return(c(e, type.I, power, obs, opt.D, opt.A, opt.E))
}

##### 3.1. D = 2: Girling and Hemming (2016) ###################################

sigma            <- sqrt(1)
E                <- c(0.1, 0.15, 0.3, 0.45, 0.75, 0.9)
rho1             <- E/(60 - 59*E)
optimal.X        <- list()
for (i in 1:length(rho1)) {
  sigma.c        <- sqrt(rho1[i]*sigma^2)
  sigma.eps      <- sqrt(sigma^2 - sigma.c^2)
  optimal.design <- admissible.ma.sw.crt(D = 2, set.T = 6, set.C = 10,
                                         set.m = 10, sigma.eps = sigma.eps,
                                         sigma.c = sigma.c, beta = 1)
  index          <- optimal.design$optimal.designs[[1]][[1]][4]
  optimal.X[[i]] <- optimal.design$switching.matrices[[6]][[10]][[index]]
}

##### 3.1. D = 2: Thompson et al. (2017) #######################################

optimal.X.t        <- list()
for (i in 1:length(rho1)) {
  sigma.c        <- sqrt(rho1[i]*sigma^2)
  sigma.eps      <- sqrt(sigma^2 - sigma.c^2)
  optimal.design.t <- admissible.ma.sw.crt(D = 2, set.T = 6, set.C = 10,
                                           set.m = 10, sigma.eps = sigma.eps,
                                           sigma.c = sigma.c, beta = 1,
                                           set.X = "equal")
  if (is.null(dim(optimal.design.t$optimal.designs[[1]][[1]]))) {
    index                 <- optimal.design.t$optimal.designs[[1]][[1]][4]
    optimal.X.t[[i]]        <-
      optimal.design.t$switching.matrices[[6]][[10]][[index]]
  } else {
    optimal.X.t[[i]]        <- list()
    for (j in 1:nrow(optimal.design.t$optimal.designs[[1]][[1]])) {
      index               <- optimal.design.t$optimal.designs[[1]][[1]][j, 4]
      optimal.X.t[[i]][[j]] <- 
        optimal.design.t$switching.matrices[[6]][[10]][[index]]
    }
  }
}

##### 3.2. D = 2: Sensitivity of optimal designs to variance parameters ########

sigma.c2         <- seq(from = 0.001, to = 0.25, length.out = 30)
sigma.eps2       <- seq(from = 0.25, to = 4, length.out = 30)
optimal.X       <- list()
optimal.indices <- matrix(0, nrow = length(sigma.c2), ncol = length(sigma.eps2))
for (i in 1:length(sigma.c2)) {
  optimal.X[[i]]   <- list()
  for (j in 1:length(sigma.eps2)) {
    optimal.design <- admissible.ma.sw.crt(D = 2, set.T = 6, set.C = 10,
                                           set.m = 10,
                                           sigma.eps = sqrt(sigma.eps2[j]),
                                           sigma.c = sqrt(sigma.c2[i]),
                                           beta = 1)
    if (!is.null(dim(optimal.design$optimal.designs[[1]][[1]]))) {
      props      <- numeric(nrow(optimal.design$optimal.designs[[1]][[1]]))
      for (k in 1:length(props)) {
        index    <- optimal.design$optimal.designs[[1]][[1]][k, 4]
        props[k] <- sum(optimal.design$switching.matrices[[6]][[10]][[index]])
      }
      optimal.indices[i, j] <-
        optimal.design$optimal.designs[[1]][[1]][which(props == max(props)), 4]
      optimal.X[[i]][[j]]   <-
        optimal.design$switching.matrices[[6]][[10]][[optimal.indices[i, j]]]
    } else {
      optimal.indices[i, j] <- optimal.design$optimal.designs[[1]][[1]][4]
      optimal.X[[i]][[j]]   <-
        optimal.design$switching.matrices[[6]][[10]][[optimal.indices[i, j]]]
    }
  }
  message("You have completed i = ", i)
}
sigmas           <- expand.grid(sigma.c2, sigma.eps2)
vector.indices   <- NULL
for (j in 1:length(sigma.eps2)) {
  vector.indices <- c(vector.indices, optimal.indices[, j])
}
data        <- data.frame(sigma.c2 = sigmas[, 1], sigma.eps2 = sigmas[, 2],
                          Design = factor(vector.indices,
                                          unique(vector.indices)))
data$Design <- mapvalues(data$Design, from = unique(vector.indices),
                         to = paste("Design", 1:length(unique(vector.indices))))
plot.opt    <- ggplot(data = data, aes(x = sigma.c2, y = sigma.eps2,
                                       colour = Design)) +
  geom_point() +
  scale_color_manual(values = c("red", "lightseagreen", "blue", "cyan",
                                "purple", "orange", "black",
                                "grey", "yellow", "hotpink",
                                "brown")) +
  xlab(expression(sigma[c]^2)) +
  ylab(expression(sigma[epsilon]^2)) +
  theme(legend.title = element_blank())
ggsave("Optimal.pdf", plot = plot.opt, units = "in", height = 5, width = 6)
indices <- c(461, 3948)
for (design in 1:length(indices)) {
  og.X  <- optimal.design$switching.matrices[[6]][[10]][[indices[design]]]
  relative.performance <- matrix(0, nrow = length(sigma.c2)*length(sigma.eps2),
                                 ncol = 3)
  for (i in 1:length(sigma.c2)) {
    for (j in 1:length(sigma.eps2)) {
      rel <- ind.design.eval(og.X, 10, 10, 6, sqrt(sigma.c2[i]),
                             sqrt(sigma.eps2[j]),
                             0, 0, "cross_sectional",
                             list(), 0.05, "none", 0)[6]/
        ind.design.eval(optimal.X[[i]][[j]], 10, 10, 6, sqrt(sigma.c2[i]),
                        sqrt(sigma.eps2[j]), 0, 0, "cross_sectional",
                        list(), 0.05, "none", 0)[6]
      relative.performance[j + (i - 1)*30, ] <- c(sigma.c2[i], sigma.eps2[j],
                                                  rel)
    }
  }
  relative.performance <- data.frame(sigma.c2 = relative.performance[, 1],
                                     sigma.eps2 = relative.performance[, 2],
                                     rel = relative.performance[, 3])
  plot.rel <- ggplot(data = relative.performance,
                     aes(x = sigma.c2, y = sigma.eps2, fill = rel)) +
    geom_tile() + scale_fill_distiller(palette = "Spectral") +
    xlab(expression(sigma[c]^2)) +
    ylab(expression(sigma[epsilon]^2)) +
    theme(legend.title = element_blank())
  ggsave(paste("Design", design, ".pdf", sep = ""), plot = plot.rel,
         units = "in", height = 5, width = 5.5)
}

##### 3.3. D = 2: Fan et al. (2018) ############################################

sigma            <- sqrt(1)
rho0             <- c(0.05, 0.1)
rho1             <- c(0.001, 0.002)
rho2             <- c(0.25, 0.5)
rho              <- as.matrix(expand.grid(rho0, rho1, rho2))
optimal.X        <- list()
for (i in 1:nrow(rho)) {
  sigma.c   <- sqrt(rho[i, 2])
  sigma.s   <- sqrt(rho[i, 3] - sigma.c^2)
  sigma.pi  <- sqrt(rho[i, 1] - sigma.c^2)
  sigma.eps <- sqrt(1 - sigma.c^2 - sigma.s^2 - sigma.pi^2)
  optimal.design <- admissible.ma.sw.crt(D = 2, set.T = 6, set.C = 10,
                                         set.m = 10, sigma.eps = sigma.eps,
                                         sigma.c = sigma.c, sigma.s = sigma.s,
                                         sigma.pi = sigma.pi, beta = 1,
                                         design = "cohort",
                                         set.X = "sw")
  if (is.null(dim(optimal.design$optimal.designs[[1]][[1]]))) {
    index                 <- optimal.design$optimal.designs[[1]][[1]][4]
    optimal.X[[i]]        <-
      optimal.design$switching.matrices[[6]][[10]][[index]]
  } else {
    optimal.X[[i]]        <- list()
    for (j in 1:nrow(optimal.design$optimal.designs[[1]][[1]])) {
      index               <- optimal.design$optimal.designs[[1]][[1]][j, 4]
      optimal.X[[i]][[j]] <- 
        optimal.design$switching.matrices[[6]][[10]][[index]]
    }
  }
  message("You have completed i = ", i)
}
X.cohort      <- optimal.X
m             <- 10
Ti            <- 6
p.theoretical <- matrix(0, nrow = 8, ncol = 5)
p.empirical   <- matrix(0, nrow = 8, ncol = 5)
for (i in 1:nrow(rho)) {
  rho0 <- rho[i, 1]; rho1 <- rho[i, 2]; rho2 <- rho[i, 3]
  psi     <- 1 + (m - 1)*rho0 - (m - 1)*rho1 - rho2
  xi      <- (m - 1)*rho1 + rho2
  gamma   <- psi + Ti*xi
  p.theoretical[i, c(1, 5)] <- (psi + 3*xi)/(2*gamma)
  p.theoretical[i, 2:4]     <- xi/gamma
  if (i %in% c(1,3)) {
    X <-X.cohort[[i]][[1]]
  } else {
    X <- X.cohort[[i]]
  }
  for (c in 1:10) {
    p.empirical[i, sum(X[c, ])] <- p.empirical[i, sum(X[c, ])] + 1
  }
}
p.empirical <- p.empirical/10

##### 3.4. D = 3: SO-HIP Study: Unconstrained ##################################

# Note: Utilised design is:
X <- matrix(c(0, 0, 0, 0, 0,0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,
              0, 0, 1, 2, 2, 2), ncol = 6, nrow = 6, byrow = T)
sigma     <- sqrt(1)
rho       <- 0.05
sigma.c   <- sqrt(rho*sigma^2)
sigma.eps <- sqrt(sigma^2 - sigma.c^2)
Sigma.inv.c.mats            <- list()
Ti                          <- 6
Sigma.inv.c.mats[[Ti]]      <- list()
m                           <- 8
Sigma.c                     <- matrix(sigma.c^2, m*Ti, m*Ti) +
                                 diag(sigma.eps^2, m*Ti, m*Ti)
Sigma.c.inv                 <- ginv(Sigma.c)
Sigma.inv.c.mats[[Ti]][[m]] <- c(Sigma.c.inv[1, 1], Sigma.c.inv[1, 2])
utilised.design             <- ind.design.eval(X, m, 6, Ti, sigma.c, sigma.eps,
                                               0, 0, "cross_sectional",
                                               Sigma.inv.c.mats, 0.05, "bonf",
                                               c(1.5*sigma,0.75*sigma))
C            <- c(2, 3, 4, 5, 6)
Ti           <- c(2, 3, 4, 5, 6)
poss.CT      <- expand.grid(Ti, C)
poss.mCT     <- NULL
for (i in 1:nrow(poss.CT)){
  poss.m     <- 2:floor(48/poss.CT[i, 1])
  for (j in 1:length(poss.m)){
    poss.mCT <- rbind(poss.mCT, c(poss.CT[i, ], poss.m[j]))
  }
}
sigma           <- sqrt(1)
rho             <- 0.05
sigma.c         <- sqrt(rho*sigma^2)
sigma.eps       <- sqrt(sigma^2 - sigma.c^2)
optimal.designs <- admissible.ma.sw.crt(poss.mCT = poss.mCT, sigma.c = sigma.c,
                                        sigma.eps = sigma.eps,
                                        power.for = "individual", beta = 0.17,
                                        delta = c(1.5*sigma, 0.75*sigma),
                                        mcc = "bonf", w = c(0, 0.5, 1 - 10^-4),
                                        cpus = 6, D = 3)

##### 3.4. D = 3: SO-HIP Study: Constrained ####################################

C            <- c(2, 3, 4, 5, 6)
Ti           <- c(4, 5, 6)
poss.CT      <- expand.grid(Ti, C)
poss.mCT        <- NULL
for (i in 1:nrow(poss.CT)){
  poss.m        <- 2:floor(48/poss.CT[i, 1])
  for (j in 1:length(poss.m)){
    poss.mCT    <- rbind(poss.mCT, c(poss.CT[i, ], poss.m[j]))
  }
}
sigma           <- sqrt(1)
rho             <- 0.05
sigma.c         <- sqrt(rho*sigma^2)
sigma.eps       <- sqrt(sigma^2 - sigma.c^2)
optimal.designs <- admissible.ma.sw.crt(poss.mCT = poss.mCT, sigma.c = sigma.c,
                                        sigma.eps = sigma.eps,
                                        power.for = "individual", beta = 0.17,
                                        delta = c(1.5*sigma, 0.75*sigma),
                                        mcc = "bonf", w = c(0, 0.5, 1 - 10^-4),
                                        cpus = 6, D = 3, set.X = "sw")

##### 3.5. D = 3: Optimal cross-sectional designs according to the value of ####
##### the cluster mean correlation #############################################
sigma             <- sqrt(1)
rho               <- seq(from = 0, to = 1, by = 0.01)
optimal.X.D       <- list()
optimal.X.A       <- list()
optimal.X.E       <- list()
for (i in 1:length(rho)) {
  sigma.c         <- sqrt(rho[i]*sigma^2)
  sigma.eps       <- sqrt(sigma^2 - sigma.c^2)
  optimal.designs <- admissible.ma.sw.crt(set.C = 6, set.T = 6, set.m = 8,
                                          sigma.c = sigma.c,
                                          sigma.eps = sigma.eps,
                                          mcc = "none", w = 0, beta = 1,
                                          cpus = 6)
  if (is.null(dim(optimal.design$optimal.designs[[1]][[1]]))) {
    index            <- optimal.design$optimal.designs[[1]][[1]][4]
    optimal.X.D[[i]] <- optimal.design$switching.matrices[[6]][[10]][[index]]
  } else {
    optimal.X.D[[i]] <- list()
    for (j in 1:nrow(optimal.design$optimal.designs[[1]][[1]])) {
      index                 <- optimal.design$optimal.designs[[1]][[1]][j, 4]
      optimal.X.D[[i]][[j]] <-
        optimal.design$switching.matrices[[6]][[10]][[index]]
    }
  }
  if (is.null(dim(optimal.design$optimal.designs[[2]][[1]]))) {
    index            <- optimal.design$optimal.designs[[2]][[1]][4]
    optimal.X.A[[i]] <- optimal.design$switching.matrices[[6]][[10]][[index]]
  } else {
    optimal.X.A[[i]] <- list()
    for (j in 1:nrow(optimal.design$optimal.designs[[2]][[1]])) {
      index                 <- optimal.design$optimal.designs[[2]][[1]][j, 4]
      optimal.X.A[[i]][[j]] <-
        optimal.design$switching.matrices[[6]][[10]][[index]]
    }
  }
  if (is.null(dim(optimal.design$optimal.designs[[3]][[1]]))) {
    index            <- optimal.design$optimal.designs[[3]][[1]][4]
    optimal.X.E[[i]] <- optimal.design$switching.matrices[[6]][[10]][[index]]
  } else {
    optimal.X.E[[i]] <- list()
    for (j in 1:nrow(optimal.design$optimal.designs[[3]][[1]])) {
      index                 <- optimal.design$optimal.designs[[3]][[1]][j, 4]
      optimal.X.E[[i]][[j]] <-
        optimal.design$switching.matrices[[6]][[10]][[index]]
    }
  }
}

##### 3.6. D = 4 ###############################################################

# Note: Assumed utilised design is:
X <- matrix(c(0, 1, 1, 2, 2, 3, 3, 3,
              0, 1, 1, 2, 2, 3, 3, 3,
              0, 0, 1, 1, 2, 2, 3, 3,
              0, 0, 1, 1, 2, 2, 3, 3,
              0, 0, 0, 1, 1, 2, 2, 3,
              0, 0, 0, 1, 1, 2, 2, 3), ncol = 8, nrow = 6, byrow = T)
sigma            <- sqrt(1)
rho              <- 0.05
sigma.c          <- sqrt(rho*sigma^2)
sigma.eps        <- sqrt(sigma^2 - sigma.c^2)
Sigma.inv.c.mats <- list()
Ti               <- 8
Sigma.inv.c.mats[[Ti]]      <- list()
m                           <- 8
Sigma.c                     <- matrix(sigma.c^2, m*Ti, m*Ti) + diag(sigma.eps^2,
                                                                    m*Ti, m*Ti)
Sigma.c.inv                 <- ginv(Sigma.c)
Sigma.inv.c.mats[[Ti]][[m]] <- c(Sigma.c.inv[1, 1], Sigma.c.inv[1, 2])
utilised.design             <-
  ind.design.eval(X, m, 6, Ti, sigma.c, sigma.eps, 0, 0, "cross_sectional",
                  Sigma.inv.c.mats, 0.05, "bonf",
                  c(1.5*sigma,0.75*sigma, 0.75*sigma))
set.seed(106)
D.optimal <- admissible.ma.sw.crt(D = 4, set.T = 8, set.C = 6, set.m = 8,
                                  sigma.c = sigma.c, sigma.eps = sigma.eps,
                                  delta = c(1.5*sigma,0.75*sigma, 0.75*sigma),
                                  optimality = "D", exhaustive = FALSE)
set.seed(406)
A.optimal <- admissible.ma.sw.crt(D = 4, set.T = 8, set.C = 6, set.m = 8,
                                  sigma.c = sigma.c, sigma.eps = sigma.eps,
                                  delta = c(1.5*sigma,0.75*sigma, 0.75*sigma),
                                  optimality = "A", exhaustive = FALSE)
set.seed(506)
E.optimal <- admissible.ma.sw.crt(D = 4, set.T = 8, set.C = 6, set.m = 8,
                                  sigma.c = sigma.c, sigma.eps = sigma.eps,
                                  delta = c(1.5*sigma,0.75*sigma, 0.75*sigma),
                                  optimality = "E", exhaustive = FALSE)
