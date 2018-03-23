################################################################################
# Function: abcd.mams()                                                        #
# Last Modified: 23/03/2018                                                    #
################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
#        K - Initial number of experimental treatment arms                     #
#        J - Maximal number of stages                                          #
#        a - Control the a-generalised type-I FWER                             #
#        b - See c below                                                       #
#        c - Power the trial to reject b out of c false null hypotheses        #
#        d - Stop the trial once d or more null hypotheses have been rejected  #
#    alpha - Control type-I FWER to alpha                                      #
#     beta - Control type-II FWER to beta                                      #
#    ratio - Allocation ratio of experimental arms to the control              #
#  delta.1 - Desired interesting treatment effect                              #
#  delta.0 - Uninteresting treatment effect                                    #
#        w - Weights used in the definition of the optimality criteria         #  
# parallel - Should parallelisation be used                                    #
#     cpus - Number of cpus to parallelise over                                #
#  summary - Print a summary of progress                                       #
#  popSize - See ?ga                                                           #
#  maxiter - See ?ga                                                           #
#      run - See ?ga                                                           #
#     seed - Number seed for reproducibility                                   #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# A list principally containing details of the parameters of the identified    #
# optimal design, and its operating characteristics                            #
################################################################################

abcd.mams <- function(K = 3, J = 2, a = 1, b = 1, c = 1, d = 1, alpha = 0.05,
                      beta = 0.10, sigma = 1, ratio = 1, delta.1 = 0.545,
                      delta.0 = 0.178, w = rep(1/3, 3),
                      parallel = TRUE, cpus = 8, summary = TRUE,
                      popSize = 100, maxiter = 500, run = maxiter,
                      seed = NULL) {
  if (is.null(seed)) {
    seed <- floor(runif(1, 1, 10^6))
    set.seed(seed)
  } else {
    set.seed(seed)
  }
  
  ##### Load Required Packages #################################################
  
  library(doParallel)
  library(doRNG)
  library(GA)
  library(iterpc)
  library(mvtnorm)
  library(parallel)
  
  ##### Input Checking #########################################################
  
  if (any(K < 1, K%%1 != 0)){
    stop("K must be an integer greater than or equal to 1.")
  }
  if (any(J < 1, J%%1 != 0)){
    stop("J must be an integer greater than or equal to 1.")
  }
  if (any(a < 1, a%%1 != 0, a > K)){
    stop("a must be an integer in [1,K].")
  }
  if (any(b < 1, b%%1 != 0, b > K)){
    stop("b must be an integer in [1,K].")
  }
  if (any(c < b, c%%1 != 0, c > K)){
    stop("c must be an integer in [b,K].")
  }
  if (any(d < 1, d%%1 != 0, d > K)){
    stop("d must be an integer in [1,K].")
  }
  if (any(alpha <= 0, alpha >= 1)){
    stop("alpha must be strictly between 0 and 1.")
  }
  if (any(beta <= 0, beta >= 1)){
    stop("beta must be strictly between 0 and 1.")
  }
  if (any(!(length(sigma) %in% c(1, 2)), any(sigma <= 0))){
    stop("sigma must be a vector of length 1 or 2 with elements greater than zero.")
  }
  if (length(sigma) == 1){
    sigma <- rep(sigma, 2)
  }
  if (ratio <= 0){
    stop("ratio must be strictly positive.")
  }
  r <- ratio
  if (delta.1 <= 0){
    stop("delta.1 must be strictly positive.")
  }
  if (delta.0 >= delta.1){
    stop("delta.0 must be less than delta.1.")
  }
  if (any(length(w) != 3, any(w < 0), w[1] + w[2] == 0)){
    stop("w must be a vector of length 3 with elements greater than or equal to 0, and with w[1] + w[2] > 0.")
  }
  if (!is.logical(parallel)){
    stop("parallel must be set to TRUE or FALSE.")
  }
  if (any(cpus < 1, cpus%%1 != 0)){
    stop("cpus must be an integer greater than or equal to 1.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }
  
  ##### Internal Function Initialisation #######################################
  
  # Initialise a function that generates possible trial outcomes and associated
  # objects which speed up computation e.g. covariance matrices,
  # boundary indices
  thetaGeneration <- function(K, J, a, b, c, d, r, sigma, delta.1, delta.0){
    perms.K.HG        <- getall(iterpc(2*J, K, replace = TRUE))
    perms.d.HG           <- perms.K.HG
    num.perms.d.HG     <- nrow(perms.d.HG)
    deg.d.HG      <- numeric(num.perms.d.HG)
    for (i in 1:num.perms.d.HG){
      deg.d.HG[i] <- nrow(getall(iterpc(table(perms.d.HG[i, ]),
                                        ordered = TRUE)))
    }
    perms.d.HG.a.fwer.indices <- rowSums(perms.d.HG%%2 == 0) >= a
    deg.d.HG.a.fwer    <- numeric(num.perms.d.HG)
    deg.d.HG.a.fwer[perms.d.HG.a.fwer.indices] <- deg.d.HG[perms.d.HG.a.fwer.indices]
    thetas.d.HG.a.fwer        <- cbind(1*(perms.d.HG%%2 == 0),
                                       ceiling(perms.d.HG/2),
                                       deg.d.HG,
                                       deg.d.HG.a.fwer,
                                       rowSums(ceiling(perms.d.HG/2)) +
                                         apply(ceiling(perms.d.HG/2), 1, max),
                                       numeric(num.perms.d.HG))
    colnames(thetas.d.HG.a.fwer) <- c(paste("psi_", 1:K, sep = ""),
                                      paste("omega_", 1:K, sep = ""),
                                      "deg_{HG}(psi,omega)",
                                      "deg_{HG,a-FWER}(psi,omega)",
                                      "N(omega)/n",
                                      "P(psi,omega|HG)")
    retain            <- rep(TRUE, nrow(thetas.d.HG.a.fwer))
    for (i in 1:nrow(thetas.d.HG.a.fwer)){
      for (j in 1:J){
        le            <- which(thetas.d.HG.a.fwer[i, (K + 1):(2*K)] <= j)
        g             <- which(thetas.d.HG.a.fwer[i, (K + 1):(2*K)] > j)
        if (all(c(length(le), length(g)) > 0)){
          if (sum(thetas.d.HG.a.fwer[i, le]) >= d){
            retain[i] <- FALSE
          }
        }
      }
    }
    thetas.d.HG.a.fwer  <- thetas.d.HG.a.fwer[retain, ]
    perms.c.delta_1     <- getall(iterpc(2*J, c, replace = TRUE))
    num.perms.c.delta_1 <- nrow(perms.c.delta_1)
    deg.c.delta_1       <- numeric(num.perms.c.delta_1)
    for (i in 1:num.perms.c.delta_1){
      deg.c.delta_1[i]  <- nrow(getall(iterpc(table(perms.c.delta_1[i, ]),
                                              ordered = TRUE)))
    }
    perms.c.delta_1.b.power.indices <- rowSums(perms.c.delta_1%%2 == 0) >= b
    deg.c.delta_1.b.power    <- numeric(num.perms.c.delta_1)
    deg.c.delta_1.b.power[perms.c.delta_1.b.power.indices] <- deg.c.delta_1[perms.c.delta_1.b.power.indices]
    if (c < K){
      perms.Kmin.c.delta_0     <- getall(iterpc(2*J, K - c, replace = TRUE))
      num.perms.Kmin.c.delta_0 <- nrow(perms.Kmin.c.delta_0)
      deg.Kmin.c.delta_0       <- numeric(num.perms.Kmin.c.delta_0)
      for (i in 1:num.perms.Kmin.c.delta_0){
        deg.Kmin.c.delta_0[i]  <- nrow(getall(iterpc(table(perms.Kmin.c.delta_0[i, ]),
                                                     ordered = TRUE)))
      }
      perms.K.d.bc.LFC         <- matrix(0, nrow = num.perms.c.delta_1*
                                         num.perms.Kmin.c.delta_0, ncol = K)
      for (i in 1:num.perms.c.delta_1){
        range    <- (1 + (i - 1)*num.perms.Kmin.c.delta_0):
                      (i*num.perms.Kmin.c.delta_0)
        i.matrix <- matrix(perms.c.delta_1[i, ],
                           nrow = num.perms.Kmin.c.delta_0, ncol = c,
                           byrow = TRUE)
        perms.K.d.bc.LFC[range, ] <- cbind(i.matrix, perms.Kmin.c.delta_0)
      }
      deg.K.d.bc.LFC         <- rep(deg.c.delta_1,
                                    each = num.perms.Kmin.c.delta_0)*
                                  rep(deg.Kmin.c.delta_0, num.perms.c.delta_1)
      deg.K.d.bc.LFC.power   <- rep(deg.c.delta_1.b.power,
                                    each = num.perms.Kmin.c.delta_0)*
                                  rep(deg.Kmin.c.delta_0,
                                      num.perms.c.delta_1)
      thetas.bc.LFC.power    <- cbind(1*(perms.K.d.bc.LFC%%2 == 0),
                                      ceiling(perms.K.d.bc.LFC/2),
                                      deg.K.d.bc.LFC,
                                      deg.K.d.bc.LFC.power,
                                      rowSums(ceiling(perms.K.d.bc.LFC/2)) +
                                        apply(ceiling(perms.K.d.bc.LFC/2), 1,
                                              max),
                                      numeric(num.perms.c.delta_1*
                                                num.perms.Kmin.c.delta_0))
    } else {
      thetas.bc.LFC.power <- cbind(1*(perms.c.delta_1%%2 == 0),
                                   ceiling(perms.c.delta_1/2),
                                   deg.c.delta_1,
                                   deg.c.delta_1.b.power,
                                   rowSums(ceiling(perms.c.delta_1/2)) +
                                     apply(ceiling(perms.c.delta_1/2), 1, max),
                                   numeric(num.perms.c.delta_1))
    }
    colnames(thetas.bc.LFC.power) <- c(paste("psi_", 1:K, sep = ""),
                                       paste("omega_", 1:K, sep = ""),
                                       "deg_{LFC}(psi,omega)",
                                       "deg_{LFC,bc-power}(psi,omega)",
                                       "N(omega)/n",
                                       "P(psi,omega|LFC)")
    retain            <- rep(TRUE, nrow(thetas.bc.LFC.power))
    for (i in 1:nrow(thetas.bc.LFC.power)){
      for (j in 1:J){
        le            <- which(thetas.bc.LFC.power[i, (K + 1):(2*K)] <= j)
        g             <- which(thetas.bc.LFC.power[i, (K + 1):(2*K)] > j)
        if (all(c(length(le), length(g)) > 0)){
          if (sum(thetas.bc.LFC.power[i, le]) >= d){
            retain[i] <- FALSE
          }
        }
      }
    }
    thetas.bc.LFC.power <- thetas.bc.LFC.power[retain, ]
    Lambda                          <- matrix(0, K*J, K*J)
    for (j in 1:J){
      Lambda[(1 + (j - 1)*K):(j*K),
             (1 + (j - 1)*K):(j*K)] <- matrix(r/(r + 1), K, K) +
                                         diag(1 - r/(r + 1), K, K)
    }
    if (J > 1){
      for (j1 in 2:J){
        for (j2 in 1:(j1 - 1)){
          Cov.Zj1Zj2 <- matrix((r/(r + 1))*sqrt(j2/j1), K, K) +
                          diag(sqrt(j2/j1) - (r/(r + 1))*sqrt(j2/j1), K, K)
          Lambda[(1 + (j2 - 1)*K):(j2*K), (1 + (j1 - 1)*K):(j1*K)] <- Cov.Zj1Zj2
          Lambda[(1 + (j1 - 1)*K):(j1*K), (1 + (j2 - 1)*K):(j2*K)] <- Cov.Zj1Zj2
        }
      }
    }
    I.div.n      <- rep((J/(sigma[1]^2 + sigma[2]^2/r))*seq_len(J)/J,
                        each = K)
    means.HG     <- list()
    l.indices.HG <- list()
    u.indices.HG <- list()
    Lambdas.HG   <- list()
    for (i in 1:nrow(thetas.d.HG.a.fwer)){
      relevant.indices   <- NULL
      psi                <- thetas.d.HG.a.fwer[i, 1:K]
      omega              <- thetas.d.HG.a.fwer[i, (K + 1):(2*K)]
      l.i                <- numeric(sum(omega))
      u.i                <- numeric(sum(omega))
      for (k in 1:K){
        relevant.indices <- c(relevant.indices,
                              seq(from = k, by = K,
                                  length.out = omega[k]))
        l.i.k            <- numeric(omega[k])
        u.i.k            <- numeric(omega[k])
        for (j in 1:omega[k]){
          if (omega[k] > j){
            l.i.k[j] <- j
          } else if (psi[k] == 0 & omega[k] == j){
            l.i.k[j] <- 2*J + 1
          } else if (psi[k] == 1 & omega[k] == j){
            l.i.k[j] <- J + j
          }
          if ((psi[k] == 0 & omega[k] == j & max(omega) > j) |
              (psi[k] == 0 & omega[k] == j & max(omega) == j &
               all(psi[omega == j] == 0)) |
              (psi[k] == 0 & omega[k] == j & max(omega) == j &
               sum(psi == 1) < d)){
            u.i.k[j] <- j
          } else if (psi[k] == 1 & omega[k] == j){
            u.i.k[j] <- 2*J + 2
          } else if (omega[k] > j | (psi[k] == 0 & omega[k] == j &
                                     max(omega) == j & sum(psi == 1) >= d)){
            u.i.k[j] <- J + j
          }
        }
        l.i[seq(from = k, by = K, length.out = omega[k])] <- l.i.k
        u.i[seq(from = k, by = K, length.out = omega[k])] <- u.i.k
      }
      means.HG[[i]]     <- numeric(length(relevant.indices))
      Lambdas.HG[[i]]   <- Lambda[sort(relevant.indices),
                                  sort(relevant.indices)]
      l.indices.HG[[i]] <- l.i[sort(relevant.indices)]
      u.indices.HG[[i]] <- u.i[sort(relevant.indices)]
    }
    I.div.n.LFC   <- list()
    l.indices.LFC <- list()
    u.indices.LFC <- list()
    Lambdas.LFC   <- list()
    delta.LFC     <- rep(c(rep(delta.1, c), rep(delta.0, K - c)), J)
    deltas.LFC    <- list()
    for (i in 1:nrow(thetas.bc.LFC.power)){
      relevant.indices   <- NULL
      psi                <- thetas.bc.LFC.power[i, 1:K]
      omega              <- thetas.bc.LFC.power[i, (K + 1):(2*K)]
      l.i                <- numeric(sum(omega))
      u.i                <- numeric(sum(omega))
      for (k in 1:K){
        relevant.indices <- c(relevant.indices,
                              seq(from = k, by = K,
                                  length.out = omega[k]))
        l.i.k            <- numeric(omega[k])
        u.i.k            <- numeric(omega[k])
        for (j in 1:omega[k]){
          if (omega[k] > j){
            l.i.k[j] <- j
          } else if (psi[k] == 0 & omega[k] == j){
            l.i.k[j] <- 2*J + 1
          } else if (psi[k] == 1 & omega[k] == j){
            l.i.k[j] <- J + j
          }
          if ((psi[k] == 0 & omega[k] == j & max(omega) > j) |
              (psi[k] == 0 & omega[k] == j & max(omega) == j &
               all(psi[omega == j] == 0)) |
              (psi[k] == 0 & omega[k] == j & max(omega) == j &
               sum(psi == 1) < d)){
            u.i.k[j] <- j
          } else if (psi[k] == 1 & omega[k] == j){
            u.i.k[j] <- 2*J + 2
          } else if (omega[k] > j | (psi[k] == 0 & omega[k] == j &
                                     max(omega) == j & sum(psi == 1) >= d)){
            u.i.k[j] <- J + j
          }
        }
        l.i[seq(from = k, by = K, length.out = omega[k])] <- l.i.k
        u.i[seq(from = k, by = K, length.out = omega[k])] <- u.i.k
      }
      deltas.LFC[[i]]    <- delta.LFC[sort(relevant.indices)]
      I.div.n.LFC[[i]]   <- I.div.n[sort(relevant.indices)]
      Lambdas.LFC[[i]]   <- Lambda[sort(relevant.indices),
                                   sort(relevant.indices)]
      l.indices.LFC[[i]] <- l.i[sort(relevant.indices)]
      u.indices.LFC[[i]] <- u.i[sort(relevant.indices)]
    }
    return(list(thetas.d.HG.a.fwer = thetas.d.HG.a.fwer,
                thetas.bc.LFC.power = thetas.bc.LFC.power,
                l.indices.HG = l.indices.HG, u.indices.HG = u.indices.HG,
                means.HG = means.HG, Lambdas.HG = Lambdas.HG,
                l.indices.LFC = l.indices.LFC, u.indices.LFC = u.indices.LFC,
                deltas.LFC = deltas.LFC, I.div.n.LFC = I.div.n.LFC,
                Lambdas.LFC = Lambdas.LFC))
  }
  
  # Initialise function that computes values of the objective function for a
  # particular design
  objFn <- function(pars, K, J, alpha, beta, r, w, N.fixed, thetas.d.HG.a.fwer,
                    thetas.bc.LFC.power, l.indices.HG, u.indices.HG, means.HG,
                    Lambdas.HG, l.indices.LFC, u.indices.LFC, deltas.LFC,
                    I.div.n.LFC, Lambdas.LFC){
    n      <- pars[1]
    bounds <- c(pars[2:(J + 1)], c(pars[2:J] + pars[(J + 2):(2*J)],
                                   pars[J + 1]), -Inf, Inf)
    for (i in 1:nrow(thetas.d.HG.a.fwer)){
      thetas.d.HG.a.fwer[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.HG[[i]]],
                                                upper = bounds[u.indices.HG[[i]]],
                                                mean = means.HG[[i]],
                                                sigma = Lambdas.HG[[i]])[1]
    }
    for (i in 1:nrow(thetas.bc.LFC.power)){
      thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                                 upper = bounds[u.indices.LFC[[i]]],
                                                 mean = deltas.LFC[[i]]*
                                                   sqrt(n*I.div.n.LFC[[i]]),
                                                 sigma = Lambdas.LFC[[i]])[1]
    }
    a.fwer   <- sum(thetas.d.HG.a.fwer[, 2*K + 2]*
                      thetas.d.HG.a.fwer[, 2*K + 4])
    bc.power <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                      thetas.bc.LFC.power[, 2*K + 4])
    ESS.HG   <- n*sum(thetas.d.HG.a.fwer[, 2*K + 1]*
                        thetas.d.HG.a.fwer[, 2*K + 3]*
                        thetas.d.HG.a.fwer[, 2*K + 4])
    ESS.LFC  <- n*sum(thetas.bc.LFC.power[, 2*K + 1]*
                        thetas.bc.LFC.power[, 2*K + 3]*
                        thetas.bc.LFC.power[, 2*K + 4])
    o        <- sum(w*c(ESS.HG, ESS.LFC, n*J*(r*K + 1))) +
                  N.fixed*(as.numeric(a.fwer > alpha)*(a.fwer - alpha)/alpha +
                             as.numeric(1 - bc.power > beta)*
                               (1 - bc.power - beta)/beta)
    return(o)
  }
  
  # Initialise a function that finds a-FWER of a single-stage design
  r1j1 <- function(pars, K, alpha, thetas.d.HG.a.fwer, l.indices.HG,
                   u.indices.HG, means.HG, Lambdas.HG){
    bounds <- c(pars, pars, -Inf, Inf)
    for (i in 1:nrow(thetas.d.HG.a.fwer)){
      thetas.d.HG.a.fwer[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.HG[[i]]],
                                                upper = bounds[u.indices.HG[[i]]],
                                                mean = means.HG[[i]],
                                                sigma = Lambdas.HG[[i]])[1]
    }
    a.fwer <- sum(thetas.d.HG.a.fwer[, 2*K + 2]*thetas.d.HG.a.fwer[, 2*K + 4])
    return(a.fwer - alpha)
  }
  
  # Initialise a function that finds power of a single-stage design
  n1j1 <- function(pars, K, beta, thetas.bc.LFC.power, l.indices.LFC,
                   u.indices.LFC, deltas.LFC, I.div.n.LFC, Lambdas.LFC, r1){
    n      <- pars
    bounds <- c(r1, r1, -Inf, Inf)
    for (i in 1:nrow(thetas.bc.LFC.power)){
      thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                                 upper = bounds[u.indices.LFC[[i]]],
                                                 mean = deltas.LFC[[i]]*
                                                   sqrt(n*I.div.n.LFC[[i]]),
                                                 sigma = Lambdas.LFC[[i]])[1]
    }
    bc.power <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                      thetas.bc.LFC.power[, 2*K + 4])
    return(bc.power - (1 - beta))
  }
  
  ##### Main Computations ######################################################
  
  # Generate required variables for a corresponding single-stage design
  theta.information   <- thetaGeneration(K, 1, a, b, c, d, r, sigma, delta.1,
                                         delta.0)
  thetas.d.HG.a.fwer  <- theta.information$thetas.d.HG.a.fwer
  thetas.bc.LFC.power <- theta.information$thetas.bc.LFC.power
  l.indices.HG        <- theta.information$l.indices.HG
  u.indices.HG        <- theta.information$u.indices.HG
  means.HG            <- theta.information$means.HG
  Lambdas.HG          <- theta.information$Lambdas.HG
  l.indices.LFC       <- theta.information$l.indices.LFC
  u.indices.LFC       <- theta.information$u.indices.LFC
  deltas.LFC          <- theta.information$deltas.LFC
  I.div.n.LFC         <- theta.information$I.div.n.LFC
  Lambdas.LFC         <- theta.information$Lambdas.LFC
  # Determine rejection boundary and sample size required by corresponding
  # single-stage design
  r1 <- uniroot(r1j1, c(-100, 100), K = K, alpha = alpha,
                thetas.d.HG.a.fwer = thetas.d.HG.a.fwer,
                l.indices.HG = l.indices.HG, u.indices.HG = u.indices.HG,
                means.HG = means.HG, Lambdas.HG = Lambdas.HG)$root
  n  <- uniroot(n1j1, c(0, 10^6), K = K, beta = beta,
                thetas.bc.LFC.power = thetas.bc.LFC.power,
                l.indices.LFC = l.indices.LFC, u.indices.LFC = u.indices.LFC,
                deltas.LFC = deltas.LFC, I.div.n.LFC = I.div.n.LFC,
                Lambdas.LFC = Lambdas.LFC, r1 = r1)$root
  # Use this to set N_fixed
  N.fixed <- n*(r*K + 1)
  # Now generate variables for desired MAMS design
  theta.information   <- thetaGeneration(K, J, a, b, c, d, r, sigma, delta.1,
                                         delta.0)
  thetas.d.HG.a.fwer  <- theta.information$thetas.d.HG.a.fwer
  thetas.bc.LFC.power <- theta.information$thetas.bc.LFC.power
  l.indices.HG        <- theta.information$l.indices.HG
  u.indices.HG        <- theta.information$u.indices.HG
  means.HG            <- theta.information$means.HG
  Lambdas.HG          <- theta.information$Lambdas.HG
  l.indices.LFC       <- theta.information$l.indices.LFC
  u.indices.LFC       <- theta.information$u.indices.LFC
  deltas.LFC          <- theta.information$deltas.LFC
  I.div.n.LFC         <- theta.information$I.div.n.LFC
  Lambdas.LFC         <- theta.information$Lambdas.LFC
  # Assign fitness function for use in optimisation
  fitness <- function(...){
    -objFn(...)
  }
  if (parallel){
    parallel <- cpus
  }
  # Use ga() to search for optimal design
  start.time <- Sys.time()
  GA         <- ga(type = "real-valued", fitness = fitness,
                   K = K, J = J, alpha = alpha, beta = beta, r = r, w = w,
                   N.fixed = N.fixed, thetas.d.HG.a.fwer = thetas.d.HG.a.fwer,
                   thetas.bc.LFC.power = thetas.bc.LFC.power,
                   l.indices.HG = l.indices.HG, u.indices.HG = u.indices.HG,
                   means.HG = means.HG, Lambdas.HG = Lambdas.HG,
                   l.indices.LFC = l.indices.LFC, u.indices.LFC = u.indices.LFC,
                   deltas.LFC = deltas.LFC, I.div.n.LFC = I.div.n.LFC,
                   Lambdas.LFC = Lambdas.LFC,
                   min = c(0, rep(-20, J), rep(0, J - 1)),
                   max = c(n, rep(20, 2*J - 1)), popSize = popSize,
                   maxiter = ceiling(maxiter), run = run,
                   parallel = parallel, seed = 10*seed)
  end.time   <- Sys.time()
  # Assign parameters of optimal design and compute its operating
  # characteristics
  n.star     <- GA@solution[1]
  n          <- ceiling(n.star)
  while (r*n %% 1 != 0){
    n        <- n + 1
  }
  bounds <- c(GA@solution[2:(J + 1)],
              c(GA@solution[2:J] + GA@solution[(J + 2):(2*J)],
                GA@solution[J + 1]), -Inf, Inf)
  a      <- GA@solution[2:(J + 1)]
  r      <- c(GA@solution[2:J] + GA@solution[(J + 2):(2*J)],
              GA@solution[J + 1])
  for (i in 1:nrow(thetas.d.HG.a.fwer)){
    thetas.d.HG.a.fwer[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.HG[[i]]],
                                              upper = bounds[u.indices.HG[[i]]],
                                              mean = means.HG[[i]],
                                              sigma = Lambdas.HG[[i]])[1]
  }
  for (i in 1:nrow(thetas.bc.LFC.power)){
    thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                               upper = bounds[u.indices.LFC[[i]]],
                                               mean = deltas.LFC[[i]]*
                                                 sqrt(n*I.div.n.LFC[[i]]),
                                               sigma = Lambdas.LFC[[i]])[1]
  }
  a.fwer   <- sum(thetas.d.HG.a.fwer[, 2*K + 2]*
                    thetas.d.HG.a.fwer[, 2*K + 4])
  bc.power <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                    thetas.bc.LFC.power[, 2*K + 4])
  ESS.HG   <- n*sum(thetas.d.HG.a.fwer[, 2*K + 1]*
                      thetas.d.HG.a.fwer[, 2*K + 3]*
                      thetas.d.HG.a.fwer[, 2*K + 4])
  ESS.LFC  <- n*sum(thetas.bc.LFC.power[, 2*K + 1]*
                      thetas.bc.LFC.power[, 2*K + 3]*
                      thetas.bc.LFC.power[, 2*K + 4])
  N.max    <- n*J*(ratio*K + 1)
  o        <- sum(w*c(ESS.HG, ESS.LFC, N.max)) +
                N.fixed*(as.numeric(a.fwer > alpha)*(a.fwer - alpha)/alpha +
                           as.numeric(1 - bc.power > beta)*
                             (1 - bc.power - beta)/beta)
  return(list(n.star = n.star, n = n, a = a, r = r,
              a.fwer = a.fwer, bc.power = bc.power, ESS.HG = ESS.HG,
              ESS.LFC = ESS.LFC, N.max = N.max, o = o,
              run.time = end.time - start.time, GA = GA))
}