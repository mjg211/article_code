################################################################################
# Function: two_stage_binary()                                                 #
# Last Modified: 13/03/2018                                                    #
################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
#        K - Initial number of experimental treatment arms                     #
#        J - Maximal number of stages                                          #
#        r - Allocation ratio of experimental arms to the control              #
#    alpha - Desired type-I FWER                                               #
#     beta - Desired type-II FWER                                              #
#    delta - Used in the definition of the alternative hypotheses              #
#    n.min - Minimal value for n to consider                                   #
#    n.max - Maximal value for n to consider                                   #
#     pESS - Response probability to use in optimisation routine               #
# parallel - Should parallelisation be used                                    #
#     cpus - Number of cpus to parallelise over                                #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# A .csv file is produced for each value of n containing details on each       #
# design's operating characteristics                                           #
################################################################################

two_stage_binary <- function(K = 2, J = 2, r = 1, alpha = 0.15, beta = 0.2,
                             delta = 0.15, n.min = 2, n.max = 54, pESS = 0.7,
                             parallel = TRUE, cpus = 7){
  
  max.a.fwer.fn <- function(p, K, J, r, n, boundaries, scenarios.a.fwer,
                            poss.control.responses){
    binomial.dens.c       <- dbinom(0:n, n, p)
    binomial.dens.e       <- dbinom(0:(r*n), r*n, p)
    stop.prob             <- list()
    stop.prob[[1]]        <- list()
    stop.prob[[2]]        <- list()
    stop.prob[[3]]        <- list()
    for (j in 1:J){
      stop.prob[[1]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[2]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[3]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      for (i in 1:nrow(poss.control.responses[[j]])){
        cum.prob                         <- matrix(0, nrow = j,
                                                   ncol = n*j + r*n*j + 1)
        poss.diff                        <- 0:(r*n) -
                                              poss.control.responses[[j]][i, 1]
        cum.prob[1, poss.diff + j*n + 1] <- binomial.dens.e
        stop.prob[[1]][[j]][i, 1]        <- sum(cum.prob[1,
                                                         1:(n*j + 1 +
                                                              boundaries[1])])
        stop.prob[[2]][[j]][i, 1]        <- sum(cum.prob[1,
                                                         (n*j + 1 + 
                                                            boundaries[2] + 1):
                                                         (n*j + r*n*j + 1)])
        stop.prob[[3]][[j]][i, 1]        <- sum(cum.prob[1,
                                                         1:(n*j + 1 +
                                                              boundaries[2])])
        if (j > 1){
          poss.diff                              <- 0:(r*n) -
                                                      poss.control.responses[[j]][i, 2]
          poss.prev                              <- (boundaries[1] + 1):
                                                      boundaries[2]
          for (l in poss.prev){
            cum.prob[j, l + poss.diff + n*j + 1] <- cum.prob[j,
                                                             l + poss.diff +
                                                               n*j + 1] +
                                                      cum.prob[j - 1,
                                                               l + n*j + 1]*
                                                        binomial.dens.e
          }
          stop.prob[[1]][[j]][i, 2] <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[3])])
          stop.prob[[2]][[j]][i, 2] <- sum(cum.prob[2, (n*j + 1 + boundaries[4] + 1):
                                                          (n*j + r*n*j + 1)])
          stop.prob[[3]][[j]][i, 2] <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[4])])
        }
      }
    }
    for (i in 1:nrow(scenarios.a.fwer)){
      max.omega <- max(scenarios.a.fwer[i, 1:K])
      for (l in 1:nrow(poss.control.responses[[max.omega]])){
        prod   <- 1
        for (j in 1:max.omega){
          prod <- prod*binomial.dens.c[1 + poss.control.responses[[max.omega]][l, j]]
        }
        for (k in 1:K){
          if (scenarios.a.fwer[i, K + k] == 1){
            prod <- prod*stop.prob[[2]][[max.omega]][l, scenarios.a.fwer[i, k]]
          } else {
            if (scenarios.a.fwer[i, k] < max.omega){
              prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.a.fwer[i, k]]
            } else {
              psi   <- scenarios.a.fwer[i, (K + 1):(2*K)]
              omega <- scenarios.a.fwer[i, 1:K]
              if (scenarios.a.fwer[i, k] == max.omega &
                    any(psi[which(omega == max.omega)] == 1)){
                prod <- prod*stop.prob[[3]][[max.omega]][l, scenarios.a.fwer[i, k]]
              } else {
                prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.a.fwer[i, k]]
              }
            }
          }
        }
        scenarios.a.fwer[i, 2*K + 3] <- scenarios.a.fwer[i, 2*K + 3] + prod
      }
    }
    print(sum(scenarios.a.fwer[, 2*K + 1]*scenarios.a.fwer[, 2*K + 3]))
    return(-sum(scenarios.a.fwer[, 2*K + 1]*scenarios.a.fwer[, 2*K + 3]))
  }
  
  min.b.power.fn <- function(p, K, J, r, n, boundaries, scenarios.b.power,
                             poss.control.responses, delta, c){
    binomial.dens.c       <- dbinom(0:n, n, p)
    binomial.dens.e.E     <- dbinom(0:(r*n), r*n, p + delta)
    binomial.dens.e.F     <- dbinom(0:(r*n), r*n, p)
    stop.prob.E             <- list()
    stop.prob.E[[1]]        <- list()
    stop.prob.E[[2]]        <- list()
    stop.prob.E[[3]]        <- list()
    stop.prob.F             <- list()
    stop.prob.F[[1]]        <- list()
    stop.prob.F[[2]]        <- list()
    stop.prob.F[[3]]        <- list()
    for (j in 1:J){
      stop.prob.E[[1]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      stop.prob.E[[2]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      stop.prob.E[[3]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      stop.prob.F[[1]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      stop.prob.F[[2]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      stop.prob.F[[3]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                      ncol = j)
      for (i in 1:nrow(poss.control.responses[[j]])){
        cum.prob.E                         <- matrix(0, nrow = j,
                                                     ncol = n*j + r*n*j + 1)
        cum.prob.F                         <- matrix(0, nrow = j,
                                                     ncol = n*j + r*n*j + 1)
        poss.diff                        <- 0:(r*n) - poss.control.responses[[j]][i, 1]
        cum.prob.E[1, poss.diff + j*n + 1] <- binomial.dens.e.E
        cum.prob.F[1, poss.diff + j*n + 1] <- binomial.dens.e.F
        stop.prob.E[[1]][[j]][i, 1]        <- sum(cum.prob.E[1, 1:(n*j + 1 + boundaries[1])])
        stop.prob.E[[2]][[j]][i, 1]        <- sum(cum.prob.E[1, (n*j + 1 + boundaries[2] + 1):
                                                               (n*j + r*n*j + 1)])
        stop.prob.E[[3]][[j]][i, 1]        <- sum(cum.prob.E[1, 1:(n*j + 1 + boundaries[2])])
        stop.prob.F[[1]][[j]][i, 1]        <- sum(cum.prob.F[1, 1:(n*j + 1 + boundaries[1])])
        stop.prob.F[[2]][[j]][i, 1]        <- sum(cum.prob.F[1, (n*j + 1 + boundaries[2] + 1):
                                                               (n*j + r*n*j + 1)])
        stop.prob.F[[3]][[j]][i, 1]        <- sum(cum.prob.F[1, 1:(n*j + 1 + boundaries[2])])
        if (j > 1){
          poss.diff                              <- 0:(r*n) - poss.control.responses[[j]][i, 2]
          poss.prev                              <- (boundaries[1] + 1):boundaries[2]
          for (l in poss.prev){
            cum.prob.E[j, l + poss.diff + n*j + 1] <- cum.prob.E[j, l + poss.diff + n*j + 1] +
                                                        cum.prob.E[j - 1, l + n*j + 1]*
                                                          binomial.dens.e.E
            cum.prob.F[j, l + poss.diff + n*j + 1] <- cum.prob.F[j, l + poss.diff + n*j + 1] +
                                                        cum.prob.F[j - 1, l + n*j + 1]*
                                                          binomial.dens.e.F
          }
          stop.prob.E[[1]][[j]][i, 2] <- sum(cum.prob.E[2, 1:(n*j + 1 + boundaries[3])])
          stop.prob.E[[2]][[j]][i, 2] <- sum(cum.prob.E[2, (n*j + 1 + boundaries[4] + 1):
                                                              (n*j + r*n*j + 1)])
          stop.prob.E[[3]][[j]][i, 2] <- sum(cum.prob.E[2, 1:(n*j + 1 + boundaries[4])])
          stop.prob.F[[1]][[j]][i, 2] <- sum(cum.prob.F[2, 1:(n*j + 1 + boundaries[3])])
          stop.prob.F[[2]][[j]][i, 2] <- sum(cum.prob.F[2, (n*j + 1 + boundaries[4] + 1):
                                                              (n*j + r*n*j + 1)])
          stop.prob.F[[3]][[j]][i, 2] <- sum(cum.prob.F[2, 1:(n*j + 1 + boundaries[4])])
        }
      }
    }
    for (i in 1:nrow(scenarios.b.power)){
      max.omega <- max(scenarios.b.power[i, 1:K])
      for (l in 1:nrow(poss.control.responses[[max.omega]])){
        prod   <- 1
        for (j in 1:max.omega){
          prod <- prod*binomial.dens.c[1 + poss.control.responses[[max.omega]][l, j]]
        }
        
        for (k in 1:c){
          if (scenarios.b.power[i, K + k] == 1){
            prod <- prod*stop.prob.E[[2]][[max.omega]][l, scenarios.b.power[i, k]]
          } else {
            if (scenarios.b.power[i, k] < max.omega){
              prod <- prod*stop.prob.E[[1]][[max.omega]][l, scenarios.b.power[i, k]]
            } else {
              psi   <- scenarios.b.power[i, (K + 1):(2*K)]
              omega <- scenarios.b.power[i, 1:K]
              if (scenarios.b.power[i, k] == max.omega &
                    any(psi[which(omega == max.omega)] == 1)){
                prod <- prod*stop.prob.E[[3]][[max.omega]][l, scenarios.b.power[i, k]]
              } else {
                prod <- prod*stop.prob.E[[1]][[max.omega]][l, scenarios.b.power[i, k]]
              }
            }
          }
        }
        if (c < K){
          for (k in (c + 1):K){
            if (scenarios.b.power[i, K + k] == 1){
              prod <- prod*stop.prob.F[[2]][[max.omega]][l, scenarios.b.power[i, k]]
            } else {
              if (scenarios.b.power[i, k] < max.omega){
                prod <- prod*stop.prob.F[[1]][[max.omega]][l, scenarios.b.power[i, k]]
              } else {
                psi   <- scenarios.b.power[i, (K + 1):(2*K)]
                omega <- scenarios.b.power[i, 1:K]
                if (scenarios.b.power[i, k] == max.omega &
                      any(psi[which(omega == max.omega)] == 1)){
                  prod <- prod*stop.prob.F[[3]][[max.omega]][l, scenarios.b.power[i, k]]
                } else {
                  prod <- prod*stop.prob.F[[1]][[max.omega]][l, scenarios.b.power[i, k]]
                }
              }
            }
          }
        }
        scenarios.b.power[i, 2*K + 3] <- scenarios.b.power[i, 2*K + 3] + prod
      }
    }
    return(sum(scenarios.b.power[, 2*K + 1]*scenarios.b.power[, 2*K + 3]))
  }
  
  max.ESS.HG.fn <- function(p, K, J, r, n, boundaries, scenarios.global,
                            poss.control.responses){
    binomial.dens.c       <- dbinom(0:n, n, p)
    binomial.dens.e       <- dbinom(0:(r*n), r*n, p)
    stop.prob             <- list()
    stop.prob[[1]]        <- list()
    stop.prob[[2]]        <- list()
    stop.prob[[3]]        <- list()
    for (j in 1:J){
      stop.prob[[1]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[2]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[3]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      for (i in 1:nrow(poss.control.responses[[j]])){
        cum.prob                         <- matrix(0, nrow = j, ncol = n*j + r*n*j + 1)
        poss.diff                        <- 0:(r*n) - poss.control.responses[[j]][i, 1]
        cum.prob[1, poss.diff + j*n + 1] <- binomial.dens.e
        stop.prob[[1]][[j]][i, 1]        <- sum(cum.prob[1, 1:(n*j + 1 + boundaries[1])])
        stop.prob[[2]][[j]][i, 1]        <- sum(cum.prob[1, (n*j + 1 + boundaries[2] + 1):
                                                              (n*j + r*n*j + 1)])
        stop.prob[[3]][[j]][i, 1]        <- sum(cum.prob[1, 1:(n*j + 1 + boundaries[2])])
        if (j > 1){
          poss.diff                              <- 0:(r*n) - poss.control.responses[[j]][i, 2]
          poss.prev                              <- (boundaries[1] + 1):boundaries[2]
          for (l in poss.prev){
            cum.prob[j, l + poss.diff + n*j + 1] <- cum.prob[j, l + poss.diff + n*j + 1] +
                                                      cum.prob[j - 1, l + n*j + 1]*
                                                        binomial.dens.e
          }
          stop.prob[[1]][[j]][i, 2]              <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[3])])
          stop.prob[[2]][[j]][i, 2]              <- sum(cum.prob[2, (n*j + 1 + boundaries[4] + 1):
                                                                   (n*j + r*n*j + 1)])
          stop.prob[[3]][[j]][i, 2]              <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[4])])
        }
      }
    }
    for (i in 1:nrow(scenarios.global)){
      max.omega <- max(scenarios.global[i, 1:K])
      for (l in 1:nrow(poss.control.responses[[max.omega]])){
        prod   <- 1
        for (j in 1:max.omega){
          prod <- prod*binomial.dens.c[1 + poss.control.responses[[max.omega]][l, j]]
        }
        for (k in 1:K){
          if (scenarios.global[i, K + k] == 1){
            prod <- prod*stop.prob[[2]][[max.omega]][l, scenarios.global[i, k]]
          } else {
            if (scenarios.global[i, k] < max.omega){
              prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.global[i, k]]
            } else {
              psi   <- scenarios.global[i, (K + 1):(2*K)]
              omega <- scenarios.global[i, 1:K]
              if (scenarios.global[i, k] == max.omega &
                    any(psi[which(omega == max.omega)] == 1)){
                prod <- prod*stop.prob[[3]][[max.omega]][l, scenarios.global[i, k]]
              } else {
                prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.global[i, k]]
              }
            }
          }
        }
        scenarios.global[i, 2*K + 3] <- scenarios.global[i, 2*K + 3] + prod
      }
    }
    return(-n*sum(scenarios.global[, 2*K + 1]*scenarios.global[, 2*K + 2]*
                    scenarios.global[, 2*K + 3]))
  }
  
  max.ESS.HA.fn <- function(p, K, J, r, n, boundaries, scenarios.global,
                            poss.control.responses, delta){
    binomial.dens.c       <- dbinom(0:n, n, p)
    binomial.dens.e       <- dbinom(0:(r*n), r*n, p + delta)
    stop.prob             <- list()
    stop.prob[[1]]        <- list()
    stop.prob[[2]]        <- list()
    stop.prob[[3]]        <- list()
    for (j in 1:J){
      stop.prob[[1]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[2]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      stop.prob[[3]][[j]] <- matrix(0, nrow = nrow(poss.control.responses[[j]]),
                                    ncol = j)
      for (i in 1:nrow(poss.control.responses[[j]])){
        cum.prob                         <- matrix(0, nrow = j, ncol = n*j + r*n*j + 1)
        poss.diff                        <- 0:(r*n) - poss.control.responses[[j]][i, 1]
        cum.prob[1, poss.diff + j*n + 1] <- binomial.dens.e
        stop.prob[[1]][[j]][i, 1]        <- sum(cum.prob[1, 1:(n*j + 1 + boundaries[1])])
        stop.prob[[2]][[j]][i, 1]        <- sum(cum.prob[1, (n*j + 1 + boundaries[2] + 1):
                                                           (n*j + r*n*j + 1)])
        stop.prob[[3]][[j]][i, 1]        <- sum(cum.prob[1, 1:(n*j + 1 + boundaries[2])])
        if (j > 1){
          poss.diff                              <- 0:(r*n) - poss.control.responses[[j]][i, 2]
          poss.prev                              <- (boundaries[1] + 1):boundaries[2]
          for (l in poss.prev){
            cum.prob[j, l + poss.diff + n*j + 1] <- cum.prob[j, l + poss.diff + n*j + 1] +
                                                      cum.prob[j - 1, l + n*j + 1]*
                                                        binomial.dens.e
          }
          stop.prob[[1]][[j]][i, 2] <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[3])])
          stop.prob[[2]][[j]][i, 2] <- sum(cum.prob[2, (n*j + 1 + boundaries[4] + 1):
                                                          (n*j + r*n*j + 1)])
          stop.prob[[3]][[j]][i, 2] <- sum(cum.prob[2, 1:(n*j + 1 + boundaries[4])])
        }
      }
    }
    for (i in 1:nrow(scenarios.global)){
      max.omega <- max(scenarios.global[i, 1:K])
      for (l in 1:nrow(poss.control.responses[[max.omega]])){
        prod   <- 1
        for (j in 1:max.omega){
          prod <- prod*binomial.dens.c[1 + poss.control.responses[[max.omega]][l, j]]
        }
        for (k in 1:K){
          if (scenarios.global[i, K + k] == 1){
            prod <- prod*stop.prob[[2]][[max.omega]][l, scenarios.global[i, k]]
          } else {
            if (scenarios.global[i, k] < max.omega){
              prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.global[i, k]]
            } else {
              psi   <- scenarios.global[i, (K + 1):(2*K)]
              omega <- scenarios.global[i, 1:K]
              if (scenarios.global[i, k] == max.omega &
                    any(psi[which(omega == max.omega)] == 1)){
                prod <- prod*stop.prob[[3]][[max.omega]][l, scenarios.global[i, k]]
              } else {
                prod <- prod*stop.prob[[1]][[max.omega]][l, scenarios.global[i, k]]
              }
            }
          }
        }
        scenarios.global[i, 2*K + 3] <- scenarios.global[i, 2*K + 3] + prod
      }
    }
    return(-n*sum(scenarios.global[, 2*K + 1]*scenarios.global[, 2*K + 2]*
                    scenarios.global[, 2*K + 3]))
  }
  
  design.performance <- function(K, J, r, n, boundaries, scenarios.global,
                                 scenarios.a.fwer, scenarios.b.power, delta, c,
                                 poss.control.responses, pESS){
    max.a.fwer   <- max.a.fwer.fn(0.5, K, J, r, n, boundaries, scenarios.a.fwer,
                                  poss.control.responses)
    if (-max.a.fwer <= alpha){
      max.a.fwer   <- optim(par = 0.5, max.a.fwer.fn, K = K, J = J, r = r, n = n,
                            boundaries = boundaries,
                            scenarios.a.fwer = scenarios.a.fwer,
                            poss.control.responses = poss.control.responses,
                            method = "Brent", lower = 0, upper = 1,
                            control = list(reltol = 1e-3))
      max.a.fwer.p <- max.a.fwer$par
      max.a.fwer.v <- -max.a.fwer$value
    } else {
      max.a.fwer.p <- 0.5
      max.a.fwer.v <- -max.a.fwer
    }
    # Determine minimal power
    if (max.a.fwer.v <= alpha){
      min.b.power <- min.b.power.fn(0.5 - delta/2, K, J, r, n, boundaries,
                                    scenarios.b.power, poss.control.responses,
                                    delta = delta, c = c)
    } else {
      min.b.power <- 0
    }
    if (min.b.power >= 1 - beta){
      min.b.power.2 <- min.b.power.fn(1 - delta, K, J, r, n, boundaries,
                                      scenarios.b.power, poss.control.responses,
                                      delta = delta, c = c)
      if (min.b.power.2 >= 1 - beta){
        min.b.power   <- optim(par = 0.5 - delta/2, min.b.power.fn, K = K,
                               J = J, r = r, n = n, boundaries = boundaries,
                               scenarios.b.power = scenarios.b.power,
                               poss.control.responses = poss.control.responses,
                               delta = delta, c = c, method = "Brent", lower = 0,
                               upper = 1 - delta, control = list(reltol = 1e-3))
        min.b.power.p <- min.b.power$par
        min.b.power.v <- min.b.power$value
      } else {
        min.b.power.p <- 1 - delta
        min.b.power.v <- min.b.power.2
      }
    } else {
      min.b.power.p <- 0.5 - delta/2
      min.b.power.v <- min.b.power
    }
    if (min.b.power.v >= 1 - beta & max.a.fwer.v <= alpha){
      max.ESS.HG   <- max.ESS.HG.fn(pESS, K, J, r, n, boundaries,
                                    scenarios.global,
                                    poss.control.responses)
      max.ESS.HG.p <- pESS
      max.ESS.HG.v <- -max.ESS.HG
      max.ESS.HA   <- max.ESS.HA.fn(pESS, K, J, r, n, boundaries,
                                    scenarios.global,
                                    poss.control.responses,
                                    delta)
      max.ESS.HA.p <- pESS
      max.ESS.HA.v <- -max.ESS.HA
    } else {
      max.ESS.HG.p <- 0
      max.ESS.HG.v <- 0
      max.ESS.HA.p <- 0
      max.ESS.HA.v <- 0
    }
    
    max.N <- J*n + J*K*r*n
    
    return(c(max.a.fwer.p, max.a.fwer.v, min.b.power.p, min.b.power.v,
             max.ESS.HG.p, max.ESS.HG.v, max.ESS.HA.p, max.ESS.HA.v,
             max.N))
  }
  
  wrapper <- function(index){
    if (J == 1){
      boundaries <- rep(designs[index, 2], 2)
    } else {
      boundaries <- c(designs[index, 2:4], designs[index, 4])
    }
    result <- design.performance(K, J, r, designs[index, 1], boundaries,
                                 scenarios.global, scenarios.a.fwer,
                                 scenarios.b.power, delta, c,
                                 poss.control.responses[[designs[index, 1]]],
                                 pESS)
    return(result)
  }
  
  ##### Main Computations #####################################################
  
  R <- 1; a <- 1; b <- 1; c <- 2
  
  all.scenarios            <- getall(iterpc(2*J, K, replace = T))
  degeneracy               <- numeric(nrow(all.scenarios))
  degeneracy.a.fwer        <- numeric(nrow(all.scenarios))
  samp.size.factor         <- numeric(nrow(all.scenarios))
  for (i in 1:nrow(all.scenarios)){
    degeneracy[i]          <- nrow(getall(iterpc(table(all.scenarios[i, ]),
                                                 ordered = T)))
    if (sum(1*is.even(all.scenarios[i, ])) >= a){
      degeneracy.a.fwer[i] <- degeneracy[i]
    }
    sum <- 0
    for (k in 1:K){
      sum <- sum + ceiling(all.scenarios[i, k]/2)*r
    }
    samp.size.factor[i] <- sum + max(ceiling(all.scenarios[i, 1:K]/2))
  }
  scenarios.global           <- cbind(ceiling(all.scenarios/2),
                                      matrix(1*is.even(all.scenarios),
                                             ncol = K),
                                      degeneracy, samp.size.factor,
                                      numeric(nrow(all.scenarios)))
  colnames(scenarios.global) <- c(paste("omega_", 1:K, sep = ""),
                                  paste("psi_", 1:K, sep = ""),
                                  "d_G(omega,psi)", "N(omega,psi)/n",
                                  "P(omega,psi|global)")
  retain            <- rep(TRUE, nrow(scenarios.global))
  for (i in 1:nrow(scenarios.global)){
    for (j in 1:J){
      le.eq         <- which(scenarios.global[i, 1:K] <= j)
      gr            <- which(scenarios.global[i, 1:K] > j)
      if (all(c(length(le.eq), length(gr)) > 0)){
        if (sum(scenarios.global[i, K + le.eq]) >= R){
          retain[i] <- FALSE
        }
      }
    }
  }
  scenarios.global  <- scenarios.global[retain, ]
  scenarios.a.fwer           <- cbind(ceiling(all.scenarios/2),
                                      matrix(1*is.even(all.scenarios), ncol = K),
                                      degeneracy.a.fwer, samp.size.factor,
                                      numeric(nrow(all.scenarios)))
  colnames(scenarios.a.fwer) <- c(paste("omega_", 1:K, sep = ""),
                                  paste("psi_", 1:K, sep = ""),
                                  "d_a-FWER(omega,psi)", "N(omega,psi)/n",
                                  "P(omega,psi|global)")
  retain            <- rep(TRUE, nrow(scenarios.a.fwer))
  for (i in 1:nrow(scenarios.a.fwer)){
    if (degeneracy.a.fwer[i] > 0){
      for (j in 1:J){
        le.eq         <- which(scenarios.a.fwer[i, 1:K] <= j)
        gr            <- which(scenarios.a.fwer[i, 1:K] > j)
        if (all(c(length(le.eq), length(gr)) > 0)){
          if (sum(scenarios.a.fwer[i, K + le.eq]) >= R){
            retain[i] <- FALSE
          }
        }
      }
    } else {
      retain[i] <- FALSE
    }
  }
  scenarios.a.fwer  <- scenarios.a.fwer[retain, ]
  all.scenarios.E          <- getall(iterpc(2*J, c, replace = T))
  degeneracy.E             <- numeric(nrow(all.scenarios.E))
  samp.size.factor.E       <- numeric(nrow(all.scenarios.E))
  for (i in 1:nrow(all.scenarios.E)){
    degeneracy             <- nrow(getall(iterpc(table(all.scenarios.E[i, ]),
                                                 ordered = T)))
    if (sum(1*is.even(all.scenarios.E[i, ])) >= b){
      degeneracy.E[i] <- degeneracy
    }
    sum <- 0
    for (k in 1:c){
      sum <- sum + ceiling(all.scenarios.E[i, k]/2)*r
    }
    if (c == K){
      samp.size.factor.E[i] <- sum + max(ceiling(all.scenarios.E[i, 1:K]/2))
    } else {
      samp.size.factor.E[i] <- sum
    }
  }
  all.scenarios.E    <- all.scenarios.E[which(degeneracy.E > 0), ]
  samp.size.factor.E <- samp.size.factor.E[which(degeneracy.E > 0)]
  degeneracy.E       <- degeneracy.E[which(degeneracy.E > 0)]
  if (c == K){
    scenarios.b.power          <- cbind(ceiling(all.scenarios.E/2),
                                        matrix(1*is.even(all.scenarios.E), ncol = K),
                                        degeneracy.E, samp.size.factor.E,
                                        numeric(nrow(all.scenarios.E)))
    colnames(scenarios.b.power) <- c(paste("omega_", 1:K, sep = ""),
                                     paste("psi_", 1:K, sep = ""),
                                     "d_b-power(omega,psi)", "N(omega,psi)/n",
                                     "P(omega,psi|global)")
  } else {
    all.scenarios.F        <- getall(iterpc(2*J, K - c, replace = T))
    degeneracy.F           <- numeric(nrow(all.scenarios.F))
    samp.size.factor.F     <- numeric(nrow(all.scenarios.F))
    for (i in 1:nrow(all.scenarios.F)){
      degeneracy.F[i]      <- nrow(getall(iterpc(table(all.scenarios.F[i, ]),
                                                 ordered = T)))
      sum <- 0
      for (k in 1:(K - c)){
        sum <- sum + ceiling(all.scenarios.F[i, k]/2)*r
      }
      samp.size.factor.F[i] <- sum
    }
    scenarios.b.power <- matrix(0, nrow = nrow(all.scenarios.E)*nrow(all.scenarios.F),
                                ncol = 2*K + 3)
    colnames(scenarios.b.power) <- c(paste("omega_", 1:K, sep = ""),
                                     paste("psi_", 1:K, sep = ""),
                                     "d_b-power(omega,psi)", "N(omega,psi)/n",
                                     "P(omega,psi|global)")
    for (i in 1:nrow(all.scenarios.E)){
      range <- (1 + (i - 1)*nrow(all.scenarios.F)):(i*nrow(all.scenarios.F))
      matrix.E <- matrix(all.scenarios.E[i, ],
                         nrow = nrow(all.scenarios.F),
                         ncol = c,
                         byrow = T)
      scenarios.b.power[range, 1:(2*K + 1)] <- cbind(ceiling(matrix.E/2),
                                                     ceiling(all.scenarios.F/2),
                                                     matrix(1*is.even(matrix.E),
                                                            ncol = c),
                                                     matrix(1*is.even(all.scenarios.F),
                                                            ncol = K - c),
                                                     degeneracy.E[i]*degeneracy.F)
      for (j in 1:nrow(all.scenarios.F)){
        scenarios.b.power[j + (i - 1)*nrow(all.scenarios.F),
                          2*K + 2] <- samp.size.factor.E[i] + samp.size.factor.F[j] +
          max(scenarios.b.power[j + (i - 1)*nrow(all.scenarios.F), 1:K])
      }
    }
  }
  retain            <- rep(TRUE, nrow(scenarios.b.power))
  for (i in 1:nrow(scenarios.b.power)){
    for (j in 1:J){
      le.eq         <- which(scenarios.b.power[i, 1:K] <= j)
      gr            <- which(scenarios.b.power[i, 1:K] > j)
      if (all(c(length(le.eq), length(gr)) > 0)){
        if (sum(scenarios.b.power[i, K + le.eq]) >= R){
          retain[i] <- FALSE
        }
      }
    }
  }
  scenarios.b.power  <- scenarios.b.power[retain, ]
  poss.control.responses               <- list()
  for (i in 2:n.max){
    poss.control.responses[[i]]        <- list()
    for (j in 1:J){
      poss.control.responses[[i]][[j]] <- getall(iterpc(i + 1, j, labels = 0:i,
                                                        ordered = T, replace = T))
    }
  }
  poss.designs        <- list()
  if (J == 1){
    poss.n <- n.min:n.max
    for (i in n.min:n.max){
      poss.designs[[i]] <- matrix(-i:(i - 1), ncol = 1)
    }
  } else {
    poss.n <- n.min:n.max
    if (r < 1){
      poss.n <- which((n.min:n.max)%%(1/r) == 0)
    } else if (r > 1) {
      poss.n <- which((r*(n.min:n.max))%%1 == 0)
    }
    for (i in poss.n){
      designs.prev      <- getall(iterpc(r*i + i, 2, labels = (-i):(r*i - 1)))
      if (J == 2){
        designs.new       <- NULL
        for (l in 1:nrow(designs.prev)){
          designs.l       <- (designs.prev[l, 1] - i + 1):(designs.prev[l, 2] + r*i - 1)
          designs.new     <- rbind(designs.new,
                                   cbind(matrix(designs.prev[l, ],
                                                nrow = length(designs.l),
                                                ncol = 2,
                                                byrow = T),
                                         designs.l))
        }
        keep <- rep(TRUE, nrow(designs.new))
        for (l in 1:nrow(designs.new)){
          if (designs.new[l, 1] < 0 & designs.new[l, 3] < 0){
            keep[l] <- FALSE
          } else if (designs.new[l, 3] <= designs.new[l, 1]){
            keep[l] <- FALSE
          }
        }
        poss.designs[[i]] <- designs.new[keep, ]
      } else {
        for (j in 2:(J - 1)){
          designs.new       <- NULL
          for (l in 1:nrow(designs.prev)){
            designs.l       <- getall(iterpc(designs.prev[l, 2*(j - 1)] -
                                               designs.prev[l, 2*(j - 1) - 1] + r*i - 1 + i,
                                             2, labels = (designs.prev[l, 2*(j - 1) - 1] - i + 1):
                                                           (designs.prev[l, 2*(j - 1)] + r*i - 1)))
            designs.new     <- rbind(designs.new,
                                     cbind(matrix(designs.prev[l, ],
                                                  nrow = nrow(designs.l),
                                                  ncol = 2*(j - 1),
                                                  byrow = T),
                                           designs.l))
          }
          keep <- rep(TRUE, nrow(designs.new))
          for (l in 1:nrow(designs.new)){
            if (all(designs.new[l, seq(from = 1, to = 2*(j - 1) - 1, by = 2)] < 0) & 
                all(designs.new[l, seq(from = 2, to = 2*(j - 1), by = 2)] >= 0) &
                designs.new[l, 2*j] < 0){
              keep[l] <- FALSE
            }
          }
          designs.prev <- designs.new[keep, ]
        }
        designs.new  <- NULL
        for (l in 1:nrow(designs.prev)){
          designs.l       <- (designs.prev[l, 2*(J - 1) - 1] - i + 1):
                               (designs.prev[l, 2*(J - 1)] + r*i - 1)
          designs.new     <- rbind(designs.new,
                                   cbind(matrix(designs.prev[l, ],
                                                nrow = length(designs.l),
                                                ncol = 2*(J - 1),
                                                byrow = T),
                                         designs.l))
        }
        keep <- rep(TRUE, nrow(designs.new))
        for (l in 1:nrow(designs.new)){
          if (all(designs.new[l, seq(from = 1, to = 2*(J - 1) - 1, by = 2)] < 0) & 
              all(designs.new[l, seq(from = 2, to = 2*(J - 1), by = 2)] >= 0) &
              designs.new[l, 2*J - 1] < 0){
            keep[l] <- FALSE
          }
        }
        poss.designs[[i]] <- designs.new[keep, ]
      }
      print(i)
    }
  }      
  
  num.designs <- 0
  for (i in poss.n){
    num.designs <- num.designs + nrow(poss.designs[[i]])
  }
  
  designs <- matrix(0, nrow = num.designs, ncol = 2*J + 9)
  current <- 1
  for (i in poss.n){
    designs[current:
              (current - 1 + nrow(poss.designs[[i]])), 1:(2*J)] <- cbind(rep(i, nrow(poss.designs[[i]])),
                                                                         poss.designs[[i]])
    current <- current + nrow(poss.designs[[i]])
  }
  
  for (i in poss.n){
    designs <- matrix(0, nrow = nrow(poss.designs[[i]]), ncol = 2*J + 9)
    designs[, 1:(2*J)] <- cbind(rep(i, nrow(poss.designs[[i]])),
                                poss.designs[[i]])
    suppressMessages(sfInit(parallel = parallel, cpus = cpus))
    sfExport("K", "J", "r", "scenarios.global", "scenarios.a.fwer",
             "scenarios.b.power", "delta", "poss.control.responses",
             "designs", "design.performance", "c",
             "max.a.fwer.fn", "min.b.power.fn", "max.ESS.HG.fn",
             "max.ESS.HA.fn", "pESS", "alpha", "beta")
    results <- sfLapply(1:nrow(poss.designs[[i]]), wrapper)
    designs[, (2*J + 1):(2*J + 9)] <- matrix(unlist(results), ncol = 9,
                                             byrow = TRUE)
    if (J == 1){
      boundaries <- c("f_1")
    } else {
      boundaries <- c("f_1", "e_1", "f_2")
    }
    colnames(designs) <- c("n", boundaries, "p.max.a.fwer", "max.a.fwer",
                           "p.min.b.power", "min.b.power", "p.max.ESS.HG",
                           "max.ESS.HG", "p.max.ESS.HA", "max.ESS.HA",
                           "max.N")
    write.csv(designs, paste("exact.designs.", i, ".csv", sep = ""))
  }
}