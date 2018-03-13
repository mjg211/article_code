################################################################################
# Function: two_stage_fisher() and fisher_wrapper()                            #
# Last Modified: 13/03/2018                                                    #
################################################################################
# Inputs (to two_stage_fisher())                                               #
# ###### #######################                                               #
#                                                                              #
#         K - Initial number of experimental treatment arms                    #
#         J - Maximal number of stages                                         #
#         r - Allocation ratio of experimental arms to the control             #
#     alpha - Desired type-I FWER                                              #
#      beta - Desired type-II FWER                                             #
#   alpha.1 - Value of alpha.1 to use in determining the design                #
#    beta.1 - Value of beta.1 to use in determining the design                 #
#     delta - Used in the definition of the alternative hypotheses             #
#     n.max - Maximal value for n to consider                                  #
#    n.init - Value of n to try initially                                      #
#      pESS - Response probability to use in optimisation routine              #
# optimized - Determines whether optimized boundaries should be found, or      #
#             those used by Jung and Sargent                                   #
#  parallel - Should parallelisation be used                                   #
#      cpus - Number of cpus to parallelise over                               #
#   summary - Determines whether a summary of progress should be printed       #
################################################################################
# Outputs (from fisher_wrapper())                                              #
# ####### #######################                                              #
#                                                                              #
# Five .csv files are produced for each value of n containing details on the   #
# design and its operating characteristics                                     #
################################################################################

two_stage_fisher <- function(K = 2, r = 1, alpha = 0.15, beta = 0.2,
                             alpha.1 = 0.075, beta.1 = 0.1, delta = 0.15,
                             n.max = 50, n.init = 40, pESS = 0.7,
                             optimized = TRUE, summary = TRUE){
  
  f1.finder <- function(p, poss.responses, K, r, n, f1, delta){
    theta.pwr <- (p + delta)*(1 - p)/(p*(1 - p - delta))
    p.vec     <- c(p, p + delta)
    n.vec     <- c(n, rep(r*n, K))
    g.z       <- numeric(1 + n + K*r*n)
    marginal.II <- 0
    for (z in 0:(n + K*r*n)){
      for (i in 1:nrow(poss.responses[[z + 1]])){
        g.z[z + 1]                          <- g.z[z + 1] +
                                                 poss.responses[[z + 1]][i, 2*K + 3]*
                                                   prod((p.vec^poss.responses[[z + 1]][i, 1:(K + 1)])*
                                                          ((1 - p.vec)^(n.vec - poss.responses[[z + 1]][i, 1:(K + 1)])))
        poss.responses[[z + 1]][i, 2*K + 3] <- poss.responses[[z + 1]][i, 2*K + 3]*
                                                 prod(theta.pwr^poss.responses[[z + 1]][i, 2:(K + 1)])
      }
      poss.responses[[z + 1]][, 2*K + 3] <- poss.responses[[z + 1]][, 2*K + 3]/
                                              sum(poss.responses[[z + 1]][, 2*K + 3])
      marginal.II                        <- marginal.II + 
                                              sum(poss.responses[[z + 1]][poss.responses[[z + 1]][, K + 3] <= f1,
                                                                          2*K + 3])*g.z[z + 1]
    }
    return(-marginal.II)
  }
  
  marginal.power.fn <- function(p, poss.responses, K, r, n, f1, e1, e2, delta){
    marginal.power.p <- 0
    theta.pwr        <- (p + delta)*(1 - p)/(p*(1 - p - delta))
    p.vec            <- c(p, p + delta)
    n.vec            <- c(n, rep(r*n, K))
    g.z              <- list()
    for (d in 1:K){
      g.z[[d]]            <- numeric(1 + n + d*r*n)
      for (z in 0:(n + d*r*n)){
        for (i in 1:nrow(poss.responses[[d]][[n]][[z + 1]])){
          g.z[[d]][z + 1] <- g.z[[d]][z + 1] +
                               poss.responses[[d]][[n]][[z + 1]][i, 2*d + 3]*
                                 prod((p.vec[1:(d + 1)]^poss.responses[[d]][[n]][[z + 1]][i, 1:(d + 1)])*
                                        ((1 - p.vec[1:(d + 1)])^(n.vec[1:(d + 1)] - poss.responses[[d]][[n]][[z + 1]][i, 1:(d + 1)])))
        }
      }
    }
    for (z1 in 0:(n + K*r*n)){
      poss.x1               <- poss.responses[[K]][[n]][[z1 + 1]]
      for (i in 1:nrow(poss.x1)){
        poss.x1[i, 2*K + 3] <- poss.x1[i, 2*K + 3]*
                                 prod(theta.pwr^poss.x1[i, 2:(K + 1)])
      }
      poss.x1[, 2*K + 3]    <- poss.x1[, 2*K + 3]/sum(poss.x1[, 2*K + 3])
      if (nrow(poss.x1) > 1){
        if (K > 1){
          marginal.power.p    <- marginal.power.p +
                                   sum(poss.x1[rowSums(poss.x1[, (K + 3):(2*K + 2)] >= e1[z1 + 1]) > 0,
                                               2*K + 3])*g.z[[K]][z1 + 1]
        } else {
          marginal.power.p    <- marginal.power.p +
                                   sum(poss.x1[poss.x1[, (K + 3):(2*K + 2)] >= e1[z1 + 1],
                                               2*K + 3])*g.z[[K]][z1 + 1]
        }
      } else {
        marginal.power.p    <- marginal.power.p +
                                 sum(poss.x1[any(poss.x1[, (K + 3):(2*K + 2)] >= e1[z1 + 1]),
                                             2*K + 3])*g.z[[K]][z1 + 1]
      }
      results.1             <- matrix(0, nrow = nrow(poss.x1), ncol = K)
      for (i in 1:nrow(poss.x1)){
        results.1[i, 1:K]   <- 1*(any(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1]) |
                                    poss.x1[i, (K + 3):(2*K + 2)] <= f1)
      }
      for (d in 1:K){
        if (nrow(poss.x1) > 1){
          if (K > 1){
            poss.x1.d <- poss.x1[apply(results.1[, 1:K], 1,
                                       function(x) length(which(x == 0)) == d), ,
                                 drop = FALSE]
          } else {
            poss.x1.d <- poss.x1[results.1[, 1] == 0, , drop = FALSE]
          }
        } else {
          if (K > 1){
            poss.x1.d <- poss.x1[length(which(results.1[, 1:K] == 0)) == d, ,
                                 drop = FALSE]
          } else {
            poss.x1.d <- poss.x1[results.1[1, 1] == 0, , drop = FALSE]
          }
        }
        if (nrow(poss.x1.d) > 0){
          for (z2 in 0:(n + d*r*n)){
            poss.z2               <- poss.responses[[d]][[n]][[z2 + 1]]
            for (i in 1:nrow(poss.z2)){
              poss.z2[i, 2*d + 3] <- poss.z2[i, 2*d + 3]*
                                       prod(theta.pwr[1:d]^poss.z2[i, 2:(d + 1)])
            }
            poss.z2[, 2*d + 3]    <- poss.z2[, 2*d + 3]/sum(poss.z2[, 2*d + 3])
            poss.x1.d.z2          <- matrix(0, nrow = nrow(poss.x1.d)*nrow(poss.z2),
                                            ncol = d + 1)
            for (x1.d.sc in 1:nrow(poss.x1.d)){
              rem.diff                   <- which(poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] < e1[z1 + 1] &
                                                    poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] > f1)
              range                      <- (1 + nrow(poss.z2)*(x1.d.sc - 1)):(nrow(poss.z2)*x1.d.sc)
              poss.x1.d.z2[range, 1:d]   <- poss.x1.d[x1.d.sc, K + 2 + rem.diff] +
                                              poss.z2[, (d + 3):(2*d + 2)]
              poss.x1.d.z2[range, d + 1] <- poss.x1.d[x1.d.sc, 2*K + 3]*poss.z2[, 2*d + 3]
            }
            if (nrow(poss.x1.d.z2) > 1 & d > 1){
              marginal.power.p <- marginal.power.p + sum(poss.x1.d.z2[rowSums(poss.x1.d.z2[, 1:d] >= e2[[z1 + 1]][d, z2 + 1]) > 0,
                                                                      d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else if (nrow(poss.x1.d.z2) == 1 & d > 1){
              marginal.power.p <- marginal.power.p + sum(poss.x1.d.z2[any(poss.x1.d.z2[, 1:d] >= e2[[z1 + 1]][d, z2 + 1]),
                                                                      d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else if (nrow(poss.x1.d.z2) > 1 & d == 1){
              marginal.power.p <- marginal.power.p + sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1],
                                                                      d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else {
              if (poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1]){
                marginal.power.p <- marginal.power.p + poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1],
                                                                    d + 1]*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
              }
            }
          }
        }
      }
    }
    return(marginal.power.p)
  }
  
  marginal.fwer.fn <- function(p, poss.responses, K, r, n, f1, e1, e2){
    marginal.fwer.p <- 0
    p.vec            <- rep(p, K + 1)
    n.vec            <- c(n, rep(r*n, K))
    g.z              <- list()
    for (d in 1:K){
      g.z[[d]]            <- numeric(1 + n + d*r*n)
      for (z in 0:(n + d*r*n)){
        for (i in 1:nrow(poss.responses[[d]][[n]][[z + 1]])){
          g.z[[d]][z + 1] <- g.z[[d]][z + 1] + poss.responses[[d]][[n]][[z + 1]][i, 2*d + 3]*
                               prod((p.vec[1:(d + 1)]^poss.responses[[d]][[n]][[z + 1]][i, 1:(d + 1)])*
                                      ((1 - p.vec[1:(d + 1)])^(n.vec[1:(d + 1)] -
                                                                 poss.responses[[d]][[n]][[z + 1]][i, 1:(d + 1)])))
        }
      }
    }
    for (z1 in 0:(n + K*r*n)){
      poss.x1               <- poss.responses[[K]][[n]][[z1 + 1]]
      if (nrow(poss.x1) > 1){
        marginal.fwer.p.new    <- marginal.fwer.p +
                                    sum(poss.x1[poss.x1[, (K + 3):(2*K + 2)] >= e1[z1 + 1],
                                                2*K + 4])*g.z[[K]][z1 + 1]
      } else {
        marginal.fwer.p.new    <- marginal.fwer.p +
                                    sum(poss.x1[any(poss.x1[, (K + 3):(2*K + 2)] >= e1[z1 + 1]),
                                                2*K + 4])*g.z[[K]][z1 + 1]
      }
      marginal.fwer.p <- marginal.fwer.p.new
      results.1            <- matrix(0, nrow = nrow(poss.x1), ncol = K)
      for (i in 1:nrow(poss.x1)){
        results.1[i, 1:K]  <- 1*(any(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1]) |
                                   poss.x1[i, (K + 3):(2*K + 2)] <= f1)
      }
      for (d in 1:K){
        if (nrow(poss.x1) > 1){
          if (K > 1){
            poss.x1.d <- poss.x1[apply(results.1[, 1:K], 1,
                                       function(x) length(which(x == 0)) == d), ,
                                 drop = FALSE]
          } else {
            poss.x1.d <- poss.x1[results.1[, 1:K] == 0, , drop = FALSE]
          }
        } else {
          poss.x1.d <- poss.x1[length(which(results.1[, 1:K] == 0)) == d, ,
                               drop = FALSE]
        }
        if (nrow(poss.x1.d) > 0){
          for (z2 in 0:(n + d*r*n)){
            poss.z2               <- poss.responses[[d]][[n]][[z2 + 1]]
            poss.x1.d.z2          <- matrix(0, nrow = nrow(poss.x1.d)*nrow(poss.z2),
                                            ncol = d + 1)
            for (x1.d.sc in 1:nrow(poss.x1.d)){
              rem.diff                   <- which(poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] < e1[z1 + 1] &
                                                    poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] > f1)
              range                      <- (1 + nrow(poss.z2)*(x1.d.sc - 1)):(nrow(poss.z2)*x1.d.sc)
              poss.x1.d.z2[range, 1:d]   <- poss.x1.d[x1.d.sc, K + 2 + rem.diff] +
                                              poss.z2[, (d + 3):(2*d + 2)]
              poss.x1.d.z2[range, d + 1] <- poss.x1.d[x1.d.sc, 2*K + 4]*poss.z2[, 2*d + 4]
            }
            if (nrow(poss.x1.d.z2) > 1 & d > 1){
              marginal.fwer.p.new <- marginal.fwer.p + sum(poss.x1.d.z2[rowSums(poss.x1.d.z2[, 1:d] >= e2[[z1 + 1]][d, z2 + 1]) > 0,
                                                                        d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else if (nrow(poss.x1.d.z2) == 1 & d > 1){
              marginal.fwer.p.new <- marginal.fwer.p + sum(poss.x1.d.z2[any(poss.x1.d.z2[, 1:d] >= e2[[z1 + 1]][d, z2 + 1]),
                                                                        d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else if (nrow(poss.x1.d.z2) > 1 & d == 1){
              marginal.fwer.p.new <- marginal.fwer.p + sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1],
                                                                        d + 1])*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
            } else {
              if (poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1]){
                marginal.fwer.p.new <- marginal.fwer.p + poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2[[z1 + 1]][d, z2 + 1],
                                                                      d + 1]*g.z[[K]][z1 + 1]*g.z[[d]][z2 + 1]
              }
            }
            marginal.fwer.p <- marginal.fwer.p.new
          }
        }
      }
    }
    return(-marginal.fwer.p)
  }
  
  ess.H0.fn <- function(p, poss.responses, K, r, n, f1, e1){
    factors        <- rep(0, K + 1)
    total.n.vec    <- c(n + K*r*n, n + K*r*n + n + n*r*(1:K))
    p.vec          <- rep(p, K + 1)
    n.vec          <- c(n, rep(r*n, K))
    g.z            <- numeric(1 + n + K*r*n)
    for (z in 0:(n + K*r*n)){
      for (i in 1:nrow(poss.responses[[K]][[n]][[z + 1]])){
        g.z[z + 1] <- g.z[z + 1] + poss.responses[[K]][[n]][[z + 1]][i, 2*K + 3]*
                        prod((p.vec[1:(K + 1)]^poss.responses[[K]][[n]][[z + 1]][i, 1:(K + 1)])*
                               ((1 - p.vec[1:(K + 1)])^(n.vec[1:(K + 1)] -
                                                          poss.responses[[K]][[n]][[z + 1]][i, 1:(K + 1)])))
      }
    }
    for (z1 in 0:(n + K*r*n)){
      poss.x1             <- poss.responses[[K]][[n]][[z1 + 1]]
      results.1           <- matrix(0, nrow = nrow(poss.x1), ncol = K)
      for (i in 1:nrow(poss.x1)){
        results.1[i, 1:K] <- 1*(any(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1]) |
                                  poss.x1[i, (K + 3):(2*K + 2)] <= f1)
        if (all(results.1[i, 1:K] == 1)){
          factors[1]      <- factors[1] + poss.x1[i, 2*K + 4]*g.z[z1 + 1]
        } else {
          factors[length(which(results.1[i, 1:K] == 0)) + 1] <- factors[length(which(results.1[i, 1:K] == 0)) + 1] +
                                                                  poss.x1[i, 2*K + 4]*g.z[z1 + 1]
        }
      }
    }
    return(-sum(factors*total.n.vec))
  }
  
  ess.H1.fn <- function(p, poss.responses, K, r, n, f1, e1, delta){
    
    factors        <- rep(0, K + 1)
    total.n.vec    <- c(n + K*r*n, n + K*r*n + n + n*r*(1:K))
    theta.pwr      <- (p + delta)*(1 - p)/(p*(1 - p - delta))
    p.vec          <- c(p, p + delta)
    n.vec          <- c(n, rep(r*n, K))
    g.z            <- numeric(1 + n + K*r*n)
    for (z in 0:(n + K*r*n)){
      for (i in 1:nrow(poss.responses[[K]][[n]][[z + 1]])){
        g.z[z + 1] <- g.z[z + 1] + poss.responses[[K]][[n]][[z + 1]][i, 2*K + 3]*
                        prod((p.vec[1:(K + 1)]^poss.responses[[K]][[n]][[z + 1]][i, 1:(K + 1)])*
                               ((1 - p.vec[1:(K + 1)])^(n.vec[1:(K + 1)] -
                                                          poss.responses[[K]][[n]][[z + 1]][i, 1:(K + 1)])))
      }
    }
    for (z1 in 0:(n + K*r*n)){
      poss.x1               <- poss.responses[[K]][[n]][[z1 + 1]]
      for (i in 1:nrow(poss.x1)){
        poss.x1[i, 2*K + 3] <- poss.x1[i, 2*K + 3]*prod(theta.pwr^poss.x1[i, 2:(K + 1)])
      }
      poss.x1[, 2*K + 3]    <- poss.x1[, 2*K + 3]/sum(poss.x1[, 2*K + 3])
      results.1             <- matrix(0, nrow = nrow(poss.x1), ncol = K)
      for (i in 1:nrow(poss.x1)){
        results.1[i, 1:K]   <- 1*(any(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1]) | 
                                    poss.x1[i, (K + 3):(2*K + 2)] <= f1)
        if (all(results.1[i, 1:K] == 1)){
          factors[1]        <- factors[1] + poss.x1[i, 2*K + 4]*g.z[z1 + 1]
        } else {
          factors[length(which(results.1[i, 1:K] == 0)) + 1] <- factors[length(which(results.1[i, 1:K] == 0)) + 1] + 
                                                                  poss.x1[i, 2*K + 4]*g.z[z1 + 1]
        }
      }
    }
    return(-sum(factors*total.n.vec))
  }
  
  try.n <- function(n, K, r, alpha, beta, alpha.1, beta.1, delta,
                    poss.responses, nchoosek.mat, optimized){
    if (optimized == TRUE){
      f1                       <- 0
      marg.type.II.stage.1     <- -optim(0.5, f1.finder, method = "Brent", 
                                         poss.responses = poss.responses[[K]][[n]],
                                         K = K, r = r, n = n, f1 = f1,
                                         delta = delta, lower = 10^-3,
                                         upper = 1 - delta[1] - 10^-3)$value
      if (marg.type.II.stage.1 < beta.1){
        while (marg.type.II.stage.1 < beta.1){
          f1                   <- f1 + 1
          marg.type.II.stage.1 <- -optim(0.5, f1.finder, method = "Brent", 
                                         poss.responses = poss.responses[[K]][[n]],
                                         K = K, r = r, n = n, f1 = f1,
                                         delta = delta, lower = 10^-3,
                                         upper = 1 - delta[1] - 10^-3)$value
          print(f1)
        }
        f1                     <- f1 - 1
      } else {
        while (marg.type.II.stage.1 > beta.1){
          f1                   <- f1 - 1
          marg.type.II.stage.1 <- -optim(0.5, f1.finder, method = "Brent", 
                                         poss.responses = poss.responses[[K]][[n]],
                                         K = K, r = r, n = n, f1 = f1,
                                         delta = delta, lower = 10^-3,
                                         upper = 1 - delta[1] - 10^-3)$value
          print(f1)
        }
      }
    } else {
      f1 <- -1
    }
    e1 <- numeric(n + K*r*n + 1)
    e2 <- list()
    for (z1 in 0:(n + K*r*n)){
      poss.x1                       <- poss.responses[[K]][[n]][[z1 + 1]]
      if (optimized == TRUE){
        if (nrow(poss.x1) > 1){
          e1.test                     <- -n
          if (K > 1){
            cond.type.I.stage.1       <- sum(poss.x1[as.logical(ceiling(rowSums(poss.x1[, (K + 3):(2*K + 2)] >= e1.test)/K)),
                                                     2*K + 4])
          } else {
            cond.type.I.stage.1       <- sum(poss.x1[poss.x1[, K + 3] >= e1.test,
                                                     2*K + 4])
          }
          while (cond.type.I.stage.1 > alpha.1){
            e1.test                   <- e1.test + 1
            if (K > 1){
              cond.type.I.stage.1     <- sum(poss.x1[as.logical(ceiling(rowSums(poss.x1[, (K + 3):(2*K + 2)] >= e1.test)/K)),
                                                     2*K + 4])
            } else {
              cond.type.I.stage.1     <- sum(poss.x1[poss.x1[, K + 3] >= e1.test,
                                                     2*K + 4])
            }
          }
          e1.test                     <- max(e1.test, f1 + 2)
        } else {
          e1.test                     <- max(max(poss.x1[, (K + 3):(2*K + 2)]) + 1,
                                             f1 + 2)
          cond.type.I.stage.1         <- 0
        }
        e1[z1 + 1]                    <- e1.test
      } else {
        if (nrow(poss.x1) == 1){
          e1.test                       <- ceiling(n*delta[1]) + 1
          cond.type.I.stage.1           <- 0
          e1[z1 + 1]                    <- e1.test
        } else {
          e1.test                       <- ceiling(n*delta[1]) + 1
          cond.type.I.stage.1           <- sum(poss.x1[poss.x1[, K + 3] >= e1.test,
                                                       2*K + 4])
          e1[z1 + 1]                    <- e1.test
        }
      }
      results.1                     <- matrix(0, nrow = nrow(poss.x1), ncol = 2*K)
      for (i in 1:nrow(poss.x1)){
        results.1[i, 1:K]           <- 1*(any(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1]) |
                                            poss.x1[i, (K + 3):(2*K + 2)] <= f1)
        results.1[i, (K + 1):(2*K)] <- 1*(poss.x1[i, (K + 3):(2*K + 2)] >= e1[z1 + 1])
      }
      e2[[z1 + 1]]  <- matrix(0, nrow = K, ncol = 1 + n + K*r*n)
      for (d in 1:K){
        if (nrow(poss.x1) > 1){
          if (K > 1){
            poss.x1.d <- poss.x1[apply(results.1[, 1:K], 1,
                                       function(x) length(which(x == 0)) == d), ,
                                 drop = FALSE]
          } else {
            poss.x1.d <- poss.x1[results.1[, 1] == 0, , drop = FALSE]
          }
          
        } else {
          if (K > 1){
            poss.x1.d <- poss.x1[length(which(results.1[, 1:K] == 0)) == d, ,
                                 drop = FALSE]
          } else {
            poss.x1.d <- poss.x1[results.1[, 1] == 0, , drop = FALSE] 
          }
        }
        if (nrow(poss.x1.d) > 0){
          for (z2 in 0:(n + d*r*n)){
            poss.z2      <- poss.responses[[d]][[n]][[z2 + 1]]
            poss.x1.d.z2 <- matrix(0, nrow = nrow(poss.x1.d)*nrow(poss.z2),
                                   ncol = d + 1)
            for (x1.d.sc in 1:nrow(poss.x1.d)){
              rem.diff   <- which(poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] < e1[z1 + 1] &
                                    poss.x1.d[x1.d.sc, (K + 3):(2*K + 2)] > f1)
              range      <- (1 + nrow(poss.z2)*(x1.d.sc - 1)):(nrow(poss.z2)*x1.d.sc)
              poss.x1.d.z2[range, 1:d]   <- poss.x1.d[x1.d.sc, K + 2 + rem.diff] +
                poss.z2[, (d + 3):(2*d + 2)]
              poss.x1.d.z2[range, d + 1] <- poss.x1.d[x1.d.sc, 2*K + 4]*poss.z2[, 2*d + 4]
            }
            if (nrow(poss.x1.d.z2) > 1){
              if (d > 1){
                e2.test                 <- 0
                cond.type.I.stage.2     <- sum(poss.x1.d.z2[as.logical(ceiling(rowSums(poss.x1.d.z2[, 1:d] >= e2.test)/d)),
                                                            d + 1])
                if (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                  while (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                    e2.test             <- e2.test + 1
                    cond.type.I.stage.2 <- sum(poss.x1.d.z2[as.logical(ceiling(rowSums(poss.x1.d.z2[, 1:d]>=e2.test)/d)),
                                                            d + 1])
                  }
                } else {
                  e2.test               <- -2*n
                  cond.type.I.stage.2   <- sum(poss.x1.d.z2[as.logical(ceiling(rowSums(poss.x1.d.z2[, 1:d] >= e2.test)/d)),
                                                            d + 1])
                  while (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                    e2.test             <- e2.test + 1
                    cond.type.I.stage.2 <- sum(poss.x1.d.z2[as.logical(ceiling(rowSums(poss.x1.d.z2[, 1:d]>=e2.test)/d)),
                                                            d + 1])
                  }
                }
              } else {
                e2.test                 <- -1
                cond.type.I.stage.2     <- sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2.test, d + 1])
                if (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                  while (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                    e2.test             <- e2.test + 1
                    cond.type.I.stage.2 <- sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2.test, d + 1])
                  }
                } else {
                  e2.test               <- -2*n
                  cond.type.I.stage.2 <- sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2.test, d + 1])
                  while (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                    e2.test             <- e2.test + 1
                    cond.type.I.stage.2 <- sum(poss.x1.d.z2[poss.x1.d.z2[, 1] >= e2.test, d + 1])
                  }
                }
              }
            } else {
              e2.test               <- -2*n
              cond.type.I.stage.2   <- sum(poss.x1.d.z2[any(poss.x1.d.z2[, 1:d] >= e2.test), d + 1])
              if (cond.type.I.stage.2 > (alpha - cond.type.I.stage.1)/K){
                e2.test             <- max(poss.x1.d.z2[, 1:d]) + 1
                cond.type.I.stage.2 <- 0
              }
            }
            e2[[z1 + 1]][d, z2 + 1] <- e2.test
          }
        } else {
          e2[[z1 + 1]][d, ]         <- NA
        }
      }
    }
    min.marginal.power   <- marginal.power.fn(0.5 - delta[1]/2, poss.responses, K, r,
                                              n, f1, e1, e2, delta)
    if (min.marginal.power >= 1 - beta){
      min.marginal.power <- optim(0.5 - delta[1]/2, marginal.power.fn,
                                  poss.responses = poss.responses,
                                  K = K, r = r, n = n, f1 = f1, e1 = e1, e2 = e2,
                                  delta = delta, method = "Brent", lower = 10^-3,
                                  upper = 1 - delta[1] - 10^-3)$value
    }
    return(list(e1 = e1, e2 = e2, f1 = f1, min.marginal.power = min.marginal.power))
  }
  
  if (summary == TRUE){
    print("Initialising all required variables...")
  }
  
  delta <- rep(delta, K)
  
  poss.n <- 1:n.max
  if (r < 1){
    poss.n <- which((1:n.max)%%(1/r) == 0)
  } else if (r > 1) {
    poss.n <- which((r*(1:n.max))%%1 == 0)
  }
  
  nchoosek.mat                 <- matrix(0, nrow = max(poss.n, r*poss.n) + 1,
                                         ncol = max(poss.n, r*poss.n))
  for (n in 1:max(poss.n, r*poss.n)){
    nchoosek.mat[1:(n + 1), n] <- choose(n, 0:n)
  }
  
  poss.responses        <- list()
  for (d in 1:K){
    poss.responses[[d]] <- list()
    for (n in poss.n){
      all.d.z.resp      <- getall(iterpc(n = max(n + 1, r*(n + 1)),
                                         r = d + 1,
                                         labels = 0:max(n + 1, r*(n + 1)),
                                         ordered = TRUE,
                                         replace = TRUE))
      keep              <- rep(FALSE, nrow(all.d.z.resp))
      for (i in 1:nrow(all.d.z.resp)){
        if (all.d.z.resp[i, 1] %in% 0:n && all(all.d.z.resp[i, 2:(d + 1)] %in% 0:(r*n))){
          keep[i]       <- TRUE
        }
      }
      all.d.z.resp      <- all.d.z.resp[keep, ]
      poss.responses[[d]][[n]] <- list()
      row.sums          <- rowSums(all.d.z.resp)
      for (z in 0:(n + d*r*n)){
        num.d.n.z       <- nrow(all.d.z.resp[which(row.sums == z), , drop = FALSE])
        poss.responses[[d]][[n]][[z + 1]]                <- matrix(0, nrow = num.d.n.z,
                                                                   ncol = (d + 1) + 1 + d + 1 + 1)
        colnames(poss.responses[[d]][[n]][[z + 1]])      <- c(paste("x_", 0:d, sep = ""),
                                                              "z", paste("diff_", 1:d, sep = ""),
                                                              "thetaless f(x|z,theta)",
                                                              "f(x|z,1)")
        poss.responses[[d]][[n]][[z + 1]][, 1:(d + 1)]   <- all.d.z.resp[which(row.sums == z), ,
                                                                         drop = FALSE]
        poss.responses[[d]][[n]][[z + 1]][, d + 2]       <- z
        for (k in 1:d){
          poss.responses[[d]][[n]][[z + 1]][, d + 2 + k] <- poss.responses[[d]][[n]][[z + 1]][, k + 1] -
                                                              poss.responses[[d]][[n]][[z + 1]][, 1]
        }
        prod   <- nchoosek.mat[1 + poss.responses[[d]][[n]][[z + 1]][, 1], n]
        for (k in 1:d){
          prod <- prod*nchoosek.mat[1 + poss.responses[[d]][[n]][[z + 1]][, k + 1], r*n]
        }
        poss.responses[[d]][[n]][[z + 1]][, 2*d + 3] <- prod
        poss.responses[[d]][[n]][[z + 1]][, 2*d + 4] <- prod/sum(prod)
      }
      print(n)
    }
  }
  
  if (summary == TRUE){
    print("Beginning search for required sample size...")
  }
  
  n.counter            <- which(poss.n == n.init)
  n                    <- poss.n[n.counter]
  
  if (summary == TRUE){
    print(paste("Trying n = ", n, "...", sep = ""))
  }
  
  all.n      <- list()
  all.n[[n]] <- try.n(n, K, r, alpha, beta, alpha.1, beta.1, delta,
                      poss.responses, nchoosek.mat, optimized)
  min.marginal.power.n <- all.n[[n]]$min.marginal.power
  if (summary == TRUE){
    print(paste("Has a minimal marginal power of ", min.marginal.power.n, sep = ""))
  }
  if (min.marginal.power.n >= 1 - beta){
    while (min.marginal.power.n >= 1 - beta & n.counter > 1){
      n.counter        <- n.counter - 1
      n                <- poss.n[n.counter]
      if (summary == TRUE){
        print(paste("Trying n = ", n, "...", sep = ""))
      }
      all.n[[n]]       <- try.n(n, K, r, alpha, beta, alpha.1, beta.1, delta,
                                poss.responses, nchoosek.mat, optimized)
      min.marginal.power.n <- all.n[[n]]$min.marginal.power
      if (summary == TRUE){
        print(paste("Has a minimal marginal power of ", min.marginal.power.n, sep = ""))
      }
    }
    if (n.counter == 1){
      n <- poss.n[n.counter]
    } else {
      n <- poss.n[n.counter + 1]
    }
  } else {
    while (min.marginal.power.n < 1 - beta & n.counter < length(poss.n)){
      n.counter        <- n.counter + 1
      n                <- poss.n[n.counter]
      if (summary == TRUE){
        print(paste("Trying n = ", n, "...", sep = ""))
      }
      all.n[[n]]       <- try.n(n, K, r, alpha, beta, alpha.1, beta.1, delta,
                                poss.responses, nchoosek.mat, optimized)
      min.marginal.power.n <- all.n[[n]]$min.marginal.power
      if (summary == TRUE){
        print(paste("Has a minimal marginal power of ", min.marginal.power.n, sep = ""))
      }
    }
    n <- poss.n[n.counter]
  }
  if (summary == TRUE){
    print(paste("Final n = ", n, sep = ""))
    print("Computing performance of final design...")
  }
  min.marginal.power <- all.n[[n]]$min.marginal.power
  max.marginal.fwer  <- -optim(0.5, marginal.fwer.fn, poss.responses = poss.responses,
                               K = K, r = r, n = n, f1 = all.n[[n]]$f1,
                               e1 = all.n[[n]]$e1, e2 = all.n[[n]]$e2,
                               lower = 10^-3, upper = 1 - 10^-3,
                               method = "Brent")$value
  ESS.H0 <- -ess.H0.fn(pESS, poss.responses = poss.responses, K = K,
                       r = r, n = n, f1 = all.n[[n]]$f1,
                       e1 = all.n[[n]]$e1)
  ESS.H1 <- -ess.H1.fn(pESS, poss.responses = poss.responses, K = K,
                       r = r, n = n, f1 = all.n[[n]]$f1, e1 = all.n[[n]]$e1,
                       delta = delta)
  max.N  <- 2*(n + r*K*n)
  score  <- c(max.marginal.fwer, min.marginal.power, ESS.H0, ESS.H1, max.N)
  return(list(alpha = alpha, alpha.1, alpha.1, beta = beta, beta.1 = beta.1,
              K = K, delta = delta, e1 = all.n[[n]]$e1, e2 = all.n[[n]]$e2,
              f1 = all.n[[n]]$f1, n = n, nchoosek.mat = nchoosek.mat,
              optimized = optimized, poss.responses = poss.responses,
              r = r, score = score))
}

fisher_wrapper <- function(index){
  result <- binary.fisher.mams(K = K, r = r, alpha = alpha, beta = beta,
                               alpha.1 = poss.designs[index, 1],
                               beta.1 = poss.designs[index, 2],
                               delta = delta, n.max = n.max, n.init = n.init,
                               pESS = pESS, optimized = optimized, summary = TRUE)
  write.csv(result$f1, paste("fisher.f1.", index, ".csv", sep = ""))
  write.csv(result$e1, paste("fisher.e1.", index, ".csv", sep = ""))
  write.csv(result$e2, paste("fisher.e2.", index, ".csv", sep = ""))
  write.csv(result$n, paste("fisher.n.", index, ".csv", sep = ""))
  write.csv(result$score, paste("fisher.score.", index, ".csv", sep = ""))
  return(result)
}
