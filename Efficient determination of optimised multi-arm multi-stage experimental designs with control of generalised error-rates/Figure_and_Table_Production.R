################################################################################
# Function: abcd_mams_analysis.R                                               #
# Last Modified: 23/03/2018                                                    #
################################################################################

##### Load Required Packages ###################################################

library(iterpc)
library(ggplot2)

##### Internal Function Initialisation #########################################

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

##### Figure 1 #################################################################

# Generate scenarios to consider
scenarios        <- expand.grid(K = c(2, 3), J = c(2, 3),
                                popSize = c(50, 100, 200))
keep             <- rep(FALSE, nrow(scenarios))
for (i in 1:nrow(scenarios)){
  if (scenarios[i, 1] == 2 & scenarios[i, 2] == 2){
    keep[i]      <- TRUE
  } else if (scenarios[i, 1] == 2 & scenarios[i, 2] == 3){
    keep[i]      <- TRUE
  } else if (scenarios[i, 1] == 3 & scenarios[i, 2] == 2){
    keep[i]      <- TRUE
  }
}
scenarios        <- scenarios[keep, ]
seed.counter     <- 1
# Run abcd.mams() 10 times for each scenario
for (i in 1:nrow(scenarios)){
  for (j in 1:10){
    result.ij    <- abcd.mams(K = scenarios[i, 1], J = scenarios[i, 2],
                              popSize = scenarios[i, 3], maxiter = 1000,
                              seed = seed.counter, cpus= 7)
    write.csv(result.ij$GA@summary, paste("result.", i, ".", j, ".csv",
                                          sep = ""))
    seed.counter <- seed.counter + 1
  }
}
# Read the results in to a data.frame() and plot using ggplot()
data             <- NULL
scenarios.index  <- rep(c(1, 2, 3), 3)
for (i in 1:nrow(scenarios)){
  for (j in 1:10){
    data         <- rbind(data, cbind(rep(paste("Scenario ", scenarios.index[i],
                                                sep = ""), 1000),
                                      rep(scenarios[i, 3], 1000),
                                      rep(paste("Replicate ", j, sep = ""),
                                          1000), 1:1000,
                                      read.csv(paste("result.", i, ".", j,
                                                     ".csv", sep = ""))[, 2]))
  }
}
df <- data.frame(Scenario = factor(data[, 1], levels = paste("Scenario", 1:3)),
                 popSize = factor(paste("popSize =", data[, 2]),
                                  levels = c("popSize = 50", "popSize = 100", "popSize = 200")),
                 Replicate = factor(data[, 3], levels = paste("Replicate", 1:10)),
                 Iteration = as.numeric(data[, 4]),
                 ObjFn = -as.numeric(data[, 5]))
p  <- ggplot(data = subset(df, df$Iteration >= 100),
             aes(Iteration, ObjFn, colour = Replicate)) +
        facet_grid(Scenario~popSize) +
        geom_line() + ylab("Best design objective function value")+
        theme(legend.position = "none")
p

##### Tables 2-7 ###############################################################

# Generate scenarios to consider
abcd.combs                    <- getall(iterpc(n = 1:3, r = 4, ordered = T,
                                               replace = T))
keep                          <- rep(TRUE, nrow(abcd.combs))
for (i in 1:nrow(abcd.combs)){
  if (abcd.combs[i, 2] > abcd.combs[i, 3]){
    keep[i]                   <- FALSE
  }
}
abcd.combs                    <- abcd.combs[keep, ]
# Run abcd.mams() 10 times for each scenario
optimal.designs               <- list()
for (i in 1:nrow(abcd.combs)){
  optimal.designs[[i]]        <- list()
  for (j in 1:10){
    optimal.designs[[i]][[j]] <- abcd.mams(a = abcd.combs[i, 1],
                                           b = abcd.combs[i, 2],
                                           c = abcd.combs[i, 3],
                                           d = abcd.combs[i, 4],
                                           cpus = 7, seed = 1000*i + j)
    write.csv(c(optimal.designs[[i]][[j]]$n,
                optimal.designs[[i]][[j]]$a,
                optimal.designs[[i]][[j]]$r,
                optimal.designs[[i]][[j]]$a.fwer,
                optimal.designs[[i]][[j]]$bc.power,
                optimal.designs[[i]][[j]]$ESS.HG,
                optimal.designs[[i]][[j]]$ESS.LFC,
                optimal.designs[[i]][[j]]$o),
              paste("optimal.designs.", i, ".", j, ".csv", sep = ""))
  }
}
# Read the optimal designs in to a matrix
results               <- list()
for (i in 1:54){
  results[[i]]        <- matrix(0, nrow = 10, ncol = 10)
  for (j in 1:10){
    results[[i]][j, ] <- read.csv(paste("optimal.designs.", i, ".", j,
                                        ".csv", sep = ""))[, 2]
  }
}
optimals              <- matrix(0, nrow = 54, ncol = 10)
for (i in 1:54){
  optimals[i, ]       <- results[[i]][which.min(results[[i]][, 10]), ]
}
# Now compute the operating characteristics of the optimal designs
set.seed(1210)
K <- 3; J <- 2; r <- 1; sigma <- c(1, 1); delta.1 <- 0.545; delta.0 <- 0.178
situations.fwer   <- 1:3
situations.power  <- expand.grid(b = 1:3, c = 1:3)
keep              <- rep(TRUE, nrow(situations.power))
for (i in 1:nrow(situations.power)){
  if (situations.power[i, 2] > situations.power[i, 3]){
    keep[i]       <- FALSE
  }
}
situations.power  <- situations.power[keep, ]
results           <- matrix(0, nrow = 54, ncol = 13)
colnames(results) <- c("1.fwer", "2.fwer", "3.fwer", "ess.hg",
                       "ess.1", "ess.2", "ess.3", "pwr.11", "pwr.12",
                       "pwr.13", "pwr.22", "pwr.23", "pwr.33")
for (z in 1:54){
  d      <- abcd.combs[z, 4]
  e      <- optimals[z, 4:5]
  f      <- optimals[z, 2:3]
  n      <- optimals[z, 1]
  bounds <- c(f, e, -Inf, Inf)
  # a = 1, b = 1, c = 1
  theta.information   <- thetaGeneration(K, J, 1, 1, 1, d, r, sigma, delta.1,
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
  results[z, 1] <- a.fwer
  results[z, 4] <- ESS.HG
  results[z, 5] <- ESS.LFC
  results[z, 8] <- bc.power
  # a = 1, b = 1, c = 2
  theta.information   <- thetaGeneration(K, J, 1, 1, 2, d, r, sigma, delta.1,
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
  
  for (i in 1:nrow(thetas.bc.LFC.power)){
    thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                               upper = bounds[u.indices.LFC[[i]]],
                                               mean = deltas.LFC[[i]]*
                                                 sqrt(n*I.div.n.LFC[[i]]),
                                               sigma = Lambdas.LFC[[i]])[1]
  }
  bc.power      <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                         thetas.bc.LFC.power[, 2*K + 4])
  results[z, 9] <- bc.power
  # a = 1, b = 1, c = 3
  theta.information   <- thetaGeneration(K, J, 1, 1, 3, d, r, sigma, delta.1,
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
  for (i in 1:nrow(thetas.bc.LFC.power)){
    thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                               upper = bounds[u.indices.LFC[[i]]],
                                               mean = deltas.LFC[[i]]*
                                                 sqrt(n*I.div.n.LFC[[i]]),
                                               sigma = Lambdas.LFC[[i]])[1]
  }
  bc.power       <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                          thetas.bc.LFC.power[, 2*K + 4])
  results[z, 10] <- bc.power
  # a = 2, b = 2, c = 2
  theta.information   <- thetaGeneration(K, J, 2, 2, 2, d, r, sigma, delta.1,
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
  results[z, 2]  <- a.fwer
  results[z, 6]  <- ESS.LFC
  results[z, 11] <- bc.power
  # a = 2, b = 2, c = 3
  theta.information   <- thetaGeneration(K, J, 2, 2, 3, d, r, sigma, delta.1,
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
  for (i in 1:nrow(thetas.bc.LFC.power)){
    thetas.bc.LFC.power[i, 2*K + 4] <- pmvnorm(lower = bounds[l.indices.LFC[[i]]],
                                               upper = bounds[u.indices.LFC[[i]]],
                                               mean = deltas.LFC[[i]]*
                                                 sqrt(n*I.div.n.LFC[[i]]),
                                               sigma = Lambdas.LFC[[i]])[1]
  }
  bc.power       <- sum(thetas.bc.LFC.power[, 2*K + 2]*
                          thetas.bc.LFC.power[, 2*K + 4])
  results[z, 12] <- bc.power
  # a = 3, b = 3, c = 3
  theta.information   <- thetaGeneration(K, J, 3, 3, 3, d, r, sigma, delta.1,
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
  results[z, 3] <- a.fwer
  results[z, 7] <- ESS.LFC
  results[z, 13] <- bc.power
}
colnames(abcd.combs) <- letters[1:4]
tables.df <- as.data.frame(cbind(abcd.combs, results))
tables.2.and.3 <- subset(tables.df, a == 2)[, -1]
tables.4.and.5 <- subset(tables.df, a == 1)[, -1]
tables.6.and.7 <- subset(tables.df, a == 3)[, -1]





