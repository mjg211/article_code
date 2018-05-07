################################################################################
# Code:   Figure_and_Table_production.R                                        #
# Author: Michael J. Grayling (mjg211@cam.ac.uk)                               #
################################################################################

# Load required packages and code
library(ggplot2)
library(lme4)
library(MASS)
library(Matrix)
library(snowfall)
library(stats)
library(dplyr)
library(tibble)
library(pbkrtest)
source("ssre.sw.R")

################################################################################
# TDS1 simulations                                                             #
################################################################################

# Initialise required variables
X           <- matrix(c(0, 0, 0, 0, 1,
                        0, 0, 0, 1, 1,
                        0, 0, 1, 1, 1,
                        0, 1, 1, 1, 1), nrow = 4, ncol = 5)
delta       <- 0.2
sigma.mat   <- matrix(c(0.26, 0.01, 0.51, 0.02, 0.77, 0.03), nrow = 3,
                      ncol = 2, byrow = TRUE)
min.max.mat <- matrix(c(1, 1000, NA, 1000, 1, NA), nrow = 3, ncol = 2,
                      byrow = T)

# Create scenarios for Tables 1, 2, 6, 7, and S. Figures 1-6 and 17-20
sigma.e2.tilde   <- c(0.26, 0.51,  0.77)
sigma.c2.tilde   <- c(0.01, 0.02, 0.03)
tau              <- c(0, delta)
t                <- 3
blinded          <- 0:1
n.min.max        <- 1:3
assume.sigma.c.0 <- 0
scenarios.1      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                                n.min.max, assume.sigma.c.0)
scenarios.1      <- cbind(scenarios.1[, 1:5],
                          min.max.mat[scenarios.1[, 6], ],
                          rep(0, nrow(scenarios.1)),
                          rep(0, nrow(scenarios.1)), scenarios.1[, 7])
scenarios.2      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, 5, 1, 1, 1)
scenarios.2      <- cbind(scenarios.2[, 1:5],
                          min.max.mat[scenarios.2[, 6], ],
                          rep(0, nrow(scenarios.2)),
                          rep(0, nrow(scenarios.2)), scenarios.2[, 7])

# Create scenarios for Tables 3, 4, and S. Figures 7-12
sigma.scen       <- 1:3
tau              <- c(0, delta)
t                <- c(2, 4)
blinded          <- 0:1
n.min.max        <- 1
assume.sigma.c.0 <- 0
scenarios.3      <- expand.grid(sigma.scen, tau, t, blinded, n.min.max,
                                assume.sigma.c.0)
scenarios.3      <- cbind(sigma.mat[scenarios.3[, 1], ], scenarios.3[, 2:4],
                          min.max.mat[scenarios.3[, 5], ],
                          rep(0, nrow(scenarios.3)),
                          rep(0, nrow(scenarios.3)), scenarios.3[, 6])

# Create scenarios for S. Table 1
sigma.e2.tilde   <- c(0.26, 0.51,  0.77)
sigma.c2.tilde   <- 0
tau              <- c(0, delta)
t                <- 3
blinded          <- 1
n.min.max        <- 1
assume.sigma.c.0 <- 1
scenarios.4      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                                n.min.max, assume.sigma.c.0)
scenarios.4      <- cbind(scenarios.4[, 1:5], min.max.mat[scenarios.4[, 6], ],
                          rep(0, nrow(scenarios.4)), rep(0, nrow(scenarios.4)),
                          scenarios.4[, 7])

# Create scenarios for Table 5 and S. Figures 13-16
sigma.e2.tilde   <- c(0.26, 0.51,  0.77)
sigma.c2.tilde   <- c(0.01, 0.02, 0.03)
tau              <- c(0, delta)
t                <- 3
blinded          <- 1
n.min.max        <- 1
assume.sigma.c.0 <- 0
scenarios.5      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                                n.min.max, assume.sigma.c.0)
scenarios.5      <- cbind(scenarios.5[, 1:5], min.max.mat[scenarios.5[, 6], ],
                          rep(0, nrow(scenarios.5)), rep(0, nrow(scenarios.5)),
                          scenarios.5[, 7])

# Create scenarios for S. Table 2
true.sigmas <- c(0.25, 0.5, 1, 2.5, 5, 10)
deltas      <- c(0.169, 0.239, 0.338, 0.534, 0.755, 1.068)
sigma.scen  <- c(0.5, 1, 1.5)
tau         <- 0:1
blinded     <- 0:1
scenarios.6 <- as.matrix(expand.grid(true.sigmas, sigma.scen, tau, blinded))

# Join scenarios 1-5 together ready for simulation, converting to matrix
colnames(scenarios.1) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.2) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.3) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.4) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.5) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
scenarios             <- as.matrix(rbind(scenarios.1, scenarios.2, scenarios.3,
                                         scenarios.4, scenarios.5))

# Store results in the following matrix
results.matrix           <- matrix(0, ncol = 40,
                                   nrow = nrow(scenarios) + nrow(scenarios.6))
colnames(results.matrix) <- c("alpha", "beta", "delta", "sigma.e.tilde",
                              "sigma.c.tilde", "n.sw", "tau", "sigma.e",
                              "sigma.c", paste("pi", 1:5, sep = ""), "n.min",
                              "n.max", "t", "blinded", "correction",
                              "assume.sigma.c.0", "REML", "kr", "replicates",
                              "Mean Reject", "Mean N", "Mean Under.P",
                              "Mean Over.P", "Mean est.sigma.c",
                              "Mean est.sigma.e", "Mean f.est.sigma.c",
                              "Mean f.est.sigma.e", "25th N", "50th N",
                              "75th N", "25th est.sigma.c",
                              "50th est.sigma.c", "75th est.sigma.c",
                              "25th est.sigma.e", "50th est.sigma.e",
                              "75th est.sigma.e")

# Simulate all above scenarios, storing interim results in .csv files
for (i in 1:nrow(results.matrix)){
  if (i <= 156){
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios[i, 1]),
                         sigma.c.tilde = sqrt(scenarios[i, 2]),
                         tau = scenarios[i, 3], t = scenarios[i, 4],
                         blinded = as.logical(scenarios[i, 5]),
                         n.min = scenarios[i, 6], kr = F,
                         n.max = scenarios[i, 7], n.components = 200,
                         var.pi = as.logical(scenarios[i, 8]),
                         sigma.pi = as.numeric(scenarios[i, 9]),
                         assume.sigma.c.0 = as.logical(scenarios[i, 10]),
                         replicates = 100000, cpus = 7, seed = i)
  } else if (i > 156 & i <= 174) {
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios[i, 1]),
                         sigma.c.tilde = sqrt(scenarios[i, 2]),
                         tau = scenarios[i, 3], t = scenarios[i, 4],
                         blinded = as.logical(scenarios[i, 5]),
                         n.min = scenarios[i, 6], kr = F, correction = delta,
                         n.max = scenarios[i, 7], n.components = 200,
                         var.pi = as.logical(scenarios[i, 8]),
                         sigma.pi = as.numeric(scenarios[i, 9]),
                         assume.sigma.c.0 = as.logical(scenarios[i, 10]),
                         replicates = 100000, cpus = 7, seed = i)
  } else {
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios.6[i - 174, 1]*
                                                       scenarios.6[i - 174, 2]),
                         sigma.c.tilde = sqrt(scenarios.6[i - 174, 1]*
                                                scenarios.6[i - 174, 2]),
                         sigma.e = sqrt(scenarios.6[i - 174, 1]),
                         sigma.c = sqrt(scenarios.6[i - 174, 1]),
                         tau = scenarios.6[i - 174, 3]*
                           deltas[which(true.sigmas ==
                                          scenarios.6[i - 174, 1])],
                         t = 3, kr = F,
                         blinded = as.logical(scenarios.6[i - 174, 4]),
                         delta = deltas[which(true.sigmas ==
                                                scenarios.6[i - 174, 1])],
                         n.min = 1, n.max = 1000, n.components = 200,
                         var.pi = F, sigma.pi = 0, assume.sigma.c.0 = F,
                         replicates = 100000, cpus = 7, seed = i)
  }
  results.matrix[i, ] <- results.i$av.results
  write.csv(results.matrix, paste("TDS1_scenario_", i, ".csv", sep = ""))
}

################################################################################
# TDS2 simulations                                                             #
################################################################################

# Initialise required variables
X           <- matrix(c(0, 1, 1, 1, 1, 1, 1, 1, 1,
                        0, 1, 1, 1, 1, 1, 1, 1, 1,
                        0, 1, 1, 1, 1, 1, 1, 1, 1,
                        0, 0, 1, 1, 1, 1, 1, 1, 1,
                        0, 0, 1, 1, 1, 1, 1, 1, 1,
                        0, 0, 1, 1, 1, 1, 1, 1, 1,
                        0, 0, 0, 1, 1, 1, 1, 1, 1,
                        0, 0, 0, 1, 1, 1, 1, 1, 1,
                        0, 0, 0, 1, 1, 1, 1, 1, 1,
                        0, 0, 0, 0, 1, 1, 1, 1, 1,
                        0, 0, 0, 0, 1, 1, 1, 1, 1,
                        0, 0, 0, 0, 1, 1, 1, 1, 1,
                        0, 0, 0, 0, 0, 1, 1, 1, 1,
                        0, 0, 0, 0, 0, 1, 1, 1, 1,
                        0, 0, 0, 0, 0, 0, 1, 1, 1,
                        0, 0, 0, 0, 0, 0, 1, 1, 1,
                        0, 0, 0, 0, 0, 0, 0, 1, 1,
                        0, 0, 0, 0, 0, 0, 0, 1, 1,
                        0, 0, 0, 0, 0, 0, 0, 0, 1,
                        0, 0, 0, 0, 0, 0, 0, 0, 1),
                      nrow = 20, ncol = 9, byrow = TRUE)
delta       <- 0.267
min.max.mat <- matrix(c(1, 1000, NA, 1000, 1, NA), nrow = 3, ncol = 2,
                      byrow = T)
sigma.mat   <- matrix(c(0.5, 1/18, 1, 1/9, 1.5, 1/6), nrow = 3, ncol = 2,
                      byrow = TRUE)

# Create scenarios for Tables 1, 2, 6, and S. Figures 1-6 and 17-20
sigma.e2.tilde   <- c(0.5, 1, 1.5)
sigma.c2.tilde   <- c(1/18, 1/9, 1/6)
tau              <- c(0, delta)
t                <- 5
blinded          <- 0:1
n.min.max        <- 1:3
assume.sigma.c.0 <- 0
scenarios.1      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                                n.min.max, assume.sigma.c.0)
scenarios.1      <- cbind(scenarios.1[, 1:5],
                          min.max.mat[scenarios.1[, 6], ],
                          rep(0, nrow(scenarios.1)),
                          rep(0, nrow(scenarios.1)), scenarios.1[, 7])
scenarios.2      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, 9, 1, 1, 1)
scenarios.2      <- cbind(scenarios.2[, 1:5],
                          min.max.mat[scenarios.2[, 6], ],
                          rep(0, nrow(scenarios.2)),
                          rep(0, nrow(scenarios.2)), scenarios.2[, 7])

# Create scenarios for Tables 3, 4, and S. Figures 7-12
sigma.scen       <- 1:3
tau              <- c(0, delta)
t                <- c(3, 7)
blinded          <- 0:1
n.min.max        <- 1
assume.sigma.c.0 <- 0
scenarios.3      <- expand.grid(sigma.scen, tau, t, blinded, n.min.max,
                                assume.sigma.c.0)
scenarios.3      <- cbind(sigma.mat[scenarios.3[, 1], ], scenarios.3[, 2:4],
                          min.max.mat[scenarios.3[, 5], ],
                          rep(0, nrow(scenarios.3)),
                          rep(0, nrow(scenarios.3)), scenarios.3[, 6])

# Create scenarios for S. Table 1
sigma.e2.tilde   <- c(0.5, 1, 1.5)
sigma.c2.tilde   <- 0
tau              <- c(0, delta)
t                <- 5
blinded          <- 1
n.min.max        <- 1
assume.sigma.c.0 <- 1
scenarios.4      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                               n.min.max, assume.sigma.c.0)
scenarios.4      <- cbind(scenarios.4[, 1:5], min.max.mat[scenarios.4[, 6], ],
                          rep(0, nrow(scenarios.4)), rep(0, nrow(scenarios.4)),
                          scenarios.4[, 7])

# Create scenarios for Table 5 and S. Fig 13-16
sigma.e2.tilde   <- c(0.5, 1, 1.5)
sigma.c2.tilde   <- c(1/18, 1/9, 1/6)
tau              <- c(0, delta)
t                <- 5
blinded          <- 1
n.min.max        <- 1
assume.sigma.c.0 <- 0
scenarios.5      <- expand.grid(sigma.e2.tilde, sigma.c2.tilde, tau, t, blinded,
                                n.min.max, assume.sigma.c.0)
scenarios.5      <- cbind(scenarios.5[, 1:5], min.max.mat[scenarios.5[, 6], ],
                          rep(0, nrow(scenarios.5)), rep(0, nrow(scenarios.5)),
                          scenarios.5[, 7])

# Create scenarios for S. Table 2
true.sigmas <- c(0.25, 0.5, 1, 2.5, 5, 10)
deltas      <- c(0.117, 0.165, 0.233, 0.368, 0.521, 0.736)
sigma.scen  <- c(0.5, 1, 1.5)
tau         <- 0:1
blinded     <- 0:1
scenarios.6 <- as.matrix(expand.grid(true.sigmas, sigma.scen, tau, blinded))

# Join all above scenarios together ready for simulation, converting to matrix
colnames(scenarios.1) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.2) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.3) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.4) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
colnames(scenarios.5) <- c("sigma.e2.tilde", "sigma.c2.tilde", "tau", "t",
                           "blinded", "n.min", "n.max", "var.pi",
                           "sigma.pi", "assume.sigma.c.0")
scenarios             <- as.matrix(rbind(scenarios.1, scenarios.2, scenarios.3,
                                         scenarios.4, scenarios.5))

# Store results in the following matrix
results.matrix           <- matrix(0, ncol = 44, nrow = nrow(scenarios) +
                                                          nrow(scenarios.6))
colnames(results.matrix) <- c("alpha", "beta", "delta",  "sigma.e.tilde",
                              "sigma.c.tilde", "n.sw", "tau", "sigma.e",
                              "sigma.c", paste("pi", 1:9, sep = ""), "n.min",
                              "n.max", "t", "blinded", "correction",
                              "assume.sigma.c.0", "REML", "kr", "replicates",
                              "Mean Reject", "Mean N", "Mean Under.P",
                              "Mean Over.P", "Mean est.sigma.c",
                              "Mean est.sigma.e", "Mean f.est.sigma.c",
                              "Mean f.est.sigma.e", "25th N", "50th N",
                              "75th N", "25th est.sigma.c",
                              "50th est.sigma.c", "75th est.sigma.c",
                              "25th est.sigma.e", "50th est.sigma.e",
                              "75th est.sigma.e")

# Simulate all scenarios, storing interim results in .csv files
for (i in 1:nrow(results.matrix)){
  if (i <= 156){
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios[i, 1]),
                         sigma.c.tilde = sqrt(scenarios[i, 2]),
                         delta = 0.267, sigma.e = sqrt(1),
                         sigma.c = sqrt(1/9), pi = rep(0, 9),
                         alpha = 0.025, beta = 0.2, n.min = scenarios[i, 6],
                         n.max = scenarios[i, 7], n.components = 100,
                         tau = scenarios[i, 3], t = scenarios[i, 4],
                         blinded = as.logical(scenarios[i, 5]),
                         var.pi = as.logical(scenarios[i, 8]),
                         sigma.pi = as.numeric(scenarios[i, 9]),
                         assume.sigma.c.0 = as.logical(scenarios[i, 10]),
                         replicates = 100000, cpus = 7, seed = 1000 + i)
  } else if (i > 156 & i <= 174) {
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios[i, 1]),
                         sigma.c.tilde = sqrt(scenarios[i, 2]),
                         delta = 0.267, sigma.e = 1, sigma.c = sqrt(1/9),
                         pi = rep(0,9), alpha = 0.025, beta = 0.2,
                         tau = scenarios[i, 3], t = scenarios[i, 4],
                         blinded = as.logical(scenarios[i, 5]),
                         n.min = scenarios[i, 6], kr = F, correction = delta,
                         n.max = scenarios[i, 7], n.components = 100,
                         var.pi = as.logical(scenarios[i, 8]),
                         sigma.pi = as.numeric(scenarios[i, 9]),
                         assume.sigma.c.0 = as.logical(scenarios[i, 10]),
                         replicates = 100000, cpus = 7, seed = 1000 + i)
  } else {
    results.i <- ssre.sw(X = X, sigma.e.tilde = sqrt(scenarios.6[i - 174, 1]*
                                                       scenarios.6[i - 174, 2]),
                         sigma.c.tilde = sqrt(scenarios.6[i - 174, 1]*
                                                scenarios.6[i - 174, 2]),
                         sigma.e = sqrt(scenarios.6[i - 174, 1]),
                         sigma.c = sqrt(scenarios.6[i - 174, 1]),
                         delta = deltas[which(true.sigmas ==
                                                scenarios.6[i - 174, 1])],
                         pi = rep(0, 9), alpha = 0.025,  beta = 0.2, n.min = 1,
                         n.max = 1000, n.components = 100,
                         tau = scenarios.6[i - 174, 3]*
                                 deltas[which(true.sigmas ==
                                                scenarios.6[i - 174, 1])],
                         t = 5, blinded = as.logical(scenarios.6[i - 174, 4]),
                         var.pi = F, sigma.pi = 0, assume.sigma.c.0 = F,
                         replicates = 100000, cpus = 7, seed = 1000 + i)
  }
  results.matrix[i, ] <- results.i$av.results
  write.csv(results.matrix, paste("TDS2_scenario_", i, ".csv", sep = ""))
}

################################################################################
# Load in final data                                                           #
################################################################################

# Load in final data (adjust working directory)
TDS1 <- read.csv("TDS1_scenario_246.csv")
TDS2 <- read.csv("TDS2_scenario_246.csv")

################################################################################
# Table 1                                                                      #
################################################################################

# TDS1 part of Table 1
sigma.c2.tilde         <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde         <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
table.1.TDS1           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 6))
table.1.TDS1[, 1:2]    <- table.1.TDS1[, c(2, 1)]
colnames(table.1.TDS1) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "ETI: Blinded", "ETI: Unblinded", "ETI: Fixed",
                            "EP: Blinded", "EP: Unblinded", "EP: Fixed")
table.1.TDS1[, 3]      <- subset(TDS1[1:36, ], tau == 0 &
                                   blinded == 1)$Mean.Reject
table.1.TDS1[, 4]      <- subset(TDS1[1:36, ], tau == 0 &
                                   blinded == 0)$Mean.Reject
table.1.TDS1[, 5]      <- subset(TDS1[109:126, ], tau == 0)$Mean.Reject
table.1.TDS1[, 6]      <- subset(TDS1[1:36, ], tau == 0.2 &
                                   blinded == 1)$Mean.Reject
table.1.TDS1[, 7]      <- subset(TDS1[1:36, ], tau == 0.2 &
                                   blinded == 0)$Mean.Reject
table.1.TDS1[, 8]      <- subset(TDS1[109:126, ], tau == 0.2)$Mean.Reject

# TDS2 part of Table 1
table.1.TDS2           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 6))
table.1.TDS2[, 1:2]    <- table.1.TDS2[, c(2, 1)]
colnames(table.1.TDS2) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "ETI: Blinded", "ETI: Unblinded", "ETI: Fixed",
                            "EP: Blinded", "EP: Unblinded", "EP: Fixed")
table.1.TDS2[, 3]      <- subset(TDS2[1:36, ], tau == 0 &
                                   blinded == 1)$Mean.Reject
table.1.TDS2[, 4]      <- subset(TDS2[1:36, ], tau == 0 &
                                   blinded == 0)$Mean.Reject
table.1.TDS2[, 5]      <- subset(TDS2[109:126, ], tau == 0)$Mean.Reject
table.1.TDS2[, 6]      <- subset(TDS2[1:36, ], tau == 0.267 &
                                   blinded == 1)$Mean.Reject
table.1.TDS2[, 7]      <- subset(TDS2[1:36, ], tau == 0.267 &
                                   blinded == 0)$Mean.Reject
table.1.TDS2[, 8]      <- subset(TDS2[109:126, ], tau == 0.267)$Mean.Reject

# Final Table 1
table.1                <- rbind(table.1.TDS1, table.1.TDS2)

################################################################################
# Table 2                                                                      #
################################################################################

# TDS1 part of Table 2
sigma.c2.tilde         <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde         <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
table.2.TDS1           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 6))
table.2.TDS1[, 1:2]    <- table.2.TDS1[, c(2, 1)]
colnames(table.2.TDS1) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "Med(N|0) Blinded", "Med(N|0) Unblinded",
                            "Med(N|0) Fixed", "Med(N|delta) Blinded",
                            "Med(N|delta) Unblinded", "Med(N|delta) Fixed")
table.2.TDS1[, 3]      <- subset(TDS1[1:36, ], tau == 0 & blinded == 1)$X50th.N
table.2.TDS1[, 4]      <- subset(TDS1[1:36, ], tau == 0 & blinded == 0)$X50th.N
table.2.TDS1[, 5]      <- subset(TDS1[109:126, ], tau == 0)$X50th.N
table.2.TDS1[, 6]      <- subset(TDS1[1:36, ], tau == 0.2 &
                                   blinded == 1)$X50th.N
table.2.TDS1[, 7]      <- subset(TDS1[1:36, ], tau == 0.2 &
                                   blinded == 0)$X50th.N
table.2.TDS1[, 8]      <- subset(TDS1[109:126, ], tau == 0.2)$X50th.N

# TDS2 part of Table 2
table.2.TDS2           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 6))
table.2.TDS2[, 1:2]    <- table.2.TDS2[, c(2, 1)]
colnames(table.2.TDS2) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "Med(N|0) Blinded", "Med(N|0) Unblinded",
                            "Med(N|0) Fixed", "Med(N|delta) Blinded",
                            "Med(N|delta) Unblinded", "Med(N|delta) Fixed")
table.2.TDS2[, 3]      <- subset(TDS2[1:36, ], tau == 0 & blinded == 1)$X50th.N
table.2.TDS2[, 4]      <- subset(TDS2[1:36, ], tau == 0 & blinded == 0)$X50th.N
table.2.TDS2[, 5]      <- subset(TDS2[109:126, ], tau == 0)$X50th.N
table.2.TDS2[, 6]      <- subset(TDS2[1:36, ], tau == 0.267 &
                                   blinded == 1)$X50th.N
table.2.TDS2[, 7]      <- subset(TDS2[1:36, ], tau == 0.267 &
                                   blinded == 0)$X50th.N
table.2.TDS2[, 8]      <- subset(TDS2[109:126, ], tau == 0.267)$X50th.N

# Final Table 2
table.2                <- rbind(table.2.TDS1, table.2.TDS2)

################################################################################
# S. Figures 1-6                                                               #
################################################################################

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = TDS1[1:36, ]$tau,
                         method = factor(rep(c("Unblinded", "Blinded"),
                                             each = 18),
                                         levels = c("Unblinded", "Blinded")),
                         sigma_c2_25th = TDS1[1:36, ]$X25th.est.sigma.c^2,
                         sigma_c2_50th = TDS1[1:36, ]$X50th.est.sigma.c^2,
                         sigma_c2_75th = TDS1[1:36, ]$X75th.est.sigma.c^2,
                         sigma_e2_25th = TDS1[1:36, ]$X25th.est.sigma.e^2,
                         sigma_e2_50th = TDS1[1:36, ]$X50th.est.sigma.e^2,
                         sigma_e2_75th = TDS1[1:36, ]$X75th.est.sigma.e^2,
                         N_25th = TDS1[1:36, ]$X25th.N,
                         N_50th = TDS1[1:36, ]$X50th.N,
                         N_75th = TDS1[1:36, ]$X75th.N)

s.fig.1.a <- ggplot(data = subset(df, tau == 0),
                    aes(x = sigmas, y = sigma_c2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[c]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_c2_25th,
                                 ymax = sigma_c2_75th), width = .25) +
               geom_hline(yintercept = 0.02, linetype = 2)
ggsave("sfig1a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.1.b <- ggplot(data = subset(df, tau == 0),
                    aes(x = sigmas, y = sigma_e2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[e]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_e2_25th,
                                 ymax = sigma_e2_75th), width = .25) +
               geom_hline(yintercept = 0.51, linetype = 2)
ggsave("sfig1b.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.2.a <- ggplot(data = subset(df, tau == 0.2),
                    aes(x = sigmas, y = sigma_c2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[c]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_c2_25th,
                                 ymax = sigma_c2_75th), width = .25) +
               geom_hline(yintercept = 0.02, linetype = 2)
ggsave("sfig2a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.2.b <- ggplot(data = subset(df, tau == 0.2),
                    aes(x = sigmas, y = sigma_e2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[e]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_e2_25th,
                                 ymax = sigma_e2_75th), width = .25) +
               geom_hline(yintercept = 0.51, linetype = 2)
ggsave("sfig2b.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.3.a <- ggplot(data = subset(df, tau == 0), aes(x = sigmas, y = N_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(italic(N)))) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig3a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.3.b <- ggplot(data = subset(df, tau == 0.2),
                    aes(x = sigmas, y = N_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(italic(N)))) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig3b.pdf", device = "pdf", width = 7, height = 5, units = "in")

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = TDS2[1:36, ]$tau,
                         method = factor(rep(c("Unblinded", "Blinded"),
                                             each = 18),
                                         levels = c("Unblinded", "Blinded")),
                         sigma_c2_25th = TDS2[1:36, ]$X25th.est.sigma.c^2,
                         sigma_c2_50th = TDS2[1:36, ]$X50th.est.sigma.c^2,
                         sigma_c2_75th = TDS2[1:36, ]$X75th.est.sigma.c^2,
                         sigma_e2_25th = TDS2[1:36, ]$X25th.est.sigma.e^2,
                         sigma_e2_50th = TDS2[1:36, ]$X50th.est.sigma.e^2,
                         sigma_e2_75th = TDS2[1:36, ]$X75th.est.sigma.e^2,
                         N_25th = TDS2[1:36, ]$X25th.N,
                         N_50th = TDS2[1:36, ]$X50th.N,
                         N_75th = TDS2[1:36, ]$X75th.N)

s.fig.4.a <- ggplot(data = subset(df, tau == 0),
                    aes(x = sigmas, y = sigma_c2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[c]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_c2_25th,
                                 ymax = sigma_c2_75th), width = .25) +
               geom_hline(yintercept = 1/9, linetype = 2)
ggsave("sfig4a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.4.b <- ggplot(data = subset(df, tau == 0),
                    aes(x = sigmas, y = sigma_e2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[e]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_e2_25th,
                                 ymax = sigma_e2_75th), width = .25) +
               geom_hline(yintercept = 1, linetype = 2)
ggsave("sfig4b.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.5.a <- ggplot(data = subset(df, tau == 0.267),
                    aes(x = sigmas, y = sigma_c2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[c]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_c2_25th,
                                 ymax = sigma_c2_75th), width = .25) +
               geom_hline(yintercept = 1/9, linetype = 2)
ggsave("sfig5a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.5.b <- ggplot(data = subset(df, tau == 0.267),
                    aes(x = sigmas, y = sigma_e2_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(sigma)[e]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_e2_25th,
                                 ymax = sigma_e2_75th), width = .25) +
               geom_hline(yintercept = 1, linetype = 2)
ggsave("sfig5b.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.6.a <- ggplot(data = subset(df, tau == 0),
                    aes(x = sigmas, y = N_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(italic(N)))) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = .25)
ggsave("sfig6a.pdf", device = "pdf", width = 7, height = 5, units = "in")

s.fig.6.b <- ggplot(data = subset(df, tau == 0.267),
                    aes(x = sigmas, y = N_50th)) +
               geom_point() + facet_grid(.~method) +
               xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                                     tilde(sigma)[e]^2, ")", sep = ""))) +
               ylab(expression(hat(italic(N)))) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = .25)
ggsave("sfig6b.pdf", device = "pdf", width = 7, height = 5, units = "in")

################################################################################
# Table 3                                                                      #
################################################################################

# TDS1 part of Table 3
table.3.TDS1           <- cbind(c(rep("Blinded", 3), rep("Unblinded", 3)),
                                rep(c("0.5*sigma.c2", "sigma.c2",
                                      "1.5*sigma.c2"), 2),
                                rep(c("0.5*sigma.e2", "sigma.e2",
                                      "1.5*sigma.e2"), 2),
                                matrix(0, nrow = 6, ncol = 6))
colnames(table.3.TDS1) <- c("Procedure", "sigma.c2.tilde", "sigma.e2.tilde",
                            "ETI: t = 2", "ETI: t = 3", "ETI: t = 4",
                            "EP: t = 2", "EP: t = 3", "EP: t = 4")
table.3.TDS1[1:3, 4]   <- subset(TDS1[127:150, ], t == 2 & tau == 0 &
                                   blinded == 1)$Mean.Reject
table.3.TDS1[4:6, 4]   <- subset(TDS1[127:150, ], t == 2 & tau == 0 &
                                   blinded == 0)$Mean.Reject
table.3.TDS1[1:3, 5]   <- table.1.TDS1[c(1, 5, 9), 3]
table.3.TDS1[4:6, 5]   <- table.1.TDS1[c(1, 5, 9), 4]
table.3.TDS1[1:3, 6]   <- subset(TDS1[127:150, ], t == 4 & tau == 0 &
                                   blinded == 1)$Mean.Reject
table.3.TDS1[4:6, 6]   <- subset(TDS1[127:150, ], t == 4 & tau == 0 &
                                   blinded == 0)$Mean.Reject
table.3.TDS1[1:3, 7]   <- subset(TDS1[127:150, ], t == 2 & tau == 0.2 &
                                   blinded == 1)$Mean.Reject
table.3.TDS1[4:6, 7]   <- subset(TDS1[127:150, ], t == 2 & tau == 0.2 &
                                   blinded == 0)$Mean.Reject
table.3.TDS1[1:3, 8]   <- table.1.TDS1[c(1, 5, 9), 6]
table.3.TDS1[4:6, 8]   <- table.1.TDS1[c(1, 5, 9), 7]
table.3.TDS1[1:3, 9]   <- subset(TDS1[127:150, ], t == 4 & tau == 0.2 &
                                   blinded == 1)$Mean.Reject
table.3.TDS1[4:6, 9]   <- subset(TDS1[127:150, ], t == 4 & tau == 0.2 &
                                   blinded == 0)$Mean.Reject
table.3.TDS1           <- as.data.frame(table.3.TDS1)

# TDS2 part of Table 3
table.3.TDS2           <- cbind(c(rep("Blinded", 3), rep("Unblinded", 3)),
                                rep(c("0.5*sigma.c2", "sigma.c2",
                                      "1.5*sigma.c2"), 2),
                                rep(c("0.5*sigma.e2", "sigma.e2",
                                      "1.5*sigma.e2"), 2),
                                matrix(0, nrow = 6, ncol = 6))
colnames(table.3.TDS2) <- c("Procedure", "sigma.c2.tilde", "sigma.e2.tilde",
                            "ETI t = 3", "ETI t = 5", "ETI t = 7",
                            "EP t = 3", "EP t = 5", "EP t = 7")
table.3.TDS2[1:3, 4]   <- subset(TDS2[127:150, ], t == 3 & tau == 0 &
                                   blinded == 1)$Mean.Reject
table.3.TDS2[4:6, 4]   <- subset(TDS2[127:150, ], t == 3 & tau == 0 &
                                   blinded == 0)$Mean.Reject
table.3.TDS2[1:3, 5]   <- table.1.TDS2[c(1, 5, 9), 3]
table.3.TDS2[4:6, 5]   <- table.1.TDS2[c(1, 5, 9), 4]
table.3.TDS2[1:3, 6]   <- subset(TDS2[127:150, ], t == 7 & tau == 0 &
                                   blinded == 1)$Mean.Reject
table.3.TDS2[4:6, 6]   <- subset(TDS2[127:150, ], t == 7 & tau == 0 &
                                   blinded == 0)$Mean.Reject
table.3.TDS2[1:3, 7]   <- subset(TDS2[127:150, ], t == 3 & tau == 0.267 &
                                   blinded == 1)$Mean.Reject
table.3.TDS2[4:6, 7]   <- subset(TDS2[127:150, ], t == 3 & tau == 0.267 &
                                   blinded == 0)$Mean.Reject
table.3.TDS2[1:3, 8]   <- table.1.TDS2[c(1, 5, 9), 6]
table.3.TDS2[4:6, 8]   <- table.1.TDS2[c(1, 5, 9), 7]
table.3.TDS2[1:3, 9]   <- subset(TDS2[127:150, ], t == 7 & tau == 0.267 &
                                   blinded == 1)$Mean.Reject
table.3.TDS2[4:6, 9]   <- subset(TDS2[127:150, ], t == 7 & tau == 0.267 &
                                   blinded == 0)$Mean.Reject
table.3.TDS2           <- as.data.frame(table.3.TDS2)

# Final Table 3 is the rbind() of the above (different column names means this
# cannot be done explicitly)

################################################################################
# Table 4                                                                      #
################################################################################

# TDS1 part of Table 4
table.4.TDS1           <- cbind(c(rep("Blinded", 3), rep("Unblinded", 3)),
                                rep(c("0.5*sigma.c2", "sigma.c2",
                                      "1.5*sigma.c2"), 2),
                                rep(c("0.5*sigma.e2", "sigma.e2",
                                      "1.5*sigma.e2"), 2),
                                matrix(0, nrow = 6, ncol = 6))
colnames(table.4.TDS1) <- c("Procedure", "sigma.c2.tilde", "sigma.e2.tilde",
                            "Med(N|0): t = 2", "Med(N|0): t = 3",
                            "Med(N|0): t = 4", "Med(N|delta): t = 2",
                            "Med(N|delta): t = 3", "Med(N|delta): t = 4")
table.4.TDS1[1:3, 4]   <- subset(TDS1[127:150, ], t == 2 & tau == 0 &
                                   blinded == 1)$X50th.N
table.4.TDS1[4:6, 4]   <- subset(TDS1[127:150, ], t == 2 & tau == 0 &
                                   blinded == 0)$X50th.N
table.4.TDS1[1:3, 5]   <- table.2.TDS1[c(1, 5, 9), 3]
table.4.TDS1[4:6, 5]   <- table.2.TDS1[c(1, 5, 9), 4]
table.4.TDS1[1:3, 6]   <- subset(TDS1[127:150, ], t == 4 & tau == 0 &
                                   blinded == 1)$X50th.N
table.4.TDS1[4:6, 6]   <- subset(TDS1[127:150, ], t == 4 & tau == 0 &
                                   blinded == 0)$X50th.N
table.4.TDS1[1:3, 7]   <- subset(TDS1[127:150, ], t == 2 & tau == 0.2 &
                                   blinded == 1)$X50th.N
table.4.TDS1[4:6, 7]   <- subset(TDS1[127:150, ], t == 2 & tau == 0.2 &
                                   blinded == 0)$X50th.N
table.4.TDS1[1:3, 8]   <- table.2.TDS1[c(1, 5, 9), 6]
table.4.TDS1[4:6, 8]   <- table.2.TDS1[c(1, 5, 9), 7]
table.4.TDS1[1:3, 9]   <- subset(TDS1[127:150, ], t == 4 & tau == 0.2 &
                                   blinded == 1)$X50th.N
table.4.TDS1[4:6, 9]   <- subset(TDS1[127:150, ], t == 4 & tau == 0.2 &
                                   blinded == 0)$X50th.N
table.4.TDS1           <- as.data.frame(table.4.TDS1)

# TDS2 part of Table 4
table.4.TDS2           <- cbind(c(rep("Blinded", 3), rep("Unblinded", 3)),
                                rep(c("0.5*sigma.c2", "sigma.c2",
                                      "1.5*sigma.c2"), 2),
                                rep(c("0.5*sigma.e2", "sigma.e2",
                                      "1.5*sigma.e2"), 2),
                                matrix(0, nrow = 6, ncol = 6))
colnames(table.4.TDS2) <- c("Procedure", "sigma.c2.tilde", "sigma.e2.tilde",
                            "Med(N|0): t = 3", "Med(N|0): t = 5",
                            "Med(N|0): t = 7", "Med(N|delta): t = 3",
                            "Med(N|delta): t = 5", "Med(N|delta): t = 7")
table.4.TDS2[1:3, 4]   <- subset(TDS2[127:150, ], t == 3 & tau == 0 &
                                   blinded == 1)$X50th.N
table.4.TDS2[4:6, 4]   <- subset(TDS2[127:150, ], t == 3 & tau == 0 &
                                   blinded == 0)$X50th.N
table.4.TDS2[1:3, 5]   <- table.2.TDS2[c(1, 5, 9), 3]
table.4.TDS2[4:6, 5]   <- table.2.TDS2[c(1, 5, 9), 4]
table.4.TDS2[1:3, 6]   <- subset(TDS2[127:150, ], t == 7 & tau == 0 &
                                   blinded == 1)$X50th.N
table.4.TDS2[4:6, 6]   <- subset(TDS2[127:150, ], t == 7 & tau == 0 &
                                   blinded == 0)$X50th.N
table.4.TDS2[1:3, 7]   <- subset(TDS2[127:150, ], t == 3 & tau == 0.267 &
                                   blinded == 1)$X50th.N
table.4.TDS2[4:6, 7]   <- subset(TDS2[127:150, ], t == 3 & tau == 0.267 &
                                   blinded == 0)$X50th.N
table.4.TDS2[1:3, 8]   <- table.2.TDS2[c(1, 5, 9), 6]
table.4.TDS2[4:6, 8]   <- table.2.TDS2[c(1, 5, 9), 7]
table.4.TDS2[1:3, 9]   <- subset(TDS2[127:150, ], t == 7 & tau == 0.267 &
                                   blinded == 1)$X50th.N
table.4.TDS2[4:6, 9]   <- subset(TDS2[127:150, ], t == 7 & tau == 0.267 &
                                   blinded == 0)$X50th.N
table.4.TDS2           <- as.data.frame(table.4.TDS2)

# Final Table 4 is the rbind() of the above (different column names means this
# cannot be done explicitly)

################################################################################
# S. Figures 7-12                                                              #
################################################################################

TDS1_sub <- TDS1[c(1, 5, 9, 10, 14, 18, 19, 23, 27, 28, 32, 36, 127:150), ]
sigmas   <- c(rep(c("0.5(sigma_c^2,sigma_e^2)", "(sigma_c^2,sigma_e^2)",
                    "1.5(sigma_c^2,sigma_e^2)"), 12))
df        <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                    tau = factor(rep(rep(c("tau = 0", "tau = delta"),
                                         each = 3), 6),
                    levels = c("tau = 0", "tau = delta")),
                    t = factor(paste("t =", TDS1_sub$t),
                               levels = c("t = 2", "t = 3", "t = 4")),
                    method = factor(c(rep(c("Unblinded", "Blinded"), each = 6),
                                      rep(c("Unblinded", "Blinded"),
                                          each = 12)),
                                    levels = c("Unblinded", "Blinded")),
                    sigma_c2_25th = TDS1_sub$X25th.est.sigma.c^2,
                    sigma_c2_50th = TDS1_sub$X50th.est.sigma.c^2,
                    sigma_c2_75th = TDS1_sub$X75th.est.sigma.c^2,
                    sigma_e2_25th = TDS1_sub$X25th.est.sigma.e^2,
                    sigma_e2_50th = TDS1_sub$X50th.est.sigma.e^2,
                    sigma_e2_75th = TDS1_sub$X75th.est.sigma.e^2,
                    N_25th = TDS1_sub$X25th.N,
                    N_50th = TDS1_sub$X50th.N,
                    N_75th = TDS1_sub$X75th.N)

s.fig.7.a <- ggplot(data = subset(df, tau == "tau = 0"),
                    aes(x = t, y = sigma_c2_50th)) +
               geom_point() + facet_grid(method~sigmas) +
               xlab(expression(italic(t))) + ylab(expression(hat(sigma)[c]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_c2_25th,
                                 ymax = sigma_c2_75th), width = .25) +
               geom_hline(yintercept = 0.02, linetype = 2)
ggsave("sfig7a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.7.b <- ggplot(data = subset(df, tau == "tau = 0"),
                    aes(x = t, y = sigma_e2_50th)) +
               geom_point() + facet_grid(method~sigmas) +
               xlab(expression(italic(t))) + ylab(expression(hat(sigma)[e]^2)) +
               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
               geom_errorbar(aes(ymin = sigma_e2_25th,
                                 ymax = sigma_e2_75th), width = .25) +
               geom_hline(yintercept = 0.51, linetype = 2)
ggsave("sfig7b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.8.a <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = sigma_c2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[c]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_c2_25th,
                                  ymax = sigma_c2_75th), width = .25) +
                geom_hline(yintercept = 0.02, linetype = 2)
ggsave("sfig8a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.8.b <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = sigma_e2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[e]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_e2_25th,
                                  ymax = sigma_e2_75th), width = .25) +
                geom_hline(yintercept = 0.51, linetype = 2)
ggsave("sfig8b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.9.a <- ggplot(data = subset(df, tau == "tau = 0"),
                     aes(x = t, y = N_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(italic(N)))) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = N_25th,
                                  ymax = N_75th), width = .25)
ggsave("sfig9a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.9.b <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = N_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(italic(N)))) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = N_25th,
                                  ymax = N_75th), width = .25)
ggsave("sfig9b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

TDS2_sub <- TDS2[c(1, 5, 9, 10, 14, 18, 19, 23, 27, 28, 32, 36, 127:150), ]
sigmas   <- c(rep(c("0.5(sigma_c^2,sigma_e^2)", "(sigma_c^2,sigma_e^2)",
                    "1.5(sigma_c^2,sigma_e^2)"), 12))
df        <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                    tau = factor(rep(rep(c("tau = 0", "tau = delta"),
                                         each = 3), 6),
                                 levels = c("tau = 0", "tau = delta")),
                    t = factor(paste("t =", TDS2_sub$t),
                               levels = c("t = 3", "t = 5", "t = 7")),
                    method = factor(c(rep(c("Unblinded", "Blinded"), each = 6),
                                      rep(c("Unblinded", "Blinded"),
                                          each = 12)),
                                    levels = c("Unblinded", "Blinded")),
                    sigma_c2_25th = TDS2_sub$X25th.est.sigma.c^2,
                    sigma_c2_50th = TDS2_sub$X50th.est.sigma.c^2,
                    sigma_c2_75th = TDS2_sub$X75th.est.sigma.e^2,
                    sigma_e2_25th = TDS2_sub$X25th.est.sigma.e^2,
                    sigma_e2_50th = TDS2_sub$X50th.est.sigma.e^2,
                    sigma_e2_75th = TDS2_sub$X75th.est.sigma.e.1^2,
                    N_25th = TDS2_sub$X25th.N,
                    N_50th = TDS2_sub$X50th.N,
                    N_75th = TDS2_sub$X75th.N)

s.fig.10.a <- ggplot(data = subset(df, tau == "tau = 0"),
                     aes(x = t, y = sigma_c2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[c]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_c2_25th,
                                  ymax = sigma_c2_75th), width = .25) +
                geom_hline(yintercept = 1/9, linetype = 2)
ggsave("sfig10a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.10.b <- ggplot(data = subset(df, tau == "tau = 0"),
                    aes(x = t, y = sigma_e2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[e]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_e2_25th,
                                  ymax = sigma_e2_75th), width = .25) +
                geom_hline(yintercept = 1, linetype = 2)
ggsave("sfig10b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.11.a <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = sigma_c2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[c]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_c2_25th,
                                  ymax = sigma_c2_75th), width = .25) +
                geom_hline(yintercept = 1/9, linetype = 2)
ggsave("sfig11a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.11.b <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = sigma_e2_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(sigma)[e]^2)) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = sigma_e2_25th,
                                  ymax = sigma_e2_75th), width = .25) +
                geom_hline(yintercept = 1, linetype = 2)
ggsave("sfig11b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.12.a <- ggplot(data = subset(df, tau == "tau = 0"),
                     aes(x = t, y = N_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(italic(N)))) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = N_25th,
                                  ymax = N_75th), width = .25)
ggsave("sfig12a.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

s.fig.12.b <- ggplot(data = subset(df, tau == "tau = delta"),
                     aes(x = t, y = N_50th)) +
                geom_point() + facet_grid(method~sigmas) +
                xlab(expression(italic(t))) +
                ylab(expression(hat(italic(N)))) +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                geom_errorbar(aes(ymin = N_25th,
                                  ymax = N_75th), width = .25)
ggsave("sfig12b.pdf", device = "pdf", width = 7, height = 4.5, units = "in")

################################################################################
# Table 5                                                                      #
################################################################################

# TDS1 part of Table 5
table.5.TDS1           <- cbind(table.1.TDS1[, 1:2], matrix(0, nrow = 9,
                                                            ncol = 4))
colnames(table.5.TDS1) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "ETI: tau_star = 0", "ETI: tau_star = delta",
                            "EP: tau_star = 0", "EP: tau_star = delta")
table.5.TDS1[, 3]      <- table.1.TDS1[, 3]
table.5.TDS1[, 5]      <- table.1.TDS1[, 6]
table.5.TDS1[, 4]      <- subset(TDS1[157:174, ], tau == 0)$Mean.Reject
table.5.TDS1[, 6]      <- subset(TDS1[157:174, ], tau == 0.2)$Mean.Reject

# TDS2 part of Table 5
table.5.TDS2           <- cbind(table.1.TDS2[, 1:2], matrix(0, nrow = 9,
                                                            ncol = 4))
colnames(table.5.TDS2) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "ETI: tau_star = 0", "ETI: tau_star = delta",
                            "EP: tau_star = 0", "EP: tau_star = delta")
table.5.TDS2[, 3]      <- table.1.TDS2[, 3]
table.5.TDS2[, 5]      <- table.1.TDS2[, 6]
table.5.TDS2[, 4]      <- subset(TDS2[157:174, ], tau == 0)$Mean.Reject
table.5.TDS2[, 6]      <- subset(TDS2[157:174, ], tau == 0.267)$Mean.Reject

# Final Table 5
table.5                <- rbind(table.5.TDS1, table.5.TDS2)

################################################################################
# S. Figures 13-16                                                             #
################################################################################

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = factor(rep(c("tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta"),
                                          each = 9),
                                      levels = c("tau = 0", "tau = delta")),
                         tau_star = factor(rep(c("tau_star = 0",
                                                 "tau_star = delta"),
                                             each = 18),
                                         levels = c("tau_star = 0",
                                                    "tau_star = delta")),
                         sigma_c2_25th = c(TDS1[19:36, ]$X25th.est.sigma.c^2,
                                           TDS1[157:174, ]$X25th.est.sigma.c^2),
                         sigma_c2_50th = c(TDS1[19:36, ]$X50th.est.sigma.c^2,
                                           TDS1[157:174, ]$X50th.est.sigma.c^2),
                         sigma_c2_75th = c(TDS1[19:36, ]$X75th.est.sigma.c^2,
                                           TDS1[157:174, ]$X75th.est.sigma.c^2),
                         sigma_e2_25th = c(TDS1[19:36, ]$X25th.est.sigma.e^2,
                                           TDS1[157:174, ]$X25th.est.sigma.e^2),
                         sigma_e2_50th = c(TDS1[19:36, ]$X50th.est.sigma.e^2,
                                           TDS1[157:174, ]$X50th.est.sigma.e^2),
                         sigma_e2_75th = c(TDS1[19:36, ]$X75th.est.sigma.e^2,
                                           TDS1[157:174, ]$X75th.est.sigma.e^2),
                         N_25th = c(TDS1[19:36, ]$X25th.N,
                                    TDS1[157:174, ]$X25th.N),
                         N_50th = c(TDS1[19:36, ]$X50th.N,
                                    TDS1[157:174, ]$X50th.N),
                         N_75th = c(TDS1[19:36, ]$X75th.N,
                                    TDS1[157:174, ]$X75th.N))

s.fig.13 <- ggplot(data = df,
                    aes(x = sigmas, y = sigma_c2_50th)) +
  geom_point() + facet_grid(tau~tau_star) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(sigma)[c]^2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_errorbar(aes(ymin = sigma_c2_25th,
                    ymax = sigma_c2_75th), width = .25) +
  geom_hline(yintercept = 0.02, linetype = 2)
ggsave("sfig13.pdf", device = "pdf", width = 7, height = 7, units = "in")

s.fig.14 <- ggplot(data = df, aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(tau~tau_star) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig14.pdf", device = "pdf", width = 7, height = 7, units = "in")

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = factor(rep(c("tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta"),
                                          each = 9),
                                      levels = c("tau = 0", "tau = delta")),
                         tau_star = factor(rep(c("tau_star = 0",
                                                 "tau_star = delta"),
                                               each = 18),
                                           levels = c("tau_star = 0",
                                                      "tau_star = delta")),
                         sigma_c2_25th = c(TDS2[19:36, ]$X25th.est.sigma.c^2,
                                           TDS2[157:174, ]$X25th.est.sigma.c^2),
                         sigma_c2_50th = c(TDS2[19:36, ]$X50th.est.sigma.c^2,
                                           TDS2[157:174, ]$X50th.est.sigma.c^2),
                         sigma_c2_75th = c(TDS2[19:36, ]$X75th.est.sigma.c^2,
                                           TDS2[157:174, ]$X75th.est.sigma.c^2),
                         sigma_e2_25th = c(TDS2[19:36, ]$X25th.est.sigma.e^2,
                                           TDS2[157:174, ]$X25th.est.sigma.e^2),
                         sigma_e2_50th = c(TDS2[19:36, ]$X50th.est.sigma.e^2,
                                           TDS2[157:174, ]$X50th.est.sigma.e^2),
                         sigma_e2_75th = c(TDS2[19:36, ]$X75th.est.sigma.e^2,
                                           TDS2[157:174, ]$X75th.est.sigma.e^2),
                         N_25th = c(TDS2[19:36, ]$X25th.N,
                                    TDS2[157:174, ]$X25th.N),
                         N_50th = c(TDS2[19:36, ]$X50th.N,
                                    TDS2[157:174, ]$X50th.N),
                         N_75th = c(TDS2[19:36, ]$X75th.N,
                                    TDS2[157:174, ]$X75th.N))

s.fig.15 <- ggplot(data = df,
                   aes(x = sigmas, y = sigma_c2_50th)) +
  geom_point() + facet_grid(tau~tau_star) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(sigma)[c]^2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_errorbar(aes(ymin = sigma_c2_25th,
                    ymax = sigma_c2_75th), width = .25) +
  geom_hline(yintercept = 1/9, linetype = 2)
ggsave("sfig15.pdf", device = "pdf", width = 7, height = 7, units = "in")

s.fig.16 <- ggplot(data = df, aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(tau~tau_star) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig16.pdf", device = "pdf", width = 7, height = 7, units = "in")

################################################################################
# Table 6                                                                      #
################################################################################

# TDS1 part of Table 6
sigma.c2.tilde         <- c("0.5*sigma.c2", "sigma.c2", "1.5*sigma.c2")
sigma.e2.tilde         <- c("0.5*sigma.e2", "sigma.e2", "1.5*sigma.e2")
table.6.TDS1           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 8))
table.6.TDS1[, 1:2]    <- table.6.TDS1[, c(2, 1)]
colnames(table.6.TDS1) <- c("sigma.c2.tilde", "sigma.e2.tilde",
                            "(n_init,1000): ETI Blinded",
                            "(n_init,1000): ETI Unblinded",
                            "(n_init,1000): EP Blinded",
                            "(n_init,1000): EP Unblinded",
                            "(1,n_init): ETI Blinded",
                            "(1,n_init): ETI Unblinded",
                            "(1,n_init): EP Blinded",
                            "(1,n_init): EP Unblinded")
table.6.TDS1[, 3]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 1
                             & n.min != 1)$Mean.Reject
table.6.TDS1[, 4]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 0
                             & n.min != 1)$Mean.Reject
table.6.TDS1[, 5]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 1
                             & n.min != 1)$Mean.Reject
table.6.TDS1[, 6]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 0
                             & n.min != 1)$Mean.Reject
table.6.TDS1[, 7]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 1
                             & n.min == 1)$Mean.Reject
table.6.TDS1[, 8]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 0
                             & n.min == 1)$Mean.Reject
table.6.TDS1[, 9]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 1
                             & n.min == 1)$Mean.Reject
table.6.TDS1[, 10] <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 0
                             & n.min == 1)$Mean.Reject

# TDS2 part of Table 6
sigma.c2.tilde         <- c("0.5*sigma.c2", "sigma.c2", "1.5*sigma.c2")
sigma.e2.tilde         <- c("0.5*sigma.e2", "sigma.e2", "1.5*sigma.e2")
table.6.TDS2           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 8))
table.6.TDS2[, 1:2]    <- table.6.TDS2[, c(2, 1)]
colnames(table.5.TDS2) <- c("sigma.c2.tilde", "sigma.e2.tilde",
                            "(n_init,1000): ETI Blinded",
                            "(n_init,1000): ETI Unblinded",
                            "(n_init,1000): EP Blinded",
                            "(n_init,1000): EP Unblinded",
                            "(1,n_init): ETI Blinded",
                            "(1,n_init): ETI Unblinded",
                            "(1,n_init): EP Blinded",
                            "(1,n_init): EP Unblinded")
table.6.TDS2[, 3]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 1 &
                               n.min != 1)$Mean.Reject
table.6.TDS2[, 4]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 0 &
                               n.min != 1)$Mean.Reject
table.6.TDS2[, 5]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 1
                             & n.min != 1)$Mean.Reject
table.6.TDS2[, 6]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 0
                             & n.min != 1)$Mean.Reject
table.6.TDS2[, 7]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 1
                             & n.min == 1)$Mean.Reject
table.6.TDS2[, 8]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 0
                             & n.min == 1)$Mean.Reject
table.6.TDS2[, 9]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 1
                             & n.min == 1)$Mean.Reject
table.6.TDS2[, 10] <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 0
                             & n.min == 1)$Mean.Reject

# Final Table 6
table.6            <- rbind(table.6.TDS1, table.6.TDS2)

################################################################################
# Table 7                                                                      #
################################################################################

# TDS1 part of Table 7
sigma.c2.tilde         <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde         <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
table.7.TDS1           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 8))
table.7.TDS1[, 1:2]    <- table.7.TDS1[, c(2, 1)]
colnames(table.7.TDS1) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "(n_init,1000): Med(N|0) Blinded",
                            "(n_init,1000): Med(N|0) Unblinded",
                            "(n_init,1000): Med(N|delta) Blinded",
                            "(n_init,1000): Med(N|delta) Unblinded",
                            "(1,n_init): Med(N|0) Blinded",
                            "(1,n_init): Med(N|0) Unblinded",
                            "(1,n_init): Med(N|delta) Blinded",
                            "(1,n_init): Med(N|delta) Unblinded")
table.7.TDS1[, 3]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 1 &
                               n.min != 1)$X50th.N
table.7.TDS1[, 4]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 0 &
                               n.min != 1)$X50th.N
table.7.TDS1[, 5]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 1 &
                               n.min != 1)$X50th.N
table.7.TDS1[, 6]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 0 &
                               n.min != 1)$X50th.N
table.7.TDS1[, 7]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 1 &
                               n.min == 1)$X50th.N
table.7.TDS1[, 8]  <- subset(TDS1[37:108, ], tau == 0 & blinded == 0 &
                               n.min == 1)$X50th.N
table.7.TDS1[, 9]  <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 1 &
                               n.min == 1)$X50th.N
table.7.TDS1[, 10] <- subset(TDS1[37:108, ], tau == 0.2 & blinded == 0 &
                               n.min == 1)$X50th.N

# TDS1 part of Table 7
sigma.c2.tilde         <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde         <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
table.7.TDS2           <- cbind(expand.grid(sigma.e2.tilde, sigma.c2.tilde),
                                matrix(0, nrow = 9, ncol = 8))
table.7.TDS2[, 1:2]    <- table.7.TDS2[, c(2, 1)]
colnames(table.6.TDS2) <- c("tilde(sigma)_c^2", "tilde(sigma)_e^2",
                            "(n_init,1000): Med(N|0) Blinded",
                            "(n_init,1000): Med(N|0) Unblinded",
                            "(n_init,1000): Med(N|delta) Blinded",
                            "(n_init,1000): Med(N|delta) Unblinded",
                            "(1,n_init): Med(N|0) Blinded",
                            "(1,n_init): Med(N|0) Unblinded",
                            "(1,n_init): Med(N|delta) Blinded",
                            "(1,n_init): Med(N|delta) Unblinded")
table.7.TDS2[, 3]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 1 &
                               n.min != 1)$X50th.N
table.7.TDS2[, 4]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 0 &
                               n.min != 1)$X50th.N
table.7.TDS2[, 5]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 1 &
                               n.min != 1)$X50th.N
table.7.TDS2[, 6]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 0 &
                               n.min != 1)$X50th.N
table.7.TDS2[, 7]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 1 &
                               n.min == 1)$X50th.N
table.7.TDS2[, 8]  <- subset(TDS2[37:108, ], tau == 0 & blinded == 0 &
                               n.min == 1)$X50th.N
table.7.TDS2[, 9]  <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 1 &
                               n.min == 1)$X50th.N
table.7.TDS2[, 10] <- subset(TDS2[37:108, ], tau == 0.267 & blinded == 0 &
                               n.min == 1)$X50th.N
# Final Table 7
table.6            <- rbind(table.6.TDS1, table.6.TDS2)

################################################################################
# S. Figures 17-20                                                             #
################################################################################

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- rbind(sigmas, sigmas, sigmas)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = factor(rep(c("tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta"),
                                          each = 9),
                                      levels = c("tau = 0", "tau = delta")),
                         method = factor(rep(c("Unblinded", "Blinded",
                                               "Unblinded", "Blinded",
                                               "Unblinded", "Blinded"),
                                               each = 18),
                                           levels = c("Unblinded", "Blinded")),
                         n.minn.max = factor(rep(c("(1,1000)",
                                                   "(n_init,1000)",
                                                   "(1,n_init)"), each = 36),
                                             levels = c("(1,1000)",
                                                        "(n_init,1000)",
                                                        "(1,n_init)")),
                         sigma_c2_25th = TDS1[1:108, ]$X25th.est.sigma.c^2,
                         sigma_c2_50th = TDS1[1:108, ]$X50th.est.sigma.c^2,
                         sigma_c2_75th = TDS1[1:108, ]$X75th.est.sigma.e^2,
                         sigma_e2_25th = TDS1[1:108, ]$X25th.est.sigma.e^2,
                         sigma_e2_50th = TDS1[1:108, ]$X50th.est.sigma.e^2,
                         sigma_e2_75th = TDS1[1:108, ]$X75th.est.sigma.e.1^2,
                         N_25th = TDS1[1:108, ]$X25th.N,
                         N_50th = TDS1[1:108, ]$X50th.N,
                         N_75th = TDS1[1:108, ]$X75th.N)

s.fig.17 <- ggplot(data = subset(df, tau == "tau = 0"),
                   aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(method~n.minn.max) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig17.pdf", device = "pdf", width = 7, height = 7, units = "in")

s.fig.18 <- ggplot(data = subset(df, tau == "tau = delta"),
                   aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(method~n.minn.max) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig18.pdf", device = "pdf", width = 7, height = 7, units = "in")

sigma.c2.tilde <- c("0.5*sigma_c^2", "sigma_c^2", "1.5*sigma_c^2")
sigma.e2.tilde <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
block          <- expand.grid(sigma.e2.tilde, sigma.c2.tilde)[, 2:1]
sigmas         <- rbind(block, block, block, block)
sigmas         <- rbind(sigmas, sigmas, sigmas)
sigmas         <- paste("(", sigmas[, 1], ", ", sigmas[, 2], ")", sep = "")
df             <- tibble(sigmas = factor(sigmas, levels = unique(sigmas)),
                         tau = factor(rep(c("tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta",
                                            "tau = 0", "tau = delta"),
                                          each = 9),
                                      levels = c("tau = 0", "tau = delta")),
                         method = factor(rep(c("Unblinded", "Blinded",
                                               "Unblinded", "Blinded",
                                               "Unblinded", "Blinded"),
                                             each = 18),
                                         levels = c("Unblinded", "Blinded")),
                         n.minn.max = factor(rep(c("(1,1000)",
                                                   "(n_init,1000)",
                                                   "(1,n_init)"), each = 36),
                                             levels = c("(1,1000)",
                                                        "(n_init,1000)",
                                                        "(1,n_init)")),
                         sigma_c2_25th = TDS2[1:108, ]$X25th.est.sigma.c^2,
                         sigma_c2_50th = TDS2[1:108, ]$X50th.est.sigma.c^2,
                         sigma_c2_75th = TDS2[1:108, ]$X75th.est.sigma.e^2,
                         sigma_e2_25th = TDS2[1:108, ]$X25th.est.sigma.e^2,
                         sigma_e2_50th = TDS2[1:108, ]$X50th.est.sigma.e^2,
                         sigma_e2_75th = TDS2[1:108, ]$X75th.est.sigma.e.1^2,
                         N_25th = TDS2[1:108, ]$X25th.N,
                         N_50th = TDS2[1:108, ]$X50th.N,
                         N_75th = TDS2[1:108, ]$X75th.N)

s.fig.19 <- ggplot(data = subset(df, tau == "tau = 0"),
                   aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(method~n.minn.max) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig19.pdf", device = "pdf", width = 7, height = 7, units = "in")

s.fig.20 <- ggplot(data = subset(df, tau == "tau = delta"),
                   aes(x = sigmas, y = N_50th)) +
  geom_point() + facet_grid(method~n.minn.max) +
  xlab(expression(paste("(", tilde(sigma)[c]^2, ",",
                        tilde(sigma)[e]^2, ")", sep = ""))) +
  ylab(expression(hat(italic(N)))) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  geom_errorbar(aes(ymin = N_25th, ymax = N_75th), width = 0.25)
ggsave("sfig20.pdf", device = "pdf", width = 7, height = 7, units = "in")

################################################################################
# S. Table 1                                                                   #
################################################################################

# TDS1 part of S. Table 1
sigma.e2.tilde           <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
s.table.1.TDS1           <- cbind(sigma.e2.tilde,
                                  matrix(0, nrow = 3, ncol = 2))
colnames(s.table.1.TDS1) <- c("sigma.e2.tilde", "ETI", "EP")
s.table.1.TDS1[, 2]      <- subset(TDS1[151:156, ], tau == 0)$Mean.Reject
s.table.1.TDS1[, 3]      <- subset(TDS1[151:156, ], tau == 0.2)$Mean.Reject

# TDS2 part of S. Table 1
sigma.e2.tilde           <- c("0.5*sigma_e^2", "sigma_e^2", "1.5*sigma_e^2")
s.table.1.TDS2           <- cbind(sigma.e2.tilde,
                                  matrix(0, nrow = 3, ncol = 2))
colnames(s.table.1.TDS2) <- c("sigma.e2.tilde", "ETI", "EP")
s.table.1.TDS2[, 2]      <- subset(TDS2[151:156, ], tau == 0)$Mean.Reject
s.table.1.TDS2[, 3]      <- subset(TDS2[151:156, ], tau == 0.267)$Mean.Reject

# Final S. Table 1
s.table.1                <- rbind(s.table.1.TDS1, s.table.1.TDS2)

################################################################################
# S. Table 2                                                                   #
################################################################################

# TDS1 part of S. Table 2
s.table.2.TDS1           <- cbind(scenarios.6[1:18, 1],
                                  scenarios.6[1:18, 1]*scenarios.6[1:18, 2],
                                  matrix(0, nrow = 18, ncol = 4))
s.table.2.TDS1[, 3]      <- subset(TDS1[175:246, ], tau == 0 &
                                   blinded == 1)$Mean.Reject
s.table.2.TDS1[, 4]      <- subset(TDS1[175:246, ], tau == 0 &
                                   blinded == 0)$Mean.Reject
s.table.2.TDS1[, 5]      <- subset(TDS1[175:246, ], tau != 0 &
                                   blinded == 1)$Mean.Reject
s.table.2.TDS1[, 6]      <- subset(TDS1[175:246, ], tau != 0 &
                                   blinded == 0)$Mean.Reject
s.table.2.TDS1           <- as_tibble(s.table.2.TDS1)
colnames(s.table.2.TDS1) <- c("sigma^2", "tilde(sigma)^2", "ETI: Blinded",
                              "ETI: Unblinded", "EP: Blinded", "EP: Unblinded")
s.table.2.TDS1           <- arrange(s.table.2.TDS1, `sigma^2`, `tilde(sigma)^2`)

# TDS2 part of S. Table 2
s.table.2.TDS2           <- cbind(scenarios.6[1:18, 1],
                                  scenarios.6[1:18, 1]*scenarios.6[1:18, 2],
                                  matrix(0, nrow = 18, ncol = 4))
s.table.2.TDS2[, 3]      <- subset(TDS2[175:246, ], tau == 0 &
                                     blinded == 1)$Mean.Reject
s.table.2.TDS2[, 4]      <- subset(TDS2[175:246, ], tau == 0 &
                                     blinded == 0)$Mean.Reject
s.table.2.TDS2[, 5]      <- subset(TDS2[175:246, ], tau != 0 &
                                     blinded == 1)$Mean.Reject
s.table.2.TDS2[, 6]      <- subset(TDS2[175:246, ], tau != 0 &
                                     blinded == 0)$Mean.Reject
s.table.2.TDS2           <- as_tibble(s.table.2.TDS2)
colnames(s.table.2.TDS2) <- c("sigma^2", "tilde(sigma)^2", "ETI: Blinded",
                              "ETI: Unblinded", "EP: Blinded", "EP: Unblinded")
s.table.2.TDS2           <- arrange(s.table.2.TDS2, `sigma^2`, `tilde(sigma)^2`)

# Final S. Table 2
s.table.2                <- rbind(s.table.2.TDS1, s.table.2.TDS2)