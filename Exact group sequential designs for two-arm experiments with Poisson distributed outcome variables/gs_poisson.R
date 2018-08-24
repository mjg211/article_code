##### Load required packages ###################################################

library(skellam)
library(mvtnorm)

##### Utility functions for normal approximation approach ######################

# Function to return information vector of normal approximation design
information_normal <- function(lambda1, lambda2, n, K) {
  return(n*(1:K)/(lambda1 + lambda2))
}

# Function to return covariance matrix of normal approximation design
covariance_normal <- function(I, K) {
  Sigma             <- diag(1, K, K)
  for (k1 in 2:K) {
    for (k2 in 1:(k1 - 1)) {
      Sigma[k1, k2] <- sqrt(I[k2]/I[k1])
      Sigma[k2, k1] <- Sigma[k1, k2]
    }
  }
  return(Sigma)
}

# Function to use to find r[k] of normal approximation design
find_rk_normal <- function(rk, r, a, Sigma, piRk) {
  Rk <- mvtnorm::pmvnorm(lower = c(a, rk), upper = c(r, Inf), sigma = Sigma)[1]
  return(Rk - piRk)
}

# Function to use to find a[k] of normal approximation design
find_ak_normal <- function(ak, r, a, I, delta, Sigma, piAk) {
  Ak <- mvtnorm::pmvnorm(lower = c(a, -Inf), upper = c(r, ak),
                         mean = delta*sqrt(I), sigma = Sigma)[1]
  return(Ak - piAk)
}

# Function to return power of normal approximation design
power_normal <- function(lambda1, lambda2, K, a, r, n) {
  I        <- information_normal(lambda1, lambda2, n, K)
  if (K > 1) {
    Sigma  <- covariance_normal(I, K)
  }
  R        <- A <- numeric(K)
  R[1]     <- stats::pnorm(r[1], mean = (lambda1 - lambda2)*sqrt(I[1]),
                           lower.tail = F)
  if (K > 1) {
    for (k in 2:K) {
      R[k] <- mvtnorm::pmvnorm(lower = c(a[1:(k - 1)], r[k]),
                               upper = c(r[1:(k - 1)], Inf),
                               mean = (lambda1 - lambda2)*sqrt(I[1:k]),
                               sigma = Sigma[1:k, 1:k])[1]
    }
  }
  return(sum(R))
}

# Function to evaluate power, a, and r for given n, in normal approximation
# design
try_n_normal <- function(n, K, piA, piR, beta, lambda0, delta, type) {
  r          <- a <- numeric(K)
  I_H0       <- information_normal(lambda0, lambda0, n, K)
  I_H1       <- information_normal(lambda0, lambda0 - delta, n, K)
  if (K > 1) {
    Sigma_H0 <- covariance_normal(I_H0, K)
    Sigma_H1 <- covariance_normal(I_H1, K)
  }
  r[1]       <- stats::qnorm(1 - piR[1])
  if (K > 1) {
    a[1]     <- stats::qnorm(piA[1], mean = delta*sqrt(I_H1[1]))
    if (a[1] > r[1]) {
      return(1)
    }
    if (K > 2) {
      for (k in 2:(K - 1)) {
        r[k] <- stats::uniroot(f = find_rk_normal, interval = c(-20, 20),
                               a = a[1:(k - 1)], r = r[1:(k - 1)],
                               Sigma = Sigma_H0[1:k, 1:k], piRk = piR[k])$root
        a[k] <- stats::uniroot(f = find_ak_normal, interval = c(-20, 20),
                               a = a[1:(k - 1)], r = r[1:(k - 1)],
                               I = I_H1[1:k], delta = delta,
                               Sigma = Sigma_H1[1:k, 1:k], piAk = piA[k])$root
        if (a[k] > r[k]) {
          return(1)
        }
      }
    }
    if (find_rk_normal(-20, r[1:(K - 1)], a[1:(K - 1)], Sigma_H0,
                       piR[K]) <= 0) {
      r[K]   <- a[K] <- -20
    } else {
      r[K]   <- a[K] <- stats::uniroot(f = find_rk_normal,
                                       interval = c(-20, 20), a = a[1:(K - 1)],
                                       r = r[1:(K - 1)], Sigma = Sigma_H0,
                                       piRk = piR[K])$root
    }
  } else {
    a        <- r
  }
  P_H1       <- power_normal(lambda0, lambda0 - delta, K, a, r, n)
  if (type == 1) {
    return(P_H1 - (1 - beta))
  } else {
    return(list(a = a, r = r))
  }
}

# Function to find required n for given spend and return the operating
# characteristics of the resulting design, in normal approximation approach
obj_fn_normal <- function(piA, piR, K, beta, lambda0, delta, n_fixed) {
  n         <- ceiling(stats::uniroot(f = try_n_normal,
                                      interval = c(10^-6, n_fixed), K = K,
                                      piA = piA, piR = piR, beta = beta,
                                      lambda0 = lambda0, delta = delta,
                                      type = 1)$root)
  des       <- try_n_normal(n, K, piA, piR, beta, lambda0, delta, 0)
  r         <- des$r
  a         <- des$a
  opchar_H0 <- int_opchar_normal(lambda0, lambda0, K, n, a, r)
  opchar_H1 <- int_opchar_normal(lambda0, lambda0 - delta, K, n, a, r)
  return(list(P_H0 = opchar_H0$P, ESS_H0 = opchar_H0$ESS, P_H1 = opchar_H1$P,
              ESS_H1 = opchar_H1$ESS, n = n, a = a, r = r))
}

# Function to return operating characteristics of normal approximation design
int_opchar_normal <- function(lambda1, lambda2, K, n, a, r) {
  n_vec    <- rep(n, K)
  I        <- information_normal(lambda1, lambda2, n, K)
  if (K > 1) {
    Sigma  <- covariance_normal(I, K)
  }
  R        <- A <- numeric(K)
  R[1]     <- stats::pnorm(r[1], mean = (lambda1 - lambda2)*sqrt(I[1]),
                           lower.tail = F)
  A[1]     <- stats::pnorm(a[1], mean = (lambda1 - lambda2)*sqrt(I[1]))
  if (K > 1) {
    for (k in 2:K) {
      R[k] <- mvtnorm::pmvnorm(lower = c(a[1:(k - 1)], r[k]),
                               upper = c(r[1:(k - 1)], Inf),
                               mean = (lambda1 - lambda2)*sqrt(I[1:k]),
                               sigma = Sigma[1:k, 1:k])[1]
      A[k] <- mvtnorm::pmvnorm(lower = c(a[1:(k - 1)], -Inf),
                               upper = c(r[1:(k - 1)], a[k]),
                               mean = (lambda1 - lambda2)*sqrt(I[1:k]),
                               sigma = Sigma[1:k, 1:k])[1]
    }
  }
  P        <- sum(R)
  cum_S    <- cumsum(S <- A + R)
  N        <- cumsum(2*n_vec)
  Med      <- ifelse(any(cum_S == 0.5), 0.5*(N[which(cum_S == 0.5)] +
                                               N[which(cum_S == 0.5) + 1]),
                     N[which(cum_S > 0.5)[1]])
  ESS      <- sum(N*S)
  VSS      <- sum(N^2*S) - ESS^2
  return(list(P = P, ESS = ESS, VSS = VSS, Med = Med, A = A, R = R, S = S,
              cum_S = cum_S))
}

##### Utility functions for exact approach #####################################

# Function to return operating characteristics of exact designs
int_opchar_exact <- function(lambda1, lambda2, K, n_vec, a, r) {
  A    <- R <- numeric(K)
  A[1] <- skellam::pskellam(a[1], n_vec[1]*lambda1, n_vec[1]*lambda2)
  R[1] <- skellam::pskellam(r[1] - 1, n_vec[1]*lambda1, n_vec[1]*lambda2,
                            lower.tail = F)
  if (K > 1) {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1),
                                          n_vec[1]*lambda1, n_vec[1]*lambda2)
    if (K > 2) {
      for (k in 2:(K - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) -
                                l, n_vec[k]*lambda1,
                              n_vec[k]*lambda2)
          A[k]       <- A[k] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(a[k] - l, n_vec[k]*lambda1,
                              n_vec[k]*lambda2)
          R[k]       <- R[k] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(r[k] - l - 1,
                              n_vec[k]*lambda1,
                              n_vec[k]*lambda2,
                              lower.tail = F)
        }
      }
    }
    for (l in (a[K - 1] + 1):(r[K - 1] - 1)) {
      A[K] <- A[K] + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(a[K] - l, n_vec[K]*lambda1,
                          n_vec[K]*lambda2)
      R[K] <- R[K] + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(r[K] - l - 1, n_vec[K]*lambda1,
                          n_vec[K]*lambda2, lower.tail = F)
    }
  }
  P     <- sum(R)
  cum_S <- cumsum(S <- A + R)
  N     <- cumsum(2*n_vec)
  Med   <- ifelse(any(cum_S == 0.5), 0.5*(N[which(cum_S == 0.5)] +
                                            N[which(cum_S == 0.5) + 1]),
                  N[which(cum_S > 0.5)[1]])
  ESS   <- sum(N*S)
  VSS   <- sum(N^2*S) - ESS^2
  return(list(P = P, ESS = ESS, VSS = VSS, Med = Med, A = A, R = R, S = S,
              cum_S = cum_S))
}

# Function to return power of exact designs
power_exact <- function(lambda1, lambda2, K, n_vec, a, r) {
  R       <- numeric(K)
  R[1]    <- skellam::pskellam(r[1] - 1, n_vec[1]*lambda1, n_vec[1]*lambda2,
                               lower.tail = F)
  if (K > 1) {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1),
                                          n_vec[1]*lambda1, n_vec[1]*lambda2)
    if (K > 2) {
      for (k in 2:(K - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) -
                                l, n_vec[k]*lambda1,
                              n_vec[k]*lambda2)
          R[k]       <- R[k] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(r[k] - l - 1,
                              n_vec[k]*lambda1,
                              n_vec[k]*lambda2,
                              lower.tail = F)
        }
      }
    }
    for (l in (a[K - 1] + 1):(r[K - 1] - 1)) {
      R[K] <- R[K] + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(r[K] - l - 1, n_vec[K]*lambda1,
                          n_vec[K]*lambda2, lower.tail = F)
    }
  }
  return(sum(R))
}

# Function to find required n for given spend and return the operating
# characteristics of the resulting design, in the exact approach
obj_fn_exact <- function(piA, piR, K, beta, lambda0, delta, n_start) {
  n           <- n_start
  power_n     <- try_n_exact(n, piA, piR, K, lambda0, delta)
  if (power_n$power < 1 - beta) {
    while (power_n$power < 1 - beta) {
      n       <- n + 1
      power_n <- try_n_exact(n, piA, piR, K, lambda0, delta)
    }
    a         <- power_n$a
    r         <- power_n$r
  } else if (power_n$power == 1 - beta) {
    a         <- power_n$a
    r         <- power_n$r
  } else {
    while (power_n$power > 1 - beta) {
      a       <- power_n$a
      r       <- power_n$r
      n       <- n - 1
      power_n <- try_n_exact(n, piA, piR, K, lambda0, delta)
    }
    n         <- n + 1
  }
  opchar_H0   <- int_opchar_exact(lambda0, lambda0, K, rep(n, K), a, r)
  opchar_H1   <- int_opchar_exact(lambda0, lambda0 - delta, K, rep(n, K), a, r)
  return(list(P_H0 = opchar_H0$P, ESS_H0 = opchar_H0$ESS, P_H1 = opchar_H1$P,
              ESS_H1 = opchar_H1$ESS, n = n, a = a, r = r))
}

# Function to evaluate power, e, and f for given n, in exact approach
try_n_exact <- function(n, piA, piR, K, lambda0, delta) {
  r            <- a <- numeric(K)
  N            <- 2*n*(1:K)
  for (k in 1:K) {
    r[k]       <- skellam::qskellam(1 - piR[k], 0.5*N[k]*lambda0,
                                    0.5*N[k]*lambda0)
    Rk         <- tI_k(lambda0, a[1:k], r[1:k], k, rep(n, k))
    if (Rk > piR[k]) {
      while (Rk > piR[k]) {
        r[k]   <- r[k] + 1
        Rk     <- tI_k(lambda0, a[1:k], r[1:k], k, rep(n, k))
      }
    } else if (Rk < piR[k]) {
      while (Rk < piR[k]) {
        r[k]   <- r[k] - 1
        Rk     <- tI_k(lambda0, a[1:k], r[1:k], k, rep(n, k))
      }
      r[k]     <- r[k] + 1
    }
    if (k < K) {
      a[k]     <- skellam::qskellam(piA[k], 0.5*N[k]*lambda0,
                                    0.5*N[k]*(lambda0 - delta))
      Ak       <- tII_k(lambda0, a[1:k], r[1:k], k, rep(n, k), delta)
      if (Ak > piA[k]) {
        while (Ak > piA[k]) {
          a[k] <- a[k] - 1
          Ak   <- tII_k(lambda0, a[1:k], r[1:k], k, rep(n, k), delta)
        }
      } else if (Ak < piA[k]) {
        while (Ak < piA[k]) {
          a[k] <- a[k] + 1
          Ak   <- tII_k(lambda0, a[1:k], r[1:k], k, rep(n, k), delta)
        }
        a[k]   <- a[k] - 1
      }
    } else {
      a[k]     <- r[k] - 1
    }
  }
  power        <- power_exact(lambda0, lambda0 - delta, K, rep(n, K), a, r)
  return(list(power = power, a = a, r = r))
}

# Function to return type-I error-rate at stage k of exact design
tI_k <- function(lambda0, a, r, stage_k, n) {
  if (stage_k == 1) {
    R_stage_k        <- skellam::pskellam(r[1] - 1, n*lambda0, n*lambda0,
                                          lower.tail = F)
  } else {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1), n[1]*lambda0,
                                          n[1]*lambda0)
    if (stage_k > 2) {
      for (k in 2:(stage_k - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) -
                                l, n[k]*lambda0,
                              n[k]*lambda0)
        }
      }
    }
    R_stage_k        <- 0
    for (l in (a[stage_k - 1] + 1):(r[stage_k - 1] - 1)) {
      R_stage_k      <- R_stage_k + pmf_C[[stage_k - 1]][l - a[stage_k - 1]]*
        skellam::pskellam(r[stage_k] - l - 1,
                          n[stage_k]*lambda0,
                          n[stage_k]*lambda0,
                          lower.tail = F)
    }
  }
  return(R_stage_k)
}

# Function to return type-II error-rate at stage k of exact design
tII_k <- function(lambda0, a, r, stage_k, n, delta) {
  lambda2            <- lambda0 - delta
  if (stage_k == 1) {
    A_stage_k        <- skellam::pskellam(a[1], n*lambda0, n*lambda2)
  } else {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1), n[1]*lambda0,
                                          n[1]*lambda2)
    if (stage_k > 2) {
      for (k in 2:(stage_k - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) -
                                l, n[k]*lambda0,
                              n[k]*lambda2)
        }
      }
    }
    A_stage_k        <- 0
    for (l in (a[stage_k - 1] + 1):(r[stage_k - 1] - 1)) {
      A_stage_k      <- A_stage_k + pmf_C[[stage_k - 1]][l - a[stage_k - 1]]*
        skellam::pskellam(a[stage_k] - l,
                          n[stage_k]*lambda0,
                          n[stage_k]*lambda2)
    }
  }
  return(A_stage_k)
}

##### Function to return normal approximation and exact designs ################

des_skellam <- function(K = 3, alpha = 0.05, beta = 0.1, lambda0 = 6.25,
                        delta = 1.1, piAk, piRk, method = "exact") {
  if (method == "exact") {
    single_stage <- obj_fn_exact(beta, alpha, 1, beta, lambda0, delta, 10)
  } else {
    single_stage <- obj_fn_normal(beta, alpha, 1, beta, lambda0, delta, 1000)
  }
  if (K > 1) {
    grid_piA   <- as.matrix(expand.grid(rep(list(piAk), K - 1)))
    grid_piA   <- grid_piA[rowSums(grid_piA) < beta, , drop = F]
    grid_piR   <- as.matrix(expand.grid(rep(list(piRk), K - 1)))
    grid_piR   <- grid_piR[rowSums(grid_piR) < alpha, , drop = F]
    designs    <- matrix(0, nrow = nrow(grid_piA)*nrow(grid_piR),
                         ncol = 5 + 4*K)
    count      <- 1
    for (k1 in 1:nrow(grid_piA)) {
      for (k2 in 1:nrow(grid_piR)) {
        if (method == "exact") {
          design_k1k2    <- obj_fn_exact(c(grid_piA[k1, ],
                                           beta - sum(grid_piA[k1, ])),
                                         c(grid_piR[k2, ],
                                           alpha - sum(grid_piR[k2, ])),
                                         K, beta, lambda0, delta,
                                         ceiling(single_stage$n/K))
        } else {
          design_k1k2    <- obj_fn_normal(c(grid_piA[k1, ],
                                            beta - sum(grid_piA[k1, ])),
                                          c(grid_piR[k2, ],
                                            alpha - sum(grid_piR[k2, ])),
                                          K, beta, lambda0, delta,
                                          single_stage$n)
        }
        designs[count, ] <- c(c(grid_piA[k1, ], beta - sum(grid_piA[k1, ])),
                              c(grid_piR[k2, ], alpha - sum(grid_piR[k2, ])),
                              design_k1k2$P_H0, design_k1k2$P_H1, design_k1k2$n,
                              design_k1k2$ESS_H0, design_k1k2$ESS_H1,
                              design_k1k2$a, design_k1k2$r)
        count            <- count + 1
      }
    }
  } else {
    designs              <- matrix(c(beta, alpha, single_stage$P_H0,
                                     single_stage$P_H1, single_stage$n,
                                     single_stage$ESS_H0, single_stage$ESS_H1,
                                     single_stage$a, single_stage$r), nrow = 1)
  }
  colnames(designs)      <- c(paste("piA", 1:K, sep = ""),
                              paste("piR", 1:K, sep = ""), "P(H0)", "P(H1)",
                              "n", "ESS(H0)", "ESS(H1)",
                              paste("a", 1:K, sep = ""),
                              paste("r", 1:K, sep = ""))
  designs                <- tibble::as_tibble(designs)
  return(designs)
}

##### Function to simulate trials given a design ###############################

sim_skellam <- function(lambda1 = 6.25, lambda2 = lambda1, K = 2, n = 45,
                        a = c(0.631,1.579), r = c(2.576, 1.579),
                        method = "normal", replicates = 100000,
                        seed = Sys.time()) {
  set.seed(seed)
  reject            <- sample_size <- numeric(100000)
  X1                <- X2 <- matrix(0, K, n)
  for (rep in 1:replicates) {
    Y1              <- Y2 <- numeric(K)
    for (k in 1:K) {
      X1[k, ]       <- rpois(n, lambda1)
      X2[k, ]       <- rpois(n, lambda2)
      Y1[k]         <- sum(X1[k, ])
      Y2[k]         <- sum(X2[k, ])
      if (method == "normal") {
        hat_lambda1 <- sum(Y1)/(n*k)
        hat_lambda2 <- sum(Y2)/(n*k)
        Tk          <- (hat_lambda1 - hat_lambda2)*
          sqrt(k*n/(hat_lambda1 + hat_lambda2))
      } else {
        Tk          <- sum(Y1) - sum(Y2)
      }
      if (Tk >= r[k]) {
        reject[rep]      <- 1
        sample_size[rep] <- 2*k*n
        break
      } else if (Tk < a[k]) {
        sample_size[rep] <- 2*k*n
        break
      }
    }
  }
  return(list(power = mean(reject), ess = mean(sample_size)))
}

##### Scenario 1: alpha = 0.05, beta = 0.1, lambda0 = 6.25, delta = 1.1 ########

# Single-stage
sc1_exact_1         <- des_skellam(K = 1)
sc1_normal_1        <- des_skellam(K = 1, method = "normal")
# Two-stage
sc1_exact_2         <- des_skellam(K = 2, piAk = seq(0.01, 0.09, 0.01),
                                   piRk = seq(0.005, 0.045, 0.005))
sc1_normal_2        <- des_skellam(K = 2, piAk = seq(0.01, 0.09, 0.01),
                                   piRk = seq(0.005, 0.045, 0.005),
                                   method = "normal")
# Three-stage
sc1_exact_3         <- des_skellam(K = 3, piAk = c(0.01, 0.03, 0.05, 0.07),
                                   piRk = c(0.01, 0.015, 0.02, 0.025, 0.03,
                                            0.035))
sc1_normal_3        <- des_skellam(K = 3, piAk = c(0.01, 0.03, 0.05, 0.07),
                                   piRk = c(0.01, 0.015, 0.02, 0.025, 0.03,
                                            0.035), method = "normal")
# Find the optimal designs
w                   <- cbind(c(1, 0, 1/2, 1/2, 0, 1/3),
                             c(0, 1, 1/2, 0, 1/2, 1/3),
                             c(0, 0, 0, 1/2, 1/2, 1/3))
sc1_optimal_exact_2 <- sc1_optimal_normal_2 <- matrix(0, nrow = 6,
                                                      ncol = 5 + 4*2)
sc1_optimal_exact_3 <- sc1_optimal_normal_3 <- matrix(0, nrow = 6,
                                                      ncol = 5 + 4*3)
for (i in 1:nrow(w)) {
  # Find the scores
  scores_exact_2    <- w[i, 1]*sc1_exact_2$`ESS(H0)` +
                         w[i, 2]*sc1_exact_2$`ESS(H1)` +
                         w[i, 3]*sc1_exact_2$n*2*2
  scores_exact_3    <- w[i, 1]*sc1_exact_3$`ESS(H0)` +
                         w[i, 2]*sc1_exact_3$`ESS(H1)` +
                         w[i, 3]*sc1_exact_3$n*2*3
  scores_normal_2   <- w[i, 1]*sc1_normal_2$`ESS(H0)` +
                         w[i, 2]*sc1_normal_2$`ESS(H1)` +
                         w[i, 3]*sc1_normal_2$n*2*2
  scores_normal_3   <- w[i, 1]*sc1_normal_3$`ESS(H0)` +
                         w[i, 2]*sc1_normal_3$`ESS(H1)` +
                         w[i, 3]*sc1_normal_3$n*2*3
  # Store which is optimal
  sc1_optimal_exact_2[i, ]     <-
    as.numeric(sc1_exact_2[which.min(scores_exact_2), ])
  sc1_optimal_exact_3[i, ]     <-
    as.numeric(sc1_exact_3[which.min(scores_exact_3), ])
  sc1_optimal_normal_2[i, ]    <-
    as.numeric(sc1_normal_2[which.min(scores_normal_2), ])
  sc1_optimal_normal_3[i, ]    <-
    as.numeric(sc1_normal_3[which.min(scores_normal_3), ])
}
sc1_optimal_exact_2            <- cbind(w, sc1_optimal_exact_2)
sc1_optimal_exact_3            <- cbind(w, sc1_optimal_exact_3)
sc1_optimal_normal_2           <- cbind(w, sc1_optimal_normal_2)
sc1_optimal_normal_3           <- cbind(w, sc1_optimal_normal_3)
colnames(sc1_optimal_exact_2)  <- c("w1", "w2", "w3", colnames(sc1_exact_2))
colnames(sc1_optimal_exact_3)  <- c("w1", "w2", "w3", colnames(sc1_exact_3))
colnames(sc1_optimal_normal_2) <- c("w1", "w2", "w3", colnames(sc1_normal_2))
colnames(sc1_optimal_normal_3) <- c("w1", "w2", "w3", colnames(sc1_normal_3))

##### Scenario 2: alpha = 0.05, beta = 0.1, lambda0 = 6.25, delta = 3.25 #######

# Single-stage
sc2_exact_1         <- des_skellam(K = 1, delta = 3.25)
sc2_normal_1        <- des_skellam(K = 1, delta = 3.25, method = "normal")
# Two-stage
sc2_exact_2         <- des_skellam(K = 2, delta = 3.25,
                                   piAk = seq(0.01, 0.08, 0.01),
                                   piRk = seq(0.005, 0.045, 0.005))
sc2_normal_2        <- des_skellam(K = 2, delta = 3.25,
                                   piAk = seq(0.01, 0.08, 0.01),
                                   piRk = seq(0.005, 0.045, 0.005),
                                   method = "normal")
# Three-stage
sc2_exact_3         <- des_skellam(K = 3, delta = 3.25,
                                   piAk = c(0.01, 0.03, 0.05, 0.07),
                                   piRk = c(0.01, 0.015, 0.02, 0.025, 0.03,
                                            0.035))
sc2_normal_3        <- des_skellam(K = 3, delta = 3.25,
                                   piAk = c(0.01, 0.03, 0.05, 0.07),
                                   piRk = c(0.01, 0.015, 0.02, 0.025, 0.03,
                                            0.035), method = "normal")
# Find the optimal designs
w                   <- cbind(c(1, 0, 1/2, 1/2, 0, 1/3),
                             c(0, 1, 1/2, 0, 1/2, 1/3),
                             c(0, 0, 0, 1/2, 1/2, 1/3))
sc2_optimal_exact_2 <- sc2_optimal_normal_2 <- matrix(0, nrow = 6,
                                                      ncol = 5 + 4*2)
sc2_optimal_exact_3 <- sc2_optimal_normal_3 <- matrix(0, nrow = 6,
                                                      ncol = 5 + 4*3)
for (i in 1:nrow(w)) {
  # Find the scores
  scores_exact_2    <- w[i, 1]*sc2_exact_2$`ESS(H0)` +
                         w[i, 2]*sc2_exact_2$`ESS(H1)` +
                         w[i, 3]*sc2_exact_2$n*2*2
  scores_exact_3    <- w[i, 1]*sc2_exact_3$`ESS(H0)` +
                         w[i, 2]*sc2_exact_3$`ESS(H1)` +
                         w[i, 3]*sc2_exact_3$n*2*3
  scores_normal_2   <- w[i, 1]*sc2_normal_2$`ESS(H0)` +
                         w[i, 2]*sc2_normal_2$`ESS(H1)` +
                         w[i, 3]*sc2_normal_2$n*2*2
  scores_normal_3   <- w[i, 1]*sc2_normal_3$`ESS(H0)` +
                         w[i, 2]*sc2_normal_3$`ESS(H1)` +
                         w[i, 3]*sc2_normal_3$n*2*3
  # Store which is optimal
  sc2_optimal_exact_2[i, ]     <-
    as.numeric(sc2_exact_2[which.min(scores_exact_2), ])
  sc2_optimal_exact_3[i, ]     <-
    as.numeric(sc2_exact_3[which.min(scores_exact_3), ])
  sc2_optimal_normal_2[i, ]    <-
    as.numeric(sc2_normal_2[which.min(scores_normal_2), ])
  sc2_optimal_normal_3[i, ]    <-
    as.numeric(sc2_normal_3[which.min(scores_normal_3), ])
}
sc2_optimal_exact_2            <- cbind(w, sc2_optimal_exact_2)
sc2_optimal_exact_3            <- cbind(w, sc2_optimal_exact_3)
sc2_optimal_normal_2           <- cbind(w, sc2_optimal_normal_2)
sc2_optimal_normal_3           <- cbind(w, sc2_optimal_normal_3)
colnames(sc2_optimal_exact_2)  <- c("w1", "w2", "w3", colnames(sc2_exact_2))
colnames(sc2_optimal_exact_3)  <- c("w1", "w2", "w3", colnames(sc2_exact_3))
colnames(sc2_optimal_normal_2) <- c("w1", "w2", "w3", colnames(sc2_normal_2))
colnames(sc2_optimal_normal_3) <- c("w1", "w2", "w3", colnames(sc2_normal_3))

##### Simulation study #########################################################

# Single-stage
sc1_normal_1_typeI   <- sim_skellam(K = 1, n = sc1_normal_1$n,
                                    a = sc1_normal_1$a1, r = sc1_normal_1$r1,
                                    seed = 1)
sc1_normal_1_power   <- sim_skellam(lambda2 = 5.15, K = 1, n = sc1_normal_1$n,
                                    a = sc1_normal_1$a1, r = sc1_normal_1$r1,
                                    seed = 2)
# Two-stage
sc1_normal_2_typeI   <- sc1_normal_2_power <- list()
sc1_optimal_normal_2 <- tibble::as_tibble(sc1_optimal_normal_2)
for (i in 1:6) {
  sc1_normal_2_typeI[[i]] <- sim_skellam(K = 2, n = sc1_optimal_normal_2$n[i],
                                         a = c(sc1_optimal_normal_2$a1[i],
                                               sc1_optimal_normal_2$a2[i]),
                                         r = c(sc1_optimal_normal_2$r1[i],
                                               sc1_optimal_normal_2$r2[i]),
                                         seed = i + 2)
  sc1_normal_2_power[[i]] <- sim_skellam(K = 2, lambda2 = 5.15,
                                         n = sc1_optimal_normal_2$n[i],
                                         a = c(sc1_optimal_normal_2$a1[i],
                                               sc1_optimal_normal_2$a2[i]),
                                         r = c(sc1_optimal_normal_2$r1[i],
                                               sc1_optimal_normal_2$r2[i]),
                                         seed = i + 8)
}
# Three-stage
sc1_normal_3_typeI        <- sc1_normal_3_power <- list()
sc1_optimal_normal_3      <- tibble::as_tibble(sc1_optimal_normal_3)
for (i in 1:5) {
  sc1_normal_3_typeI[[i]] <- sim_skellam(K = 3, n = sc1_optimal_normal_3$n[i],
                                         a = c(sc1_optimal_normal_3$a1[i],
                                               sc1_optimal_normal_3$a2[i],
                                               sc1_optimal_normal_3$a3[i]),
                                         r = c(sc1_optimal_normal_3$r1[i],
                                               sc1_optimal_normal_3$r2[i],
                                               sc1_optimal_normal_3$r3[i]),
                                         seed = i + 14)
  sc1_normal_3_power[[i]] <- sim_skellam(K = 3, lambda2 = 5.15,
                                         n = sc1_optimal_normal_3$n[i],
                                         a = c(sc1_optimal_normal_3$a1[i],
                                               sc1_optimal_normal_3$a2[i],
                                               sc1_optimal_normal_3$a3[i]),
                                         r = c(sc1_optimal_normal_3$r1[i],
                                               sc1_optimal_normal_3$r2[i],
                                               sc1_optimal_normal_3$r3[i]),
                                         seed = i + 19)
}
sc2_normal_1_typeI   <- sim_skellam(K = 1, n = sc2_normal_1$n,
                                    a = sc2_normal_1$a1, r = sc2_normal_1$r1,
                                    seed = 25)
sc2_normal_1_power   <- sim_skellam(lambdaE = 5.15, K = 1, n = sc2_normal_1$n,
                                    a = sc2_normal_1$a1, r = sc2_normal_1$r1,
                                    seed = 26)
# Two-stage
sc2_normal_2_typeI   <- sc2_normal_2_power <- list()
sc2_optimal_normal_2 <- tibble::as_tibble(sc2_optimal_normal_2)
for (i in 1:3) {
  sc2_normal_2_typeI[[i]] <- sim_skellam(K = 2, n = sc2_optimal_normal_2$n[i],
                                         a = c(sc2_optimal_normal_2$a1[i],
                                               sc2_optimal_normal_2$a2[i]),
                                         r = c(sc2_optimal_normal_2$r1[i],
                                               sc2_optimal_normal_2$r2[i]),
                                         seed = i + 26)
  sc2_normal_2_power[[i]] <- sim_skellam(K = 2, lambda2 = 5.15,
                                         n = sc2_optimal_normal_2$n[i],
                                         a = c(sc2_optimal_normal_2$a1[i],
                                               sc2_optimal_normal_2$a2[i]),
                                         r = c(sc2_optimal_normal_2$r1[i],
                                               sc2_optimal_normal_2$r2[i]),
                                         seed = i + 29)
}
# Three-stage
sc2_normal_3_typeI        <- sc2_normal_3_power <- list()
sc2_optimal_normal_3      <- tibble::as_tibble(sc2_optimal_normal_3)
for (i in 1:3) {
  sc2_normal_3_typeI[[i]] <- sim_skellam(K = 3, n = sc2_optimal_normal_3$n[i],
                                         a = c(sc2_optimal_normal_3$a1[i],
                                               sc2_optimal_normal_3$a2[i],
                                               sc2_optimal_normal_3$a3[i]),
                                         r = c(sc2_optimal_normal_3$r1[i],
                                               sc2_optimal_normal_3$r2[i],
                                               sc2_optimal_normal_3$r3[i]),
                                         seed = i + 31)
  sc2_normal_3_power[[i]] <- sim_skellam(K = 3, lambdaE = 5.15,
                                         n = sc2_optimal_normal_3$n[i],
                                         a = c(sc2_optimal_normal_3$a1[i],
                                               sc2_optimal_normal_3$a2[i],
                                               sc2_optimal_normal_3$a3[i]),
                                         r = c(sc2_optimal_normal_3$r1[i],
                                               sc2_optimal_normal_3$r2[i],
                                               sc2_optimal_normal_3$r3[i]),
                                         seed = i + 34)
}
sc1_typeI     <- rbind(c(sc1_normal_1_typeI$power, sc1_normal_1_typeI$ess),
                       cbind(unlist(sc1_normal_2_typeI)[seq(1, 11, 2)],
                             unlist(sc1_normal_2_typeI)[seq(2, 12, 2)]),
                       cbind(unlist(sc1_normal_3_typeI)[seq(1, 9, 2)],
                             unlist(sc1_normal_3_typeI)[seq(2, 10, 2)]))
sc1_power     <- rbind(c(sc1_normal_1_power$power, sc1_normal_1_power$ess),
                       cbind(unlist(sc1_normal_2_power)[seq(1, 11, 2)],
                             unlist(sc1_normal_2_power)[seq(2, 12, 2)]),
                       cbind(unlist(sc1_normal_3_power)[seq(1, 9, 2)],
                             unlist(sc1_normal_3_power)[seq(2, 10, 2)]))
sc1           <- cbind(sc1_typeI, sc1_power)
rownames(sc1) <- NULL
colnames(sc1) <- c("P(H0)", "ESS(H0)", "P(H1)", "ESS(H1)")
sc2_typeI     <- rbind(c(sc2_normal_1_typeI$power, sc2_normal_1_typeI$ess),
                       cbind(unlist(sc2_normal_2_typeI)[seq(1, 5, 2)],
                             unlist(sc2_normal_2_typeI)[seq(2, 6, 2)]),
                       cbind(unlist(sc2_normal_3_typeI)[seq(1, 5, 2)],
                             unlist(sc2_normal_3_typeI)[seq(2, 6, 2)]))
sc2_power     <- rbind(c(sc2_normal_1_power$power, sc2_normal_1_power$ess),
                       cbind(unlist(sc2_normal_2_power)[seq(1, 5, 2)],
                             unlist(sc2_normal_2_power)[seq(2, 6, 2)]),
                       cbind(unlist(sc2_normal_3_power)[seq(1, 5, 2)],
                             unlist(sc2_normal_3_power)[seq(2, 6, 2)]))
sc2           <- cbind(sc2_typeI, sc2_power)
rownames(sc2) <- NULL
colnames(sc2) <- c("P(H0)", "ESS(H0)", "P(H1)", "ESS(H1)")
