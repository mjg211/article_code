##### Load required packages ###################################################

library(ggplot2)
library(mvtnorm)
library(skellam)
library(snowfall)
library(tibble)
library(tidyr)

##### Utility functions for normal approximation approach ######################

# Function to return covariance matrix of normal approximation design
normal_covariance  <- function(sqrt_I) {
  Sigma            <- diag(length(sqrt_I))
  for (k1 in seq_len(length(sqrt_I) - 1) + 1L) {
    vec            <- 1:(k1 - 1)
    Sigma[k1, vec] <- Sigma[vec, k1] <- sqrt_I[vec]/sqrt_I[k1]
  }
  Sigma
}

# Function to use to find a[k] of normal approximation design
normal_find_ak     <- function(ak, piAk, mean, Sigma, l, r) {
  mvtnorm::pmvnorm(l, c(r, ak), mean, sigma = Sigma)[1] - piAk
}

# Function to use to find r[k] of normal approximation design
normal_find_rk     <- function(rk, piRk, Sigma, a, u) {
  mvtnorm::pmvnorm(c(a, rk), u, sigma = Sigma)[1] - piRk
}

# Function to return information vector of normal approximation design
normal_information <- function(lambda1, lambda2, nv) {
  nv/(lambda1 + lambda2)
}

# Function to find required n for given spend and return the operating
# characteristics of the resulting design, in normal approximation approach
normal_obj_fn      <- function(piA, piR, K, beta, delta, min_Lambda0,
                               max_Lambda0, min_Lambda1, max_Lambda1,
                               lambda_ESS, n_fixed) {
  n            <- ceiling(stats::uniroot(normal_try_n, c(1e-10, n_fixed),
                                         K           = K,
                                         beta        = beta,
                                         delta       = delta,
                                         max_Lambda1 = max_Lambda1,
                                         piA         = piA,
                                         piR         = piR,
                                         type        = "root")$root)
  nv           <- n*(1:K)
  des          <- normal_try_n(n, K, beta, delta, max_Lambda1, piA, piR,
                               "design")
  Sigma        <- normal_covariance(nv)
  power_H0     <- normal_power(max_Lambda0, max_Lambda0, nv, des$a, des$r)
  power_H1_min <- normal_power(min_Lambda1, min_Lambda1 - delta, nv, des$a,
                               des$r)
  power_H1_max <- normal_power(max_Lambda1, max_Lambda1 - delta, nv, des$a,
                               des$r)
  ess0         <- normal_opchar(lambda_ESS, lambda_ESS, nv, des$a,
                                des$r)$ESS
  ess1         <- normal_opchar(lambda_ESS, lambda_ESS - delta, nv, des$a,
                                des$r)$ESS
  list(P_H0_min = power_H0,
       P_H0_max = power_H0,
       P_H1_min = power_H1_min,
       P_H1_max = power_H1_max,
       ESS0     = ess0,
       ESS1     = ess1,
       n        = n,
       a        = des$a,
       r        = des$r)
}

# Function to return operating characteristics of normal approximation design
normal_opchar      <- function(lambda1, lambda2, nv, a, r) {
  K             <- length(nv)
  ndiff         <- c(nv[1], nv[-1] - nv[-length(nv)])
  mean          <- (lambda1 - lambda2)*sqrt(normal_information(lambda1, lambda2,
                                                               nv))
  Sigma         <- normal_covariance(sqrt(nv))
  A             <- c(stats::pnorm(a[1], mean[1]), numeric(K - 1))
  R             <- c(stats::pnorm(r[1], mean[1], lower.tail = F),
                     numeric(K - 1))
  if (K > 1) {
    vec_old     <- 1
    if (K > 2) {
      for (k in seq_len(K - 2) + 1L) {
        vec_new <- c(vec_old, k)
        R[k]    <- mvtnorm::pmvnorm(c(a[vec_old], r[k]), c(r[vec_old], Inf),
                                    mean[vec_new],
                                    sigma = Sigma[vec_new, vec_new])[1]
        A[k]    <- mvtnorm::pmvnorm(c(a[vec_old], -Inf), c(r[vec_old], a[k]),
                                    mean[vec_new],
                                    sigma = Sigma[vec_new, vec_new])[1]
        vec_old <- vec_new
      }
    }
    R[K]        <- mvtnorm::pmvnorm(c(a[vec_old], r[K]), c(r[vec_old], Inf),
                                    mean, sigma = Sigma)[1]
    A[K]        <- 1 - sum(R) - sum(A[vec_old])
  }
  cum_S         <- cumsum(S <- A + R)
  N             <- 2*nv
  Med           <- ifelse(any(cum_S == 0.5), nv[which(cum_S == 0.5)] +
                            nv[which(cum_S == 0.5) + 1],
                          N[which(cum_S > 0.5)[1]])
  ESS           <- sum(N*S)
  list(P     = sum(R),
       ESS   = ESS,
       VSS   = sum(N^2*S) - ESS^2,
       Med   = Med,
       A     = A,
       R     = R,
       S     = S,
       cum_S = cum_S)
}

# Function to return power of normal approximation design
normal_power       <- function(lambda1, lambda2, nv, a, r, Sigma) {
  if (missing(Sigma)) {
    mean    <- sqrt(normal_information(lambda1, lambda2, nv))
    Sigma   <- normal_covariance(mean)
    mean    <- (lambda1 - lambda2)*mean
  } else {
    mean    <- (lambda1 - lambda2)*sqrt(normal_information(lambda1, lambda2,
                                                           nv))
  }
  P         <- stats::pnorm(r[1], mean[1], lower.tail = F)
  vec_old   <- 1L
  for (k in seq_len(length(nv) - 1) + 1L) {
    vec_new <- c(vec_old, k)
    P       <- P + mvtnorm::pmvnorm(c(a[vec_old], r[k]), c(r[vec_old], Inf),
                                    mean[vec_new],
                                    sigma = Sigma[vec_new, vec_new])[1]
    vec_old <- vec_new
  }
  P
}

# Function to evaluate power, a, and r for given n, in normal approximation
# design
normal_try_n       <- function(n, K, beta, delta, max_Lambda1, piA, piR, type) {
  nv          <- n*(1:K)
  a           <- r <- mean_H0 <- numeric(K)
  mean_H1     <- delta*sqrt(normal_information(max_Lambda1, max_Lambda1 - delta,
                                               nv))
  Sigma       <- normal_covariance(sqrt(nv))
  r[1]        <- stats::qnorm(1 - piR[1])
  if (K > 1) {
    a[1]      <- stats::qnorm(piA[1], mean_H1[1])
    if (a[1] > r[1]) {
      return(1)
    }
    vec_old   <- 1
    for (k in seq_len(K - 2) + 1L) {
      vec_new <- c(vec_old, k)
      a[k]    <- stats::uniroot(normal_find_ak, c(-20, 20),
                                piAk  = piA[k],
                                mean  = mean_H1[vec_new],
                                Sigma = Sigma[vec_new, vec_new],
                                l     = c(a[vec_old], -Inf),
                                r     = r[vec_old])$root
      r[k]    <- stats::uniroot(normal_find_rk, c(-20, 20),
                                piRk  = piR[k],
                                Sigma = Sigma[vec_new, vec_new],
                                a     = a[vec_old],
                                u     = c(r[vec_old], Inf))$root
      vec_old <- vec_new
      if (a[k] > r[k]) {
        return(1)
      }
    }
    if (normal_find_rk(-20, piR[K], Sigma, a[vec_old],
                       c(r[vec_old], Inf)) <= 0) {
      r[K]    <- a[K] <- -20
    } else {
      r[K]    <- a[K] <- stats::uniroot(normal_find_rk, c(-20, 20),
                                        piRk  = piR[K],
                                        Sigma = Sigma,
                                        a     = a[vec_old],
                                        u     = c(r[vec_old], Inf))$root
    }
  } else {
    a         <- r
  }
  P_H1        <- normal_power(max_Lambda1, max_Lambda1 - delta, nv, a, r, Sigma)
  if (type == "root") {
    P_H1 - (1 - beta)
  } else {
    list(a = a, r = r)
  }
}

##### Utility functions for exact approach #####################################

# Function to find required n for given spend and return the operating
# characteristics of the resulting design, in the exact approach
exact_obj_fn   <- function(piA, piR, K, beta, delta, min_Lambda0, max_Lambda0,
                           min_Lambda1, max_Lambda1, lambda_ESS, n_start) {
  n             <- n_start
  power_n       <- exact_try_n(n, K, delta, min_Lambda0, max_Lambda0,
                               min_Lambda1, max_Lambda1, piA, piR)
  if (power_n$power < 1 - beta) {
    while (power_n$power < 1 - beta) {
      n         <- n + 1
      power_n   <- exact_try_n(n, K, delta, min_Lambda0, max_Lambda0,
                               min_Lambda1, max_Lambda1, piA, piR)
    }
  } else if (power_n$power > 1 - beta) {
    while (power_n$power > 1 - beta) {
      n         <- n - 1
      power_n   <- exact_try_n(n, K, delta, min_Lambda0, max_Lambda0,
                               min_Lambda1, max_Lambda1, piA, piR)
    }
    n           <- n + 1
  }
  power_n       <- exact_try_n(n, K, delta, min_Lambda0, max_Lambda0,
                               min_Lambda1, max_Lambda1, piA, piR)
  power_H0_min  <- exact_power(min_Lambda0, min_Lambda0, n*(1:K), power_n$a,
                               power_n$r)
  power_H0_max  <- exact_power(max_Lambda0, max_Lambda0, n*(1:K), power_n$a,
                               power_n$r)
  power_H1_min  <- exact_power(min_Lambda1, min_Lambda1 - delta, n*(1:K),
                               power_n$a, power_n$r)
  power_H1_max  <- exact_power(max_Lambda1, max_Lambda1 - delta, n*(1:K),
                               power_n$a, power_n$r)
  ess0          <- exact_opchar(lambda_ESS, lambda_ESS, n*(1:K), power_n$a,
                                power_n$r)$ESS
  ess1          <- exact_opchar(lambda_ESS, lambda_ESS - delta, n*(1:K),
                                power_n$a, power_n$r)$ESS
  list(P_H0_min = power_H0_min,
       P_H0_max = power_H0_max,
       P_H1_min = power_H1_min,
       P_H1_max = power_H1_max,
       ESS0     = ess0,
       ESS1     = ess1,
       n        = n,
       a        = power_n$a,
       r        = power_n$r)
}

# Function to return operating characteristics of exact designs
exact_opchar   <- function(lambda1, lambda2, nv, a, r) {
  K                  <- length(nv)
  ndiff              <- c(nv[1], nv[-1] - nv[-K])
  mean1              <- ndiff*lambda1
  mean2              <- ndiff*lambda2
  A                  <- c(skellam::pskellam(a[1], mean1[1], mean2[1]),
                          numeric(K - 1))
  R                  <- c(skellam::pskellam(r[1] - 1, mean1[1], mean2[1],
                                            lower.tail = F), numeric(K - 1))
  if (K > 1) {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1),
                                          mean1[1], mean2[1])
    if (K > 2) {
      for (k in 2:(K - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] +
            pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) - l,
                              mean1[k], mean2[k])
          A[k]       <- A[k] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(a[k] - l, mean1[k], mean2[k])
          R[k]       <- R[k] + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(r[k] - l - 1, mean1[k],
                              mean2[k], lower.tail = F)
        }
      }
    }
    for (l in (a[K - 1] + 1):(r[K - 1] - 1)) {
      A[K]           <- A[K] + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(a[K] - l, mean1[K], mean2[K])
      R[K]           <- R[K] + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(r[K] - l - 1, mean1[K],
                          mean2[K], lower.tail = F)
    }
  }
  cum_S              <- cumsum(S <- A + R)
  N                  <- 2*nv
  Med                <- ifelse(any(cum_S == 0.5), nv[which(cum_S == 0.5)] +
                                 nv[which(cum_S == 0.5) + 1],
                               N[which(cum_S > 0.5)[1]])
  ESS                <- sum(N*S)
  list(P     = sum(R),
       ESS   = ESS,
       VSS   = sum(N^2*S) - ESS^2,
       Med   = Med,
       A     = A,
       R     = R,
       S     = S,
       cum_S = cum_S)
}

# Function to return power of exact designs
exact_power    <- function(lambda1, lambda2, nv, a, r) {
  K                  <- length(nv)
  ndiff              <- c(nv[1], nv[-1] - nv[-K])
  mean1              <- ndiff*lambda1
  mean2              <- ndiff*lambda2
  P                  <- skellam::pskellam(r[1] - 1, mean1[1], mean2[1],
                                          lower.tail = F)
  if (K > 1) {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1), mean1[1],
                                          mean2[1])
    if (K > 2) {
      for (k in 2:(K - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] +
            pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) - l,
                              mean1[k], mean2[k])
          P          <- P + pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::pskellam(r[k] - l - 1, mean1[k],
                              mean2[k], lower.tail = F)
        }
      }
    }
    for (l in (a[K - 1] + 1):(r[K - 1] - 1)) {
      P              <- P + pmf_C[[K - 1]][l - a[K - 1]]*
        skellam::pskellam(r[K] - l - 1, mean1[K],
                          mean2[K], lower.tail = F)
    }
  }
  P
}

# Function to evaluate power, e, and f for given n, in exact approach
exact_try_n    <- function(n, K, delta, min_Lambda0, max_Lambda0, min_Lambda1,
                           max_Lambda1, piA, piR) {
  a            <- r <- numeric(K)
  nv           <- n*(1:K)
  ndiff        <- c(nv[1], nv[-1] - nv[-length(nv)])
  mean_H0_min  <- ndiff*min_Lambda0
  mean_H0_max  <- ndiff*max_Lambda0
  mean1_H1_min <- ndiff*min_Lambda1
  mean2_H1_min <- ndiff*(min_Lambda1 - delta)
  mean1_H1_max <- ndiff*max_Lambda1
  mean2_H1_max <- ndiff*(max_Lambda1 - delta)
  for (k in 1:K) {
    vec        <- 1:k
    r[k]       <- skellam::qskellam(1 - piR[k], mean_H0_max[k], mean_H0_max[k])
    Rk         <- max(exact_typeI_k(mean_H0_min[vec], a[vec], r[vec]),
                      exact_typeI_k(mean_H0_max[vec], a[vec], r[vec]))
    if (Rk > piR[k]) {
      while (Rk > piR[k]) {
        r[k]   <- r[k] + 1
        Rk     <- max(exact_typeI_k(mean_H0_min[vec], a[vec], r[vec]),
                      exact_typeI_k(mean_H0_max[vec], a[vec], r[vec]))
      }
    } else if (Rk < piR[k]) {
      while (Rk < piR[k]) {
        r[k]   <- r[k] - 1
        Rk     <- max(exact_typeI_k(mean_H0_min[vec], a[vec], r[vec]),
                      exact_typeI_k(mean_H0_max[vec], a[vec], r[vec]))
      }
      r[k]     <- r[k] + 1
    }
    if (k < K) {
      a[k]     <- skellam::qskellam(piA[k], mean1_H1_max[k], mean2_H1_max[k])
      Ak       <- max(exact_typeII_k(mean1_H1_min[vec], mean2_H1_min[vec],
                                     a[1:k], r[1:k]),
                      exact_typeII_k(mean1_H1_max[vec], mean2_H1_max[vec],
                                     a[1:k], r[1:k]))
      if (Ak > piA[k]) {
        while (Ak > piA[k]) {
          a[k] <- a[k] - 1
          Ak   <- max(exact_typeII_k(mean1_H1_min[vec], mean2_H1_min[vec],
                                     a[1:k], r[1:k]),
                      exact_typeII_k(mean1_H1_max[vec], mean2_H1_max[vec],
                                     a[1:k], r[1:k]))
        }
      } else if (Ak < piA[k]) {
        while (Ak < piA[k]) {
          a[k] <- a[k] + 1
          Ak   <- max(exact_typeII_k(mean1_H1_min[vec], mean2_H1_min[vec],
                                     a[1:k], r[1:k]),
                      exact_typeII_k(mean1_H1_max[vec], mean2_H1_max[vec],
                                     a[1:k], r[1:k]))
        }
        a[k]   <- a[k] - 1
      }
    } else {
      a[K]     <- r[K] - 1
    }
  }
  list(power = min(exact_power(min_Lambda1, min_Lambda1 - delta, nv, a, r),
                   exact_power(max_Lambda1, max_Lambda1 - delta, nv, a, r)),
       a     = a,
       r     = r)
}

# Function to return type-I error-rate at stage k of exact design
exact_typeI_k  <- function(mean_H0, a, r) {
  stage_k            <- length(mean_H0)
  if (stage_k == 1) {
    R_stage_k        <- skellam::pskellam(r - 1, mean_H0, mean_H0,
                                          lower.tail = F)
  } else {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1), mean_H0[1],
                                          mean_H0[1])
    if (stage_k > 2) {
      for (k in 2:(stage_k - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] +
            pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) -  l,
                              mean_H0[k], mean_H0[k])
        }
      }
    }
    R_stage_k        <- 0
    for (l in (a[stage_k - 1] + 1):(r[stage_k - 1] - 1)) {
      R_stage_k      <- R_stage_k +
        pmf_C[[stage_k - 1]][l - a[stage_k - 1]]*
        skellam::pskellam(r[stage_k] - l - 1,
                          mean_H0[stage_k],
                          mean_H0[stage_k], lower.tail = F)
    }
  }
  R_stage_k
}

# Function to return type-II error-rate at stage k of exact design
exact_typeII_k <- function(mean1_H1, mean2_H1, a, r) {
  stage_k            <- length(mean1_H1)
  if (stage_k == 1) {
    A_stage_k        <- skellam::pskellam(a[1], mean1_H1, mean2_H1)
  } else {
    pmf_C            <- list()
    pmf_C[[1]]       <- skellam::dskellam((a[1] + 1):(r[1] - 1), mean1_H1[1],
                                          mean2_H1[1])
    if (stage_k > 2) {
      for (k in 2:(stage_k - 1)) {
        pmf_C[[k]]   <- numeric(r[k] - a[k] - 1)
        for (l in (a[k - 1] + 1):(r[k - 1] - 1)) {
          pmf_C[[k]] <- pmf_C[[k]] +
            pmf_C[[k - 1]][l - a[k - 1]]*
            skellam::dskellam((a[k] + 1):(r[k] - 1) - l,
                              mean1_H1[k], mean2_H1[k])
        }
      }
    }
    A_stage_k        <- 0
    for (l in (a[stage_k - 1] + 1):(r[stage_k - 1] - 1)) {
      A_stage_k      <- A_stage_k + pmf_C[[stage_k - 1]][l - a[stage_k - 1]]*
        skellam::pskellam(a[stage_k] - l,
                          mean1_H1[stage_k],
                          mean2_H1[stage_k])
    }
  }
  A_stage_k
}

##### Function to return normal approximation and exact designs ################

des_skellam <- function(K = 2L, alpha = 0.05, beta = 0.2, delta = 2.25,
                        min_Lambda0 = 15, max_Lambda0 = 30,
                        min_Lambda1 = min_Lambda0, max_Lambda1 = max_Lambda0,
                        lambda_ESS = 15, piAk, piRk, method = "exact") {
  if (method == "exact") {
    single_stage         <- exact_obj_fn(beta, alpha, 1, beta, delta,
                                         min_Lambda0, max_Lambda0, min_Lambda1,
                                         max_Lambda1, lambda_ESS, 10)
  } else {
    single_stage         <- normal_obj_fn(beta, alpha, 1, beta, delta,
                                          min_Lambda0, max_Lambda0, min_Lambda1,
                                          max_Lambda1, lambda_ESS, 1e10)
  }
  if (K > 1) {
    grid_piA             <- as.matrix(expand.grid(rep(list(piAk), K - 1)))
    grid_piA             <- grid_piA[rowSums(grid_piA) < beta, , drop = F]
    grid_piR             <- as.matrix(expand.grid(rep(list(piRk), K - 1)))
    grid_piR             <- grid_piR[rowSums(grid_piR) < alpha, , drop = F]
    designs              <- matrix(0, nrow = nrow(grid_piA)*nrow(grid_piR),
                                   ncol = 7 + 4*K)
    count                <- 1
    for (k1 in 1:nrow(grid_piA)) {
      for (k2 in 1:nrow(grid_piR)) {
        if (method == "exact") {
          design_k1k2    <- exact_obj_fn(c(grid_piA[k1, ],
                                           beta - sum(grid_piA[k1, ])),
                                         c(grid_piR[k2, ],
                                           alpha - sum(grid_piR[k2, ])),
                                         K, beta, delta, min_Lambda0,
                                         max_Lambda0, min_Lambda1, max_Lambda1,
                                         lambda_ESS,
                                         ceiling(1.1*single_stage$n/K))
        } else {
          design_k1k2    <- normal_obj_fn(c(grid_piA[k1, ],
                                            beta - sum(grid_piA[k1, ])),
                                          c(grid_piR[k2, ],
                                            alpha - sum(grid_piR[k2, ])),
                                          K, beta, delta, min_Lambda0,
                                          max_Lambda0, min_Lambda0, max_Lambda1,
                                          lambda_ESS, 1.1*single_stage$n)
        }
        designs[count, ] <- c(c(grid_piA[k1, ], beta - sum(grid_piA[k1, ])),
                              c(grid_piR[k2, ], alpha - sum(grid_piR[k2, ])),
                              design_k1k2$P_H0_min, design_k1k2$P_H0_max,
                              design_k1k2$P_H1_min, design_k1k2$P_H1_max,
                              design_k1k2$n, design_k1k2$ESS0,
                              design_k1k2$ESS1, design_k1k2$a, design_k1k2$r)
        count            <- count + 1
        message("k2 = ", k2)
      }
      message("k1 = ", k1)
    }
  } else {
    designs              <- matrix(c(beta, alpha, single_stage$P_H0_min,
                                     single_stage$P_H0_max,
                                     single_stage$P_H1_min,
                                     single_stage$P_H1_max, single_stage$n,
                                     single_stage$ESS0, single_stage$ESS1,
                                     single_stage$a, single_stage$r), 1)
  }
  colnames(designs)      <- c(paste("piA", 1:K, sep = ""),
                              paste("piR", 1:K, sep = ""),
                              "P(min_Lambda0,min_Lambda0)",
                              "P(max_Lambda0,max_Lambda0)",
                              "P(min_Lambda0,min_Lambda0 - delta)",
                              "P(max_Lambda0,max_Lambda0 - delta)", "n",
                              "ESS0", "ESS1",
                              paste("a", 1:K, sep = ""),
                              paste("r", 1:K, sep = ""))
  tibble::as_tibble(designs)
}

##### Function to simulate trials given a design ###############################

sim_skellam <- function(lambda1 = 15, lambda2 = lambda1, K = 2L, n = 41L,
                        a = c(0.61, 1.74), r = c(1.96, 1.74),
                        method = "normal", replicates = 100000,
                        seed = Sys.time()) {
  set.seed(seed)
  nv                     <- n*(1:K)
  reject                 <- sample_size <- integer(replicates)
  Tk                     <- numeric(replicates)
  for (rep in 1:replicates) {
    Y1k                  <- Y2k <- 0
    for (k in 1:K) {
      Y1k                <- Y1k + rpois(1, n*lambda1)
      Y2k                <- Y2k + rpois(1, n*lambda2)
      if (method == "normal") {
        hat_lambda1      <- Y1k/nv[k]
        hat_lambda2      <- Y2k/nv[k]
        Tk[rep]          <- (hat_lambda1 - hat_lambda2)*
                              sqrt(nv[k]/(hat_lambda1 + hat_lambda2))
        if (is.na(Tk[rep])) {
          Tk[rep]        <- 0
        }
      } else {
        Tk[rep]          <- Y1k - Y2k
      }
      if (Tk[rep] >= r[k]) {
        reject[rep]      <- 1
        sample_size[rep] <- 2*k*n
        break
      } else if (Tk[rep] < a[k]) {
        sample_size[rep] <- 2*k*n
        break
      }
    }
  }
  list(power = mean(reject),
       ess   = mean(sample_size),
       Tk    = Tk)
}

##### Example power curves (Supplementary Figure 1) ############################

typeI      <- power <- numeric(1000)
lambda1    <- seq(15, 30, length.out = 1000)
# WARNING: The following loop has a long execution time
for (i in 1:1000) {
  typeI[i] <- exact_power(lambda1[i], lambda1[i], 42*(1:2), c(40, 111),
                          c(118, 112))
  power[i] <- exact_power(lambda1[i], lambda1[i] - 2.25, 42*(1:2), c(40, 111),
                          c(118, 112))
}
data       <- tibble::tibble(lambda1 = rep(lambda1, 2),
                             P       = c(typeI, 1 - power),
                             type    = factor(rep(c("  Type-I error-rate  ",
                                                    "  Type-II error-rate  "),
                                                  each = 1000)))
s_fig_1    <- ggplot2::ggplot(data,
                              ggplot2::aes(x      = lambda1,
                                           y      = P,
                                           colour = type)) +
  ggplot2::geom_line() +
  ggplot2::xlab(expression(lambda[1])) +
  ggplot2::ylab("Probability") +
  ggplot2::ylim(c(0, 0.2)) +
  ggplot2::theme_grey() +
  ggplot2::theme(legend.title    = ggplot2::element_blank(),
                 legend.position = "bottom")
ggplot2::ggsave("s_fig_1.pdf", plot = s_fig_1, device = "pdf", width = 5,
                height = 4, units = "in")

##### Grid search for maximal type-I and type-II error-rates (see ##############
##### Supplementary Material) ##################################################

typeI               <- function(lambda1, nv, a, r) {
  exact_power(lambda1, lambda1, nv, a, r)
}
power               <- function(lambda1, nv, a, r) {
  exact_power(lambda1, lambda1 - 2.25, nv, a, r)
}
max_typeI           <- function(nv, a, r) {
  maxima  <- stats::optimize(typeI, c(15, 30), maximum = T,
                             nv = nv, a = a, r = r)
  lambdas <- c(15, maxima$maximum, 30)
  typeI   <- c(exact_power(15, 15, nv, a, r),
               maxima$objective,
               exact_power(30, 30, nv, a, r))
  c(lambdas[which.max(typeI)], max(typeI))
}
min_power           <- function(nv, a, r) {
  minima  <- stats::optimize(power, c(15, 30),
                             nv = nv, a = a, r  = r)
  lambdas <- c(15, minima$minimum, 30)
  power   <- c(exact_power(15, 15 - 2.25, nv, a, r),
               minima$objective,
               exact_power(30, 30 - 2.25, nv, a, r))
  c(lambdas[which.min(power)], min(power))
}
wrapper             <- function(n) {
  stage_1_bounds          <- expand.grid(a1 = (-2*n):(2*n - 2),
                                         r1 = (-2*n + 2):(2*n))
  stage_1_bounds          <-
    stage_1_bounds[which(stage_1_bounds[, 1] < stage_1_bounds[, 2] - 1), ]
  bounds                  <- expand.grid(s1 = 1:nrow(stage_1_bounds),
                                         r2 = (-2*n):(2*n))
  bounds                  <- cbind(stage_1_bounds[bounds[, 1], ], bounds[, 2])
  location_max_typeI      <- location_min_power <- numeric(nrow(bounds))
  for (i in 1:nrow(bounds)) {
    location_max_typeI[i] <- max_typeI(n*(1:2),
                                       c(bounds[i, 1], bounds[i, 3] - 1),
                                       c(bounds[i, 2], bounds[i, 3]))[1]
    location_min_power[i] <- min_power(n*(1:2),
                                       c(bounds[i, 1], bounds[i, 3] - 1),
                                       c(bounds[i, 2], bounds[i, 3]))[1]
  }
  results                 <- cbind(location_max_typeI, location_min_power)
  write.csv(cbind(location_max_typeI, location_min_power),
            paste("location_max_min_AHIeg_", n, ".csv", sep = ""))
  list(location_max_typeI = location_max_typeI,
       location_min_power = location_min_power)
}
sfInit(parallel = T, cpus = 6)
sfLibrary(skellam)
sfLibrary(stats)
sfExport("exact_power", "max_typeI", "min_power", "power", "typeI")
# WARNING: The following line has a long execution time
results             <- sfLapply(1:20, wrapper)
sfStop()
num_designs         <- 0
for (n in 1:20) {
  num_designs       <- num_designs + length(results[[n]]$location_max_typeI)
}
data_mat            <- matrix(0, num_designs, 2)
counter             <- 1
for (n in 1:5) {
  range             <- counter:(counter +
                                  length(results[[n]]$location_max_typeI) - 1)
  data_mat[range, ] <- cbind(results[[n]]$location_max_typeI,
                             results[[n]]$location_min_power)
  counter           <- counter + length(results[[n]]$location_max_typeI)
}
data                <- tibble::tibble(lambda_max_typeI = data_mat[, 1],
                                      lambda_min_power = data_mat[, 2])
length(which(data$lambda_max_typeI == 15))
length(which(!((data$lambda_max_typeI == 15) |
                 (data$lambda_max_typeI == 30))))
length(which(data$lambda_max_typeI == 30))
length(which(data$lambda_min_power == 15))
length(which(!((data$lambda_min_power == 15) |
                 (data$lambda_min_power == 30))))
length(which(data$lambda_min_power == 30))

##### Scenario 1: alpha = 0.05, beta = 0.2, delta = 2.25, ######################
##### Lambda0 = Lambda1 = [15, 30] #############################################

# Single-stage
sc1_exact_1                    <- des_skellam(K      = 1L)
sc1_normal_1                   <- des_skellam(K      = 1L,
                                              method = "normal")
# Two-stage (error-spending)
sc1_exact_2                    <- des_skellam(piAk   = seq(0.02, 0.18, 0.02),
                                              piRk   = seq(0.005, 0.045, 0.005))
sc1_normal_2                   <- des_skellam(piAk   = seq(0.02, 0.18, 0.02),
                                              piRk   = seq(0.005, 0.045, 0.005),
                                              method = "normal")
# Two-stage (globally optimal)
sc1_global_2                   <-
  tibble::as_tibble(exact_des_two_stage_cpp(0.05, 0.2, 2.25, 1,
                                            ceiling(0.75*sc1_exact_1$n),
                                            15, 30, 15, 30, 15, 1))
colnames(sc1_global_2)         <- c("n", "a1", "r1", "r2", "typeI_min",
                                    "typeI_max", "1 - typeII_min",
                                    "1 - typeII_max", "ESS0", "ESS1")
# Three-stage
# WARNING: The following line has a long execution time
sc1_exact_3                    <-
  des_skellam(K      = 3L,
              piAk   = c(0.03, 0.06, 0.09, 0.12),
              piRk   = c(0.01, 0.015, 0.02, 0.025, 0.03, 0.035))
sc1_normal_3                   <-
  des_skellam(K      = 3L,
              piAk   = c(0.03, 0.06, 0.09, 0.12),
              piRk   = c(0.01, 0.015, 0.02, 0.025, 0.03, 0.035),
              method = "normal")
# Find the optimal designs
w                              <- cbind(c(1, 0, 1/2, 1/2,   0, 1/3),
                                        c(0, 1, 1/2,   0, 1/2, 1/3),
                                        c(0, 0,   0, 1/2, 1/2, 1/3))
sc1_optimal_exact_2            <- sc1_optimal_normal_2 <- matrix(0, 6, 7 + 4*2)
sc1_optimal_global_2           <- matrix(0, 6, 10)
sc1_optimal_exact_3            <- sc1_optimal_normal_3 <- matrix(0, 6, 7 + 4*3)
sc1_optimal_exact_2            <- cbind(sc1_optimal_exact_2, 0)
for (i in 1:nrow(w)) {
  # Find the scores
  scores_exact_2               <- w[i, 1]*sc1_exact_2$ESS0 +
    w[i, 2]*sc1_exact_2$ESS1 +
    w[i, 3]*sc1_exact_2$n*2*2
  scores_exact_3               <- w[i, 1]*sc1_exact_3$ESS0 +
    w[i, 2]*sc1_exact_3$ESS1 +
    w[i, 3]*sc1_exact_3$n*2*3
  scores_global_2              <- w[i, 1]*sc1_global_2$ESS0 +
    w[i, 2]*sc1_global_2$ESS1 +
    w[i, 3]*sc1_global_2$n*2*2
  scores_normal_2              <- w[i, 1]*sc1_normal_2$ESS0 +
    w[i, 2]*sc1_normal_2$ESS1 +
    w[i, 3]*sc1_normal_2$n*2*2
  scores_normal_3              <- w[i, 1]*sc1_normal_3$ESS0 +
    w[i, 2]*sc1_normal_3$ESS1 +
    w[i, 3]*sc1_normal_3$n*2*3
  # Store which is optimal
  sc1_optimal_exact_2[i, ]     <-
    c(as.numeric(sc1_exact_2[which.min(scores_exact_2), ]),
      min(scores_exact_2)/min(scores_global_2))
  sc1_optimal_exact_3[i, ]     <-
    as.numeric(sc1_exact_3[which.min(scores_exact_3), ])
  sc1_optimal_global_2[i, ]    <-
    as.numeric(sc1_global_2[which.min(scores_global_2), ])
  sc1_optimal_normal_2[i, ]    <-
    as.numeric(sc1_normal_2[which.min(scores_normal_2), ])
  sc1_optimal_normal_3[i, ]    <-
    as.numeric(sc1_normal_3[which.min(scores_normal_3), ])
}
sc1_optimal_exact_2            <-
  tibble::as_tibble(cbind(w, sc1_optimal_exact_2))
sc1_optimal_global_2           <-
  tibble::as_tibble(cbind(w, sc1_optimal_global_2))
sc1_optimal_exact_3            <-
  tibble::as_tibble(cbind(w, sc1_optimal_exact_3))
sc1_optimal_normal_2           <-
  tibble::as_tibble(cbind(w, sc1_optimal_normal_2))
sc1_optimal_normal_3           <-
  tibble::as_tibble(cbind(w, sc1_optimal_normal_3))
colnames(sc1_optimal_exact_2)  <- c("w1", "w2", "w3", colnames(sc1_exact_2),
                                    "O_NO/O_G")
colnames(sc1_optimal_global_2) <- c("w1", "w2", "w3", colnames(sc1_global_2))
colnames(sc1_optimal_exact_3)  <- c("w1", "w2", "w3", colnames(sc1_exact_3))
colnames(sc1_optimal_normal_2) <- c("w1", "w2", "w3", colnames(sc1_normal_2))
colnames(sc1_optimal_normal_3) <- c("w1", "w2", "w3", colnames(sc1_normal_3))
# Perform a quick confirmation that each of the error-rates of the exact designs
# is actually controlled
max_typeI_exact_2              <- min_power_exact_2 <- max_typeI_exact_3 <-
  min_power_exact_3 <- max_typeI_global_2 <-
  min_power_global_2 <- numeric(6)
for (i in 1:6) {
  max_typeI_exact_2[i]         <- max_typeI(sc1_optimal_exact_2$n[i]*(1:2),
                                            c(sc1_optimal_exact_2$a1[i],
                                              sc1_optimal_exact_2$a2[i]),
                                            c(sc1_optimal_exact_2$r1[i],
                                              sc1_optimal_exact_2$r2[i]))[2]
  min_power_exact_2[i]         <- min_power(sc1_optimal_exact_2$n[i]*(1:2),
                                            c(sc1_optimal_exact_2$a1[i],
                                              sc1_optimal_exact_2$a2[i]),
                                            c(sc1_optimal_exact_2$r1[i],
                                              sc1_optimal_exact_2$r2[i]))[2]
  max_typeI_global_2[i]        <- max_typeI(sc1_optimal_global_2$n[i]*(1:2),
                                            c(sc1_optimal_global_2$a1[i],
                                              sc1_optimal_global_2$a2[i]),
                                            c(sc1_optimal_global_2$r1[i],
                                              sc1_optimal_global_2$r2[i]))[2]
  min_power_global_2[i]        <- min_power(sc1_optimal_global_2$n[i]*(1:2),
                                            c(sc1_optimal_global_2$a1[i],
                                              sc1_optimal_global_2$a2[i]),
                                            c(sc1_optimal_global_2$r1[i],
                                              sc1_optimal_global_2$r2[i]))[2]
  max_typeI_exact_3[i]         <- max_typeI(sc1_optimal_exact_3$n[i]*(1:3),
                                            c(sc1_optimal_exact_3$a1[i],
                                              sc1_optimal_exact_3$a2[i],
                                              sc1_optimal_exact_3$a3[i]),
                                            c(sc1_optimal_exact_3$r1[i],
                                              sc1_optimal_exact_3$r2[i],
                                              sc1_optimal_exact_3$r3[i]))[2]
  min_power_exact_3[i]         <- min_power(sc1_optimal_exact_3$n[i]*(1:3),
                                            c(sc1_optimal_exact_3$a1[i],
                                              sc1_optimal_exact_3$a2[i],
                                              sc1_optimal_exact_3$a3[i]),
                                            c(sc1_optimal_exact_3$r1[i],
                                              sc1_optimal_exact_3$r2[i],
                                              sc1_optimal_exact_3$r3[i]))[2]
  print(i)
}
max(c(max_typeI_exact_2, max_typeI_exact_3, max_typeI_global_2)) # < 0.05
min(c(min_power_exact_2, min_power_exact_3, min_power_global_2)) # > 0.8

##### Scenario 2: alpha = 0.05, beta = 0.1, delta = 1.1, #######################
##### Lambda0 = Lambda1 = 6.25 #################################################

# Single-stage
sc2_exact_1                    <- des_skellam(K           = 1L,
                                              beta        = 0.1,
                                              delta       = 1.1,
                                              min_Lambda0 = 6.25,
                                              max_Lambda0 = 6.25,
                                              lambda_ESS  = 6.25)
sc2_normal_1                   <- des_skellam(K           = 1L,
                                              beta        = 0.1,
                                              delta       = 1.1,
                                              min_Lambda0 = 6.25,
                                              max_Lambda0 = 6.25,
                                              lambda_ESS  = 6.25,
                                              method      = "normal")
# Two-stage
sc2_exact_2                    <- des_skellam(beta        = 0.1,
                                              delta       = 1.1,
                                              min_Lambda0 = 6.25,
                                              max_Lambda0 = 6.25,
                                              lambda_ESS  = 6.25,
                                              piAk   = seq(0.01, 0.09, 0.01),
                                              piRk   = seq(0.005, 0.045, 0.005))
sc2_normal_2                   <- des_skellam(beta        = 0.1,
                                              delta       = 1.1,
                                              min_Lambda0 = 6.25,
                                              max_Lambda0 = 6.25,
                                              lambda_ESS  = 6.25,
                                              piAk   = seq(0.01, 0.09, 0.01),
                                              piRk   = seq(0.005, 0.045, 0.005),
                                              method = "normal")
# Three-stage
# WARNING: The following line has a long execution time
sc2_exact_3                    <-
  des_skellam(K           = 3L,
              beta        = 0.1,
              delta       = 1.1,
              min_Lambda0 = 6.25,
              max_Lambda0 = 6.25,
              lambda_ESS  = 6.25,
              piAk        = c(0.01, 0.03, 0.05, 0.07),
              piRk        = c(0.01, 0.015, 0.02, 0.025, 0.03, 0.035))
sc2_normal_3                   <-
  des_skellam(K           = 3L,
              beta        = 0.1,
              delta       = 1.1,
              min_Lambda0 = 6.25,
              max_Lambda0 = 6.25,
              lambda_ESS  = 6.25,
              piAk        = c(0.01, 0.03, 0.05, 0.07),
              piRk        = c(0.01, 0.015, 0.02, 0.025, 0.03, 0.035),
              method      = "normal")
# Find the optimal designs
w                              <- cbind(c(1, 0, 1/2, 1/2,   0, 1/3),
                                        c(0, 1, 1/2,   0, 1/2, 1/3),
                                        c(0, 0,   0, 1/2, 1/2, 1/3))
sc2_optimal_exact_2            <- sc2_optimal_normal_2 <- matrix(0, 6, 7 + 4*2)
sc2_optimal_exact_3            <- sc2_optimal_normal_3 <- matrix(0, 6, 7 + 4*3)
for (i in 1:nrow(w)) {
  # Find the scores
  scores_exact_2               <- w[i, 1]*sc2_exact_2$ESS0 +
    w[i, 2]*sc2_exact_2$ESS1 +
    w[i, 3]*sc2_exact_2$n*2*2
  scores_exact_3               <- w[i, 1]*sc2_exact_3$ESS0 +
    w[i, 2]*sc2_exact_3$ESS1 +
    w[i, 3]*sc2_exact_3$n*2*3
  scores_normal_2              <- w[i, 1]*sc2_normal_2$ESS0 +
    w[i, 2]*sc2_normal_2$ESS1 +
    w[i, 3]*sc2_normal_2$n*2*2
  scores_normal_3              <- w[i, 1]*sc2_normal_3$ESS0 +
    w[i, 2]*sc2_normal_3$ESS1 +
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
sc2_optimal_exact_2            <-
  tibble::as_tibble(cbind(w, sc2_optimal_exact_2))
sc2_optimal_exact_3            <-
  tibble::as_tibble(cbind(w, sc2_optimal_exact_3))
sc2_optimal_normal_2           <-
  tibble::as_tibble(cbind(w, sc2_optimal_normal_2))
sc2_optimal_normal_3           <-
  tibble::as_tibble(cbind(w, sc2_optimal_normal_3))
colnames(sc2_optimal_exact_2)  <- c("w1", "w2", "w3", colnames(sc2_exact_2))
colnames(sc2_optimal_exact_3)  <- c("w1", "w2", "w3", colnames(sc2_exact_3))
colnames(sc2_optimal_normal_2) <- c("w1", "w2", "w3", colnames(sc2_normal_2))
colnames(sc2_optimal_normal_3) <- c("w1", "w2", "w3", colnames(sc2_normal_3))

##### Simulation study #########################################################

wrapper <- function(i) {
  design    <- des_skellam(delta       = lambda_frac[i, 1]*lambda_frac[i, 2],
                           min_Lambda0 = lambda_frac[i, 1],
                           max_Lambda0 = lambda_frac[i, 1],
                           lambda_ESS  = lambda_frac[i, 1],
                           piAk        = 0.1,
                           piRk        = 0.025,
                           method      = "normal")
  sim_typeI <- sim_skellam(lambda1    = lambda_frac[i, 1],
                           lambda2    = lambda_frac[i, 1],
                           K          = 2L,
                           n          = design$n,
                           a          = c(design$a1, design$a2),
                           r          = c(design$r1, design$r2),
                           method     = "normal",
                           replicates = 100000,
                           seed       = i)
  sim_power <- sim_skellam(lambda1    = lambda_frac[i, 1],
                           lambda2    = lambda_frac[i, 1]*
                             (1 - lambda_frac[i, 2]),
                           K          = 2L,
                           n          = design$n,
                           a          = c(design$a1, design$a2),
                           r          = c(design$r1, design$r2),
                           method     = "normal",
                           replicates = 100000,
                           seed       = i)
  print(i)
  return(c(design$n, design$`P(max_Lambda0,max_Lambda0)`,
           design$`P(max_Lambda0,max_Lambda0 - delta)`,
           design$ESS0, design$ESS1, sim_typeI$power, sim_power$power,
           sim_typeI$ess, sim_power$ess))
}

lambda_frac  <- expand.grid(lambda = seq(1, 10, 0.5),
                            frac   = seq(0.25, 0.75, 0.025))
delta        <- lambda_frac[, 1]*lambda_frac[, 2]
suppressMessages(sfInit(parallel = T, cpus = 4))
suppressMessages(sfLibrary(mvtnorm))
suppressMessages(sfLibrary(stats))
sfExport("lambda_frac", "des_skellam", "sim_skellam", "normal_covariance",
         "normal_find_ak", "normal_find_rk", "normal_information",
         "normal_obj_fn", "normal_opchar", "normal_power", "normal_try_n",
         "max_typeI", "typeI", "exact_power")
simulations_new  <- lapply(1:nrow(lambda_frac), wrapper)
simulations2 <- matrix(unlist(simulations_new), nrow(lambda_frac), 9,
                       byrow = T)
sfStop()
est_vs_sim   <- tibble::tibble(lambda           = rep(lambda_frac[, 1], 4),
                               `delta %`        =
                                 rep(100*lambda_frac[, 2], 4),
                               n                = rep(simulations2[, 1], 4),
                               estimate         = c(simulations2[, 2],
                                                    1 - simulations2[, 3],
                                                    simulations2[, 4],
                                                    simulations2[, 5]),
                               empirical        = c(simulations2[, 6],
                                                    1 - simulations2[, 7],
                                                    simulations2[, 8],
                                                    simulations2[, 9]),
                               difference       = estimate - empirical,
                               type             =
                                 factor(rep(c("type-I", "type-II", "ESS0",
                                              "ESS1"),
                                            each = nrow(lambda_frac))))

fig_1a <- ggplot2::ggplot(dplyr::filter(est_vs_sim, type == "type-I"),
                          ggplot2::aes(x = type, y = difference)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x  = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::ylab("Difference: theoretical and empirical type-I error-rates")
ggplot2::ggsave("fig_1a.pdf", plot = fig_1a, device = "pdf", width = 2.5,
                height = 4, units = "in")

fig_1b <- ggplot2::ggplot(dplyr::filter(est_vs_sim, type == "type-II"),
                          ggplot2::aes(x = type, y = difference)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x  = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::ylab("Difference: theoretical and empirical type-II error-rates")
ggplot2::ggsave("fig_1b.pdf", plot = fig_1b, device = "pdf", width = 2.5,
                height = 4, units = "in")



