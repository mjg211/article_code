#install.packages("ggplot")
#install.packages("ggthemes")
#install.packages("mvtnorm")
##install.packages("patchwork")
#install.packages("snowfall")
#install.packages("tibble")
#install.packages("tidyr")
library(ggplot2)
library(ggthemes)
library(patchwork)
library(snowfall)
library(tibble)
library(tidyr)

##### Futility stopping only ###################################################

# Helper functions
root_f2  <- function(f2, f1, CovZ, alpha) {
  pmvnorm(c(f1, f2), c(Inf, Inf), c(0, 0), sigma = CovZ) - alpha
}

root_n02 <- function(n02, f1, rho, sigma2, ratio, alpha, power, delta, CovZ) {
  n01        <- rho*n02
  sqrtI      <- sqrt(1/c(sigma2[1]/n01 + sigma2[2]/(ratio*n01),
                         sigma2[1]/n02 + sigma2[2]/(ratio*n02)))
  CovZ[1, 2] <- CovZ[2, 1] <- sqrtI[1]/sqrtI[2]
  f2         <- uniroot(root_f2, c(-10, 10), f1 = f1, CovZ = CovZ,
                        alpha = alpha)$root
  pmvnorm(c(f1, f2), c(Inf, Inf), delta*sqrtI, sigma = CovZ) - power
}

opchar   <- function(pars, sigma2, ratio, alpha, power, delta, CovZ) {
  f1      <- pars[1]
  rho     <- pars[2]
  n02     <- uniroot(root_n02, c(1, 1e5), f1 = f1, rho = rho,
                     sigma2 = sigma2, ratio = ratio, alpha = alpha,
                     power = power, delta = delta, CovZ = CovZ)$root
  n01     <- rho*n02
  n11     <- ratio*n01
  n       <- c(n01 + n11, n02 + ratio*n02)
  S0      <- pnorm(f1)
  S0      <- c(S0, 1 - S0)
  S1      <- pnorm(f1 - delta/sqrt(sigma2[1]/n01 + sigma2[2]/n11))
  S1      <- c(S1, 1 - S1)
  ess0    <- sum(n*S0)
  ess1    <- sum(n*S1)
  if (S0[1] < 0.5) {
    med0  <- n[2]
  } else if (S0[1] == 0.5) {
    med0  <- n[1] + 0.5*n[2]
  } else {
    med0  <- n[1]
  }
  if (S1[1] < 0.5) {
    med1  <- n[2]
  } else if (S1[1] == 0.5) {
    med1  <- n[1] + 0.5*n[2]
  } else {
    med1  <- n[1]
  }
  c(ess0, ess1, med0, med1, sqrt(sum(n^2*S0) - ess0^2),
    sqrt(sum(n^2*S1) - ess1^2), n[2], S1[1])
}

wrapper  <- function(i) {
  c(designs[i, 1:2], opchar(designs[i, 1:2], sigma2, ratio, alpha, power,
                            delta, CovZ))
}

# Example parameters
sigma2                     <- c(1, 1)
delta                      <- 0.3
ratio                      <- 1
alpha                      <- 0.025
power                      <- 0.9
CovZ                       <- diag(1, 2, 2)
# Combinations of f1 and rho to consider
designs                    <-
  as.matrix(expand.grid(f1  = seq(-3.000, 1.950, 0.0025),
                        rho = seq( 0.010, 0.999, 0.0025)))
designs                    <- rbind(designs, c(-Inf, 1)) # Fixed-sample
# Evaluate them using parallelisation
sfInit(parallel = TRUE, cpus = 7)
sfLibrary(mvtnorm)
sfLibrary(stats)
sfExport("designs", "sigma2", "delta", "ratio", "alpha", "power", "CovZ",
         "opchar", "root_f2", "root_n02")
output                     <- sfLapply(1:nrow(designs), wrapper)
designs                    <- matrix(unlist(output), nrow(designs), 10,
                                     byrow = TRUE)
designs[nrow(designs), 1]  <- NA_real_
designs[nrow(designs), 10] <- 0
designs                    <- tibble(f1            = designs[, 1],
                                     rho           = designs[, 2],
                                     `ESS(0)`      = designs[, 3],
                                     `ESS(delta)`  = designs[, 4],
                                     `MSS(0)`      = designs[, 5],
                                     `MSS(delta)`  = designs[, 6],
                                     `SDSS(0)`     = designs[, 7],
                                     `SDSS(delta)` = designs[, 8],
                                     `max(n)`      = designs[, 9],
                                     `PIE(0)`      = 0,
                                     `PIE(delta)`  = designs[, 10])
# Search for what is optimal for each w under H0 and H1
optimal0_SDSS              <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Null",
                                     optimality = "SDSS")
optimal1_SDSS              <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Alternative",
                                     optimality = "SDSS")
optimal0_MSS               <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Null",
                                     optimality = "MSS")
optimal1_MSS               <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Alternative",
                                     optimality = "MSS")
optimal0_PIE               <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Null",
                                     optimality = "PIE")
optimal1_PIE               <- tibble(w          = seq(0, 1, 0.01),
                                     f1         = NA_real_,
                                     rho        = NA_real_,
                                     scenario   = "Alternative",
                                     optimality = "PIE")
for (i in 1:nrow(optimal0_SDSS)) {
  score_0                  <- optimal0_SDSS$w[i]*designs$`ESS(0)` +
    (1 - optimal0_SDSS$w[i])*designs$`SDSS(0)`
  optimal_0                <- which.min(score_0)
  optimal0_SDSS$f1[i]      <- designs$f1[optimal_0]
  optimal0_SDSS$rho[i]     <- designs$rho[optimal_0]
  score_1                  <- optimal1_SDSS$w[i]*designs$`ESS(delta)` +
    (1 - optimal1_SDSS$w[i])*designs$`SDSS(delta)`
  optimal_1                <- which.min(score_1)
  optimal1_SDSS$f1[i]      <- designs$f1[optimal_1]
  optimal1_SDSS$rho[i]     <- designs$rho[optimal_1]

  score_0                  <- optimal0_MSS$w[i]*designs$`ESS(0)` +
    (1 - optimal0_MSS$w[i])*designs$`MSS(0)`
  optimal_0                <- which.min(score_0)
  optimal0_MSS$f1[i]       <- designs$f1[optimal_0]
  optimal0_MSS$rho[i]      <- designs$rho[optimal_0]
  score_1                  <- optimal1_MSS$w[i]*designs$`ESS(delta)` +
    (1 - optimal1_MSS$w[i])*designs$`MSS(delta)`
  optimal_1                <- which.min(score_1)
  optimal1_MSS$f1[i]       <- designs$f1[optimal_1]
  optimal1_MSS$rho[i]      <- designs$rho[optimal_1]

  score_0                  <-
    optimal0_PIE$w[i]*designs$`ESS(0)`/max(designs$`ESS(0)`)
  optimal_0                <- which.min(score_0)
  optimal0_PIE$f1[i]       <- designs$f1[optimal_0]
  optimal0_PIE$rho[i]      <- designs$rho[optimal_0]
  score_1                  <-
    optimal1_PIE$w[i]*designs$`ESS(delta)`/max(designs$`ESS(delta)`) +
    (1 - optimal1_PIE$w[i])*designs$`PIE(delta)`/max(designs$`PIE(delta)`)
  optimal_1                <- which.min(score_1)
  optimal1_PIE$f1[i]       <- designs$f1[optimal_1]
  optimal1_PIE$rho[i]      <- designs$rho[optimal_1]

  if (i%%100 == 0) {
    message(i)
  }
}
optimal                    <- rbind(optimal0_SDSS, optimal1_SDSS, optimal0_MSS,
                                   optimal1_MSS, optimal0_PIE, optimal1_PIE)
optimal                    <- gather(optimal, key = "key", value = "value",
                                     f1:rho)
optimal$scenario[which(optimal$scenario == "Null")]        <-
  paste0("theta == 0")
optimal$scenario[which(optimal$scenario == "Alternative")] <-
  paste0("theta == delta")
optimal$optimality[which(optimal$optimality == "SDSS")]    <-
  paste0("italic(SDSS)(theta)")
optimal$optimality[which(optimal$optimality == "MSS")]     <-
  paste0("italic(MSS)(theta)")
optimal$optimality[which(optimal$optimality == "PIE")]     <-
  paste0("italic(PIE)(theta)")
optimal$optimality                                         <-
  factor(optimal$optimality, c("italic(SDSS)(theta)", "italic(MSS)(theta)",
                               "italic(PIE)(theta)"))
figure3a                    <- ggplot() +
  geom_line(data = dplyr::filter(optimal, key == "f1"),
            mapping = aes(w, value, colour = key)) +
  ylab(expression(Optimal~~italic(f)[1])) +
  facet_grid(optimality~scenario, labeller = label_parsed) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  scale_colour_colorblind(breaks = c("f1", "rho"),
                          labels = list(bquote(italic(f)[1]),
                                        bquote(italic(t)[1]))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w)))
figure3b                    <- ggplot() +
  geom_line(data = dplyr::filter(optimal, key == "rho"),
            mapping = aes(w, value, colour = key)) +
  ylab(expression(Optimal~~italic(t)[1])) +
  facet_grid(optimality~scenario, labeller = label_parsed) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  scale_colour_colorblind(breaks = c("f1", "rho"),
                          labels = list(bquote(italic(f)[1]),
                                        bquote(italic(t)[1]))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w)))
figure3 <- figure3a + figure3b
ggsave("figure3.pdf", figure3, device = "pdf", width = 6.3, height = 6.3)

##### Wang-Tsiatis #############################################################

# Helper functions
covariance    <- function(sqrt_I) {
  CovZ            <- diag(length(sqrt_I))
  for (k1 in seq_len(length(sqrt_I) - 1) + 1L) {
    vec           <- 1:(k1 - 1)
    CovZ[k1, vec] <- CovZ[vec, k1] <- sqrt_I[vec]/sqrt_I[k1]
  }
  CovZ
}

power         <- function(tau, e, fu, sqrt_I, CovZ, J, seqs) {
  means <- tau*sqrt_I
  P     <- pnorm(e[1], means[1], lower.tail = FALSE)
  P2    <- pmvnorm(c(fu[1], e[2]), c(e[1], Inf), means[seqs[[2]]],
                   sigma = CovZ[seqs[[2]], seqs[[2]]])[1]
  if (!is.nan(P2)) {
    P   <- P + P2
  }
  if (J == 3) {
    P3  <- pmvnorm(c(fu[seqs[[2]]], e[3]), c(e[seqs[[2]]], Inf),
                   means[seqs[[3]]], sigma = CovZ[seqs[[3]], seqs[[3]]])[1]
    if (!is.nan(P3)) {
      P <- P + P3
    }
  }
  if (is.nan(P)) {
    0
  } else {
    P
  }
}

root_C        <- function(C, alpha, bound_factor, Delta, CovZ, J, Jm1, seqs,
                          numeric_J) {
  e <- C*bound_factor
  alpha - power(0, e, c(-e[seqs[[Jm1]]], C), numeric_J, CovZ, J, seqs)
}

root_n0J      <- function(n0J, delta, sigma2, ratio, CovZ, rho, J, seqs, e, fu,
                          beta) {
  n0 <- n0J*rho
  (1 - beta) - power(delta, e, fu, sqrt(n0/(sigma2[1] + sigma2[2]/ratio)), CovZ,
                     J, seqs)
}

opchar_int    <- function(tau, e, f, sqrt_I, CovZ, n, J, seqs) {
  len_tau          <- length(tau)
  means            <- matrix(tau, ncol = 1)%*%matrix(sqrt_I, 1)
  opchar           <- matrix(0, len_tau, 6 + 4*J)
  E                <- Fu <- numeric(J)
  for (t in 1:len_tau) {
    Fu[1]          <- pnorm(f[1], mean = means[t, 1])
    E[1]           <- pnorm(e[1], mean = means[t, 1], lower.tail = FALSE)
    if (J == 3) {
      for (j in 2:(J - 1)) {
        jm1        <- j - 1
        Fu[j]      <-
          pmvnorm(c(f[seqs[[jm1]]], -Inf), c(e[seqs[[jm1]]], f[j]),
                  means[t, seqs[[j]]], sigma = CovZ[seqs[[j]], seqs[[j]]],
                  algorithm = GenzBretz(abseps = 0.0000001,
                                        maxpts = 1000000))[1]
        E[j]       <-
          pmvnorm(c(f[seqs[[jm1]]], e[j]), c(e[seqs[[jm1]]], Inf),
                  means[t, seqs[[j]]], sigma = CovZ[seqs[[j]], seqs[[j]]],
                  algorithm = GenzBretz(abseps = 0.0000001,
                                        maxpts = 1000000))[1]
      }
    }
    Jm1            <- J - 1
    E[J]           <-
      pmvnorm(c(f[seqs[[Jm1]]], e[J]), c(e[seqs[[Jm1]]], Inf),
              means[t, ], sigma = CovZ,
              algorithm = GenzBretz(abseps = 0.0000001, maxpts = 1000000))[1]
    Fu[J]          <- 1 - sum(E) - sum(Fu[seqs[[Jm1]]])
    Fu[is.nan(Fu)] <- 0
    E[is.nan(E)]   <- 0
    cum_S          <- cumsum(S <- E + Fu)
    MSS            <- ifelse(any(cum_S == 0.5),
                             0.5*(n[which(cum_S == 0.5)] +
                                    n[which(cum_S == 0.5) + 1]),
                             n[which(cum_S > 0.5)[1]])
    if (sum(n^2*S) - sum(n*S)^2 < 0) {
      SDSS         <- 0
    } else {
      SDSS         <- sqrt(sum(n^2*S) - sum(n*S)^2)
    }
    opchar[t, ]    <- c(tau[t], sum(E), sum(n*S),
                        SDSS, MSS, E, Fu, S, cum_S,
                        n[J])
  }
  colnames(opchar) <- c("tau", "P(tau)", "ESS(tau)", "SDSS(tau)",
                        "MSS(tau)", paste(rep(c("E", "F", "S"), each = J),
                                          rep(seqs[[J]], 3), "(tau)", sep = ""),
                        paste("cum{S", seqs[[J]], "(tau)}", sep = ""), "max(n)")
  as_tibble(opchar)
}

opchar        <- function(pars, sigma2, ratio, alpha, beta, delta) {
  Delta        <- pars[1]
  rho          <- c(pars[-1], 1)
  J            <- length(rho)
  seqs         <- list()
  seqs[[J]]    <- 1:J
  for (j in seqs[[J]][-J]) {
    seqs[[j]]  <- 1:j
  }
  bound_factor <- rho^(Delta - 0.5)
  Jm1          <- J - 1
  CovZ         <- covariance(rho/sqrt(sigma2[1] + sigma2[2]/ratio))
  C            <- uniroot(root_C, c(1e-6, 1e6), alpha = alpha,
                          bound_factor = bound_factor, Delta = Delta,
                          CovZ = CovZ, J = J, Jm1 = Jm1, seqs = seqs,
                          numeric_J = numeric(J), tol = .Machine$double.eps^0.5,
                          maxiter = 10000)$root
  e            <- C*bound_factor
  fu           <- c(-e[seqs[[Jm1]]], C)
  n0J          <- uniroot(root_n0J, c(1e-6, 1e6), delta = delta,
                          sigma2 = sigma2, ratio = ratio, CovZ = CovZ,
                          rho = rho, J = J, seqs = seqs, e = e, fu = fu,
                          beta = beta, tol = .Machine$double.eps^0.5,
                          maxiter = 10000)$root
  n0           <- n0J*rho
  n            <- n0*(1 + ratio)
  sqrt_I       <- sqrt(n0/(sigma2[1] + sigma2[2]/ratio))
  opchar_H0    <- opchar_int(0, e, fu, sqrt_I, CovZ, n, J, seqs)
  opchar_H1    <- opchar_int(delta, e, fu, sqrt_I, CovZ, n, J, seqs)
  as.numeric(c(opchar_H0[3:5], opchar_H1[3:5], n[J],
               sum(opchar_H0[6:(4 + J)]), sum(opchar_H1[(6 + J):(4 + 2*J)])))
}

##### Example parameters
sigma2 <- c(1, 1)
delta  <- 0.3
ratio  <- 1
alpha  <- 0.025
beta   <- 0.1

### 2-stage

CovZ                       <- diag(1, 2, 2)

wrapper                    <- function(i) {
  c(designs[i, 1:2], opchar(designs[i, 1:2], sigma2, ratio, alpha, beta, delta))
}

# Combinations of Delta and rho to consider
designs                       <-
  as.matrix(expand.grid(Delta = seq(0, 0.5, 0.1),
                        rho   = seq(0.0100, 0.9999, 0.0001)))
designs                       <- rbind(designs, c(0.5, 1))
# Evaluate them using parallelisation
sfInit(parallel = TRUE, cpus = 7)
sfLibrary(mvtnorm)
sfLibrary(stats)
sfLibrary(tibble)
sfExport("designs", "sigma2", "delta", "ratio", "alpha", "beta", "power",
         "CovZ", "opchar", "root_f2", "root_n02", "covariance", "root_C",
         "root_n0J", "opchar_int")
output2                       <- sfLapply(1:nrow(designs), wrapper)
designs                       <- matrix(unlist(output2), nrow(designs), 11,
                                        byrow = TRUE)
designs[nrow(designs), 1]     <- NA_real_
designs[nrow(designs), 10:11] <- 0
# Easy way to add in how fixed-sample design works
designs                       <- tibble(Delta         = designs[, 1],
                                        rho           = designs[, 2],
                                        `ESS(0)`      = designs[, 3],
                                        `SDSS(0)`     = designs[, 4],
                                        `MSS(0)`      = designs[, 5],
                                        `ESS(delta)`  = designs[, 6],
                                        `SDSS(delta)` = designs[, 7],
                                        `MSS(delta)`  = designs[, 8],
                                        `max(n)`      = designs[, 9],
                                        `PIE(0)`      = designs[, 10],
                                        `PIE(delta)`  = designs[, 10])
# Search for what is optimal for each w under H1
optimal1_SDSS                 <- tibble(w          = rep(seq(0, 1, 0.0005), 1),
                                        rho        = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "SDSS")
optimal1_SDSS                 <- rbind(optimal1_SDSS, optimal1_SDSS,
                                       optimal1_SDSS, optimal1_SDSS,
                                       optimal1_SDSS, optimal1_SDSS)
optimal1_SDSS$Delta           <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_SDSS)/6)
optimal1_MSS                  <- tibble(w          = rep(seq(0, 1, 0.0005), 1),
                                        rho        = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "MSS")
optimal1_MSS                  <- rbind(optimal1_MSS, optimal1_MSS, optimal1_MSS,
                                       optimal1_MSS, optimal1_MSS, optimal1_MSS)
optimal1_MSS$Delta            <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_MSS)/6)
optimal1_PIE                  <- tibble(w          = rep(seq(0, 1, 0.0005), 1),
                                        rho        = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "PIE")
optimal1_PIE                  <- rbind(optimal1_PIE, optimal1_PIE, optimal1_PIE,
                                       optimal1_PIE, optimal1_PIE, optimal1_PIE)
optimal1_PIE$Delta            <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_PIE)/6)
for (i in 1:nrow(optimal1_SDSS)) {
  designs_sub                 <-
    dplyr::filter(designs, Delta %in% c(optimal1_SDSS$Delta[i], NA))
  score_1                     <- optimal1_SDSS$w[i]*designs_sub$`ESS(delta)` +
    (1 - optimal1_SDSS$w[i])*designs_sub$`SDSS(delta)`
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_SDSS$rho[i]        <- designs_sub$rho[which_min[length(which_min)]]
  score_1                     <- optimal1_MSS$w[i]*designs_sub$`ESS(delta)` +
    (1 - optimal1_MSS$w[i])*designs_sub$`MSS(delta)`
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_MSS$rho[i]         <- designs_sub$rho[which_min[length(which_min)]]
  score_1                     <-
    optimal1_PIE$w[i]*designs_sub$`ESS(delta)`/max(designs_sub$`ESS(delta)`) +
    (1 - optimal1_PIE$w[i])*designs_sub$`PIE(delta)`/
    max(designs_sub$`PIE(delta)`)
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_PIE$rho[i]         <- designs_sub$rho[which_min[length(which_min)]]
  if (i%%100 == 0) {
    print(i)
  }
}
optimal                       <- rbind(optimal1_SDSS, optimal1_MSS,
                                       optimal1_PIE)
optimal$optimality[which(optimal$optimality == "SDSS")] <-
  "italic(SDSS)(delta)"
optimal$optimality[which(optimal$optimality == "MSS")]  <-
  "italic(MSS)(delta)"
optimal$optimality[which(optimal$optimality == "PIE")]  <-
  "italic(PIE)(delta)"
optimal$Delta                 <- as.character(optimal$Delta)
optimal$optimality            <- factor(optimal$optimality,
                                        unique(optimal$optimality))
figure2                       <- ggplot() +
  geom_line(data = optimal, mapping = aes(w, rho, colour = Delta)) +
  facet_grid(optimality~., labeller = label_parsed) +
  ylab("Optimal value") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w))) +
  ylab(expression(italic(t)[1])) +
  scale_colour_colorblind() +
  labs(color = expression(Delta)) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figure2.pdf", figure2, device = "pdf", width = 6.3, height = 6.3)

### 3-stage

CovZ                          <- diag(1, 3, 3)

wrapper                       <- function(i) {
  c(designs[i, 1:3], opchar(designs[i, 1:3], sigma2, ratio, alpha, beta, delta))
}

designs                       <-
  as.matrix(expand.grid(Delta = seq(0, 0.5, 0.1),
                        rho1  = seq(0.001, 0.999, by = 0.002),
                        rho2  = seq(0.001, 0.999, by = 0.002)))
designs                       <- designs[which(designs[, 2] < designs[, 3]), ]
designs                       <- rbind(designs, c(0.5, 1, 1))
sfInit(parallel = TRUE, cpus = 7)
sfLibrary(mvtnorm)
sfLibrary(stats)
sfLibrary(tibble)
sfExport("designs", "sigma2", "delta", "ratio", "alpha", "beta", "CovZ",
         "opchar", "opchar_int", "root_C", "root_n0J", "covariance", "power")
output3                       <- sfLapply(1:nrow(designs), wrapper)
designs                       <- matrix(unlist(output3), nrow(designs), 12,
                                        byrow = TRUE)
designs[nrow(designs), 1]     <- NA_real_
designs[nrow(designs), 11:12] <- 0
designs                       <- tibble(Delta         = designs[, 1],
                                        rho1          = designs[, 2],
                                        rho2          = designs[, 3],
                                        `ESS(0)`      = designs[, 4],
                                        `SDSS(0)`     = designs[, 5],
                                        `MSS(0)`      = designs[, 6],
                                        `ESS(delta)`  = designs[, 7],
                                        `SDSS(delta)` = designs[, 8],
                                        `MSS(delta)`  = designs[, 9],
                                        `max(n)`      = designs[, 10],
                                        `PIE(0)`      = designs[, 10],
                                        `PIE(delta)`  = designs[, 11])
optimal1_SDSS                 <- tibble(w          = rep(seq(0, 1, 0.05), 1),
                                        rho1       = NA_real_,
                                        rho2       = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "SDSS")
optimal1_SDSS                 <- rbind(optimal1_SDSS, optimal1_SDSS,
                                       optimal1_SDSS, optimal1_SDSS,
                                       optimal1_SDSS,optimal1_SDSS)
optimal1_SDSS$Delta           <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_SDSS)/6)
optimal1_MSS                  <- tibble(w          = rep(seq(0, 1, 0.05), 1),
                                        rho1       = NA_real_,
                                        rho2       = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "MSS")
optimal1_MSS                  <- rbind(optimal1_MSS, optimal1_MSS, optimal1_MSS,
                                       optimal1_MSS, optimal1_MSS, optimal1_MSS)
optimal1_MSS$Delta            <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_MSS)/6)
optimal1_PIE                  <- tibble(w          = rep(seq(0, 1, 0.05), 1),
                                        rho1       = NA_real_,
                                        rho2       = NA_real_,
                                        scenario   = "Alternative",
                                        optimality = "PIE")
optimal1_PIE                  <- rbind(optimal1_PIE, optimal1_PIE, optimal1_PIE,
                                       optimal1_PIE, optimal1_PIE, optimal1_PIE)
optimal1_PIE$Delta            <- rep(seq(0, 0.5, 0.1),
                                     each = nrow(optimal1_PIE)/6)
for (i in 1:nrow(optimal1_SDSS)) {
  designs_sub                 <-
    dplyr::filter(designs, Delta %in% c(optimal1_SDSS$Delta[i], NA))
  score_1                     <- optimal1_SDSS$w[i]*designs_sub$`ESS(delta)` +
    (1 - optimal1_SDSS$w[i])*designs_sub$`SDSS(delta)`
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_SDSS$rho1[i]       <- designs_sub$rho1[which_min[length(which_min)]]
  optimal1_SDSS$rho2[i]       <- designs_sub$rho2[which_min[length(which_min)]]
  score_1                     <- optimal1_MSS$w[i]*designs_sub$`ESS(delta)` +
    (1 - optimal1_MSS$w[i])*designs_sub$`MSS(delta)`
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_MSS$rho1[i]        <- designs_sub$rho1[which_min[length(which_min)]]
  optimal1_MSS$rho2[i]        <- designs_sub$rho2[which_min[length(which_min)]]
  score_1                     <- optimal1_PIE$w[i]*designs_sub$`ESS(delta)`/
    max(designs_sub$`ESS(delta)`) +
    (1 - optimal1_PIE$w[i])*designs_sub$`PIE(delta)`/
    max(designs_sub$`PIE(delta)`)
  which_min                   <- which(score_1 <= min(score_1) + 1e-6)
  optimal1_PIE$rho1[i]        <- designs_sub$rho1[which_min[length(which_min)]]
  optimal1_PIE$rho2[i]        <- designs_sub$rho2[which_min[length(which_min)]]
  if (i%%100 == 0) {
    print(i)
  }
}
optimal                       <- rbind(optimal1_SDSS, optimal1_MSS,
                                       optimal1_PIE)
optimal$optimality[which(optimal$optimality == "SDSS")] <-
  paste0("italic(SDSS)(delta)")
optimal$optimality[which(optimal$optimality == "MSS")]  <-
  paste0("italic(MSS)(delta)")
optimal$optimality[which(optimal$optimality == "PIE")]  <-
  paste0("italic(PIE)(delta)")
optimal                       <- tidyr::gather(optimal, key = "key",
                                               value = "value", rho1:rho2)
optimal$Delta                 <- as.character(optimal$Delta)
optimal$key[optimal$key == "rho1"] <- "italic(t)[1]"
optimal$key[optimal$key == "rho2"] <- "italic(t)[2]"
optimal$optimality            <- factor(optimal$optimality,
                                        unique(optimal$optimality))
sfigure1                      <- ggplot() +
  geom_line(data = optimal, mapping = aes(w, value, colour = Delta)) +
  geom_point(data = optimal, mapping = aes(w, value, colour = Delta),
             size = 0.5) +
  facet_grid(optimality~key, labeller = label_parsed) +
  ylab("Optimal value") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w))) +
  scale_colour_colorblind() +
  labs(color = expression(Delta)) +
  theme_bw() +
  theme(legend.position = "bottom")
sfigure1
ggsave("sfigure1.pdf", sfigure1, device = "pdf", width = 6.3, height = 6.3)

##### Figure 1 #################################################################

opchar2       <- function(pars, sigma2, ratio, alpha, beta, delta, tau) {
  Delta        <- pars[1]
  rho          <- c(pars[-1], 1)
  J            <- length(rho)
  seqs         <- list()
  for (j in 1:J) {
    seqs[[j]]  <- 1:j
  }
  bound_factor <- rho^(Delta - 0.5)
  Jm1          <- J - 1
  CovZ         <- covariance(1/sqrt(sigma2[1]/rho + sigma2[2]/(ratio*rho)))
  C            <- uniroot(root_C, c(1e-6, 1e6), alpha = alpha,
                          bound_factor = bound_factor, Delta = Delta,
                          CovZ = CovZ, J = J, Jm1 = Jm1, seqs = seqs,
                          numeric_J = numeric(J), tol = .Machine$double.eps^0.5,
                          maxiter = 10000)$root
  e            <- C*bound_factor
  fu           <- c(-e[seqs[[Jm1]]], C)
  n0J          <- uniroot(root_n0J, c(1e-6, 1e6), delta = delta,
                          sigma2 = sigma2, ratio = ratio, CovZ = CovZ,
                          rho = rho, J = J, seqs = seqs, e = e, fu = fu,
                          beta = beta, tol = .Machine$double.eps^0.5,
                          maxiter = 10000)$root
  n0           <- n0J*rho
  n1           <- n0*ratio
  n            <- n0 + n1
  sqrt_I       <- 1/sqrt(sigma2[1]/n0 + sigma2[2]/n1)
  ess          <- sdss <- mss <- pie <- numeric(length(tau))
  for (i in 1:length(tau)) {
    opchars    <- opchar_int(tau[i], e, fu, sqrt_I, CovZ, n, J, seqs)
    ess[i]     <- opchars[3]
    sdss[i]    <- opchars[4]
    mss[i]     <- opchars[5]
    if (tau[i] <= 0) {
      pie[i]   <- sum(opchars[6:(4 + J)])
    } else {
      pie[i]   <- sum(opchars[(6 + J):(4 + 2*J)])
    }
  }
  cbind(ess, sdss, mss, pie)
}

designs       <- as.matrix(expand.grid(Delta  = seq(0, 0.5, by = 0.1),
                                       rho    = c(0.25, 0.5, 0.75)))
results       <- NULL
tau           <- seq(from = -2*delta, to = 3*delta, length.out = 500)
for (i in 1:nrow(designs)) {
  results     <- rbind(results,
                       cbind(matrix(designs[i, ], length(tau), 2, byrow = TRUE),
                             tau,
                             opchar2(designs[i, ], sigma2, ratio, alpha, beta,
                                     delta, tau)))
  message(i)
}
results       <- tibble::tibble(Delta = unlist(results[, 1]),
                                rho   = unlist(results[, 2]),
                                tau   = unlist(results[, 3]),
                                ess   = unlist(results[, 4]),
                                sdss  = unlist(results[, 5]),
                                mss   = unlist(results[, 6]),
                                pie   = unlist(results[, 7]))
results       <- tidyr::gather(results, key = "key", value = "value", ess:pie)
results$key[results$key == "ess"]  <- "italic(ESS)(tau)"
results$key[results$key == "sdss"] <- "italic(SDSS)(tau)"
results$key[results$key == "mss"]  <- "italic(MSS)(tau)"
results$key[results$key == "pie"]  <- "italic(PIE)(tau)"
results$key <- factor(results$key, levels = c("italic(ESS)(tau)",
                                              "italic(SDSS)(tau)",
                                              "italic(MSS)(tau)",
                                              "italic(PIE)(tau)"))
results$Delta <- as.character(results$Delta)
results$rho   <- as.character(paste0("italic(t)[1] == ", results$rho))
figure1       <- ggplot(results, aes(tau, value, colour = Delta)) +
  geom_line() +
  facet_grid(key~rho, labeller = label_parsed, scales = "free_y") +
  scale_colour_colorblind() +
  labs(color = expression(Delta)) +
  theme_bw() +
  xlab(expression(tau)) +
  ylab("Value") +
  theme(legend.position = "bottom")
ggsave("figure1.pdf", figure1, device = "pdf", width = 6.3, height = 6.3)
