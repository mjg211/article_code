##### Futility stopping only ###################################################

# Helper functions
root_f2    <- function(f2, f1, CovZ, alpha) {
  pmvnorm(c(f1, f2), sigma = CovZ)[1] - alpha
}

root_n02   <- function(n02, f1, rho, sigma2, ratio, alpha, power, delta, CovZ) {
  n01        <- rho*n02
  sqrtI      <- sqrt(1/c(sigma2[1]/n01 + sigma2[2]/(ratio*n01),
                         sigma2[1]/n02 + sigma2[2]/(ratio*n02)))
  CovZ[1, 2] <- CovZ[2, 1] <- sqrtI[1]/sqrtI[2]
  f2         <- uniroot(root_f2, c(-1e3, 1e3), f1 = f1, CovZ = CovZ,
                        alpha = alpha)$root
  pmvnorm(c(f1, f2), c(Inf, Inf), delta*sqrtI, sigma = CovZ)[1] - power
}

opchar     <- function(pars, sigma2, ratio, alpha, power, delta, CovZ) {
  f1   <- pars[1]
  rho  <- pars[2]
  n02  <- uniroot(root_n02, c(1, 1e6), f1 = f1, rho = rho,
                  sigma2 = sigma2, ratio = ratio, alpha = alpha,
                  power = power, delta = delta, CovZ = CovZ)$root
  n01  <- rho*n02
  n11  <- ratio*n01
  n    <- c(n01 + n11, n02 + ratio*n02)
  S0   <- pnorm(f1)
  S0   <- c(S0, 1 - S0)
  S1   <- pnorm(f1 - delta/sqrt(sigma2[1]/n01 + sigma2[2]/n11))
  S1   <- c(S1, 1 - S1)
  ess0 <- sum(n*S0)
  ess1 <- sum(n*S1)
  if (S0[1] < 0.5) {
    med0 <- n[2]
  } else if (S0[1] == 0.5) {
    med0 <- n[1] + 0.5*n[2]
  } else {
    med0 <- n[1]
  }
  if (S1[1] < 0.5) {
    med1 <- n[2]
  } else if (S1[1] == 0.5) {
    med1 <- n[1] + 0.5*n[2]
  } else {
    med1 <- n[1]
  }
  c(ess0, ess1, med0, med1, sqrt(sum(n^2*S0) - ess0^2), sqrt(sum(n^2*S1) - ess1^2), n[2])
}

wrapper    <- function(i) {
  c(designs[i, 1:2], opchar(designs[i, 1:2], sigma2, ratio, alpha, power,
                            delta, CovZ))
}

# Example parameters
sigma2             <- c(1, 1)
delta              <- 1
ratio              <- 1
alpha              <- 0.025
power              <- 0.9
CovZ               <- diag(1, 2, 2)
# Combinations of f1 and rho to consider
designs            <- as.matrix(expand.grid(f1  = seq(-0.25, 1.95, 0.005),
                                            rho = seq(0.01, 0.99, 0.005)))
designs            <- cbind(designs, matrix(0, nrow(designs), 7))
# Evaluate them using parallelisation
sfInit(parallel = T, cpus = 7)
sfLibrary(stats)
sfLibrary(mvtnorm)
sfExport("designs", "sigma2", "delta", "ratio", "alpha", "power", "CovZ",
         "opchar", "root_f2", "root_n02")
output             <- sfLapply(1:nrow(designs), wrapper)
designs            <- matrix(unlist(output), nrow(designs), 9, byrow = T)
# Easy way to add in how fixed-sample design works
min_n              <- min(designs[, 9]) - 1e-6
designs            <- rbind(designs,
                            c(NA, 1, min_n, min_n, min_n, min_n, 0, 0, min_n))
designs            <- tibble::tibble(f1            = designs[, 1],
                                     rho           = designs[, 2],
                                     `ESS(0)`      = designs[, 3],
                                     `ESS(delta)`  = designs[, 4],
                                     `MSS(0)`      = designs[, 5],
                                     `MSS(delta)`  = designs[, 6],
                                     `SDSS(0)`     = designs[, 7],
                                     `SDSS(delta)` = designs[, 8],
                                     `max(n)`      = designs[, 9])
# Search for what is optimal for each w under H0 and H1
optimal0           <- tibble::tibble(w        = rep(seq(0, 1, 0.0005), 1),
                                     f1       = NA_real_,
                                     rho      = NA_real_,
                                     scenario = "Null")
optimal1           <- tibble::tibble(w        = rep(seq(0, 1, 0.0005), 1),
                                     f1       = NA_real_,
                                     rho      = NA_real_,
                                     scenario = "Alternative")
for (i in 1:nrow(optimal0)) {
  score_0          <- optimal0$w[i]*designs$`ESS(0)` +
    (1 - optimal0$w[i])*designs$`SDSS(0)`
  optimal_0        <- which.min(score_0)
  optimal0$f1[i]   <- designs$f1[optimal_0]
  optimal0$rho[i]  <- designs$rho[optimal_0]
  score_d          <- optimal1$w[i]*designs$`ESS(delta)` +
    (1 - optimal1$w[i])*designs$`SDSS(delta)`
  optimal_d        <- which.min(score_d)
  optimal1$f1[i]   <- designs$f1[optimal_d]
  optimal1$rho[i]  <- designs$rho[optimal_d]
  print(i)
}
optimal            <- rbind(optimal0, optimal1)
optimal            <- tidyr::gather(optimal, key = "key", value = "value",
                                    f1:rho)
optimal$scenario[which(optimal$scenario == "Null")]        <-
  paste0("theta == 0")
optimal$scenario[which(optimal$scenario == "Alternative")] <-
  paste0("theta == delta")
# Plot the findings
figure3            <- ggplot() +
  geom_line(data = optimal, mapping = aes(w, value, colour = key)) +
  ylab("Optimal value(s)") +
  facet_grid(~scenario, labeller = label_parsed) +
  ggplot2::theme_grey(base_size = 11) + ggplot2::theme(
    panel.background =
      ggplot2::element_rect(fill = "white",
                            colour = NA),
    panel.border = ggplot2::element_rect(fill = NA,
                                         colour = "grey70",
                                         size = 0.5),
    panel.grid.major =
      ggplot2::element_line(colour = "grey87", size = 0.25),
    panel.grid.minor =
      ggplot2::element_line(colour = "grey87", size = 0.125),
    axis.ticks = ggplot2::element_line(colour = "grey70",
                                       size = 0.25),
    legend.key = ggplot2::element_rect(fill = "white",
                                       colour = NA),
    strip.background = ggplot2::element_rect(fill = "grey70",
                                             colour = NA),
    strip.text =
      ggplot2::element_text(colour = "white",
                            size = ggplot2::rel(0.8)),
    legend.position = "bottom",
    plot.margin = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3), "cm"),
    complete = TRUE) +
  scale_colour_colorblind(breaks = c("f1","rho"),
                          labels = list(bquote(italic(f)[1]),
                                        bquote(italic(t)[1]))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w)))
ggsave("figure3.pdf", figure3, device = "pdf", width = 6, height = 4)
# Note: Get Supplementary Figure 3 by replacing SDSS with MSS in search for
# optimal designs

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
  P     <- stats::pnorm(e[1], mean = means[1], lower.tail = F)
  for (j in seqs[[J]][-1]) {
    jm1 <- j - 1
    P   <- P + pmvnorm(c(fu[seqs[[jm1]]], e[j]),
                       c(e[seqs[[jm1]]], Inf),
                       means[seqs[[j]]],
                       sigma = CovZ[seqs[[j]], seqs[[j]]])[1]
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
  n1 <- n0*ratio
  (1 - beta) - power(delta, e, fu, 1/sqrt(sigma2[1]/n0 + sigma2[2]/n1), CovZ, J,
                     seqs)
}

opchar_int    <- function(tau, e, f, sqrt_I, CovZ, n) {
  J                <- length(e)
  len_tau          <- length(tau)
  means            <- matrix(tau, ncol = 1)%*%matrix(sqrt_I, 1)
  opchar           <- matrix(0, len_tau, 6 + 4*J)
  E                <- Fu <- numeric(J)
  for (t in 1:len_tau) {
    Fu[1]          <- stats::pnorm(f[1], mean = means[t, 1])
    E[1]           <- stats::pnorm(e[1], mean = means[t, 1], lower.tail = F)
    if (J > 2) {
      for (j in 2:(J - 1)) {
        Fu[j]      <-
          mvtnorm::pmvnorm(c(f[1:(j - 1)], -Inf), c(e[1:(j - 1)], f[j]),
                           means[t, 1:j], sigma = CovZ[1:j, 1:j],
                           algorithm = GenzBretz(abseps = 0.0000001,
                                                 maxpts = 1000000))[1]
        E[j]       <-
          mvtnorm::pmvnorm(c(f[1:(j - 1)], e[j]), c(e[1:(j - 1)], Inf),
                           means[t, 1:j], sigma = CovZ[1:j, 1:j],
                           algorithm = GenzBretz(abseps = 0.0000001,
                                                 maxpts = 1000000))[1]
      }
    }
    E[J]           <-
      mvtnorm::pmvnorm(c(f[1:(J - 1)], e[J]), c(e[1:(J - 1)], Inf),
                       means[t, ], sigma = CovZ,
                       algorithm = GenzBretz(abseps = 0.0000001,
                                             maxpts = 1000000))[1]
    Fu[J]          <- 1 - sum(E) - sum(Fu[1:(J - 1)])
    cum_S          <- cumsum(S <- E + Fu)
    MSS            <- ifelse(any(cum_S == 0.5),
                             0.5*(n[which(cum_S == 0.5)] +
                                    n[which(cum_S == 0.5) + 1]),
                             n[which(cum_S > 0.5)[1]])
    opchar[t, ]    <- c(tau[t], sum(E), sum(n*S),
                        sqrt(sum(n^2*S) - sum(n*S)^2), MSS, E, Fu, S, cum_S,
                        n[J])
  }
  colnames(opchar) <- c("tau", "P(tau)", "ESS(tau)", "SDSS(tau)",
                        "MSS(tau)", paste(rep(c("E", "F", "S"), each = J),
                                          rep(1:J, 3), "(tau)", sep = ""),
                        paste("cum{S", 1:J, "(tau)}", sep = ""), "max(n)")
  tibble::as_tibble(opchar)
}

opchar        <- function(pars, sigma2, ratio, alpha, beta, delta) {
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
  opchar_H0    <- opchar_int(0, e, fu, sqrt_I, CovZ, n)
  opchar_H1    <- opchar_int(delta, e, fu, sqrt_I, CovZ, n)
  maxESS       <- 0 # Removed for speed
    #-optim(par = 0.5*delta, fn = minus_ess, e = e, f = fu, sqrt_I = sqrt_I,
    #       lower = -delta, upper = delta, CovZ = CovZ, n = n,
    #       method = "Brent")$value
  as.numeric(c(opchar_H0[3:5], opchar_H1[3:5], maxESS, n[J]))
}

minus_ess     <- function(tau, e, f, sqrt_I, CovZ, n) {
  J        <- length(e)
  means    <- tau*sqrt_I
  S        <- c(stats::pnorm(f[1], mean = means[1]) +
                  stats::pnorm(e[1], mean = means[1], lower.tail = F),
                numeric(J - 1))
  if (J > 2) {
    for (j in 2:(J - 1)) {
      S[j] <-
        mvtnorm::pmvnorm(c(f[1:(j - 1)], -Inf), c(e[1:(j - 1)], f[j]),
                         means[1:j], sigma = CovZ[1:j, 1:j])[1] +
        mvtnorm::pmvnorm(c(f[1:(j - 1)], e[j]), c(e[1:(j - 1)], Inf),
                         means[1:j], sigma = CovZ[1:j, 1:j])[1]
    }
  }
  S[J]     <- 1 - sum(S[1:(J - 1)])
  -sum(n*S)
}

wrapper       <- function(i) {
  c(designs[i, 1:2], opchar(designs[i, 1:2], sigma2, ratio, alpha, beta, delta))
}

# Example parameters
sigma2             <- c(1, 1)
delta              <- 1
ratio              <- 1
alpha              <- 0.025
beta               <- 0.1
delta              <- 1
CovZ               <- diag(1, 2, 2)
# Combinations of Delta and rho to consider
designs            <- as.matrix(expand.grid(Delta  = seq(0, 0.5, by = 0.1),
                                            rho    = seq(0.0001, 0.9999,
                                                         0.0001)))
# Evaluate them using parallelisation
sfInit(parallel = T, cpus = 7)
sfLibrary(stats)
sfLibrary(mvtnorm)
sfExport("designs", "sigma2", "delta", "ratio", "alpha", "beta", "power",
         "CovZ", "opchar", "root_f2", "root_n02", "covariance", "root_C",
         "root_n0J", "opchar_int")
output2            <- sfLapply(1:nrow(designs), wrapper)
designs            <- matrix(unlist(output2), nrow(designs), 10, byrow = T)
# Easy way to add in how fixed-sample design works
min_n              <- min(designs[, 10]) - 1e-6
designs            <- rbind(designs,
                            c(NA, 1, min_n, 0, min_n, min_n, 0, min_n, min_n,
                              min_n))
designs            <- tibble::tibble(Delta         = designs[, 1],
                                     rho           = designs[, 2],
                                     `ESS(0)`      = designs[, 3],
                                     `SDSS(0)`     = designs[, 4],
                                     `MSS(0)`      = designs[, 5],
                                     `ESS(delta)`  = designs[, 6],
                                     `SDSS(delta)` = designs[, 7],
                                     `MSS(delta)`  = designs[, 8],
                                     `maxESS`      = designs[, 9],
                                     `max(n)`      = designs[, 10])
# Search for what is optimal for each w under H1
optimal1           <- tibble::tibble(w        = rep(seq(0, 1, 0.0005), 1),
                                     rho      = NA_real_,
                                     scenario = "Alternative")
optimal1           <- rbind(optimal1, optimal1, optimal1, optimal1, optimal1,
                            optimal1)
optimal1$Delta     <- rep(seq(0, 0.5, by = 0.1), each = nrow(optimal1)/6)
for (i in 1:nrow(optimal1)) {
  designs_sub      <- dplyr::filter(designs,
                                    Delta %in% c(optimal1$Delta[i], NA) &
                                      rho > 0.001)
  score_d          <- optimal1$w[i]*designs_sub$`ESS(delta)` +
    (1 - optimal1$w[i])*designs_sub$`SDSS(delta)`
  optimal_d        <- which.min(score_d)
  optimal1$rho[i]  <- designs_sub$rho[optimal_d]
  if (i%%100 == 0) {
    print(i)
  }
}
optimal            <- optimal1
optimal$Delta      <- as.character(optimal$Delta)
optimal$rho[which(optimal$rho < 0.01)] <- 1
figure2            <- ggplot() +
  geom_line(data = optimal, mapping = aes(w, rho, colour = Delta)) +
  ylab("Optimal value") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab(expression(italic(w))) +
  ylab(expression(italic(t)[1])) +
  scale_colour_colorblind() +
  labs(color = expression(Delta)) +
  ggplot2::theme_grey(base_size = 11) + ggplot2::theme(
    panel.background =
      ggplot2::element_rect(fill = "white",
                            colour = NA),
    panel.border = ggplot2::element_rect(fill = NA,
                                         colour = "grey70",
                                         size = 0.5),
    panel.grid.major =
      ggplot2::element_line(colour = "grey87", size = 0.25),
    panel.grid.minor =
      ggplot2::element_line(colour = "grey87", size = 0.125),
    axis.ticks = ggplot2::element_line(colour = "grey70",
                                       size = 0.25),
    legend.key = ggplot2::element_rect(fill = "white",
                                       colour = NA),
    strip.background = ggplot2::element_rect(fill = "grey70",
                                             colour = NA),
    strip.text =
      ggplot2::element_text(colour = "white",
                            size = ggplot2::rel(0.8)),
    legend.position = "bottom",
    plot.margin = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3), "cm"),
    complete = TRUE)
ggsave("figure2.pdf", figure2, device = "pdf", width = 6, height = 5)
