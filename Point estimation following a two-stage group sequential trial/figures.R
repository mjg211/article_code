##### Load required packages ###################################################

#install.packages("mvtnorm")
#install.packages("patchwork")
#install.packages("tibble")
#install.packages("tidyverse")
library(mvtnorm)
library(patchwork)
library(snowfall)
library(tidyverse)

##### Load required local functions ############################################

bias_mle            <- function(theta, I1, I2, l, u) {
  ((I2 - I1)/(I2*sqrt(I1)))*(dnorm(u, theta*sqrt(I1), 1) -
                               dnorm(l, theta*sqrt(I1), 1))
}

cmue_optim          <- function(theta, k, z, I1, I2, l, u) {
  if (k == 1) {
    S1       <- pnorm(l, theta*sqrt(I1), 1) + pnorm(u, theta*sqrt(I1), 1,
                                                    lower.tail = FALSE)
    if (z >= u) {
      num    <- pnorm(l, theta*sqrt(I1)) + pnorm(z, theta*sqrt(I1)) -
        pnorm(u, theta*sqrt(I1))
    } else {
      num    <- pnorm(z, theta*sqrt(I1))
    }
    int      <- num/S1
  } else {
    if (is.infinite(u)) {
      log_S2 <- pnorm(theta*sqrt(I1) - l, log.p = TRUE)
    } else {
      log_S2 <- log(pnorm(u, theta*sqrt(I1)) - pnorm(l, theta*sqrt(I1)))
    }
    z        <- seq(qnorm(0.0001, z), z, length.out = 1000)
    cf2z     <- exp(f2(z, theta, I1, I2, l, u, log = TRUE) - log_S2)
    int      <- int.simpson2(z, cf2z)
  }
  (int - 0.5)^2
}

est_gs              <- function(theta, I1, I2, l, u, theta0) {
  L1                <- pnorm(l, theta*sqrt(I1), 1)
  U1                <- pnorm(u, theta*sqrt(I1), 1, lower.tail = FALSE)
  S1                <- L1 + U1
  E11               <- E12 <- E21 <- E22 <- matrix(0, length(theta), 9)
  z1                <- seq(qnorm(1e-8, theta[1]*sqrt(I1)),
                           qnorm(1 - 1e-8, theta[length(theta)]*sqrt(I1)),
                           length.out = 10000)
  z2                <- seq(qnorm(1e-8, theta[1]*sqrt(I2)),
                           qnorm(1 - 1e-8, theta[length(theta)]*sqrt(I2)),
                           length.out = 10000)
  est_mle1          <- z1/sqrt(I1)
  c                 <- (est_mle1 - theta0)^2/((est_mle1 - theta0)^2 + I1^-1)
  theta_star        <- c*est_mle1 + (1 - c)*theta0
  est_mle2          <- z2/sqrt(I2)
  est_mae11         <- est_mae12  <- est_cmle1  <- est_cmue1 <- est_cwmae1 <-
    est_mae21  <- est_mae22  <- est_cmle2 <- est_cmue2  <-
    est_mue2   <- est_umvue2 <- numeric(10000)
  for (i in 1:10000) {
    if (z1[i] < l | z1[i] > u) {
      est_mae11[i]  <- uniroot(mae_root, c(-1e10, 1e10), mle = est_mle1[i],
                               I1 = I1, I2 = I2, l = l, u = u)$root
      est_mae12[i]  <- est_mle1[i] - bias_mle(est_mle1[i], I1, I2, l, u)
      est_cmle1[i]  <- optim(est_mle1[i], min_clog_likelihood, k = 1,
                             z = z1[i], I1 = I1, I2 = I2, l = l, u = u,
                             method = "BFGS")$par
      est_cmue1[i]  <- optim(est_mle1[i], cmue_optim, k = 1, z = z1[i],
                             I1 = I1, I2 = I2, l = l, u = u,
                             method = "BFGS")$par
      est_cwmae1[i] <- est_mle1[i] - bias_mle(theta_star[i], I1, I2, l, u)
    }
    est_mae21[i]    <- uniroot(mae_root, c(-1e10, 1e10), mle = est_mle2[i],
                               I1 = I1, I2 = I2, l = l, u = u)$root
    est_mae22[i]    <- est_mle2[i] - bias_mle(est_mle2[i], I1, I2, l, u)
    est_cmle2[i]    <- try(optim(est_mle2[i], min_clog_likelihood, k = 2,
                                 z = z2[i], I1 = I1, I2 = I2, l = l, u = u,
                                 method = "BFGS")$par, silent = TRUE)
    est_cmue2[i]    <- try(optim(est_mle2[i], cmue_optim, k = 2, z = z2[i],
                                 I1 = I1, I2 = I2, l = l, u = u,
                                 method = "BFGS")$par, silent = TRUE)
    est_umvue2[i]   <- umvue_2(z2[i], I1, I2, l, u)
    est_mue2[i]     <- uniroot(mue_2_root, c(-1e10, 1e10), z = z2[i],
                               I1 = I1, I2 = I2, l = l, u = u)$root
    if (i%%100 == 0) {
      message(i, "/10000 estimates computed")
    }
  }
  est_cmle2         <- as.numeric(est_cmle2)
  est_cmue2         <- as.numeric(est_cmue2)
  failed_cmle2      <- which(is.na(est_cmle2))
  failed_cmue2      <- which(is.na(est_cmue2))
  if (length(failed_cmle2) > 0) {
    for (i in 1:length(failed_cmle2)) {
      if (failed_cmle2[i] > 1 & failed_cmle2[i] < 10000) {
        est_cmle2[failed_cmle2[i]] <- 0.5*(est_cmle2[failed_cmle2[i] - 1] +
                                             est_cmle2[failed_cmle2[i] + 1])
      } else if (failed_cmle2[i] == 1) {
        est_cmle2[failed_cmle2[i]] <- est_cmle2[failed_cmle2[i] + 1]
      } else if (failed_cmle2[i] == 10000) {
        est_cmle2[failed_cmle2[i]] <- est_cmle2[failed_cmle2[i] - 1]
      }
    }
  }
  if (length(failed_cmue2) > 0) {
    for (i in 1:length(failed_cmue2)) {
      if (failed_cmue2[i] > 1 & failed_cmue2[i] < 10000) {
        est_cmue2[failed_cmue2[i]] <- 0.5*(est_cmue2[failed_cmue2[i] - 1] +
                                             est_cmue2[failed_cmue2[i] + 1])
      } else if (failed_cmue2[i] == 1) {
        est_cmue2[failed_cmue2[i]] <- est_cmue2[failed_cmue2[i] + 1]
      } else if (failed_cmue2[i] == 10000) {
        est_cmue2[failed_cmue2[i]] <- est_cmue2[failed_cmue2[i] - 1]
      }
    }
  }
  est_cmle2[is.na(est_cmle2)]      <- 0
  est_cmue2[is.na(est_cmue2)]      <- 0
  est_cumvue2       <- (z2*sqrt(I2) - I1*est_umvue2)/(I2 - I1)
  est11             <- rbind(est_mae11, est_mae12, est_cmle1, est_cmue1,
                             est_mle1, est_cwmae1)
  est12             <- est11*est11
  est21             <- rbind(est_mae21, est_mae22, est_cmle2, est_cmue2,
                             est_cumvue2, est_mle2, est_mue2, est_umvue2)
  est22             <- est21*est21
  for (i in 1:length(theta)) {
    f1z1            <- f1(z1, theta[i], I1, l, u)
    f2z2            <- f2(z2, theta[i], I1, I2, l, u)
    integrand11     <- matrix(f1z1, 6, 10000, byrow = TRUE)*est11
    E11[i, ]        <- c(int.simpson2(z1, integrand11[1, ]),
                         int.simpson2(z1, integrand11[2, ]),
                         int.simpson2(z1, integrand11[3, ]),
                         int.simpson2(z1, integrand11[4, ]),
                         rep(int.simpson2(z1, integrand11[5, ]), 4),
                         int.simpson2(z1, integrand11[6, ]))/S1[i]
    integrand12     <- matrix(f1z1, 6, 10000, byrow = TRUE)*est12
    E12[i, ]        <- c(int.simpson2(z1, integrand12[1, ]),
                         int.simpson2(z1, integrand12[2, ]),
                         int.simpson2(z1, integrand12[3, ]),
                         int.simpson2(z1, integrand12[4, ]),
                         rep(int.simpson2(z1, integrand12[5, ]), 4),
                         int.simpson2(z1, integrand12[6, ]))/S1[i]
    integrand21     <- matrix(f2z2, 8, 10000, byrow = TRUE)*est21
    temp            <- int.simpson2(z2, integrand21[6, ])
    E21[i, ]        <- c(int.simpson2(z2, integrand21[1, ]),
                         int.simpson2(z2, integrand21[2, ]),
                         int.simpson2(z2, integrand21[3, ]),
                         int.simpson2(z2, integrand21[4, ]),
                         int.simpson2(z2, integrand21[5, ]),
                         temp,
                         int.simpson2(z2, integrand21[7, ]),
                         int.simpson2(z2, integrand21[8, ]),
                         temp)/(1 - S1[i])
    integrand22     <- matrix(f2z2, 8, 10000, byrow = TRUE)*est22
    temp            <- int.simpson2(z2, integrand22[6, ])
    E22[i, ]        <- c(int.simpson2(z2, integrand22[1, ]),
                         int.simpson2(z2, integrand22[2, ]),
                         int.simpson2(z2, integrand22[3, ]),
                         int.simpson2(z2, integrand22[4, ]),
                         int.simpson2(z2, integrand22[5, ]),
                         temp,
                         int.simpson2(z2, integrand22[7, ]),
                         int.simpson2(z2, integrand22[8, ]),
                         temp)/(1 - S1[i])
  }
  list(E11 = E11, E12 = E12, E21 = E21, E22 = E22, S1 = S1,
       failed = length(failed_cmle2) + length(failed_cmue2))
}

f1                  <- function(z, theta, I1, l, u, log = FALSE) {
  f                       <- dnorm(z, theta*sqrt(I1), 1, log = log)
  f[which(z > l & z < u)] <- 0
  f
}

f2                  <- function(z, theta, I1, I2, l, u, log = FALSE) {
  if (!log) {
    exp(-0.5*(z - theta*sqrt(I2))^2)*
      (pnorm(u, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)) -
         pnorm(l, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)))/sqrt(2*pi)
  } else {
    if (is.infinite(u)) {
      -0.5*(z - theta*sqrt(I2))^2 +
        pnorm((z*sqrt(I1/I2) - l)/sqrt((I2 - I1)/I2), log.p = TRUE) -
        log(sqrt(2*pi))
    } else {
      -0.5*(z - theta*sqrt(I2))^2 +
        log((pnorm(u, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)) -
               pnorm(l, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)))) - log(sqrt(2*pi))
    }
  }
}

mae_root            <- function(theta, mle, I1, I2, l, u) {
  mle - theta - bias_mle(theta, I1, I2, l, u)
}

min_clog_likelihood <- function(theta, k, z, I1, I2, l, u) {
  if (k == 1) {
    if (is.infinite(u)) {
      -f1(z, theta, I1, l, u, log = TRUE) +
        pnorm(l, theta*sqrt(I1), 1, log.p = TRUE)
    } else {
      -f1(z, theta, I1, l, u, log = TRUE) +
        log(pnorm(u, theta*sqrt(I1), 1, lower.tail = FALSE) +
              pnorm(l, theta*sqrt(I1), 1))
    }
  } else {
    if (is.infinite(u)) {
      -f2(z, theta, I1, I2, l, u, log = TRUE) +
        pnorm(theta*sqrt(I1) - l, log.p = TRUE)
    } else {
      -f2(z, theta, I1, I2, l, u, log = TRUE) +
        log(pnorm(u, theta*sqrt(I1)) - pnorm(l, theta*sqrt(I1)))
    }
  }
}

mue_2_root          <- function(theta, z, I1, I2, l, u) {
  1 - pnorm(u, theta*sqrt(I1), 1) +
    pmvnorm(lower = c(l, z), upper = c(u, Inf),
            mean  = theta*sqrt(c(I1, I2)),
            sigma = rbind(c(1, sqrt(I1/I2)), c(sqrt(I1/I2), 1)))[1] - 0.5
}

produce_figure      <- function(output, theta, save_as) {
  len_theta  <- length(theta)
  E11        <- output$E11
  E12        <- output$E12
  E21        <- output$E21
  E22        <- output$E22
  S1         <- output$S1
  data       <- 
    tibble::tibble(theta                      = rep(theta, 9),
                   estimator                  = rep(c("MAE1", "MAE2", "CMLE",
                                                      "CMUE", "CUMVUE", "MLE",
                                                      "MUE", "UMVUE", "CWMAE"),
                                                    each = len_theta),
                   S1                         = rep(S1, 9),
                   `E(hat(theta)|theta,1)`    = as.vector(E11),
                   `E(hat(theta)^2|theta,1)`  = as.vector(E12),
                   `E(hat(theta)|theta,2)`    = as.vector(E21),
                   `E(hat(theta)^2|theta,2)`  = as.vector(E22),
                   `Var(hat(theta)|theta,1)`  =
                     `E(hat(theta)^2|theta,1)` - `E(hat(theta)|theta,1)`^2,
                   `Var(hat(theta)|theta,2)`  =
                     `E(hat(theta)^2|theta,2)` - `E(hat(theta)|theta,2)`^2,
                   `Bias(hat(theta)|theta,1)` = `E(hat(theta)|theta,1)` - theta,
                   `Bias(hat(theta)|theta,2)` = `E(hat(theta)|theta,2)` - theta,
                   `RMSE(hat(theta)|theta,1)` =
                     sqrt(`Var(hat(theta)|theta,1)` +
                            `Bias(hat(theta)|theta,1)`^2),
                   `RMSE(hat(theta)|theta,2)` =
                     sqrt(`Var(hat(theta)|theta,2)` +
                            `Bias(hat(theta)|theta,2)`^2),
                   `Bias(hat(theta)|theta)`   =
                     S1*`Bias(hat(theta)|theta,1)` +
                     (1 - S1)*`Bias(hat(theta)|theta,2)`,
                   `RMSE(hat(theta)|theta)`   =
                     S1*`RMSE(hat(theta)|theta,1)` +
                     (1 - S1)*`RMSE(hat(theta)|theta,2)`)
  A          <- ggplot(data, aes(x = theta, y = `Bias(hat(theta)|theta,1)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(Bias), "(", hat(theta), "|", theta, ",1)"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(-25, 20)
  B          <- ggplot(data, aes(x = theta, y = `RMSE(hat(theta)|theta,1)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(RMSE), "(", hat(theta), "|", theta, ",1)"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(0, 25)
  C          <- ggplot(data, aes(x = theta, y = `Bias(hat(theta)|theta,2)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(Bias), "(", hat(theta), "|", theta, ",2)"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(-1, 25)
  D          <- ggplot(data, aes(x = theta, y = `RMSE(hat(theta)|theta,2)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(RMSE), "(", hat(theta), "|", theta, ",2)"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(0, 25)
  E          <- ggplot(data, aes(x = theta, y = `Bias(hat(theta)|theta)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(Bias), "(", hat(theta), "|", theta, ")"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(-2, 1)
  Fp         <- ggplot(data, aes(x = theta, y = `RMSE(hat(theta)|theta)`,
                                 colour = estimator)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab(expression(theta)) +
    ylab(expression(paste(italic(RMSE), "(", hat(theta), "|", theta, ")"))) +
    scale_color_manual(values = c("#000000", "#999999", "#E69F00", "#56B4E9",
                                  "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                  "#CC79A7")) #+ ylim(0, 7.5)
  figure     <- ((A + B) / (C + D) / (E + Fp)) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
  ggsave(paste0(save_as, ".pdf"), figure, device = "pdf", width = 6.3,
         height = 9.7, units = "in")
}

umvue_2             <- function(z, I1, I2, l, u) {
  est   <- z/sqrt(I2) - ((I2 - I1)/(I2*sqrt(I1)))*
    (dnorm(u, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)) -
       dnorm(l, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)))/
    (pnorm(u, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)) -
       pnorm(l, z*sqrt(I1/I2), sqrt((I2 - I1)/I2)))
  if (is.infinite(est)) {
    est <- z/sqrt(I2)
  }
  est
}

##### Produce figure for Example 1 (Figure 1) ##################################

output1 <- est_gs(theta = seq(-2, 2, 0.05), I1 = 5.291785, I2 = 10.58357,
                  u = (sqrt(2/(1:2))*1.977)[1], l = -u, theta0 = 1)
produce_figure(output1, theta = seq(-2, 2, 0.05), save_as = "figure1")

##### Produce figures for Example 2 (Supplementary Figures 1-2) ################

output2  <- est_gs(theta = seq(-20, 20, 0.5), I1 = 0.0575, I2 = 0.086875, l = 0,
                   u = Inf, theta0 = 10)
# Note: Figure 2 is found by modifying axis limit of Supplementary Figure 1
produce_figure(output2, theta = seq(-20, 20, 0.5), save_as = "sfigure1")
# For better insight into the CMLE in this case, consider the conditional
# log-likelihood when (k,z) = (1,-0.01)
z        <- -0.01
I1       <- 0.0575
l        <- 0
data     <- tibble::tibble(theta    = seq(-10, 1000, length.out = 1e4),
                           clog_lik = dnorm(z - theta*sqrt(I1), log = TRUE) -
                             pnorm(l - theta*sqrt(I1), log.p = TRUE))
sfigure2 <- ggplot(data, aes(theta, clog_lik)) + geom_line() + theme_bw() +
  xlab(expression(theta)) + ylab("Conditional log-likelihood")
ggsave("sfigure2.pdf", sfigure2, device = "pdf", width = 6.3, height = 5,
       units = "in")

##### Produce figure for Example 3 (Figure 3) ##################################

output3 <- est_gs(theta = seq(-20, 20, 0.5), I1 = 0.03185288, I2 = 0.06370576,
                  l = 0.2694422, u = 2.72955, theta0 = 10)
produce_figure(output3, theta = seq(-20, 20, 0.5), save_as = "figure3")
