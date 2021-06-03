################################################################################
# Author(s):     Michael J Grayling (michael.grayling@newcastle.ac.uk)         #
# Description:   Reproduces the figures from 'Planning a randomised trial from #
#                an external pilot: When can we gain from Bayesian-frequentist #
#                methods?'                                                     #
# Filename:      figures.R                                                     #
# Last modified: 03 Jun 2021                                                   #
################################################################################

##### Functions ################################################################

theme_externalpilot             <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_family = base_family,
                      base_size   = base_size) +
    ggplot2::theme(axis.ticks       = ggplot2::element_line(colour = "grey70",
                                                            size   = 0.25),
                   complete         = T,
                   legend.key       = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   legend.position  = "bottom",
                   #legend.title     = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   panel.border     = ggplot2::element_rect(fill   = NA,
                                                            colour = "grey70",
                                                            size   = 0.5),
                   panel.grid.major = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.25),
                   panel.grid.minor = ggplot2::element_line(colour = "grey87",
                                                            size   = 0.125),
                   plot.margin      = ggplot2::unit(c(0.3, 0.5, 0.3, 0.3),
                                                    "cm"),
                   plot.title       = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_rect(fill   = "grey70",
                                                            colour = NA),
                   strip.text       =
                     ggplot2::element_text(colour = "white",
                                           size   = ggplot2::rel(0.8)))
}

initialise_data                 <- function(designs, density, samples) {
  len_designs              <- length(designs)
  opchar                   <- tibble::tibble(Design            = character(),
                                             nEP               = numeric(),
                                             alpha             = numeric(),
                                             beta              = numeric(),
                                             delta             = numeric(),
                                             sigma             = numeric(),
                                             method            = character(),
                                             q                 = numeric(),
                                             a                 = numeric(),
                                             b                 = numeric(),
                                             `design scenario` = character())
  for (i in 1:len_designs) {
    opchar                 <-
      tibble::add_case(opchar,
                       Design            = designs[[i]]$Design,
                       nEP               = designs[[i]]$nEP,
                       alpha             = designs[[i]]$alpha,
                       beta              = designs[[i]]$beta,
                       delta             = designs[[i]]$delta,
                       sigma             = designs[[i]]$sigma,
                       method            = designs[[i]]$method,
                       q                 = designs[[i]]$q,
                       a                 = designs[[i]]$a,
                       b                 = designs[[i]]$b,
                       `design scenario` = designs[[i]]$`design scenario`)
  }
  columns_to_add           <- c("Optimal n", "E(N)", "SD(N)", "E(P)", "SD(P)",
                                "E(CI width)", "SD(CI width)", "E(CrI width)",
                                "SD(CrI width)", "Linear", "Quadratic", "Mixed")
  opchar[, columns_to_add] <- NA_real_
  dist_n                   <-
    tibble::tibble(Design            = rep(opchar$Design, each = density),
                   `design scenario` = rep(opchar$`design scenario`,
                                           each = density),
                   n                 = NA_real_,
                   pdf               = NA_real_,
                   cdf               = NA_real_)
  dist_p                   <-
    tibble::tibble(Design            = rep(opchar$Design, each = density),
                   `design scenario` = rep(opchar$`design scenario`,
                                           each = density),
                   p                 = rep(seq(0, 1, length.out = density),
                                           len_designs),
                   pdf               = NA_real_,
                   cdf = NA_real_)
  if (is.null(samples)) {
    list(dist_n = dist_n,
         dist_p = dist_p,
         opchar = opchar)
  } else {
    list(dist_n   = dist_n,
         dist_p   = dist_p,
         opchar   = opchar,
         rsamples = tibble::tibble(Design     = rep(opchar$Design,
                                                    each = samples),
                                   sigma2_hat = NA_real_,
                                   n          = NA_real_,
                                   p          = NA_real_))
  }
}

int_internal                    <- function(designs, density, progress) {
  len_designs <- length(designs)
  int_data    <- initialise_data(designs, density, NULL)
  for (i in 1:len_designs) {
    if (all(progress, i == 1)) {
      message("  Beginning integral based evaluation for design ", i, " of ",
              len_designs, "..")
    } else if (all(progress, i < len_designs)) {
      message("..beginning integral based evaluation for design ", i, " of ",
              len_designs, "..")
    } else if (progress) {
      message("..beginning integral based evaluation for design ", i, " of ",
              len_designs, ".")
    }
    int_data  <- int_design(i, designs, int_data, density)
  }
  int_data
}

sim_internal                    <- function(designs, density, samples,
                                            progress) {
  len_designs <- length(designs)
  sim_data    <- initialise_data(designs, density, samples)
  for (i in 1:len_designs) {
    if (all(progress, i == 1)) {
      message("  Beginning simulation based evaluation for design ", i,
              " of ", len_designs, "..")
    } else if (all(progress, i < len_designs)) {
      message("..beginning simulation based evaluation for design ", i,
              " of ", len_designs, "..")
    } else if (progress) {
      message("..beginning simulation based evaluation for design ", i,
              " of ", len_designs, ".")
    }
    sim_data  <- sim_design(i, designs, density, samples, sim_data)
  }
  sim_data
}

int_cdf_n_internal              <- function(n, design) {
  seq_n <- seq(design$n_min, n, length.out = 100)
  pdf_n <- int_pdf_n(seq_n, design)
  cdf_n <- fda.usc::int.simpson2(seq_n, pdf_n, method = "TRAPZ")[1]
  ifelse(cdf_n > 1, 1, cdf_n)
}
int_cdf_n                       <- Vectorize(int_cdf_n_internal, "n")

int_cdf_power_internal          <- function(p, design) {
  seq_p <- seq(0, p, length.out = 100)
  pdf_p <- int_pdf_power(seq_p, design)
  cdf_p <- fda.usc::int.simpson2(seq_p, pdf_p, method = "TRAPZ")[1]
  ifelse(cdf_p > 1, 1, cdf_p)
}
int_cdf_power                   <- Vectorize(int_cdf_power_internal, "p")

ci_width                        <- function(sigma2_hat, design, exponent) {
  (design$nu*sigma2_hat/stats::qgamma(0.025, shape = 0.5*design$nu,
                                      rate = 0.5) -
     design$nu*sigma2_hat/stats::qgamma(0.975, shape = 0.5*design$nu,
                                        rate = 0.5))^exponent
}

ci_width_integrand              <- function(sigma2_hat, design, exponent) {
  ci_width(sigma2_hat, design, exponent)*int_pdf_sigma2_hat(sigma2_hat, design)
}

cri_width                       <- function(sigma2_hat, design, exponent) {
  b_posterior <- design$nu*sigma2_hat + design$b
  (invgamma::qinvgamma(0.975, shape = design$a_posterior, rate = b_posterior) -
      invgamma::qinvgamma(0.025, shape = design$a_posterior,
                          rate = b_posterior))^exponent
}

cri_width_integrand             <- function(sigma2_hat, design, exponent) {
  cri_width(sigma2_hat, design, exponent)*int_pdf_sigma2_hat(sigma2_hat, design)
}

int_df_n                        <- function(n, design) {
  len_n     <- length(n)
  df        <- numeric(len_n)
  for (i in 1:len_n) {
    if (any(n[i] < design$n_min, n[i] > design$n_max)) {
      df[i] <- 0
    } else {
      df[i] <- int_df_n_internal(n[i], design)
    }
  }
  df
}

int_df_n_internal               <- function(n, design) {
  int_pdf_sigma2_hat(
    stats::uniroot(root_int_df_n, c(design$sigma_min, design$sigma_max),
                   n       = n,
                   design  = design,
                   maxiter = 1e6L)$root, design)
}

int_df_power                    <- function(p, design) {
  int_pdf_n(((stats::qnorm(p) + design$c)/design$power_factor)^2, design)
}

e_ci_width                      <- function(design, exponent) {
  stats::integrate(ci_width_integrand, 0, Inf,
                   design   = design,
                   exponent = exponent)$value
}

e_cri_width                     <- function(design, exponent) {
  stats::integrate(cri_width_integrand, 0, Inf,
                   design   = design,
                   exponent = exponent)$value
}

e_n                             <- function(design, exponent) {
  seq_n <- seq(design$n_min, design$n_max, length.out = 100)
  fda.usc::int.simpson2(seq_n, e_n_integrand(seq_n, design, exponent),
                        method = "TRAPZ")[1]
}

e_n_integrand                   <- function(n, design, exponent) {
  n^exponent*int_pdf_n(n, design)
}

e_penalty                       <- function(design) {
  seq_n <- seq(design$n_min, design$n_max, length.out = 100)
  lin   <- penalty_integrand(seq_n, "lin",  design)
  quad  <- penalty_integrand(seq_n, "quad", design)
  mix   <- penalty_integrand(seq_n, "mix",  design)
  c(fda.usc::int.simpson2(seq_n, lin,  method = "TRAPZ")[1],
    fda.usc::int.simpson2(seq_n, quad, method = "TRAPZ")[1],
    fda.usc::int.simpson2(seq_n, mix,  method = "TRAPZ")[1])
}

e_power                         <- function(design, exponent) {
  seq_n <- seq(design$n_min, design$n_max, length.out = 100)
  fda.usc::int.simpson2(seq_n, power_integrand(seq_n, design, exponent),
                        method = "TRAPZ")[1]
}

n_f                             <- function(sigma2_hat, design) {
  if (design$method == "bayes_post") {
    n          <- n_bayes_posterior(sigma2_hat, design)
  } else if (design$method == "freq_t") {
    n          <- n_freq_t(sigma2_hat, design)
  } else {
    if (design$method == "bayes_mean") {
      sigma2_n <- (design$nu*sigma2_hat + design$b)/(design$a_posterior - 1)
    } else if (design$method == "bayes_mode") {
      sigma2_n <- (design$nu*sigma2_hat + design$b)/(design$a_posterior + 1)
    } else if (design$method == "bayes_q") {
      sigma2_n <- invgamma::qinvgamma(design$q,
                                      shape = design$a_posterior,
                                      rate  = design$nu*sigma2_hat + design$b)
    } else if (design$method == "freq_q") {
      sigma2_n <- design$nu*sigma2_hat/
        stats::qgamma(1 - design$q, shape = design$nu/2, rate = 0.5)
    } else if (design$method == "freq_z") {
      sigma2_n <- sigma2_hat
    }
    n          <- n_freq_z(sigma2_n, design)
  }
  n[which(n <= 2)] <- 2
  n
}

n_bayes_posterior_integrand     <- function(sigma2, sigma2_hat, design) {
  posterior_dist(sigma2, sigma2_hat, design)*n_freq_z(sigma2, design)
}

n_bayes_posterior_internal      <- function(sigma2_hat, design) {
  stats::integrate(n_bayes_posterior_integrand, 0, Inf,
                   sigma2_hat = sigma2_hat,
                   design     = design)$value
}
n_bayes_posterior               <- Vectorize(n_bayes_posterior_internal,
                                             "sigma2_hat")

n_freq_z                        <- function(sigma2, design) {
  design$n_factor*sigma2
}

n_freq_t                        <- function(sigma2, design) {
  stats::uniroot(root_n_freq_t, interval = c(2, .Machine$double.xmax),
                 sigma2  = sigma2,
                 design  = design,
                 maxiter = 1e5)$root
}

root_n_freq_t                   <- function(n, sigma2, design) {
  ncp   <- suppressWarnings(stats::qt(1 - design$alpha, n - 2))
  if (is.nan(ncp)) {
    ncp <- .Machine$double.xmax
  }
  o     <- n - 4*stats::qt(1 - design$beta, design$nEP - 2,
                           ncp)^2*sigma2/design$delta^2
  if (is.infinite(o)) {
    -.Machine$double.xmax
  } else {
    o
  }
}

int_pdf_n                       <- function(n, design) {
  int_df_n(n, design)/design$normalisation_n
}

int_pdf_power                   <- function(p, design) {
  int_df_power(p, design)/design$normalisation_power
}

int_pdf_sigma2_hat              <- function(sigma2_hat, design) {
  stats::dgamma(sigma2_hat,
                shape = design$nu,
                rate  = design$nu/design$sigma2_true)
}

penalty                         <- function(n, penalty, design) {
  if (penalty == "lin") {
    abs(n - design$n_optimal)
  } else if (penalty == "quad") {
    (n - design$n_optimal)^2
  } else if (penalty == "mix") {
    f           <- n
    which_le    <- which(n <= design$n_optimal)
    which_g     <- which(n > design$n_optimal)
    f[which_le] <- (f[which_le] - design$n_optimal)^2
    f[which_g]  <- abs(f[which_g] - design$n_optimal)
    f
  }
}

penalty_integrand               <- function(n, penalty, design) {
  penalty(n, penalty, design)*int_pdf_n(n, design)
}

posterior_dist                  <- function(sigma2, sigma2_hat, design) {
  invgamma::dinvgamma(sigma2,
                      shape = design$a_posterior,
                      rate  = design$nu*sigma2_hat + design$b)
}

power                           <- function(n, design) {
  stats::pnorm(sqrt(n)*design$power_factor - design$c)
}

power_integrand                 <- function(n, design, exponent) {
  power(n, design)^exponent*int_pdf_n(n, design)
}

prior                           <- function(sigma2, design) {
  invgamma::dinvgamma(sigma2, shape = design$a, rate = design$b)
}

quantile_int_cdf_n_internal     <- function(q, design) {
  stats::uniroot(root_quantile_int_cdf_n, c(design$n_min, design$n_max),
                 q      = q,
                 design = design)$root
}
quantile_int_cdf_n              <- Vectorize(quantile_int_cdf_n_internal, "q")

quantile_int_cdf_power_internal <- function(q, design) {
  stats::uniroot(root_quantile_int_cdf_power, c(0, 1),
                 q      = q,
                 design = design)$root
}
quantile_int_cdf_power          <- Vectorize(quantile_int_cdf_power_internal,
                                             "q")

root_int_df_n                   <- function(sigma2_hat, n, design) {
  n_f(sigma2_hat, design) - n
}

root_quantile_int_cdf_n         <- function(n, q, design) {
  int_cdf_n(n, design) - q
}

root_quantile_int_cdf_power     <- function(p, q, design) {
  int_cdf_power(p, design) - q
}

update_designs                  <- function(designs) {
  for (i in 1:length(designs)) {
    designs[[i]]$c            <- stats::qnorm(1 - designs[[i]]$alpha)
    designs[[i]]$sigma2_true  <- designs[[i]]$sigma^2
    designs[[i]]$n_factor     <-
      4*(designs[[i]]$c + stats::qnorm(1 - designs[[i]]$beta))^2/
      designs[[i]]$delta^2
    designs[[i]]$power_factor <-
      designs[[i]]$delta/sqrt(4*designs[[i]]$sigma2_true)
    designs[[i]]$n_optimal    <- n_freq_z(designs[[i]]$sigma2_true,
                                          designs[[i]])
    designs[[i]]$nu           <- 0.5*(designs[[i]]$nEP - 2)
    designs[[i]]$a_posterior  <- designs[[i]]$a + designs[[i]]$nu
  }
  designs
}

int_design                      <- function(i, designs, int_data, density) {
  design                             <- designs[[i]]
  design$sigma_min                   <-
    stats::qgamma(1e-16, shape = design$nu, rate = design$nu/design$sigma2_true)
  design$sigma_max                   <-
    stats::qgamma(1 - 1e-16,
                  shape = design$nu,
                  rate  = design$nu/design$sigma2_true)
  design$n_min                       <- n_f(design$sigma_min, design)
  design$n_max                       <- n_f(design$sigma_max, design)
  if (is.infinite(design$n_max)) {
    design$n_max                     <- .Machine$double.xmax
  }
  if (design$n_max > design$n_min) {
    design$normalisation_n             <-
      stats::integrate(int_df_n, design$n_min, design$n_max,
                       design = design)$value
    #design$normalisation_power         <-
    #  stats::integrate(int_df_power, 0, 1, design = design)$value
    int_data$opchar[i, 12]             <- design$n_optimal
    int_data$opchar[i, 13]             <- e_n(design, exponent = 1)
    int_data$opchar[i, 14]             <- sqrt(e_n(design, exponent = 2) -
                                                 int_data$opchar[i, 13]^2)
    int_data$opchar[i, 15]             <- e_power(design, exponent = 1)
    int_data$opchar[i, 16]             <- sqrt(e_power(design, exponent = 2) -
                                                 int_data$opchar[i, 15]^2)
    if (design$method %in% c("freq_q", "freq_t", "freq_z")) {
      int_data$opchar[i, 17]           <- e_ci_width(design, exponent = 1)
      int_data$opchar[i, 18]           <- sqrt(e_ci_width(design, exponent = 2) -
                                                 int_data$opchar[i, 17]^2)
    } else {
      int_data$opchar[i, 19]           <- e_cri_width(design, exponent = 1)
      int_data$opchar[i, 20]           <- sqrt(e_cri_width(design, exponent = 2) -
                                                 int_data$opchar[i, 19]^2)

    }
    penalties                          <- e_penalty(design)
    int_data$opchar[i, 21]             <- penalties[1]
    int_data$opchar[i, 22]             <- penalties[2]
    int_data$opchar[i, 23]             <- penalties[3]
    range_dist                         <- (1 + density*(i - 1)):(density*i)
    int_data$dist_n[range_dist, 3]     <- seq(design$n_min, design$n_max,
                                              length.out = density)
    p                                  <- int_data$dist_p$p[range_dist]
    n                                  <- int_data$dist_n$n[range_dist]
    int_data$dist_n[range_dist, 4]     <- int_pdf_n(n, design)
    int_data$dist_n[range_dist, 5]     <- int_cdf_n(n, design)
    #int_data$dist_p[range_dist, 4]     <- int_pdf_power(p, design)
    #int_data$dist_p[range_dist, 5]     <- int_cdf_power(p, design)
  } else {
    int_data$opchar[i, 12]             <- design$n_optimal
    int_data$opchar[i, 13]             <- design$n_min
    int_data$opchar[i, 14]             <- 0
    int_data$opchar[i, 15]             <- power(design$n_min, design)
    int_data$opchar[i, 16]             <- 0
    if (design$method %in% c("freq_q", "freq_t", "freq_z")) {
      int_data$opchar[i, 17]           <- e_ci_width(design, exponent = 1)
      int_data$opchar[i, 18]           <- sqrt(e_ci_width(design, exponent = 2) -
                                                 int_data$opchar[i, 17]^2)
    } else {
      int_data$opchar[i, 19]           <- e_cri_width(design, exponent = 1)
      int_data$opchar[i, 20]           <- sqrt(e_cri_width(design, exponent = 2) -
                                                 int_data$opchar[i, 19]^2)

    }
    int_data$opchar[i, 21]             <- penalty(design$n_min, "lin", design)
    int_data$opchar[i, 22]             <- penalty(design$n_min, "quad", design)
    int_data$opchar[i, 23]             <- penalty(design$n_min, "mix", design)
    #range_dist                         <- (1 + density*(i - 1)):(density*i)
    #int_data$dist_n[range_dist, 3]     <- seq(design$n_min, design$n_max,
    #                                          length.out = density)
    #p                                  <- int_data$dist_p$p[range_dist]
    #n                                  <- int_data$dist_n$n[range_dist]
    #int_data$dist_n[range_dist, 4]     <- int_pdf_n(n, design)
    #int_data$dist_n[range_dist, 5]     <- int_cdf_n(n, design)
    #int_data$dist_p[range_dist, 4]     <- int_pdf_power(p, design)
    #int_data$dist_p[range_dist, 5]     <- int_cdf_power(p, design)
  }
  int_data
}

sim_design                      <- function(i, designs, density, samples,
                                            sim_data) {
  design                                      <- designs[[i]]
  range_samples                               <-
    (1 + samples*(i - 1)):(samples*i)
  sigma2_hat                                  <-
    stats::rgamma(samples,
                  shape = design$nu,
                  rate  = design$nu/design$sigma2_true)
  sim_data$rsamples$sigma2_hat[range_samples] <- sigma2_hat
  sim_data$rsamples$n[range_samples]          <- n_f(sigma2_hat, design)
  sim_data$rsamples$p[range_samples]          <-
    power(sim_data$rsamples$n[range_samples], design)
  sim_data$opchar[i, 12]                      <- design$n_optimal
  sim_data$opchar[i, 13]                      <-
    mean(sim_data$rsamples$n[range_samples])
  sim_data$opchar[i, 14]                      <-
    sqrt(mean(sim_data$rsamples$n[range_samples]^2) - sim_data$opchar[i, 13]^2)
  sim_data$opchar[i, 15]                      <-
    mean(sim_data$rsamples$p[range_samples])
  sim_data$opchar[i, 16]                      <-
    sqrt(mean(sim_data$rsamples$p[range_samples]^2) - sim_data$opchar[i, 15]^2)
  if (design$method %in% c("freq_q", "freq_t", "freq_z")) {
    sim_data$opchar[i, 17]                    <-
      mean(ci_width(sigma2_hat, design, exponent = 1))
    sim_data$opchar[i, 18]                    <-
      sqrt(mean(ci_width(sigma2_hat, design, exponent = 2)) -
             sim_data$opchar[i, 17]^2)
  } else {
    sim_data$opchar[i, 19]                    <-
      mean(cri_width(sigma2_hat, design, exponent = 1))
    sim_data$opchar[i, 20]                    <-
      sqrt(mean(cri_width(sigma2_hat, design, exponent = 2)) -
             sim_data$opchar[i, 19]^2)

  }
  sim_data$opchar[i, 21]                      <-
    mean(penalty(sim_data$rsamples$n[range_samples], "lin",  design))
  sim_data$opchar[i, 22]                      <-
    mean(penalty(sim_data$rsamples$n[range_samples], "quad", design))
  sim_data$opchar[i, 23]                      <-
    mean(penalty(sim_data$rsamples$n[range_samples], "mix",  design))
  range_dist                                  <-
    (1 + density*(i - 1)):(density*i)
  sim_data$dist_n[range_dist, 3]              <-
    seq(min(sim_data$rsamples$n[range_samples]),
        max(sim_data$rsamples$n[range_samples]), length.out = density)
  for (j in 1:density) {
    sim_data$dist_n[range_dist[j], 5]         <-
      length(which(sim_data$rsamples$n[range_samples] <=
                     sim_data$dist_n[range_dist[j], 2]))
    sim_data$dist_p[range_dist[j], 5]         <-
      length(which(sim_data$rsamples$p[range_samples] <=
                     sim_data$dist_p[range_dist[j], 2]))
  }
  sim_data$dist_n[range_dist, 5]              <-
    sim_data$dist_n[range_dist, 3]/samples
  sim_data$dist_p[range_dist, 5]              <-
    sim_data$dist_p[range_dist, 3]/samples
  sim_data
}

sim <- function(designs, samples = 10000L, cpus = 1, progress = F) {

  ##### Input checks ###########################################################

  check_designs(designs)
  check_samples(samples)
  check_cpus(cpus)
  check_progress(progress)

  ##### Main computations ######################################################

  outputs          <- list()
  outputs$sim_data <- sim_internal(designs, samples, cpus, progress)
  outputs$plots    <- plot_sim(outputs$sim_data)

  ##### Outputting #############################################################

  outputs

}

int <- function(designs, cpus = 1, progress = F) {

  ##### Input checks ###########################################################

  check_designs(designs)
  check_cpus(cpus)
  check_progress(progress)

  ##### Main computations ######################################################

  designs          <- update_designs(designs)
  outputs          <- list()
  outputs$int_data <- int_internal(designs, cpus, progress)
  outputs$plots    <- plot_int(outputs$int_data)

  ##### Outputting #############################################################

  outputs

}

des <- function(nEP = 70, alpha = 0.05, beta = 0.2, delta = 0.2, sigma = 1,
                method = "freq_z", q, a, b, ab) {

  ##### Input checks ###########################################################

  check_nEP(nEP)
  check_alpha(alpha)
  check_beta(beta)
  check_delta(delta)
  check_sigma(sigma)
  check_method(method)
  check_q(q, method)
  #check_ab(a, b, ab, method)

  ##### Main computations ######################################################

  outputs                                  <- list()
  counter                                  <- 1L
  for (i in 1:length(nEP)) {
    for (j in 1:length(alpha)) {
      for (k in 1:length(beta)) {
        for (l in 1:length(delta)) {
          for (m in 1:length(sigma)) {
            for (n in 1:length(method)) {
              if (method[n] == "freq_q") {
                for (o in 1:length(q)) {
                  outputs[[counter]]       <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = q[o],
                                                   a      = NA_real_,
                                                   b      = NA_real_)
                  counter                  <- counter + 1L
                }
              } else if (method[n] %in% c("bayes_mean", "bayes_mode",
                                          "bayes_post")) {
                if (missing(ab)) {
                  for (p in 1:length(a)) {
                    if (any(method[n] != "bayes_mean",
                            all(method[n] == "bayes_mean", a[p] > 1))) {
                      for (r in 1:length(b)) {
                        outputs[[counter]] <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = NA,
                                                   a      = a[p],
                                                   b      = b[r])
                        counter            <- counter + 1L
                      }
                    }
                  }
                } else {
                  for (p in 1:nrow(ab)) {
                    if (any(method[n] != "bayes_mean",
                            all(method[n] == "bayes_mean", ab[p, 1] > 1))) {
                      outputs[[counter]]   <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = NA,
                                                   a      = ab[p, 1],
                                                   b      = ab[p, 2])
                      counter              <- counter + 1L
                    }
                  }
                }
              } else if (method[n] == "bayes_q") {
                for (o in 1:length(q)) {
                  if (missing(ab)) {
                    for (p in 1:length(a)) {
                      for (r in 1:length(b)) {
                        outputs[[counter]] <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = q[o],
                                                   a      = a[p],
                                                   b      = b[r])
                        counter            <- counter + 1L
                      }
                    }
                  } else {
                    for (p in 1:nrow(ab)) {
                      outputs[[counter]]   <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = q[o],
                                                   a      = ab[p, 1],
                                                   b      = ab[p, 2])
                      counter              <- counter + 1L
                    }
                  }
                }
              } else {
                outputs[[counter]]         <- list(nEP    = nEP[i],
                                                   alpha  = alpha[j],
                                                   beta   = beta[k],
                                                   delta  = delta[l],
                                                   sigma  = sigma[m],
                                                   method = method[n],
                                                   q      = NA,
                                                   a     = NA,
                                                   b     = NA)
                counter                    <- counter + 1L
              }
            }
          }
        }
      }
    }
  }
  method_names                             <- c("bayes_mean" = "Bayes mean",
                                                "bayes_mode" = "Bayes mode",
                                                "bayes_post" = "Bayes post",
                                                "bayes_q"    = "Bayes q",
                                                "freq_q"     = "Freq q",
                                                "freq_t"     = "Freq t",
                                                "freq_z"     = "Freq z")
  len_outputs                              <- length(outputs)
  if (length(unique(nEP)) > 1) {
    for (i in 1:len_outputs) {
      outputs[[i]]$Design                  <- paste("nEP =", outputs[[i]]$nEP)
    }
    prefix                                 <- ", "
  } else {
    prefix                                 <- ""
  }
  if (length(unique(method)) > 1) {
    for (i in 1:len_outputs) {
      outputs[[i]]$Design                  <-
        paste0(outputs[[i]]$Design, prefix, "method = ",
               method_names[outputs[[i]]$method])
    }
    prefix                                 <- ", "
  }
  for (i in 1:len_outputs) {
    if (outputs[[i]]$method %in% c("bayes_q", "freq_q")) {
      outputs[[i]]$Design                  <-
        paste0(outputs[[i]]$Design, prefix, "q = ", outputs[[i]]$q)
    }
    if (outputs[[i]]$method %in% c("bayes_mean", "bayes_mode", "bayes_post",
                                   "bayes_q")) {
      outputs[[i]]$Design                  <-
        paste0(outputs[[i]]$Design, prefix, "a = ", outputs[[i]]$a, ", b = ",
               outputs[[i]]$b)
    }
    outputs[[i]]$`design scenario`         <-
      paste0("alpha = ", outputs[[i]]$alpha, ", beta = ", outputs[[i]]$beta,
             ", delta = ", outputs[[i]]$delta, ", sigma = ", outputs[[i]]$sigma)
  }

  ##### Outputting #############################################################

  outputs

}

comp <- function(designs, samples = 10000L, cpus = 1, progress = F) {

  ##### Input checks ###########################################################

  check_designs(designs)
  check_samples(samples)
  check_cpus(cpus)
  check_progress(progress)

  ##### Main computations ######################################################

  designs          <- update_designs(designs)
  outputs          <- list()
  outputs$int_data <- int_internal(designs, cpus, progress)
  outputs$sim_data <- sim_internal(designs, samples, cpus, progress)
  outputs$plots    <- plot_comp(outputs$int_data, outputs$sim_data)

  ##### Outputting #############################################################

  outputs

}

check_belong       <- function(input, input_name, allowed) {
  if (!all(input %in% allowed)) {
    stop(input_name, " must contain values in ",
         paste(allowed, collapse = ", "), " only")
  }
}

check_duplicates   <- function(input, input_name) {
  if (length(unique(input)) < length(input)) {
    warning(input_name, " contains duplicated values")
  }
}

check_finite       <- function(input, input_name) {
  if (any(is.infinite(input))) {
    if (length(input) == 1) {
      stop(input_name, " must be finite")
    } else {
      stop(input_name, " must contain finite values only")
    }
  }
}

check_g            <- function(input, input_name, g = 0) {
  if (any(input <= g)) {
    if (length(input) == 1) {
      stop(input_name, " must be greater than ", g)
    } else {
      stop(input_name, " must contain values greater than ", g, " only")
    }
  }
}

check_ge           <- function(input, input_name, ge = 0) {
  if (any(input < ge)) {
    if (length(input) == 1) {
      stop(input_name, " must be greater than or equal to ", ge)
    } else {
      stop(input_name, " must contain values greater than or equal to ", ge,
           " only")
    }
  }
}

check_integer      <- function(input, input_name) {
  if (any(input%%1 != 0)) {
    if (length(input) == 1) {
      stop(input_name, " must be a whole number")
    } else {
      stop(input_name, " must contain whole numbers only")
    }
  }
}

check_l            <- function(input, input_name, l = 1) {
  if (any(input >= l)) {
    if (length(input) == 1) {
      stop(input_name, " must be less than ", l)
    } else {
      stop(input_name, " must contain values less than ", l, " only")
    }
  }
}

check_le           <- function(input, input_name, le = 1) {
  if (any(input > le)) {
    if (length(input) == 1) {
      stop(input_name, " must be less than or equal to ", le)
    } else {
      stop(input_name, " must contain values less than or equal to ", le,
           " only")
    }
  }
}

check_length       <- function(input, input_name, len = 1) {
  if (length(input) != len) {
    stop(input_name, " must be of length ", len)
  }
}

check_logical      <- function(input, input_name) {
  if (!is.logical(input)) {
    stop(input_name, " must be logical")
  }
}

check_not_missing  <- function(input, input_name, condition, condition_name) {
  if (condition) {
    warning(input_name, " has been specified but this will have no affect ",
            "given the value of ", condition_name)
  }
}

check_numeric      <- function(input, input_name) {
  if (!is.numeric(input)) {
    stop(input_name, " must be numeric")
  }
}

check_vector       <- function(input, input_name) {
  if (!is.vector(input)) {
    stop(input_name, " must be a vector")
  }
}

check_a0           <- function(a0, method) {
  if (missing(a0)) {
    if (any(c("bayes_q", "freq_q") %in% method)) {
      stop("a0 must be specified when method contains bayes_q or freq_q")
    }
  } else {
    if (any(c("bayes_q", "bayes_mean", "bayes_mode",
              "bayes_post") %in% method)) {
      check_numeric(a0, "a0")
      check_vector(a0, "a0")
      check_g(a0, "a0")
      check_finite(a0, "a0")
      check_duplicates(a0, "a0")
    } else {
      check_not_missing(a0, "a0", T, "method")
    }
  }
}

check_alpha        <- function(alpha) {
  check_numeric(alpha, "alpha")
  check_vector(alpha, "alpha")
  check_g(alpha, "alpha")
  check_l(alpha, "alpha")
  check_duplicates(alpha, "alpha")
}

check_b0           <- function(b0, method) {
  if (missing(b0)) {
    if (any(c("bayes_q", "freq_q") %in% method)) {
      stop("b0 must be specified when method contains bayes_q or freq_q")
    }
  } else {
    if (any(c("bayes_q", "bayes_mean", "bayes_mode",
              "bayes_post") %in% method)) {
      check_numeric(b0, "b0")
      check_vector(b0, "b0")
      check_g(b0, "b0")
      check_finite(b0, "b0")
      check_duplicates(b0, "b0")
    } else {
      check_not_missing(b0, "b0", T, "method")
    }
  }
}

check_beta         <- function(beta) {
  check_numeric(beta, "beta")
  check_vector(beta, "beta")
  check_g(beta, "beta")
  check_l(beta, "beta")
  check_duplicates(beta, "beta")
}

check_cpus         <- function(cpus) {
  check_numeric(cpus, "cpus")
  check_vector(cpus, "cpus")
  check_length(cpus, "cpus")
  check_ge(cpus, "cpus", 1)
  check_finite(cpus, "cpus")
  check_integer(cpus, "cpus")
}

check_delta        <- function(delta) {
  check_numeric(delta, "delta")
  check_vector(delta, "delta")
  check_g(delta, "delta")
  check_finite(delta, "delta")
  check_duplicates(delta, "delta")
}

check_designs      <- function(designs) {
  if (!is.list(designs)) {
    stop("designs must be a list")
  }
  for (i in 1:length(designs)) {
    if (any(is.null(designs[[i]]$nEP), is.null(designs[[i]]$alpha),
            is.null(designs[[i]]$beta), is.null(designs[[i]]$delta),
            is.null(designs[[i]]$sigma), is.null(designs[[i]]$method),
            is.null(designs[[i]]$q), is.null(designs[[i]]$a0),
            is.null(designs[[i]]$b0))) {
      stop("designs[[", i, "]] is missing required design components")
    }
  }
}

check_include_freq <- function(include_freq) {
  check_logical(include_freq, "include_freq")
}

check_method       <- function(method) {
  check_belong(method, "method", c("bayes_mean", "bayes_mode", "bayes_post",
                                   "bayes_q", "freq_q", "freq_t", "freq_z"))
  check_duplicates(method, "method")
}

check_nEP          <- function(nEP) {
  check_numeric(nEP, "nEP")
  check_vector(nEP, "nEP")
  check_g(nEP, "nEP")
  check_finite(nEP, "nEP")
}

check_progress     <- function(progress) {
  check_logical(progress, "progress")
}

check_q            <- function(q, method) {
  if (missing(q)) {
    if (any(c("bayes_q", "freq_q") %in% method)) {
      stop("q must be specified when method contains bayes_q or freq_q")
    }
  } else {
    if (any(c("bayes_q", "freq_q") %in% method)) {
      check_numeric(q, "q")
      check_vector(q, "q")
      check_g(q, "q")
      check_l(q, "q")
      check_duplicates(q, "q")
    } else {
      check_not_missing(q, "q", T, "method")
    }
  }
}

check_samples      <- function(samples) {
  check_numeric(samples, "samples")
  check_vector(samples, "samples")
  check_numeric(samples, "samples")
  check_vector(samples, "samples")
  check_length(samples, "samples")
  check_ge(samples, "samples", 1)
  check_finite(samples, "samples")
  check_integer(samples, "samples")
}

check_sigma        <- function(sigma) {
  check_numeric(sigma, "sigma")
  check_vector(sigma, "sigma")
  check_g(sigma, "sigma")
  check_finite(sigma, "sigma")
  check_duplicates(sigma, "sigma")
}

check_sigma_hat    <- function(sigma_hat) {
  check_numeric(sigma_hat, "sigma_hat")
  check_vector(sigma_hat, "sigma_hat")
  check_g(sigma_hat, "sigma_hat")
  check_finite(sigma_hat, "sigma_hat")
  check_duplicates(sigma_hat, "sigma_hat")
}

library(ggplot2)

##### Figure 1 #################################################################

designs                     <- des(nEP    = 70,
                                   alpha  = 0.025,
                                   beta   = 0.1,
                                   delta  = 0.35,
                                   sigma  = 1,
                                   method = c("bayes_mean", "bayes_mode",
                                              "bayes_post", "bayes_q", "freq_q",
                                              "freq_t", "freq_z"),
                                   q      = 0.5,
                                   ab     = matrix(c(41, 40), 1, 2))
designs                     <- update_designs(designs)
int_data                    <- int_internal(designs, 1000L, T)
names_des                   <-
  c("Posterior mean", "Posterior mode",
    "Posterior distribution",
    "Posterior quantile (*q* = 0.5)",
    "Quantile of the pooled sample variance (*q* = 0.5)",
    "Pooled sample variance with t inflation",
    "Pooled sample variance")
names_levels                <-
  c("Pooled sample variance", "Quantile of the pooled sample variance (*q* = 0.5)",
    "Pooled sample variance with t inflation", "Posterior mean",
    "Posterior mode", "Posterior quantile (*q* = 0.5)",
    "Posterior distribution")
int_data$dist_n$plot_design <- factor(rep(names_des, each = 1000),
                                      levels = names_levels)
int_data$dist_n$type        <- factor(c(rep("Bayesian-frequentist", 4000),
                                        rep("Frequentist", 3000)),
                                      c("Bayesian-frequentist", "Frequentist"))
int_data$dist_n             <- dplyr::filter(int_data$dist_n, n <= 750)
figure1                     <-
  ggplot2::ggplot(int_data$dist_n,
                  ggplot2::aes(n, cdf, colour = plot_design,
                               linetype = type)) +
  ggplot2::geom_vline(data = int_data$opchar,
                      ggplot2::aes(xintercept = `Optimal n`), linetype = 2,
                      size = 1.05) +
  ggplot2::geom_line(size = 1.05) +
  ggplot2::xlab(expression(italic(n))) +
  ggplot2::ylab(expression(paste(italic(F)[italic(N)], "(", italic(n), ")",
                                 sep = ""))) +
  theme_externalpilot() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2),
                  linetype = ggplot2::guide_legend(ncol = 1)) +
  ggplot2::theme(legend.title = ggplot2::element_blank()) +
  ggplot2::theme(legend.text = ggtext::element_markdown(size = 7)) +
  ggplot2::scale_colour_manual(values = as.vector(khroma::colour("muted")(7)),
                               labels = names_levels)
ggplot2::ggsave("figure1.pdf", figure1, device = "pdf", height = 7,
                width = 6.3, units = "in")

##### Supplementary Figures 1-2 ################################################

designs         <- des(nEP    = seq(20, 250, 10),
                       alpha  = 0.025,
                       beta   = 0.1,
                       delta  = 0.35,
                       sigma  = sqrt(c(0.5, 1, 2)),
                       method = "bayes_mean",
                       ab     = rbind(c(1.01, 0.01),
                                      c(2, 1),
                                      c(6, 5),
                                      c(11, 10),
                                      c(21, 20),
                                      c(31, 30),
                                      c(41, 40)))
designs         <- update_designs(designs)
int_data        <- int_internal(designs, 10L, T)
int_data$opchar <- dplyr::arrange(int_data$opchar, sigma, a, b, nEP)
omit            <- seq(from = 1, to = 504 - 23, by = 24)
seq1            <- (1:nrow(int_data$opchar))[-omit]
seq2            <- seq1 - 1
plot_data       <-
  tibble::tibble(nEP                             = int_data$opchar$nEP[-omit],
                 a                               = int_data$opchar$a[-omit],
                 b                               = int_data$opchar$b[-omit],
                 sigma2                          =
                   paste("sigma^2 ==", int_data$opchar$sigma[-omit]^2),
                 `italic(RECrIW(n[EP],sigma^2))` =
                   100*(1 - int_data$opchar$`E(CrI width)`[seq1]/
                          int_data$opchar$`E(CrI width)`[seq2]),
                 `italic(ECrIW(n[EP],sigma^2))`  =
                   int_data$opchar$`E(CrI width)`[-omit])
plot_data$ab    <- paste0("a = ", plot_data$a, ", b = ", plot_data$b)
plot_data$ab    <- factor(plot_data$ab, levels = c("a = 1.01, b = 0.01",
                                                   "a = 2, b = 1",
                                                   "a = 6, b = 5",
                                                   "a = 11, b = 10",
                                                   "a = 21, b = 20",
                                                   "a = 31, b = 30",
                                                   "a = 41, b = 40"))
plot_data       <-
  tidyr::gather(plot_data, key = "key", value = "value",
                `italic(RECrIW(n[EP],sigma^2))`:`italic(ECrIW(n[EP],sigma^2))`)
supfigure1         <-
  ggplot2::ggplot(plot_data, ggplot2::aes(nEP, value, colour = ab)) +
  ggplot2::geom_line(size = 0.75) + ggplot2::geom_point(size = 0.75) +
  ggplot2::facet_grid(key~sigma2,
                      scales   = "free_y",
                      labeller = ggplot2::label_parsed) +
  ggplot2::scale_colour_viridis_d(labels = c("*a* = 1.01, *b* = 0.01",
                                             "*a* = 2, *b* = 1",
                                             "*a* = 6, *b* = 5",
                                             "*a* = 11, *b* = 10",
                                             "*a* = 21, *b* = 20",
                                             "*a* = 31, *b* = 30",
                                             "*a* = 41, *b* = 40")) +
  ggplot2::xlab(expression(italic(n[EP]))) +
  ggplot2::ylab("Value") +
  theme_externalpilot() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::theme(legend.title = ggplot2::element_blank()) +
  ggplot2::theme(legend.text = ggtext::element_markdown())
ggplot2::ggsave("supfigure1.pdf", supfigure1, device = "pdf", height = 7,
                width = 6.3, units = "in")

hellinger <- function(a, b, nEP, sigma2_hat) {
  a1 <- a + 0.5*(nEP - 2)
  b1 <- b + 0.5*(nEP - 2)*sigma2_hat
  a2 <- a + 0.5*(nEP - 10 - 2)
  b2 <- b + 0.5*(nEP - 10 - 2)*sigma2_hat
  1 - exp(0.5*a1*log(b1) + 0.5*a2*log(b2) + lgamma(0.5*(a1 + a2)) -
            0.5*lgamma(a1) - 0.5*lgamma(a2) - 0.5*(a1 + a2)*log(0.5*(b1 + b2)))
}

ab           <- rbind(c(1.01, 0.01),
                      c(2, 1),
                      c(6, 5),
                      c(11, 10),
                      c(21, 20),
                      c(31, 30),
                      c(41, 40))
nEP          <- seq(20, 250, 10)
sigma2_hat   <- seq(from       = qgamma(0.01, 0.5*(20 - 2), 0.5*(20 - 2)),
                    to         = qgamma(0.99, 0.5*(20 - 2), 0.5*(20 - 2)),
                    length.out = 1000)
data         <- expand.grid(nEP, sigma2_hat, 1:7)
data         <- as.matrix(cbind(data[, 1:2], ab[data[, 3], ], 0))
for (i in 1:nrow(data)) {
  data[i, 5] <- hellinger(data[i, 3], data[i, 4], data[i, 1], data[i, 2])
}
data         <- tibble::tibble(nEP            = data[, 1],
                               `hat(sigma)^2` = data[, 2],
                               a              = data[, 3],
                               b              = data[, 4],
                               H              = data[, 5])
data$ab      <- paste("italic(a) == ", data$a, "~~italic(b) == ", data$b,
                      sep = "")
data$ab      <- factor(data$ab,
                       levels = c("italic(a) == 1.01~~italic(b) == 0.01",
                                  "italic(a) == 2~~italic(b) == 1",
                                  "italic(a) == 6~~italic(b) == 5",
                                  "italic(a) == 11~~italic(b) == 10",
                                  "italic(a) == 21~~italic(b) == 20",
                                  "italic(a) == 31~~italic(b) == 30",
                                  "italic(a) == 41~~italic(b) == 40"))
omit         <- seq(from = 1, to = 168000 - 23, by = 24)
seq1         <- (1:nrow(data))[-omit]
seq2         <- seq1 - 1
plot_data    <- tibble::tibble(nEP            = data$nEP[-omit],
                               `hat(sigma)^2` = data$`hat(sigma)^2`[-omit],
                               a              = data$a[-omit],
                               b              = data$b[-omit],
                               ab             = data$ab[-omit],
                               `rel H`        =
                                 100*(1 - data$H[seq1]/data$H[seq2]))
supfigure2   <-
  ggplot2::ggplot(plot_data, ggplot2::aes(nEP, `hat(sigma)^2`, z = `rel H`)) +
  ggplot2::geom_contour_filled() +
  ggplot2::facet_wrap(~ab, labeller = ggplot2::label_parsed) +
  ggplot2::scale_fill_viridis_d() +
  ggplot2::xlab(expression(italic(n[EP]))) +
  ggplot2::ylab(expression(hat(sigma)^2)) +
  theme_externalpilot() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::theme(legend.title = ggplot2::element_blank())
ggplot2::ggsave("supfigure2.pdf", supfigure2, device = "pdf", height = 7,
                width = 6.3, units = "in")

##### Figure 2 and Supplementary Figures 3-4 ###################################

designs   <- des(nEP    = 70,
                 alpha  = 0.025,
                 beta   = 0.1,
                 delta  = 0.35,
                 sigma  = sqrt(c(emdbook::lseq(0.1, 1, ceiling(101/2)),
                                 emdbook::lseq(1, 10, ceiling(101/2))[-1])),
                 method = c("bayes_mean", "bayes_mode", "bayes_post",
                            "bayes_q", "freq_q", "freq_t", "freq_z"),
                 q      = c(0.5, 0.9),
                 ab     = rbind(c(1.01, 0.01),
                                c(2, 1),
                                c(6, 5),
                                c(11, 10),
                                c(21, 20),
                                c(31, 30),
                                c(41, 40)))
designs   <- update_designs(designs)
st_time   <- Sys.time()
int_data  <- int_internal(designs, 1, T)
end_time  <- Sys.time()
data_wide <- int_data$opchar
designs   <- des(nEP    = 70,
                 alpha  = 0.025,
                 beta   = 0.1,
                 delta  = 0.35,
                 sigma  = sqrt(c(emdbook::lseq(0.5, 1, ceiling(101/2)),
                                 emdbook::lseq(1, 2, ceiling(101/2))[-1])),
                 method = c("bayes_mean", "bayes_mode", "bayes_post",
                            "bayes_q", "freq_q", "freq_t", "freq_z"),
                 q      = c(0.5, 0.9),
                 ab     = rbind(c(1.01, 0.01),
                                c(2, 1),
                                c(6, 5),
                                c(11, 10),
                                c(21, 20),
                                c(31, 30),
                                c(41, 40)))
designs   <- update_designs(designs)
st_time   <- Sys.time()
int_data  <- int_internal(designs, 1, T)
end_time  <- Sys.time()
data_zoom <- int_data$opchar
data      <- dplyr::bind_rows(data_wide, data_zoom)
data_wide <- data_wide_3
data_zoom <- data_zoom_3
data_wide_freq         <- dplyr::filter(data_wide, method %in% c("freq_q", "freq_t", "freq_z"))
data_wide_freq$bayes   <- "Frequentist"
data_wide_freq$method2 <- NA
data_wide_freq$method2[which(data_wide_freq$method == "freq_q" & data_wide_freq$q == 0.5)] <- "Quantile of the pooled sample variance (q = 0.5)"
data_wide_freq$method2[which(data_wide_freq$method == "freq_q" & data_wide_freq$q == 0.9)] <- "Quantile of the pooled sample variance (q = 0.9)"
data_wide_freq$method2[which(data_wide_freq$method == "freq_t")] <- "Pooled sample variance with t inflation"
data_wide_freq$method2[which(data_wide_freq$method == "freq_z")] <- "Pooled sample variance"
data_wide_bayes         <- dplyr::filter(data_wide, !(method %in% c("freq_q", "freq_t", "freq_z")))
data_wide_bayes$bayes   <- "Bayesian-Frequentist"
data_wide_bayes$method2 <- NA
data_wide_bayes$method2[which(data_wide_bayes$method == "bayes_mean")] <- "Posterior mean"
data_wide_bayes$method2[which(data_wide_bayes$method == "bayes_mode")] <- "Posterior mode"
data_wide_bayes$method2[which(data_wide_bayes$method == "bayes_post")] <- "Posterior distribution"
data_wide_bayes$method2[which(data_wide_bayes$method == "bayes_q" & data_wide_bayes$q == 0.5)] <- "Posterior quantile (q = 0.5)"
data_wide_bayes$method2[which(data_wide_bayes$method == "bayes_q" & data_wide_bayes$q == 0.9)] <- "Posterior quantile (q = 0.9)"
data_wide_bayes$ab <- paste("italic(a) == ", data_wide_bayes$a, "~~italic(b) == ", data_wide_bayes$b,
                                            sep = "")
data_wide_freq_1 <- data_wide_freq
data_wide_freq_1$ab <- "italic(a) == 1.01~~italic(b) == 0.01"
data_wide_freq_2 <- data_wide_freq
data_wide_freq_2$ab <- "italic(a) == 2~~italic(b) == 1"
data_wide_freq_3 <- data_wide_freq
data_wide_freq_3$ab <- "italic(a) == 6~~italic(b) == 5"
data_wide_freq_4 <- data_wide_freq
data_wide_freq_4$ab <- "italic(a) == 11~~italic(b) == 10"
data_wide_freq_5 <- data_wide_freq
data_wide_freq_5$ab <- "italic(a) == 21~~italic(b) == 20"
data_wide_freq_6 <- data_wide_freq
data_wide_freq_6$ab <- "italic(a) == 31~~italic(b) == 30"
data_wide_freq_7 <- data_wide_freq
data_wide_freq_7$ab <- "italic(a) == 41~~italic(b) == 40"
data_wide <- rbind(data_wide_bayes,
                   data_wide_freq_1,
                   data_wide_freq_2,
                   data_wide_freq_3,
                   data_wide_freq_4,
                   data_wide_freq_5,
                   data_wide_freq_6,
                   data_wide_freq_7)

data_wide$sigma2 <- data_wide$sigma^2
data_wide        <- tidyr::gather(data_wide, key = "key", value = "value",
                                  Linear:Mixed)
data_wide$key    <- factor(data_wide$key, levels = c("Linear", "Quadratic", "Mixed"))
data_wide$ab     <- factor(data_wide$ab,
                           levels = c("italic(a) == 1.01~~italic(b) == 0.01",
                                      "italic(a) == 2~~italic(b) == 1",
                                      "italic(a) == 6~~italic(b) == 5",
                                      "italic(a) == 11~~italic(b) == 10",
                                      "italic(a) == 21~~italic(b) == 20",
                                      "italic(a) == 31~~italic(b) == 30",
                                      "italic(a) == 41~~italic(b) == 40"))
data_wide$method2 <- factor(data_wide$method2,
                            levels = c("Pooled sample variance",
                                       "Pooled sample variance with t inflation",
                                       "Quantile of the pooled sample variance (q = 0.5)",
                                       "Quantile of the pooled sample variance (q = 0.9)",
                                       "Posterior mean",
                                       "Posterior mode",
                                       "Posterior quantile (q = 0.5)",
                                       "Posterior quantile (q = 0.9)",
                                       "Posterior distribution"))

supfigure4 <- ggplot(data_wide, aes(sigma2, as.numeric(value), colour = method2,
                                 linetype = bayes)) +
  geom_line() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  ph2rand:::theme_ph2rand() + xlab(expression(sigma^2)) + ylab("Loss") +
  facet_grid(key~ab, labeller=label_parsed,
                                                   scales = "free_y") +
  ggplot2::scale_colour_manual(values = as.vector(khroma::colour("muted")(9))) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.text=element_text(size=5))

ggsave("supfigure4.pdf", supfigure3, device = "pdf", height = 6.3,
       width = 9.7, units = "in")

data_zoom_freq         <- dplyr::filter(data_zoom, method %in% c("freq_q", "freq_t", "freq_z"))
data_zoom_freq$bayes   <- "Frequentist"
data_zoom_freq$method2 <- NA
data_zoom_freq$method2[which(data_zoom_freq$method == "freq_q" & data_zoom_freq$q == 0.5)] <- "Quantile of the pooled sample variance (q = 0.5)"
data_zoom_freq$method2[which(data_zoom_freq$method == "freq_q" & data_zoom_freq$q == 0.9)] <- "Quantile of the pooled sample variance (q = 0.9)"
data_zoom_freq$method2[which(data_zoom_freq$method == "freq_t")] <- "Pooled sample variance with t inflation"
data_zoom_freq$method2[which(data_zoom_freq$method == "freq_z")] <- "Pooled sample variance"


data_zoom_bayes         <- dplyr::filter(data_zoom, !(method %in% c("freq_q", "freq_t", "freq_z")))
data_zoom_bayes$bayes   <- "Bayesian-Frequentist"
data_zoom_bayes$method2 <- NA
data_zoom_bayes$method2[which(data_zoom_bayes$method == "bayes_mean")] <- "Posterior mean"
data_zoom_bayes$method2[which(data_zoom_bayes$method == "bayes_mode")] <- "Posterior mode"
data_zoom_bayes$method2[which(data_zoom_bayes$method == "bayes_post")] <- "Posterior distribution"
data_zoom_bayes$method2[which(data_zoom_bayes$method == "bayes_q" & data_zoom_bayes$q == 0.5)] <- "Posterior quantile (q = 0.5)"
data_zoom_bayes$method2[which(data_zoom_bayes$method == "bayes_q" & data_zoom_bayes$q == 0.9)] <- "Posterior quantile (q = 0.9)"
data_zoom_bayes$ab <- paste("italic(a) == ", data_zoom_bayes$a, "~~italic(b) == ", data_zoom_bayes$b,
                            sep = "")

data_zoom_freq_1    <- data_zoom_freq
data_zoom_freq_1$ab <- "italic(a) == 1.01~~italic(b) == 0.01"
data_zoom_freq_2    <- data_zoom_freq
data_zoom_freq_2$ab <- "italic(a) == 2~~italic(b) == 1"
data_zoom_freq_3    <- data_zoom_freq
data_zoom_freq_3$ab <- "italic(a) == 6~~italic(b) == 5"
data_zoom_freq_4    <- data_zoom_freq
data_zoom_freq_4$ab <- "italic(a) == 11~~italic(b) == 10"
data_zoom_freq_5    <- data_zoom_freq
data_zoom_freq_5$ab <- "italic(a) == 21~~italic(b) == 20"
data_zoom_freq_6    <- data_zoom_freq
data_zoom_freq_6$ab <- "italic(a) == 31~~italic(b) == 30"
data_zoom_freq_7    <- data_zoom_freq
data_zoom_freq_7$ab <- "italic(a) == 41~~italic(b) == 40"

data_zoom <- rbind(data_zoom_bayes,
                   data_zoom_freq_1,
                   data_zoom_freq_2,
                   data_zoom_freq_3,
                   data_zoom_freq_4,
                   data_zoom_freq_5,
                   data_zoom_freq_6,
                   data_zoom_freq_7)

data_zoom$sigma2 <- data_zoom$sigma^2
data_zoom        <- tidyr::gather(data_zoom, key = "key", value = "value", Linear:Mixed)
data_zoom$key    <- factor(data_zoom$key, levels = c("Linear", "Quadratic", "Mixed"))
data_zoom$ab     <- factor(data_zoom$ab,
                           levels = c("italic(a) == 1.01~~italic(b) == 0.01",
                                      "italic(a) == 2~~italic(b) == 1",
                                      "italic(a) == 6~~italic(b) == 5",
                                      "italic(a) == 11~~italic(b) == 10",
                                      "italic(a) == 21~~italic(b) == 20",
                                      "italic(a) == 31~~italic(b) == 30",
                                      "italic(a) == 41~~italic(b) == 40"))
data_zoom$method2 <- factor(data_zoom$method2,
                            levels = c("Pooled sample variance",
                                       "Pooled sample variance with t inflation",
                                       "Quantile of the pooled sample variance (q = 0.5)",
                                       "Quantile of the pooled sample variance (q = 0.9)",
                                       "Posterior mean",
                                       "Posterior mode",
                                       "Posterior quantile (q = 0.5)",
                                       "Posterior quantile (q = 0.9)",
                                       "Posterior distribution"))
figure2 <- ggplot(dplyr::filter(data_zoom, key %in% c("Quadratic", "Mixed") &
                                     ab %in% c("italic(a) == 2~~italic(b) == 1",
                                               "italic(a) == 6~~italic(b) == 5",
                                               "italic(a) == 11~~italic(b) == 10",
                                               "italic(a) == 21~~italic(b) == 20",
                                               "italic(a) == 31~~italic(b) == 30")),
                     aes(sigma2, as.numeric(value), colour = method2,
                         linetype = bayes)) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  ph2rand:::theme_ph2rand() + xlab(expression(sigma^2)) + ylab("Loss") +
  facet_grid(key~ab, labeller=label_parsed,
             scales = "free_y") +
  ggplot2::scale_colour_manual(values = as.vector(khroma::colour("muted")(9))) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.text=element_text(size=5))
ggsave("figure2.pdf", figure2, device = "pdf", height = 6.3,
       width = 9.7, units = "in")
supfigure3 <- ggplot(data_zoom, aes(sigma2, as.numeric(value), colour = method2,
                                    linetype = bayes)) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() +
  ph2rand:::theme_ph2rand() + xlab(expression(sigma^2)) + ylab("Loss") +
  facet_grid(key~ab, labeller=label_parsed,
             scales = "free_y") +
  ggplot2::scale_colour_manual(values = as.vector(khroma::colour("muted")(9))) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.text=element_text(size=5), axis.text.x = element_text(angle=45))
ggsave("supfigure3.pdf", supfigure3, device = "pdf", height = 6.3,
       width = 9.7, units = "in")

##### Figure 3 #################################################################

designs                 <- des(nEP    = 70,
                               alpha  = 0.025,
                               beta   = 0.1,
                               delta  = 0.35,
                               sigma  = sqrt(c(0.5, 0.75, 1, 1.25, 1.5, 1.75,
                                               2)),
                               method = c("bayes_q", "freq_z"),
                               q      = 0.5,
                               ab     = matrix(c(11, 10), 1, 2))
designs                 <- update_designs(designs)
int_data                <- int_internal(designs, 1000L, T)
int_data$dist_n$sigma2  <- factor(rep(c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2),
                                      each = 2000))
int_data_fig            <- rbind(int_data$dist_n)
int_data_fig$sigma2     <- rep(int_data$opchar$sigma^2, each = 1000)
int_data_fig$sigma2     <- paste("sigma^2 ==", int_data_fig$sigma2)
int_data_fig$Design     <- rep(rep(c("Posterior quantile (q = 0.5)",
                                     "Pooled sample variance"), each = 1000), 7)
int_data_fig            <- dplyr::filter(int_data_fig, n < 1250)
figure3                 <-
  ggplot2::ggplot(int_data_fig,
                  ggplot2::aes(n, pdf, colour = Design)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~sigma2, labeller = ggplot2::label_parsed) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 0.5"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(0.5))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 0.75"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(0.75))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 1"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(1))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 1.25"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(1.25))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 1.5"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(1.5))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 1.75"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(1.75))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::geom_vline(data = dplyr::filter(int_data_fig, sigma2 == "sigma^2 == 2"),
                      ggplot2::aes(xintercept = dplyr::filter(int_data2$opchar, sigma == sqrt(2))$`Optimal n`[1]),
                      linetype = 2) +
  ggplot2::scale_colour_manual(values = as.vector(khroma::colour("muted")(2))) +
  ggplot2::xlab(expression(italic(n))) +
  ggplot2::ylab(expression(paste(italic(f)[italic(N)], "(", italic(n), ")",
                                 sep = ""))) +
  theme_externalpilot() +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
  ggplot2::theme(legend.title = ggplot2::element_blank()) +
  ggplot2::theme(legend.text = ggtext::element_markdown())
ggsave("figure3.pdf", figure3, device = "pdf", height = 5, width = 6.3,
       units = "in")
