################################################################################
# File:        psa_response.R                                                  #
# Description: Replicates the analysis reported in "Grayling MJ, McMenamin M,  #
#              Chandler R, Heer R, Wason JMS (2021) Improving the power of PSA #
#              response analyses"                                              #
################################################################################

##### Required packages ########################################################

library(boot)
library(car)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(MASS)
library(patchwork)
library(readr)
library(tibble)

##### Required functions #######################################################

clopper_pearson <- function(s, m, alpha = 0.05) {
  if (s == 0) {
    Clow <- 0
    Cupp <- 1 - (alpha/2)^(1/m)
  } else if (s == m) {
    Clow <- (alpha/2)^(1/m)
    Cupp <- 1
  } else {
    Clow <- stats::qbeta(alpha/2, s, m - s + 1)
    Cupp <- stats::qbeta(1 - alpha/2, s + 1, m - s)
  }
  c(Clow, Cupp)
}

cp_root         <- function(m_new, p, aug_width) {
  cp <- 100*clopper_pearson(p*m_new/100, m_new)
  (cp[2] - cp[1]) - aug_width
}

deltamethod     <- function(lm1, d) {
  prob            <- meanprob(lm1,d)
  templm1         <- lm1
  templm1$coef[1] <- templm1$coef[1] + 0.0001
  derivative      <- (meanprob(templm1, d) - prob)/0.0001
  var.meanprob    <- (summary(lm1)$coef[2]^2)*derivative^2
  inv.logit(c(prob, prob - 1.96*sqrt(var.meanprob),
              prob + 1.96*sqrt(var.meanprob)))
}

meanprob        <- function(lm1, d) {
  mean  <- lm1$coef[1]
  sigma <- summary(lm1)$sigma
  prob  <- 1 - pnorm(d, mean, sigma)
  logit(mean(prob))
}

plot_theme      <- function(base_size = 11, base_family = "") {
  ggplot2::theme_grey(base_family = base_family,
                      base_size   = base_size) +
    ggplot2::theme(axis.ticks       = ggplot2::element_line(colour = "grey70",
                                                            size   = 0.25),
                   complete         = T,
                   legend.key       = ggplot2::element_rect(fill   = "white",
                                                            colour = NA),
                   legend.position  = "bottom",
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

##### Analysis #################################################################

# Make sure to have the data.csv and included.csv files in the same working
# directory
screened          <- readr::read_csv("/Users/michaelgrayling/Downloads/PSA response - Standardised (1).csv")
included          <- dplyr::filter(screened, `Final inclusion?` == "Yes")
included_mcrpc    <- dplyr::filter(included, `Type of PC` == "mCRPC")
data              <- readr::read_csv("data.csv")
data$Arm          <- as.factor(data$Arm)
data$Article      <- as.factor(data$Article)
# To change to all data, rather than just mcrpc, comment out the following row
included_analysis <- included_mcrpc
# and replace it with the following
#included_analysis <- included
results           <- matrix(NA, nrow = nrow(included_analysis), ncol = 9)
sample_sizes      <- numeric(nrow(included_analysis))
for (i in 1:nrow(included_analysis)) {
  data_i          <-
    dplyr::filter(data, Article == included_analysis$Filename[i] &
                                  Arm == included_analysis$Arm[i])
  sample_sizes[i] <- nrow(data_i)
  psadata         <-
    list(psa = -data_i$Outcome,
         d   =
           -as.numeric(included_analysis$`Dichotomisation threshold(s) (%)`[i]))
  truncationpoint <- 0.01
  y               <- (100 - psadata$psa)/100
  d               <- (100 - psadata$d)/100
  y               <- replace(y, which(y < truncationpoint), truncationpoint)
  lm_boxcox       <- lm(y1~1, data.frame(y1 = y))
  lambda          <- boxCox(lm_boxcox)$x[which.max(boxCox(lm_boxcox)$y)]
  y               <- (y^lambda - 1)/lambda
  d               <- (d^lambda - 1)/lambda
  lm_augbin       <- lm(y~1)
  S               <- ifelse(y < d, 1, 0)
  results_i       <- c(binconf(x = sum(S), n = length(S)),
                       deltamethod(lm_augbin, d))
  results_i[4:6]  <- 1 - results_i[4:6]
  results_i[5:6]  <- results_i[6:5]
  results[i, ]    <-
    100*c(results_i[1:3], results_i[3] - results_i[2], results_i[4:6],
          results_i[6] - results_i[5],
          ((results_i[3] - results_i[2]) - (results_i[6] - results_i[5]))/
            (results_i[3] - results_i[2]))
  print(i)
}

included_analysis$`hat(p) (Standard)`  <- results[, 1]
included_analysis$`LCI(p) (Standard)`  <- results[, 2]
included_analysis$`UCI(p) (Standard)`  <- results[, 3]
included_analysis$`l(p) (Standard)`    <- results[, 4]
included_analysis$`hat(p) (Augmented)` <- results[, 5]
included_analysis$`LCI(p) (Augmented)` <- results[, 6]
included_analysis$`UCI(p) (Augmented)` <- results[, 7]
included_analysis$`l(p) (Augmented)`   <- results[, 8]
included_analysis$`Efficiency gain`    <- results[, 9]

p1 <- ggplot(included_analysis, aes(`hat(p) (Standard)`, `hat(p) (Augmented)`,
                           colour = `hat(p) (Standard)`)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "gray75") +
  geom_point(size = 0.5) +
  coord_fixed() +
  plot_theme() +
  xlim(0, 100) +
  ylim(0, 100) +
  xlab("Standard point estimate (%)") +
  ylab("Augmented point estimate (%)") +
  scale_colour_viridis_c(limits = c(0,100))+
  theme(legend.position = "bottom") +
  labs(colour = "Standard point\nestimate (%)")

p2 <- ggplot(included_analysis, aes(`l(p) (Standard)`, `l(p) (Augmented)`,
                           colour = `hat(p) (Standard)`)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "gray75") +
  geom_point(size = 0.5) +
  coord_fixed() +
  plot_theme() +
  xlim(0, 65) +
  ylim(0, 65) +
  xlab("Standard confidence\ninterval length (%)") +
  ylab("Augmented confidence\ninterval length (%)") +
  scale_colour_viridis_c(limits = c(0, 100)) +
  theme(legend.position = "bottom") +
  labs(colour = "Standard point\nestimate (%)")

p3 <- ggplot(included_analysis, aes(1, `Efficiency gain`)) +
  geom_boxplot(colour = "gray75", alpha = 0.05, outlier.shape = NA) +
  geom_jitter(aes(colour = `hat(p) (Standard)`), size = 0.5) +
  plot_theme() +
  ylim(0, 100) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Efficiency gain: Confidence interval width reduction (%)") +
  scale_colour_viridis_c(limits = c(0, 100)) +
  theme(legend.position = "bottom") +
  labs(colour = "Standard point\nestimate (%)")

mod_n      <- numeric(nrow(included_analysis))
for (i in 1:nrow(included_analysis)) {
  mod_n[i] <- stats::uniroot(cp_root, c(1e-6, 1e6),
                             p = included_analysis$`hat(p) (Standard)`[i],
                             aug_width =
                               included_analysis$`l(p) (Augmented)`[i])$root
  print(i)
}
results    <- cbind(results, mod_n)

included_analysis$`Augmented implied n` <- mod_n

included_analysis$`Augmented implied %` <-
  100*(included_analysis$`Augmented implied n`/
         as.numeric(included_analysis$`Number of patients assumed in computing point estimate`) - 1)

p4 <- ggplot(included_analysis, aes(1, `Augmented implied %`)) +
  geom_boxplot(colour = "gray75", alpha = 0.05, outlier.shape = NA) +
  geom_jitter(aes(colour = `hat(p) (Standard)`), size = 0.5) +
  plot_theme() +
  coord_cartesian(ylim = c(0, 500)) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Efficiency gain: Implied increase in sample size (%)") +
  scale_colour_viridis_c(limits = c(0, 100)) +
  theme(legend.position = "bottom") +
  labs(colour = "Standard point\nestimate (%)")
((p1 / p2) | p3 | p4) +
  plot_layout(widths = c(2.325, 1, 1), height = c(1, 1, 2, 2),
              guides = 'collect') +
  plot_annotation(tag_levels = "A") & theme(legend.position = 'bottom')
ggsave("figure2.pdf", device = "pdf", width = 6.25, height = 6.75, units = "in")

sensitivity            <- NULL
for (i in 1:nrow(included)) {
  sensitivity_i        <- matrix(0, 100, 9)
  data_i               <- dplyr::filter(data, Article == included$Filename[i] &
                                          Arm == included$Arm[i])
  sample_sizes[i]      <- nrow(data_i)
  data_i_og            <-
    list(psa = -data_i$Outcome,
         d   = -as.numeric(included$`Dichotomisation threshold(s) (%)`[i]))
  for (j in 1:100) {
    data_ij            <- data_i_og
    data_ij$psa        <- data_ij$psa + runif(length(data_ij$psa), -1, 1)

    psadata            <- data_ij
    truncationpoint    <- 0.01
    y                  <- (100 - psadata$psa)/100
    d                  <- (100  -psadata$d)/100
    y                  <- replace(y, which(y < truncationpoint), truncationpoint)
    lm_boxcox          <- lm(y1~1, data.frame(y1 = y))
    lambda             <- boxCox(lm_boxcox)$x[which.max(boxCox(lm_boxcox)$y)]
    y                  <- (y^lambda - 1)/lambda
    d                  <- (d^lambda - 1)/lambda
    lm_augbin          <- lm(y~1)
    S                  <- ifelse(y < d, 1, 0)
    results_ij         <- c(binconf(x = sum(S), n = length(S)),
                            deltamethod(lm_augbin, d))
    sensitivity_i[j, ] <- 100*c(results_ij[1:3], results_ij[3] - results_ij[2],
                                results_ij[4:6], results_ij[6] - results_ij[5],
                                ((results_ij[3] - results_ij[2]) -
                                   (results_ij[6] - results_ij[5]))/
                                  (results_ij[3] - results_ij[2]))
    message("i = ", i, ", j = ", j)
  }
  sensitivity          <- rbind(sensitivity, sensitivity_i)
}
sensitivity_mod_n      <- numeric(nrow(sensitivity))
for (i in 1:nrow(sensitivity)) {
  sensitivity_mod_n[i] <- stats::uniroot(cp_root, c(1e-6, 1e6),
                                         p = sensitivity[i, 1],
                                         aug_width = sensitivity[i, 8])$root
  print(i)
}
sensitivity            <-
  cbind(sensitivity,
        100*(sensitivity_mod_n/rep(as.numeric(included$`Number of patients assumed in computing point estimate`),
                                   each = 100) - 1))

sensitivity            <- tibble::as_tibble(sensitivity)
sensitivity$Trial      <- factor(rep(1:121, each = 100))
labels_sens            <- character(121)
labels_sens[c(1, 25, 50, 75, 100, 121)] <- c("1", "25", "50", "75", "100",
                                             "121")
p1 <- ggplot(sensitivity, aes(Trial, V9)) +
  geom_boxplot(outlier.size = 0.5) +
  plot_theme() + xlab("Trial") +
  ylab("Efficiency gain: Confidence interval\nwidth reduction (%)") +
  scale_x_discrete(labels = labels_sens)
p2 <- ggplot(sensitivity, aes(Trial, sensitivity_mod_n)) +
  geom_boxplot(outlier.size = 0.5) +
  plot_theme() + xlab("Trial") +
  ylab("Efficiency gain: Implied increase\nin sample size (%)") +
  scale_x_discrete(labels = labels_sens) + ylim(0, 500)
(p1 / p2) + plot_annotation(tag_levels = "A")
ggsave("sfigure1.pdf", device = "pdf", width = 6.25, height = 9, units = "in")

medians_ci      <- sd_ci <- medians_n <- sd_n <-numeric(nrow(included))
for (i in 1:nrow(included)) {
  medians_ci[i] <- quantile(sensitivity[(1 + 100*(i - 1)):(100*i), 9])
  medians_n[i]  <- quantile(sensitivity[(1 + 100*(i - 1)):(100*i), 10])
  sd_ci[i]      <- sd(sensitivity[(1 + 100*(i - 1)):(100*i), 9])
  sd_n[i]       <- sd(sensitivity[(1 + 100*(i - 1)):(100*i), 10])
  print(i)
}
included$`Number of bars clipped` <-
  as.numeric(included$`Number of bars clipped`)
included$`Clip point (%)`[22]     <- 100
results_clip                      <- matrix(NA, nrow = nrow(included),
                                            ncol = 11)
multipliers                       <- seq(1.5, 10, 0.5)
medians_clip                      <- lqrs_clip <- uqrs_clip <-
  medians_clip_n <- lqrs_clip_n <- uqrs_clip_n <- numeric(length(multipliers))
for (z in 1:length(multipliers)) {
  for (i in 1:nrow(included)) {
    if (included$`Number of bars clipped`[i] > 0) {
      data_i          <- dplyr::filter(data, Article == included$Filename[i] &
                                         Arm == included$Arm[i])
      data_i          <-
        list(psa = -data_i$Outcome,
             d   = -as.numeric(included$`Dichotomisation threshold(s) (%)`[i]))
      data_i$psa[data_i$psa == -included$`Clip point (%)`[i]] <-
        multipliers[z]*data_i$psa[data_i$psa == -included$`Clip point (%)`[i]]
      psadata         <- data_i
      truncationpoint <- 0.01
      y               <- (100 - psadata$psa)/100
      d               <- (100 - psadata$d)/100
      y               <- replace(y, which(y < truncationpoint), truncationpoint)
      lm_boxcox       <- lm(y1~1, data.frame(y1 = y))
      lambda          <- boxCox(lm_boxcox)$x[which.max(boxCox(lm_boxcox)$y)]
      y               <- (y^lambda - 1)/lambda
      d               <- (d^lambda - 1)/lambda
      lm_augbin       <- lm(y~1)
      S               <- ifelse(y < d, 1, 0)
      results_i       <- c(binconf(x = sum(S), n = length(S)),
                           deltamethod(lm_augbin, d))
      results_clip[i, 1:9] <- 100*c(results_i[1:3], results_i[3] - results_i[2],
                                    results_i[4:6], results_i[6] - results_i[5],
                                    ((results_i[3] - results_i[2]) -
                                       (results_i[6] - results_i[5]))/
                                      (results_i[3] - results_i[2]))

      results_clip[i, 10]  <- stats::uniroot(cp_root, c(1e-6, 1e6),
                                             p = results_clip[i, 1],
                                             aug_width =
                                               results_clip[i, 8])$root
    }
    print(i)
  }
  results_clip[, 11] <-
    100*(results_clip[, 10]/
           as.numeric(included$`Number of patients assumed in computing point estimate`) - 1)
  medians_clip[z]    <- quantile(results_clip[!is.na(results_clip[, 9]), 9],
                                 0.5)
  lqrs_clip[z]       <- quantile(results_clip[!is.na(results_clip[, 9]), 9],
                                 0.25)
  uqrs_clip[z]       <- quantile(results_clip[!is.na(results_clip[, 9]), 9],
                                 0.75)
  medians_clip_n[z]  <- quantile(results_clip[!is.na(results_clip[, 9]), 11],
                                 0.5)
  lqrs_clip_n[z]     <- quantile(results_clip[!is.na(results_clip[, 9]), 11],
                                 0.25)
  uqrs_clip_n[z]     <- quantile(results_clip[!is.na(results_clip[, 9]), 11],
                                 0.75)
  message("***z = ", z, "***")
}
data_clip <-
  tibble::tibble(multiplier = c(1, multipliers),
                 medians    = c(quantile(results[!is.na(results_clip[, 9]), 9],
                                         0.5), medians_clip),
                 lqr        = c(quantile(results[!is.na(results_clip[, 9]), 9],
                                         0.25), lqrs_clip),
                 uqr        = c(quantile(results[!is.na(results_clip[, 9]), 9],
                                         0.75), uqrs_clip),
                 medians_n  = c(quantile(results[!is.na(results_clip[, 9]), 11],
                                         0.5), medians_clip_n),
                 lqr_n      = c(quantile(results[!is.na(results_clip[, 9]), 11],
                                         0.25), lqrs_clip_n),
                 uqr_n      = c(quantile(results[!is.na(results_clip[, 9]), 11],
                                         0.75), uqrs_clip_n))
p1 <- ggplot(data_clip, aes(multiplier, medians)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = lqr, ymax = uqr), width = 0.1) +
  xlab("Clipped bars multiplied by") +
  ylab("Efficiency gain: Confidence interval\nwidth reduction (%), Median [IQR]") +
  plot_theme() + ylim(0, 50)
p2 <- ggplot(data_clip, aes(multiplier, medians_n)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = lqr_n, ymax = uqr_n), width = 0.1) +
  xlab("Clipped bars multiplied by") +
  ylab("Efficiency gain: Implied increase\nin sample size (%), Median [IQR]") +
  plot_theme()
(p1 / p2) + plot_annotation(tag_levels = "A")
ggsave("sfigure2.pdf", device = "pdf", width = 6.25, height = 9, units = "in")
