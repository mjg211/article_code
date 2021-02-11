################################################################################
# R code to reproduce the results from "Response adaptive intervention         #
# allocation in stepped-wedge cluster randomised trials" by Grayling MJ, Wason #
# JMS, and Villar SS                                                           #
################################################################################

library(ggplot2)
library(iterpc)
library(readr)
library(Rfast)
library(snowfall)

################################################################################
##### Simulate #################################################################
################################################################################

components_ra_swcrt  <- function(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                 sigma_s2, sigma_e2, replicates, prefix,
                                 parallel, cpus) {
  C                                         <- length(S)
  Z                                         <- matrix(0, nrow = C*Ti,
                                                      ncol = 2*C + C*Ti)
  for (i in 1:C) {
    for (j in 1:Ti) {
      Z[(i - 1)*Ti + j, c(i, C + (i - 1)*Ti + j, C + C*Ti + i)] <- 1
    }
  }
  G                                         <-
    diag(c(rep(sigma_c2, C), rep(sigma_pi2, C*Ti), rep(sigma_s2/m, C)))
  R                                         <- diag(rep(sigma_e2/m, C*Ti))
  Sigma                                     <- Z%*%G%*%t(Z) + R
  chol_Sigma                                <- chol(Sigma)
  responses                                 <- matrix(0, C*Ti, replicates)
  for (i in 1:replicates) {
    responses[, i]                          <-
      as.numeric(stats::rnorm(C*Ti)%*%chol_Sigma)
  }
  vecs_t                                    <- list()
  for (t in c(set_T, Ti)) {
    vecs_t_1                                <- vecs_t[[t]] <- 1:t
    for (c in 2:C) {
      vecs_t[[t]]                           <- c(vecs_t[[t]],
                                                 Ti*(c - 1) + vecs_t_1)
    }
  }
  X_id                                      <- gls_factors <- list()
  for (t in c(set_T, Ti)) {
    S_combs                                 <-
      iterpc::getall(iterpc::iterpc(n = t + 1, r = C, labels = 1:(t + 1),
                                    ordered = F, replace = T))
    keep                                    <- rep(T, nrow(S_combs))
    init_S                                  <- S
    init_S[which(init_S > set_T[1])]        <- set_T[1] + 1
    for (i in 1:nrow(S_combs)) {
      if (length(unique(S_combs[i, ])) < 2) {
        keep[i]                             <- F
      } else {
        int_S                               <- S_combs[i, ]
        int_S[which(int_S > set_T[1])]      <- set_T[1] + 1
        if (!all(int_S == init_S)) {
          keep[i]                           <- F
        }
      }
    }
    S_combs                                 <- S_combs[keep, , drop = F]
    X_id[[t]]                               <- S_combs
    gls_factors[[t]]                        <- matrix(0, sum(keep), C*t)
    D_t                                     <- matrix(0L, C*t, 1 + t)
    D_t[, 1]                                <- 1L
    for (i in 1:C) {
      for (j in 1:t) {
        if (j > 1) {
          D_t[(i - 1)*t + j, j]             <- 1L
        }
      }
    }
    Z                                       <- matrix(0, nrow = C*t,
                                                      ncol = 2*C + C*t)
    for (i in 1:C) {
      for (j in 1:t) {
        Z[(i - 1)*t + j, c(i, C + (i - 1)*t + j, C + C*t + i)] <- 1L
      }
    }
    G                                       <-
      diag(c(rep(sigma_c2, C), rep(sigma_pi2, C*t), rep(sigma_s2/m, C)))
    R                                       <- diag(rep(sigma_e2/m, C*t))
    Sigma                                   <- Z%*%G%*%t(Z) + R
    Sigma_inv_t                             <- MASS::ginv(Sigma)
    for (z in 1:nrow(S_combs)) {
      new_X                                 <- matrix(0, C, t)
      for (j in 1:C) {
        if (S_combs[z, j] <= t) {
          new_X[j, S_combs[z, j]:t]         <- 1
        }
      }
      D_t[, t + 1]                          <- 0L
      for (i in 1:C) {
        for (j in 1:t) {
          if (new_X[i, j] == 1) {
            D_t[(i - 1)*t + j, t + 1]       <- 1L
          }
        }
      }
      gls_factors[[t]][z, ]                 <-
        (MASS::ginv(t(D_t)%*%Sigma_inv_t%*%D_t)%*%t(D_t)%*%Sigma_inv_t)[t + 1, ]
      message("t = ", t, ", des% = ", round(100*z/nrow(S_combs), 3))
    }
  }
  gls_factors_t1                            <-
    Rfast::colsums(gls_factors[[set_T[1]]]%*%responses[vecs_t[[set_T[1]]], ])
  sfInit(parallel = T, cpus = 7)
  sfLibrary(Rfast)
  sfLibrary(utils)
  sfExport("gls_factors", "vecs_t", "responses", "set_T", "Ti", "prefix")
  results                                   <- sfLapply(1:replicates,
                                                        wrapper_gls)
  sfStop()
  remove(results)
  remove(responses)
  nrow_poss_combs                           <- poss_combs <- names_poss_combs <-
                                               list()
  for (t in set_T) {
    nrow_poss_combs[[t]]                    <- poss_combs[[t]]       <-
                                               names_poss_combs[[t]] <- list()
    for (c in 1:(C - 1)) {
      poss_combs[[t]][[c]]                  <-
        iterpc::getall(iterpc::iterpc(n = Ti + 1 - t, r = c,
                                      labels  = (t + 1):(Ti + 1), ordered = F,
                                      replace = T))
      nrow_poss_combs[[t]][[c]]             <- nrow(poss_combs[[t]][[c]])
      names_poss_combs[[t]][[c]]            <-
        character(nrow_poss_combs[[t]][[c]])
      for (i in 1:nrow_poss_combs[[t]][[c]]) {
        names_poss_combs[[t]][[c]][i]       <- paste(poss_combs[[t]][[c]][i, ],
                                                     collapse = ",")
      }
    }
  }
  sigma2                                    <-
    sigma_e2 + sigma_c2 + sigma_pi2 + sigma_s2
  alpha0                                    <- (sigma_c2 + sigma_pi2)/sigma2
  alpha1                                    <- sigma_c2/sigma2
  alpha2                                    <- (sigma_c2 + sigma_s2)/sigma2
  psi                                       <-
    1 + (m - 1)*alpha0 - (m - 1)*alpha1 - alpha2
  xi                                        <- (m - 1)*alpha1 + alpha2
  int_S                                     <- S
  int_S[which(int_S > set_T[1])]            <- set_T[1] + 1L
  factor                                    <- set_T[1] + 1 - int_S
  U                                         <- sum(factor)
  V                                         <- sum(factor^2)
  W                                         <-
    sum(sapply(1:set_T[1], function(t) sum(int_S <= t))^2)
  gamma                                     <- psi + set_T[1]*xi
  I1                                        <-
    ((C*U - W)*gamma + (U^2 - C*V)*xi)*m/(C*sigma2*gamma*psi)
  S_matrix                                  <-
    iterpc::getall(iterpc::iterpc(n = Ti - set_T[1] + 1,
                                  r = C - length(which(S <= set_T[1])),
                                  labels  = (set_T[1] + 1):(Ti + 1),
                                  ordered = F, replace = T))
  S_matrix                                  <-
    cbind(matrix(S[which(S <= set_T[1])], nrow(S_matrix),
                 length(which(S <= set_T[1])), byrow = T), S_matrix)
  I_vector                                  <- numeric(nrow(S_matrix))
  I_names                                   <- character(nrow(S_matrix))
  factor                                    <- Ti + 1 - S_matrix
  U                                         <- rowSums(factor)
  V                                         <- rowSums(factor^2)
  gamma                                     <- psi + Ti*xi
  for (i in 1:nrow(S_matrix)) {
    W                                       <-
      sum(sapply(1:Ti, function(t) sum(S_matrix[i, ] <= t))^2)
    I_vector[i]                             <-
      ((C*U[i] - W)*gamma + (U[i]^2 - C*V[i])*xi)*m/(C*sigma2*gamma*psi)
    I_names[i]                              <- paste(S_matrix[i, ],
                                                     collapse = ",")
  }
  names(I_vector)                           <- I_names
  I_t1                                      <-
    I_vector[paste(paste(S[S <= set_T[1]], collapse = ","),
                   names_poss_combs[[set_T[1]]][[sum(S > set_T[1])]],
                   sep = ",")]
  which_0_t1                                <- (S > set_T[1])
  num_0_t1                                  <- sum(which_0_t1)
  permissible_S_t1                          <-
    poss_combs[[set_T[1]]][[num_0_t1]]
  indices_t1                                <-
    Rfast::rowsums(Ti + 1 - permissible_S_t1) + 1
  list(alpha            = alpha,
       alpha0           = alpha0,
       alpha1           = alpha1,
       alpha2           = alpha2,
       C                = C,
       clusters         = rep(rep(1:C, each = m), Ti),
       cpus             = cpus,
       gls_factors_t1   = gls_factors_t1,
       I_vector         = I_vector,
       indices_t1       = indices_t1,
       m                = m,
       mCt              = m*C*(1:Ti),
       names_poss_combs = names_poss_combs,
       nrow_poss_combs  = nrow_poss_combs,
       num_0_t1         = num_0_t1,
       parallel         = parallel,
       permissible_S_t1 = permissible_S_t1,
       poss_combs       = poss_combs,
       prefix           = prefix,
       psi              = psi,
       rem_t1           = Ti - set_T[1],
       replicates       = replicates,
       S                = S,
       scaled_info_t1   = I_t1/max(I_t1),
       seq_C            = 1:C,
       seq_Ti           = 1:Ti,
       set_T            = set_T,
       sigma_c2         = sigma_c2,
       sigma_e2         = sigma_e2,
       sigma_pi2        = sigma_pi2,
       sigma_s2         = sigma_s2,
       sigma2           = sigma2,
       sqrt_I1          = sqrt(I1),
       Ti               = Ti,
       Tip1             = Ti + 1,
       which_0_t1       = which_0_t1,
       X_id             = X_id,
       xi               = xi,
       zero_Cn          = numeric(C*m))
}

sim_ra_swcrt         <- function(sc, prefix, design, components) {
  sfInit(parallel = components$parallel, cpus = components$cpus)
  sfExport("design", "components")
  sfLibrary(readr)
  sfLibrary(Rfast)
  sfLibrary(stats)
  results <- sfLapply(1:components$replicates, wrapper_sim_ra_swcrt)
  sfStop()
  data    <- matrix(unlist(results), components$replicates,
                    (5 + length(components$S))*length(design$theta), byrow = T)
  write.csv(data, paste0(prefix, "_", sc, ".csv"))
}

sim_swcrt            <- function(prefix, design, components) {
  sfInit(parallel = components$parallel, cpus = components$cpus)
  sfExport("design", "components")
  sfLibrary(readr)
  sfLibrary(Rfast)
  sfLibrary(stats)
  results <- sfLapply(1:components$replicates, wrapper_sim_swcrt)
  sfStop()
  data    <- matrix(unlist(results), components$replicates,
                    4*length(design$theta), byrow = T)
  write.csv(matrix(colmeans(data), length(design$theta), 4, byrow = T),
            paste0(prefix, ".csv"))
}

wrapper_gls          <- function(i) {
  for (t in c(set_T[-1], Ti)) {
    utils::write.table(Rfast::rowsums(gls_factors[[t]]%*%
                                        responses[vecs_t[[t]], i]),
                       paste0(prefix, "_analysis_", t, "_", i, ".csv"),
                       row.names = F, col.names = F)
  }
  return(T)
}

wrapper_sim_ra_swcrt <- function(rep) {
  len                             <- 5 + length(components$S)
  obj_fn_reps                     <- numeric(len*length(design$theta))
  gls_factors                     <- list()
  for (t in c(components$set_T[-1], components$Ti)) {
    gls_factors[[t]]              <-
      read_csv(paste0(components$prefix, "_analysis_", t, "_", rep, ".csv"),
               col_names = F, col_types = "d")
  }
  for (th in 1:length(design$theta)) {
    curr_S                        <- components$S
    Z_t                           <-
      (components$gls_factors_t1[rep] + design$theta[th])*components$sqrt_I1
    eff                           <-
      dbinom(0:(components$num_0_t1*components$rem_t1),
             components$num_0_t1*components$rem_t1,
             pnorm((Z_t - design$eta)/
                     (design$gamma*(1 - t/components$Ti))))[components$indices_t1]
    curr_S                        <-
      c(curr_S[!components$which_0_t1],
        components$permissible_S_t1[which.max(
          design$w*components$scaled_info_t1 + (1 - design$w)*eff/max(eff)), ])
    count                         <- 1L
    if (length(components$set_T) > 1) {
      for (t in components$set_T[-1]) {
        if (curr_S[components$C] > t) {
          rem_t                   <- components$Ti - t
          count                   <- count + 1L
          tp1                     <- t + 1L
          int_S                   <- curr_S
          int_S[which(int_S > t)] <- tp1
          curr_X_id               <-
            which(apply(components$X_id[[t]], 1,
                        function(x) { all(x == int_S) }))
          factor                  <- tp1 - int_S
          U                       <- sum(factor)
          V                       <- sum(factor^2)
          W                       <-
            sum(sapply(1:t, function(t) sum(int_S <= t))^2)
          gamma                   <- components$psi + t*components$xi
          I                       <-
            ((components$C*U - W)*gamma + (U^2 - components$C*V)*components$xi)*
            components$m/(components$C*components$sigma2*gamma*components$psi)
          Z_t                     <-
            (as.numeric(gls_factors[[t]][curr_X_id, 1]) +
               design$theta[th])*sqrt(I)
          which_0                   <- (curr_S > t)
          num_0                     <- sum(which_0)
          permissible_S             <- components$poss_combs[[t]][[num_0]]
          info                    <-
            components$I_vector[paste0(
              paste(curr_S[!which_0], collapse = ","), ",",
              components$names_poss_combs[[t]][[num_0]])]
          if (components$g == 1) {
            if (Z_t > components$nu) {
              rep_tp1               <- t(rep(tp1, num_0))
              eff                   <- apply(permissible_S[, which_0], 1,
                                             function(x) all(x == rep_tp1))
            } else {
              rep_Tip1              <- t(rep(components$Tip1, num_0))
              eff                   <- apply(permissible_S[, which_0], 1,
                                             function(x) all(x == rep_Tip1))
            }
          } else {
            eff                     <-
              dbinom(0:(num_0*rem_t), num_0*rem_t,
                     pnorm((Z_t - design$eta)/
                             (design$gamma*(1 - t/components$Ti)))
                     )[rowsums(components$Tip1 - permissible_S) + 1]
          }
          curr_S                    <-
            c(curr_S[!which_0],
              permissible_S[which.max(design$w*info/max(info) +
                                        (1 - design$w)*eff/max(eff)), ])
        }
      }
    }
    curr_X_id                       <-
      which(apply(components$X_id[[components$Ti]], 1,
                  function(x) { all(x == curr_S) }))
    factor                          <- components$Tip1 - curr_S
    U                               <- sum(factor)
    V                               <- sum(factor^2)
    W                               <-
      sum(sapply(components$seq_Ti, function(t) sum(curr_S <= t))^2)
    gamma                           <-
      components$psi + components$Ti*components$xi
    I                               <-
      ((components$C*U - W)*gamma + (U^2 - components$C*V)*components$xi)*
      components$m/(components$C*components$sigma2*gamma*components$psi)
    Z_t                             <-
      (as.numeric(gls_factors[[components$Ti]][curr_X_id, 1]) +
         design$theta[th])*sqrt(I)
    obj_fn_reps[(1 + len*(th - 1)):(len*th)] <-
      c(Z_t/sqrt(I),
        Z_t,
        as.numeric(Z_t >= qnorm(1 - components$alpha)),
        sum(components$Tip1 - curr_S)/(components$C*components$Ti),
        count, curr_S)
  }
  obj_fn_reps
}

wrapper_sim_swcrt    <- function(rep) {
  obj_fn_reps                            <- numeric(4*length(design$theta))
  gls_factors                            <- list()
  gls_factors[[components$Ti]]           <-
    read_csv(paste0(components$prefix, "_analysis_", components$Ti, "_", rep,
                    ".csv"), col_names = F, col_types = "d")
  curr_S                                 <- components$S
  curr_X_id                              <-
    which(apply(components$X_id[[components$Ti]], 1,
                function(x) { all(x == curr_S) }))
  factor                                 <- components$Tip1 - curr_S
  U                                      <- sum(factor)
  V                                      <- sum(factor^2)
  W                                      <-
    sum(sapply(components$seq_Ti, function(t) sum(curr_S <= t))^2)
  gamma                                  <-
    components$psi + components$Ti*components$xi
  I                                      <-
    ((components$C*U - W)*gamma + (U^2 - components$C*V)*components$xi)*
    components$m/(components$C*components$sigma2*gamma*components$psi)
  for (th in 1:length(design$theta)) {
    Z_t                                  <-
      (as.numeric(gls_factors[[components$Ti]][curr_X_id, 1]) +
         design$theta[th])*sqrt(I)
    obj_fn_reps[(1 + 4*(th - 1)):(4*th)] <-
      c(Z_t/sqrt(I),
        Z_t,
        as.numeric(Z_t >= qnorm(1 - components$alpha)),
        sum(components$Tip1 - curr_S)/(components$C*components$Ti))
  }
  obj_fn_reps
}

##### All ######################################################################

replicates <- 100000
mu         <- 0
parallel   <- T
cpus       <- 11
scenarios  <- as.matrix(expand.grid(w     = c(1/1000, 1/4, 1/3, 1/2, 2/3, 3/4,
                                              999/1000),
                                    eta   = -1:4,
                                    gamma = c(1, 2.5, 5)))

##### TDS1 #####################################################################

X         <- rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 1, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 0, 1),
                   c(0, 0, 0, 0, 0, 0, 0, 0, 1))
C         <- nrow(X)
Ti        <- ncol(X)
m         <- 7
S         <- c(rep(2:5, each = 3), rep(6:9, each = 2))
alpha     <- 0.05
beta      <- 0.2
delta     <- 0.24
theta     <- delta*seq(-1, 2, by = 0.5)
sigma_c2  <- 1/9
sigma_e2  <- 1
sigma_pi2 <- sigma_s2 <- 0
pi        <- numeric(Ti)

# TDS1: {3}

set.seed(13)
set_T      <- 3
prefix     <- "sc1_3"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

# TDS1: {4}

set.seed(14)
set_T      <- 4
prefix     <- "sc1_4"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

# TDS1: {5}

set.seed(15)
set_T      <- 5
prefix     <- "sc1_5"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

# TDS1: {3, 6}

set.seed(136)
set_T      <- c(3, 6)
prefix     <- "sc1_36"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

##### TDS2 #####################################################################

X         <- rbind(c(0, 1, 1, 1, 1),
                   c(0, 0, 1, 1, 1),
                   c(0, 0, 0, 1, 1),
                   c(0, 0, 0, 0, 1))
C         <- nrow(X)
Ti        <- ncol(X)
m         <- 70
S         <- 2:5
alpha     <- 0.05
beta      <- 0.1
delta     <- 0.2
theta     <- delta*seq(-1, 2, by = 0.5)
sigma_c2  <- 0.02
sigma_e2  <- 0.51
sigma_pi2 <- sigma_s2 <- 0
pi        <- numeric(Ti)

# TDS2: {3}

set.seed(23)
set_T      <- 3
prefix     <- "sc2_3"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

# TDS2: {2, 3, 4}

set.seed(2234)
set_T      <- 2:4
prefix     <- "sc2_234"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

##### TDS3 #####################################################################

X         <- rbind(c(0, 1, 1, 1),
                   c(0, 1, 1, 1),
                   c(0, 1, 1, 1),
                   c(0, 1, 1, 1),
                   c(0, 0, 1, 1),
                   c(0, 0, 1, 1),
                   c(0, 0, 1, 1),
                   c(0, 0, 1, 1),
                   c(0, 0, 0, 1),
                   c(0, 0, 0, 1),
                   c(0, 0, 0, 1),
                   c(0, 0, 0, 1))
C         <- nrow(X)
Ti        <- ncol(X)
m         <- 10
S         <- rep(2:4, each = 4)
alpha     <- 0.025
beta      <- 0.11
delta     <- 2
theta     <- delta*seq(-1, 2, by = 0.5)
rho       <- 0.33
tau       <- 0.7
pi        <- 0.9
sigma     <- 5
sigma_c2  <- pi*rho*sigma^2
sigma_pi2 <- rho*sigma^2*(1 - pi)
sigma_s2  <- (1 - rho)*tau*sigma^2
sigma_e2  <- (1 - rho)*sigma^2*(1 - tau)
pi        <- numeric(Ti)

# TDS3: {2}

set.seed(32)
set_T      <- 2
prefix     <- "sc3_2"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

# TDS3: {2, 3}

set.seed(323)
set_T      <- 2:3
prefix     <- "sc3_23"
components <- components_ra_swcrt(S, Ti, m, set_T, alpha, sigma_c2, sigma_pi2,
                                  sigma_s2, sigma_e2, replicates, prefix,
                                  parallel, cpus)
for (sc in 1:nrow(scenarios)) {
  scenario <- sim_ra_swcrt(sc, prefix, list(w     = scenarios[sc, 1],
                                            eta   = scenarios[sc, 2],
                                            gamma = scenarios[sc, 3],
                                            theta = theta), components)
}

################################################################################
##### Analyse ##################################################################
################################################################################

analyse_ra_sw_data <- function(set_T, X, delta, alpha, beta, eta_range,
                               gamma_range, eta_focus, gamma_focus, w_focus,
                               prefix) {
  scenarios                          <-
    as.matrix(expand.grid(w      =
                            paste("w =", c("1/1000", "1/4", "1/3", "1/2", "2/3",
                                           "3/4", "999/1000")),
                          eta    = paste("eta ==", -1:4),
                          gamma  = paste("gamma ==", c(1, 2.5, 5))))
  nrow_scenarios                     <- nrow(scenarios)
  setwd(paste0("/Users/michaelgrayling/Documents/Work/Papers/Response adaptive",
               " intervention allocation in stepped-wedge cluster randomised ",
               "trials/RA SW-CRT/", prefix))
  C                                  <- nrow(X)
  Ti                                 <- ncol(X)
  theta                              <- delta*seq(-1, 2, 0.5)
  len_theta                          <- length(theta)
  tds_sum                            <-
    tibble::tibble(w         = rep(scenarios[, 1], each = len_theta),
                   eta       = rep(scenarios[, 2], each = len_theta),
                   gamma     = rep(scenarios[, 3], each = len_theta),
                   theta     = rep(theta, nrow_scenarios),
                   theta_hat = numeric(len_theta*nrow_scenarios),
                   power     = numeric(len_theta*nrow_scenarios),
                   prop_I    = numeric(len_theta*nrow_scenarios))
  tds_avX                            <-
    tibble::tibble(theta     = rep(theta, each = C*Ti),
                   C         = rep(rep(1:C, each = Ti), len_theta),
                   Ti        = rep(1:Ti, length(theta)/Ti),
                   Xbar_ij   = numeric(length(Ti)))
  for (i in 1:nrow(scenarios)) {
    scenario_i                       <-
      suppressMessages(readr::read_csv(paste0(prefix, "_", i, ".csv")))
    tds_sum[(1 + (i - 1)*len_theta):(i*len_theta), 5:7]      <-
      cbind(colMeans(scenario_i[, seq(from = 2, by = 5 + C,
                                      length.out = len_theta)]),
            colMeans(scenario_i[, seq(from = 4, by = 5 + C,
                                      length.out = len_theta)]),
            colMeans(scenario_i[, seq(from = 5, by = 5 + C,
                                      length.out = len_theta)]))
    if (all(scenarios[i, 1] == paste("w =", w_focus),
            scenarios[i, 2] == paste("eta ==", eta_focus),
            scenarios[i, 3] == paste("gamma ==", gamma_focus))) {
      for (j in 1:len_theta) {
        X_i                          <- matrix(0, C, Ti)
        for (k in 1:100000) {
          S                          <-
            as.numeric(scenario_i[k, (7 + (j - 1)*(5 + C)):
                                    (7 + (j - 1)*(5 + C) + C - 1)])
          for (c in 1:C) {
            if (S[c] <= Ti) {
              X_i[c, S[c]:Ti]        <- X_i[c, S[c]:Ti] + 1
            }
          }
        }
        tds_avX[(1 + (j - 1)*C*Ti):(C*Ti + (j - 1)*C*Ti), 4] <-
          as.vector(t(X_i))/100000
        message("j = ", j)
      }
    }
    remove(scenario_i)
    message(i)
  }
  tds_sum$w                          <-
    factor(tds_sum$w, levels = paste("w =", c("1/1000", "1/4", "1/3", "1/2",
                                              "2/3", "3/4", "999/1000")))
  tds_sum$eta                        <- as.factor(tds_sum$eta)
  tds_sum$gamma                      <- as.factor(tds_sum$gamma)
  plots                              <- list()
  min_X                              <- X
  min_X[which(X[, set_T[1]] == 0), ] <- 0
  min_prop                           <- sum(min_X)/(C*Ti)
  max_X                              <- X
  max_X[, (set_T[1] + 1):Ti]         <- 1
  max_prop                           <- sum(max_X)/(C*Ti)
  cbbPalette                         <-
    c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
      "#D55E00", "#CC79A7")
  plots$prop_I                       <-
    ggplot(tds_sum, aes(theta, prop_I, colour = w)) +
    facet_grid(eta~gamma, labeller = label_parsed) +
    geom_hline(yintercept = min_prop, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_hline(yintercept = sum(X)/(C*Ti), colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_hline(yintercept = max_prop, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = 0, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_line(size = 1/3) +
    geom_point(size = 1/3) +
    ph2rand:::theme_ph2rand() +
    scale_colour_manual(values = cbbPalette) +
    ylab(paste("Emp. av. prop. of cl.-per. in interv. con. (EACP)")) +
    xlab(expression(theta)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
          text        = element_text(size = 8))
  ggsave(paste0(prefix, "_prop_I.pdf"), plot = plots$prop_I, width = 5,
         height = 5, units = "in", device = "pdf")
  plots$power                        <-
    ggplot(tds_sum, aes(theta, power, colour = w)) +
    facet_grid(eta~gamma, labeller = label_parsed) +
    geom_hline(yintercept = alpha, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = 0, colour = "gray75", linetype = 2, size = 1/3) +
    geom_hline(yintercept = 1 - beta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_line(size = 1/3) +
    geom_point(size = 1/3) +
    ph2rand:::theme_ph2rand() +
    scale_colour_manual(values = cbbPalette) +
    ylab("Empirical rejection probability (ERP)") +
    xlab(expression(theta)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
          text        = element_text(size = 8))
  ggsave(paste0(prefix, "_power.pdf"), plot = plots$power, width = 5,
         height = 5, units = "in", device = "pdf")
  plots$prop_I_range                 <-
    ggplot(dplyr::filter(tds_sum, eta %in% paste("eta ==", eta_range) &
                           gamma %in% paste("gamma ==", gamma_range)),
           aes(theta, prop_I, colour = w)) +
    facet_grid(eta~gamma, labeller = label_parsed) +
    geom_hline(yintercept = min_prop, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_hline(yintercept = sum(X)/(C*Ti), colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_hline(yintercept = max_prop, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = 0, colour = "gray75", linetype = 2, size = 1/3) +
    geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_line(size = 1/3) +
    geom_point(size = 1/3) +
    ph2rand:::theme_ph2rand() +
    scale_colour_manual(values = cbbPalette) +
    ylab(paste("Emp. av. prop. of cl.-per. in interv. con. (EACP)")) +
    xlab(expression(theta)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
          text        = element_text(size = 8))
  ggsave(paste0(prefix, "_prop_I_range.pdf"), plot = plots$prop_I_range,
         width = 5, height = 5, units = "in", device = "pdf")
  plots$power_range                  <-
    ggplot(dplyr::filter(tds_sum, eta %in% paste("eta ==", eta_range) &
                           gamma %in% paste("gamma ==", gamma_range)),
           aes(theta, power, colour = w)) +
    facet_grid(eta~gamma, labeller = label_parsed) +
    geom_hline(yintercept = alpha, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_hline(yintercept = 1 - beta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_vline(xintercept = 0, colour = "gray75", linetype = 2, size = 1/3) +
    geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_line(size = 1/3) +
    geom_point(size = 1/3) +
    ph2rand:::theme_ph2rand() +
    scale_colour_manual(values = cbbPalette) +
    ylab("Empirical rejection probability (ERP)") +
    xlab(expression(theta)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
          text        = element_text(size = 8))
  ggsave(paste0(prefix, "_power_range.pdf"), plot = plots$power_range,
         width = 5, height = 5, units = "in", device = "pdf")
  data_focus                         <-
    dplyr::filter(tds_sum, eta == paste("eta ==", eta_focus) &
                    gamma == paste("gamma ==", gamma_focus))
  data_focus                         <-
    tidyr::gather(data_focus, key = "key", value = "value", power:prop_I)
  data_hline                         <-
    tibble::tibble(key = factor(rep(c("power", "prop_I"), each = 2)),
                   value = c(alpha, 1 - beta, min_prop, max_prop))
  data_hline                         <- rbind(data_hline,
                                              c("prop_I", sum(X)/(C*Ti)))
  data_hline$value                   <- as.numeric(data_hline$value)
  facet_labs                         <-
    c("Empirical rejection probability (ERP)",
      paste("Emp. av. prop. of cl.-per. in interv. con. (EACP)"))
  names(facet_labs)                  <- c("power", "prop_I")
  plots$focus                        <-
    ggplot(data_focus, aes(theta, value, colour = w)) +
    facet_wrap(~key, nrow = 2, scales = "free_y",
               labeller = labeller(key = facet_labs)) +
    geom_hline(data = data_hline, aes(yintercept = value), linetype = 2,
               colour = "gray75", size = 1/3) +
    geom_vline(xintercept = 0, colour = "gray75", linetype = 2, size = 1/3) +
    geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
               size = 1/3) +
    geom_line(size = 1/3) + geom_point(size = 1/3) +
    ph2rand:::theme_ph2rand() + scale_colour_manual(values = cbbPalette) +
    ylab("Value") + xlab(expression(theta)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
          text        = element_text(size = 8))
  ggsave(paste0(prefix, "_focus.pdf"), plot = plots$focus, width = 5,
         height = 5, units = "in", device = "pdf")
  tds_avX$theta_fac                  <-
    factor(paste("theta ==", tds_avX$theta),
           levels = paste("theta ==", theta))
  plots$avX                          <-
    ggplot(tds_avX, aes(Ti, rev(C), fill = Xbar_ij)) + geom_tile() +
    facet_wrap(~theta_fac, nrow = 2, labeller = label_parsed) +
    scale_fill_viridis_c(name = expression(bar(italic(X))[italic(Tij)])) +
    xlab("Time period") + ylab("Cluster") +
    geom_vline(xintercept = set_T + 0.5, linetype = 2, colour = "white",
               size = 1/3) +
    ph2rand:::theme_ph2rand() +
    theme(text = element_text(size = 8))
  if (prefix %in% c("sc1_3", "sc1_4", "sc1_5", "sc1_36")) {
    plots$avX                        <- plots$avX +
      scale_y_continuous(breaks = rev(c(1, 5, 10, 15, 20)),
                         labels = c(1, 5, 10, 15, 20)) +
      scale_x_continuous(breaks = c(1, 5, 9))
  } else if (prefix %in% c("sc2_3", "sc2_234")) {
    plots$avX                        <- plots$avX +
      scale_y_continuous(breaks = rev(c(1, 2, 3, 4)), labels = 1:4) +
      scale_x_continuous(breaks = c(1, 3, 5))
  } else {
    plots$avX                        <- plots$avX +
      scale_y_continuous(breaks = rev(c(1, 4, 8, 12)),
                         labels = c(1, 4, 8, 12)) +
      scale_x_continuous(breaks = c(1, 2, 3, 4))
  }
  ggsave(paste0(prefix, "_avX.pdf"), plot = plots$avX, width = 5, height = 5,
         units = "in", device = "pdf")
  remove(tds_avX)
  return(list(plots = plots, data = tds_sum))
}

##### TDS1 #####################################################################

X       <- rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 1, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 1, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 1, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 1, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 0, 1, 1, 1, 1),
                 c(0, 0, 0, 0, 0, 0, 1, 1, 1),
                 c(0, 0, 0, 0, 0, 0, 1, 1, 1),
                 c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                 c(0, 0, 0, 0, 0, 0, 0, 1, 1),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 1),
                 c(0, 0, 0, 0, 0, 0, 0, 0, 1))
alpha   <- 0.05
beta    <- 0.2
delta   <- 0.24
tds1_3  <- analyse_ra_sw_data(3, X, delta, alpha, beta, -1:1, c(1, 2.5), 0, 2.5,
                              "1/2", "sc1_3")
tds1_4  <- analyse_ra_sw_data(4, X, delta, alpha, beta, -1:1, c(1, 2.5), 0, 2.5,
                              "1/2", "sc1_4")
tds1_5  <- analyse_ra_sw_data(5, X, delta, alpha, beta, -1:1, c(1, 2.5), 0, 2.5,
                              "1/2", "sc1_5")
tds1_36 <- analyse_ra_sw_data(c(3, 6), X, delta, alpha, beta, -1:1, c(1, 2.5),
                              0, 2.5, "1/2", "sc1_36")

tds1_3_focus         <- dplyr::filter(tds1_3$data, eta == "eta == 0" &
                                        gamma == "gamma == 2.5")
tds1_3_focus$set_ti  <- "{3}"
tds1_4_focus         <- dplyr::filter(tds1_4$data, eta == "eta == 0" &
                                        gamma == "gamma == 2.5")
tds1_4_focus$set_ti  <- "{4}"
tds1_5_focus         <- dplyr::filter(tds1_5$data, eta == "eta == 0" &
                                        gamma == "gamma == 2.5")
tds1_5_focus$set_ti  <- "{5}"
tds1_36_focus        <- dplyr::filter(tds1_36$data, eta == "eta == 0" &
                                        gamma == "gamma == 2.5")
tds1_36_focus$set_ti <- "{3,6}"
tds1_focus           <- dplyr::bind_rows(tds1_3_focus, tds1_4_focus,
                                         tds1_5_focus, tds1_36_focus)
tds1_focus           <- tidyr::gather(tds1_focus, key = "key", value = "value",
                                      power:prop_I)
tds1_focus$set_ti    <- factor(tds1_focus$set_ti,
                               levels = c("{3,6}", "{3}", "{4}", "{5}"))
facet_labs           <-
  c("Empirical rejection probability (ERP)",
    paste("Emp. av. prop. of cl.-per. in interv. con. (EACP)"))
names(facet_labs)    <- c("power", "prop_I")
data_hline           <- tibble::tibble(key   = factor(rep("power", 2)),
                                       value = c(alpha, 1 - beta))
data_hline$value     <- as.numeric(data_hline$value)
plot_set_ti          <-
  ggplot(tds1_focus, aes(theta, value, colour = w)) +
  facet_grid(key~set_ti, scales = "free_y",
             labeller = labeller(key = facet_labs)) +
  geom_hline(data = data_hline, aes(yintercept = value), linetype = 2,
             colour = "gray75", size = 1/3) +
  geom_vline(xintercept = 0, colour = "gray75", linetype = 2, size = 1/3) +
  geom_vline(xintercept = delta, colour = "gray75", linetype = 2,
             size = 1/3) +
  geom_line(size = 1/3) + geom_point(size = 1/3) +
  ph2rand:::theme_ph2rand() + scale_colour_manual(values = cbbPalette) +
  ylab("Value") + xlab(expression(theta)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.75),
        text        = element_text(size = 8))
ggsave("tds1_set_ti.pdf", plot = plot_set_ti, width = 5, height = 5,
       units = "in", device = "pdf")

##### TDS2 #####################################################################
X        <- rbind(c(0, 1, 1, 1, 1),
                  c(0, 0, 1, 1, 1),
                  c(0, 0, 0, 1, 1),
                  c(0, 0, 0, 0, 1))
alpha    <- 0.05
beta     <- 0.1
delta    <- 0.2
tds2_3   <- analyse_ra_sw_data(3, X, delta, alpha, beta, -1:1, c(1, 2.5), 0,
                               2.5, "1/2", "sc2_3")
tds2_234 <- analyse_ra_sw_data(2:4, X, delta, alpha, beta, -1:1, c(1, 2.5), 0,
                               2.5, "1/2", "sc2_234")

##### TDS3 #####################################################################
X       <- rbind(c(0, 1, 1, 1),
                 c(0, 1, 1, 1),
                 c(0, 1, 1, 1),
                 c(0, 1, 1, 1),
                 c(0, 0, 1, 1),
                 c(0, 0, 1, 1),
                 c(0, 0, 1, 1),
                 c(0, 0, 1, 1),
                 c(0, 0, 0, 1),
                 c(0, 0, 0, 1),
                 c(0, 0, 0, 1),
                 c(0, 0, 0, 1))
alpha   <- 0.025
beta    <- 0.11
delta   <- 2
tds3_2  <- analyse_ra_sw_data(2, X, delta, alpha, beta, -1:1, c(1, 2.5), 0,
                              2.5, "1/2", "sc3_2")
tds3_23 <- analyse_ra_sw_data(2:3, X, delta, alpha, beta, -1:1, c(1, 2.5), 0,
                              2.5, "1/2", "sc3_23")
