################################################################################
# Function: ssre.sw()                                                          #
# Last Modified: 01/05/2018                                                    #
################################################################################
# Inputs                                                                       #
# ######                                                                       #
#                                                                              #
#                X - Matrix of binary treatment indicators                     #
#    sigma.e.tilde - Assumed residual s.d.                                     #
#          sigma.e - True residual s.d.                                        #
#    sigma.c.tilde - Assumed between cluster s.d.                              #
#          sigma.c - True between cluster s.d.                                 #
#            alpha - Desired type-I error-rate                                 #
#             beta - Desired type-II error-rate                                #
#            delta - Desired treatment effect                                  #
#              tau - True treatment effect                                     #
#               pi - Period effects                                            #
#            n.min - Minimal group size allowed post re-estimation             #
#            n.max - Maximal group size allowed post re-estimation             #
#                t - Re-estimation point                                       #
#          blinded - Logical. Indicates where blinded SSRE is used             #
#       correction - The value of tau_star if using blinded SSRE               #
# assume.sigma.c.0 - Logical. Indicates whether to assum sigma_c is 0          #
#             REML - Logical. Indicates whether to use REML estimation         #
#     n.components - The components required for SSRE are stored for values of #
#                    n up to this number. The larger this is the faster the    #
#                    SSRE can be performed, but it uses more memory            #
#               kr - Logical. Indicates where to use Kenward-Roger in the      #
#                    final hypothesis test                                     #
#           var.pi - Logical. Indicates whether the period effects should be   #
#                    treated as random variables                               #
#         sigma.pi - s.d. of the period effects if var.pi=T                    #
#         parallel - Logical. Should parallelisation be used                   #
#             cpus - Number of cpus to parallelise over if parallel=T          #
#       replicates - Number of replicates to use in estimating the average     #
#                    performance                                               #
#          summary - Logical. Should a summary of progress be printed          #
#             seed - Number seed for reproducibility                           #
################################################################################
# Outputs                                                                      #
# #######                                                                      #
#                                                                              #
# A list containing the inputs variables along with matrices summarising the   #
# performance of the considered re-estimation procedure.                       #
################################################################################

ssre.sw <- function(X = matrix(c(0, 0, 0, 0, 1,
                                 0, 0, 0, 1, 1,
                                 0, 0, 1, 1, 1,
                                 0, 1, 1, 1, 1), nrow = 4, ncol = 5),
                    sigma.e.tilde = sqrt(0.51), sigma.e = sqrt(0.51),
                    sigma.c.tilde = sqrt(0.02), sigma.c = sqrt(0.02),
                    alpha = 0.05, beta = 0.1, delta = 0.2, tau = 0,
                    pi = rep(0, 5), n.min = NA, n.max = 200, t = 1,
                    blinded = TRUE, correction = 0, assume.sigma.c.0 = FALSE,
                    REML = TRUE, n.components = 100, kr = FALSE,
                    var.pi = FALSE, sigma.pi = 0, parallel = TRUE,
                    cpus = 4, replicates = 100000, summary = TRUE,
                    seed = Sys.time()){
  
  # Set seeds for each replicate
  set.seed(seed)
  seeds <- sample(.Machine$integer.max, replicates)
  
  ##### ERROR CHECKING #########################################################
  
  if (!is.matrix(X) | ncol(X) == 1 | nrow(X) == 1 | !all(X %in% c(0, 1))){
    stop("X must be an indicator matrix with at least 2 rows and 2 columns.")
  }
  if (sigma.e.tilde <= 0){
    stop("Assumed within person standard deviation sigma.e must be strictly
         positive.")
  }
  if (sigma.e <= 0){
    stop("Within person standard deviation sigma.e must be strictly positive.")
  }
  if (sigma.c.tilde < 0){
    stop("Assumed between cluster standard deviation sigma.e must be strictly
         positive.")
  }
  if (sigma.c <= 0){
    stop("Assumed between cluster standard deviation sigma.e must be strictly
         positive.")
  }
  if ((alpha <= 0) | (alpha >= 1)){
    stop("Type-I error rate alpha must be strictly between 0 and 1.")
  }
  if ((beta <= 0) | (beta >= 1)){
    stop("Type-II error rate beta must be strictly between 0 and 1.")
  }
  if (delta <= 0){
    stop("Clinically relevant difference delta to power for must be strictly
         positive.")
  }
  if (!is.numeric(tau)){
    stop("True treatment effect tau must be numeric.")
  }
  if (!is.vector(pi) | (length(pi) != ncol(X))){
    stop("pi must be a vector of length ncol(X).")
  }
  if (!is.na(n.min)){
    if ((n.min%%1 != 0) | (n.min < 1)){
      stop("n.min must be either NA or a whole number greater than or equal to
           1.")
    }
  }
  if (!is.na(n.max)){
    if ((n.max%%1 != 0) | (n.max < 2)){
      stop("n.max must be either NA or a whole number greater than or equal to
           2.")
    }
  }
  if (!(t %in% 1:ncol(X))){
    stop("t must be in 1:ncol(X).")
  }
  if (!is.logical(blinded)){
    stop("blinded must be set to TRUE or FALSE.")
  }
  if (!is.numeric(correction)){
    stop("correction must be numeric.")
  }
  if (!is.logical(assume.sigma.c.0)){
    stop("assume.sigma.c.0 must be logical.")
  }
  if (!is.logical(REML)){
    stop("REML must be logical.")
  }
  if (!is.logical(var.pi)){
    stop("var.pi must be set to TRUE or FALSE.")
  }
  if (sigma.pi < 0){
    stop("Variance of period effects sigma.pi must be greater than or equal to
         0.")
  }
  if (!is.logical(parallel)){
    stop("parallel must be set to TRUE or FALSE.")
  }
  if ((cpus%%1 != 0) | (cpus < 1)){
    stop("cpus must be a whole number greater than or equal to 1.")
  }
  if ((replicates%%1 != 0) | (replicates < 1)){
    stop("replicates must be a whole number greater than or equal to 1.")
  }
  if (!is.logical(summary)){
    stop("summary must be set to TRUE or FALSE.")
  }
  
  ##### FUNCTION INITIALISATION ###############################################
  
  # Initialise a function which can assess power for given design parameters.
  # This is used to re-estimate the required value of n
  reestnFinder <- function(n.new, n.old, C, Ti, t, delta, alpha, beta, est.sigma.e,
                           est.sigma.c, X.des.mats, diag.mats, one.mats,
                           n.components, X){
    # Use the pre-stored design and covariance matrices to compute the
    # covariance matrix of the fixed effects for the inputs. That is, compute
    # Cov(\hat{beta}, \hat{beta})
    N.T              <- n.old*C*t + n.new*C*(Ti - t)
    # If n.new is larger than n.components we need to generate the components
    if (n.new > n.components){
      X.des <- list()
      n.vec <- c(rep(n.old, t), rep(n.new, Ti - t))
      pers  <- rbind(rep(0, Ti - 1), diag(1, Ti - 1, Ti - 1))
      for (c in 1:C){
        count      <- 1
        X.des[[c]] <- matrix(0, nrow = n.old*t + n.new*(Ti - t), ncol = Ti + 1)
        for (k in 1:Ti){
          for (j in 1:n.vec[k]){
            X.des[[c]][count, ] <- c(1, pers[k, ], X[c, k])
            count <- count + 1
          }
        }
      }
      X.des.mats.ret <- X.des
      # Now generate matrices used in creating Sigma later
      diag.mats.ret  <- Diagonal(n.old*t + n.new*(Ti - t))
      one.mats.ret   <- as(matrix(1, nrow = n.old*t + n.new*(Ti - t),
                                  ncol = n.old*t + n.new*(Ti - t)), "dspMatrix")
      X.des.mats[[n.new]] <- X.des.mats.ret
      diag.mats[[n.new]]  <- diag.mats.ret
      one.mats[[n.new]]   <- one.mats.ret
    }
    # Compute what the covariance matrix of theta would be for this design
    X.des            <- X.des.mats[[n.new]]
    Sigma.c.inv      <- diag.mats[[n.new]]/(est.sigma.e^2) -
                          ((est.sigma.c^2)/
                             (est.sigma.e^4 +
                                (est.sigma.c*est.sigma.e)^2*(N.T/C)))*
                                  one.mats[[n.new]]
    XT.Sigma.inv.X   <- matrix(0, nrow = Ti + 1, ncol = Ti + 1)
    for (i in 1:C){
      XT.Sigma.inv.X <- XT.Sigma.inv.X + t(X.des[[i]])%*%Sigma.c.inv%*%X.des[[i]]
    }
    Covariance       <- ginv(as.matrix(XT.Sigma.inv.X))
    # Use this to assess power
    PR               <- pt(q = qt(1 - alpha, df = N.T - C - Ti),
                           ncp = delta*sqrt(1/Covariance[Ti + 1, Ti + 1]),
                           df = N.T - C - Ti, lower.tail = FALSE)
    Score            <- PR - (1 - beta)
    return(Score)
  }
  
  # Initialise a function which evaluates the information level of a
  # conventional SW-CRT design, using Hussey & Hughes' formula
  informationSWCRT <- function(n, X, sigma.e, sigma.c){
    C      <- nrow(X)
    Ti     <- ncol(X)
    sigma2 <- sigma.e^2/n
    U      <- sum(X)
    V      <- sum(rowSums(X)^2)
    W      <- sum(colSums(X)^2)
    # Hussey and Hughes
    I      <- ((C*U - W)*sigma2 + (U^2 + C*Ti*U - Ti*W - C*V)*sigma.c^2)/
                (C*sigma2*(sigma2 + Ti*sigma.c^2))
    return(I)
  }
  
  # Initialise a function which evaluates the power of a conventional SW-CRT
  # design for a given per-cluster per-period sample size. Used in determining
  # the initial guess at the required value of n
  nSWCRT <- function(n, X, delta, alpha, beta, sigma.e, sigma.c){
    # Determine information level for input n
    I     <- informationSWCRT(n, X, sigma.e, sigma.c)
    # Evaluate 1 - power
    PnotR <- pt(q = qt(1 - alpha, df = n*nrow(X)*ncol(X) - nrow(X) - ncol(X)),
                ncp = delta*sqrt(I),
                df = n*nrow(X)*ncol(X) - nrow(X) - ncol(X))
    # Compute how far we are from desired type-II error-rate
    Score <- (beta - PnotR)^2
    return(Score)
  }
  
  # Initialise a function which simulates the outcome of a single trial for all
  # chosen design parameters
  singleTrial <- function(rep){
    # Set seed for this replicate
    set.seed(seeds[rep])
    # Assign/initialise design variables for before the interim re-assessment
    response  <- numeric(C*t*n)
    period    <- periods[[1]]
    treatment <- treatments[[1]]
    cluster   <- clusters[[1]]
    # If we're considering the effect under variable period effects, update the
    # value of mean
    if (var.pi == TRUE){
      new.pi                            <- numeric(Ti)
      for (i in 1:Ti){
        new.pi[i]                       <- rnorm(1, pi[i], sigma.pi)
      }
      mean                              <- numeric(C*t*n)
      for (c in 1:C){
        mean.c                          <- NULL
        for (j in 1:t){
          mean.c[(1 + (j - 1)*n):(j*n)] <- rep(tau*X[c, j] + new.pi[j], n)
        }
        mean[(1 + (c - 1)*t*n):(c*t*n)] <- mean.c
      }
    } else {
      mean                              <- means[[1]]
    }
    # Simulate the response for each trial participant in periods 1,...,t
    for (c in 1:C){
      response[(1 + (c - 1)*t*n):(c*t*n)] <- as.numeric(rnorm(n*t)%*%
                                                          chol.Sigma.22) +
                                               mean[(1 + (c - 1)*t*n):(c*t*n)]
    }
    # If t < T then determine the re-estimated value of n
    if (t < Ti){
      # If using the blinded method
      if (blinded == TRUE){
        # Compute S_{Ct}^2 and bar(S)_{Ct-t}^2
        if (assume.sigma.c.0 == FALSE){
          S.Ct.2      <- 0
          for (i in 1:C){
            for (j in 1:t){
              S.Ct.2  <- S.Ct.2 + sum((response[(1 +
                                                   (i - 1)*t*n +
                                                   (j - 1)*n):((i - 1)*t*n +
                                                                 j*n)] -
                                         mean(response[(1 +
                                                          (i - 1)*t*n +
                                                          (j - 1)*n):
                                                         ((i - 1)*t*n +
                                                            j*n)]))^2)
            }
          }
          S.Ct.2      <- (1/(n*C*t - C*t))*S.Ct.2
          S.mean.j.2 <- 0
          for (i in 1:C){
            for (j in 1:t){
              S.mean.j.2 <- S.mean.j.2 + sum((mean(response[cluster == i &
                                                              period == j]) -
                                                mean(response[period == j]))^2)
            }
          }
          S.mean.j.2 <- (n/(C*t - t))*S.mean.j.2
          # Use these to assign \hat{\sigma}_e^2 and \hat{\sigma}_c^2
          est.sigma.e  <- sqrt(S.Ct.2)
          if (S.mean.j.2 >= est.sigma.e^2 + factor){
            est.sigma.c  <- sqrt((S.mean.j.2 - est.sigma.e^2 - factor)/n)
          } else if (S.mean.j.2 >= est.sigma.e^2) {
            est.sigma.c  <- sqrt((S.mean.j.2 - est.sigma.e^2)/n)
          } else {
            est.sigma.c  <- 0
          }
        } else { # If assuming sigma.c = 0 we can jump to here
          est.sigma.c <- 0
          est.sigma.e <- sqrt((1/(n*C*t - 1))*sum((response -
                                                     mean(response))^2))
        }
      } else { # If using the unblinded method
        df.analysis <- data.frame(Period = factor(period,
                                                  levels = unique(period)),
                                  Cluster = factor(cluster,
                                                   levels = unique(cluster)),
                                  Treatment = factor(treatment,
                                                     levels = unique(treatment)),
                                  Response = response)
        # Fit the appropriate LMM
        if (t == 1 && sum(X[, 1:t]) == 0){
          interim.analysis <- lmer(Response ~ (1 | Cluster), data = df.analysis,
                                   REML = REML)
        } else if (t == 1 && sum(X[, 1:t]) != 0){
          interim.analysis <- lmer(Response ~ Treatment + (1 | Cluster),
                                   data = df.analysis, REML = REML)
        } else if (t != 1 && sum(X[, 1:t]) == 0){
          interim.analysis <- lmer(Response ~ Period + (1 | Cluster),
                                   data = df.analysis, REML = REML)
        } else {
          interim.analysis <- lmer(Response ~ Period + Treatment +
                                     (1 | Cluster),
                                   data = df.analysis, REML = REML)
        }
        if (assume.sigma.c.0 == FALSE) { # If assuming sigma.c = 0 we can jump
                                         # to here
          est.sigma.c      <- sqrt(VarCorr(interim.analysis)[[1]][1])
        } else {
          est.sigma.c      <- 0
        }
        est.sigma.e        <- sigma(interim.analysis)
      }
      # Now use reestnFinder() to determine the re-estimated value of n for the
      # remainder of the trial
      Under.P   <- 0
      Over.P    <- 0
      # Evaluate at n.min to check if this would provide enough power
      f.n.min   <- reestnFinder(n.min, n, C, Ti, t, delta, alpha, beta,
                                est.sigma.e, est.sigma.c, X.des.mats, diag.mats,
                                one.mats, n.components, X)
      if (f.n.min >= 0){
        reest.n <- n.min
        Over.P  <- 1
      } else {
        # If it doesn't check if n.components provides enough power
        f.n.components <- reestnFinder(n.components, n, C, Ti, t, delta, alpha,
                                       beta, est.sigma.e, est.sigma.c,
                                       X.des.mats, diag.mats, one.mats,
                                       n.components, X)
        if (f.n.components > 0) {
          # If it does then iteratively search until we find the minimal value
          # of n which provides the desired power
          a     <- n.min
          b     <- n.components
          fa    <- f.n.min
          fb    <- f.n.components
          check <- 0
          while ((b - a) > 1){
            c       <- ceiling((a + b)/2)
            fc      <- reestnFinder(c, n, C, Ti, t, delta, alpha, beta,
                                    est.sigma.e,  est.sigma.c, X.des.mats,
                                    diag.mats, one.mats, n.components, X)
            if (fc == 0){
              check <- 1
              break
            } else {
              if (((fa < 0) && (fc < 0)) || ((fa > 0) & (fc > 0))) {
                a   <- c
                fa  <- fc
              } else {
                b   <- c
                fb  <- fc
              }
            }
          }
          if (check == 0){
            reest.n <- b
          } else {
            reest.n <- c
          }
        } else { # If it doesn't check whether n.max provides enough power
          f.n.max <- reestnFinder(n.max, n, C, Ti, t, delta, alpha, beta,
                                  est.sigma.e, est.sigma.c, X.des.mats, diag.mats,
                                  one.mats, n.components, X)
          if (f.n.max > 0){
            # If it does then iteratively search until we find the minimal value
            # of n which provides the desired power
            a     <- n.components
            b     <- n.max
            fa    <- f.n.components
            fb    <- f.n.max
            check <- 0
            while ((b - a) > 1){
              c       <- ceiling((a + b)/2)
              fc      <- reestnFinder(c, n, C, Ti, t, delta, alpha, beta,
                                      est.sigma.e,  est.sigma.c, X.des.mats,
                                      diag.mats, one.mats, n.components, X)
              if (fc == 0){
                check <- 1
                break
              } else {
                if (((fa < 0) && (fc < 0)) || ((fa > 0) & (fc > 0))) {
                  a   <- c
                  fa  <- fc
                } else {
                  b   <- c
                  fb  <- fc
                }
              }
            }
            if (check == 0){
              reest.n <- b
            } else {
              reest.n <- c
            }
          } else { # If n.max doesn't provide enough power then it is the best
                   # we can do
            reest.n   <- n.max
            Under.P   <- 1
          }
        }
      }
      # If n.reest > n.components we have to generate the required components
      if (reest.n > n.components){
        period.c                          <- numeric((Ti - t)*reest.n)
        for (j in 1:(Ti - t)){
          period.c[(1 + (j - 1)*reest.n):(j*reest.n)] <- rep(j + t, reest.n)
        }
        periods.ret   <- c(periods[[1]], rep(period.c, C))
        treatment.2   <- numeric(C*(Ti - t)*reest.n)
        cluster.2     <- numeric(C*(Ti - t)*reest.n)
        mean.2        <- list()
        for (c in 1:C){
          treatment.c <- numeric(reest.n*(Ti - t))
          mean.c      <- numeric(reest.n*(Ti - t))
          for (j in 1:(Ti - t)){
            treatment.c[(1 + (j - 1)*reest.n):(j*reest.n)] <-
              rep(X[c, j + t], reest.n)
            mean.c[(1 + (j - 1)*reest.n):(j*reest.n)]      <-
              rep(tau*X[c, j + t] + pi[j + t], reest.n)
          }
          treatment.2[(1 + (c - 1)*(Ti - t)*reest.n):
                        (c*(Ti - t)*reest.n)] <- treatment.c
          cluster.2[(1 + (c - 1)*(Ti - t)*reest.n):
                      (c*(Ti - t)*reest.n)]   <- rep(c, reest.n*(Ti - t))
          mean.2[[c]]  <- mean.c
        }
        treatments.ret <- c(treatments[[1]], treatment.2)
        clusters.ret   <- c(clusters[[1]], cluster.2)
        means.ret      <- mean.2
        periods[[2]][[reest.n]]    <- periods.ret
        treatments[[2]][[reest.n]] <- treatments.ret
        clusters[[2]][[reest.n]]   <- clusters.ret
        means[[2]][[reest.n]]      <- means.ret
        Sigma.11         <- as(matrix(sigma.c^2, reest.n*(Ti - t), reest.n*(Ti - t)),
                               "dspMatrix") + Diagonal(x = rep(sigma.e^2,
                                                               reest.n*(Ti - t)))
        Sigma.12         <- as(matrix(sigma.c^2, ncol = n*t,
                                        nrow = reest.n*(Ti - t)), "dgeMatrix")
        Sigma.22.inv     <- Diagonal(n*t)/(sigma.e^2) -
                              ((sigma.c^2)/(sigma.e^4 +
                                (sigma.c*sigma.e)^2*(n*t)))*
                                  as(matrix(1, nrow = n*t, ncol = n*t), "dspMatrix")
        Sigma.fact       <- Sigma.12%*%Sigma.22.inv
        Sigma.fact.mats[[reest.n]]       <- Sigma.fact
        Sigma.c.bar      <- as(Sigma.11 - Sigma.fact%*%t(Sigma.12), "dspMatrix")
        chol.Sigma.c.bar <- matrix(0, nrow = reest.n*(Ti - t), ncol = reest.n*(Ti - t))
        chol.Sigma.c.bar[upper.tri(chol.Sigma.c.bar,
                                   diag = TRUE)] <- chol(Sigma.c.bar)@x
        
        chol.Sigma.c.bar.mats[[reest.n]] <- chol.Sigma.c.bar
      }
      # Simulate patient responses for the remainder of the trial. First
      # assign/initialise required variables
      response.2      <- numeric(C*(Ti - t)*reest.n)
      period          <- periods[[2]][[reest.n]]
      treatment       <- treatments[[2]][[reest.n]]
      cluster         <- clusters[[2]][[reest.n]]
      # If we're considering the effect under variable period effects, mean.2
      # needs to be generated
      if (var.pi == TRUE){
        mean.2        <- list()
        for (c in 1:C){
          mean.c      <- numeric(reest.n*(Ti - t))
          for (j in 1:(Ti - t)){
            mean.c[(1 + (j - 1)*reest.n):(j*reest.n)] <- rep(tau*X[c, j + t] +
                                                               new.pi[j + t],
                                                             reest.n)
          }
          mean.2[[c]] <- mean.c
        }
      } else { # Otherwise mean.2 can be assigned using pre-stored variables
        mean.2        <- means[[2]][[reest.n]]
      }
      # Simulate the response for each trial participant in periods t+1,...,T
      for (c in 1:C){
        response.2[(1 + (c - 1)*(Ti - t)*reest.n):
                     (c*(Ti - t)*reest.n)] <-
          as.numeric(rnorm((Ti - t)*reest.n)%*%
                       chol.Sigma.c.bar.mats[[reest.n]] + 
                       mean.2[[c]] + t(Sigma.fact.mats[[reest.n]]%*%
                                         (response[(1 + (c - 1)*n*t):(n*t*c)] - 
                                            mean[(1 + (c - 1)*n*t):(n*t*c)])))
      }
      response     <- c(response, response.2)
    } else {
      reest.n      <- 0
      Under.P      <- 0
      Over.P       <- 0
    }
    # Perform the final analysis using the Hussey & Hughes model
    df.analysis <- data.frame(Period = factor(period, levels = unique(period)),
                              Cluster = factor(cluster,
                                               levels = unique(cluster)),
                              Treatment = factor(treatment,
                                                 levels = unique(treatment)),
                              Response = response)
    final.analysis <- lmer(Response ~ Period + Treatment + (1 | Cluster),
                           data = df.analysis, REML = REML)
    f.est.sigma.c  <- sqrt(VarCorr(final.analysis)[[1]][1])
    f.est.sigma.e  <- sigma(final.analysis)
    if (t == Ti){
      est.sigma.c  <- f.est.sigma.c
      est.sigma.e  <- f.est.sigma.e
    }
    N          <- n*C*t + reest.n*C*(Ti - t)
    if (kr == FALSE){ # If we're not using KR we're nearly done
      tr.eff   <- fixef(final.analysis)[names(fixef(final.analysis)) ==
                                          "Treatment1"]
      pos      <- nrow(vcov(final.analysis))
      I.tr     <- 1/vcov(final.analysis)[pos, pos]
      TS       <- tr.eff*sqrt(I.tr)
      # Now compare this to the critical value
      e        <- qt(1 - alpha, df = N - C - Ti)
      if (TS > e){
        Reject <- 1
      } else {
        Reject <- 0
      }
    } else { # Otherwise we have to use pbkrtest
      reduced <- update(final.analysis, . ~ . - Treatment)
      p.value <- KRmodcomp(final.analysis, reduced)$test[1, "p.value"]/2
      if (p.value < alpha){
        Reject <- 1
      } else {
        Reject <- 0
      }
    }
    # First compute the test statistic
    # Return information on the otucome of the trial
    output   <- c(Reject, N, Under.P, Over.P, est.sigma.c, est.sigma.e,
                  f.est.sigma.c, f.est.sigma.e)
    return(output)
  }
  
  # Initialise a wrapper to use in parallelisation
  wrapper <- function(rep){
    result <- singleTrial(rep)
    return(result)
  }
  
  # Initialise a function which generates design/covariance matrices for the 
  # outcome variables
  componentsN <- function(n){
    # First generate X
    X.des <- list()
    n.vec <- c(rep(n.sw, t), rep(n, Ti - t))
    pers  <- rbind(rep(0, Ti - 1), diag(1, Ti - 1, Ti - 1))
    for (c in 1:C){
      count      <- 1
      X.des[[c]] <- matrix(0, nrow = n.sw*t + n*(Ti - t), ncol = Ti + 1)
      for (k in 1:Ti){
        for (j in 1:n.vec[k]){
          X.des[[c]][count, ] <- c(1, pers[k, ], X[c, k])
          count <- count + 1
        }
      }
    }
    X.des.mats.ret <- X.des
    # Now generate matrices used in creating Sigma later
    diag.mats.ret  <- Diagonal(n.sw*t + n*(Ti - t))
    one.mats.ret   <- as(matrix(1, nrow = n.sw*t + n*(Ti - t),
                                ncol = n.sw*t + n*(Ti - t)), "dspMatrix")
    # Now assign vectors which store information regarding period/treatment/
    # cluster to use in fitting the LMM
    period.c                          <- numeric((Ti - t)*n)
    for (j in 1:(Ti - t)){
      period.c[(1 + (j - 1)*n):(j*n)] <- rep(j + t, n)
    }
    periods.ret   <- c(periods[[1]], rep(period.c, C))
    treatment.2   <- numeric(C*(Ti - t)*n)
    cluster.2     <- numeric(C*(Ti - t)*n)
    mean.2        <- list()
    for (c in 1:C){
      treatment.c <- numeric(n*(Ti - t))
      mean.c      <- numeric(n*(Ti - t))
      for (j in 1:(Ti - t)){
        treatment.c[(1 + (j - 1)*n):(j*n)] <- rep(X[c, j + t], n)
        mean.c[(1 + (j - 1)*n):(j*n)]      <- rep(tau*X[c, j + t] + pi[j + t],
                                                  n)
      }
      treatment.2[(1 + (c - 1)*(Ti - t)*n):(c*(Ti - t)*n)] <- treatment.c
      cluster.2[(1 + (c - 1)*(Ti - t)*n):(c*(Ti - t)*n)]   <- rep(c, n*(Ti - t))
      mean.2[[c]]  <- mean.c
    }
    treatments.ret <- c(treatments[[1]], treatment.2)
    clusters.ret   <- c(clusters[[1]], cluster.2)
    means.ret      <- mean.2
    # Return outputs
    output      <- list()
    output[[1]] <- X.des.mats.ret
    output[[2]] <- diag.mats.ret
    output[[3]] <- one.mats.ret
    output[[4]] <- periods.ret
    output[[5]] <- treatments.ret
    output[[6]] <- clusters.ret
    output[[7]] <- means.ret
    return(output)
  }
  
  # Initialise an additional wrapper to use in parallelisation
  wrapper.2 <- function(n){
    result <- componentsN(n)
    return(result)
  }
  
  ##### MAIN COMPUTATIONS ######################################################
  
  # Print summary of progress if desired
  if (summary == TRUE){
    print("Initialising all required variables...")
  }
  
  # Find required sample size of a conventional SW-CRT design using assumed
  # variance values
  C       <- nrow(X)
  Ti      <- ncol(X)
  n.sw    <- suppressWarnings(optim(par = 10, fn = nSWCRT, X = X, delta = delta,
                                    alpha = alpha, beta = beta,
                                    sigma.e = sigma.e.tilde,
                                    sigma.c = sigma.c.tilde)$par)
  n.sw    <- ceiling(n.sw)
  # If n.min was not specified, then it is set to n.sw
  if (is.na(n.min)){
    n.min <- n.sw
  }
  if (is.na(n.max)){
    n.max <- n.sw
  }
  # Assign variables that are used in re-estimated required sample size, and in
  # fitting LMMs
  X.des.mats            <- list()
  diag.mats             <- list()
  one.mats              <- list()
  Sigma.fact.mats       <- list()
  chol.Sigma.c.bar.mats <- list()
  periods               <- list()
  treatments            <- list()
  clusters              <- list()
  means                 <- list()
  period.c              <- numeric(t*n.sw)
  for (j in 1:t){
    period.c[(1 + (j - 1)*n.sw):(j*n.sw)] <- rep(j, n.sw)
  }
  periods[[1]]  <- rep(period.c, C)
  treatment     <- numeric(C*t*n.sw)
  cluster       <- numeric(C*t*n.sw)
  mean          <- numeric(C*t*n.sw)
  for (c in 1:C){
    treatment.c <- numeric(t*n.sw)
    mean.c      <- NULL
    for (j in 1:t){
      treatment.c[(1 + (j - 1)*n.sw):(j*n.sw)] <- rep(X[c, j], n.sw)
      mean.c[(1 + (j - 1)*n.sw):(j*n.sw)]      <- rep(tau*X[c, j] + pi[j], n.sw)
    }
    treatment[(1 + (c - 1)*t*n.sw):(c*t*n.sw)] <- treatment.c
    cluster[(1 + (c - 1)*t*n.sw):(c*t*n.sw)]   <- rep(c, n.sw*t)
    mean[(1 + (c - 1)*t*n.sw):(c*t*n.sw)]      <- mean.c
  }
  treatments[[1]]  <- treatment
  clusters[[1]]    <- cluster
  means[[1]]       <- mean
  periods[[2]]     <- list()
  treatments[[2]]  <- list()
  clusters[[2]]    <- list()
  means[[2]]       <- list()
  Sigma.11         <- as(matrix(sigma.c^2, n.components*(Ti - t),
                                n.components*(Ti - t)),
                         "dspMatrix") + Diagonal(x = rep(sigma.e^2,
                                                         n.components*(Ti - t)))
  Sigma.12         <- as(matrix(sigma.c^2, ncol = n.sw*t,
                                nrow = n.components*(Ti - t)), "dgeMatrix")
  Sigma.22.inv     <- Diagonal(n.sw*t)/(sigma.e^2) -
                        ((sigma.c^2)/(sigma.e^4 +
                                        (sigma.c*sigma.e)^2*(n.sw*t)))*
                        as(matrix(1, nrow = n.sw*t, ncol = n.sw*t), "dspMatrix")
  Sigma.fact       <- Sigma.12%*%Sigma.22.inv
  Sigma.c.bar      <- as(Sigma.11 - Sigma.fact%*%t(Sigma.12), "dspMatrix")
  chol.Sigma.c.bar <- matrix(0, nrow = n.components*(Ti - t),
                             ncol = n.components*(Ti - t))
  chol.Sigma.c.bar[upper.tri(chol.Sigma.c.bar,
                             diag = TRUE)] <- chol(Sigma.c.bar)@x
  chol.Sigma.c.bar <- as(chol.Sigma.c.bar, "dtpMatrix")
  remove(Sigma.11)
  remove(Sigma.12)
  remove(Sigma.22.inv)
  remove(Sigma.c.bar)
  # If the design will re-estimate the required sample size, additional matrices
  # and vectors need to be generated
  if (t < Ti){
    sink("NULL")
    # Parallelise for speed
    suppressMessages(sfInit(parallel = parallel, cpus = cpus))
    sfExport("n.sw", "t", "Ti", "C", "X", "tau", "pi", "sigma.c",
             "sigma.e", "periods", "treatments", "clusters", "means",
             "componentsN")
    suppressMessages(sfLibrary(Matrix))
    results <- sfLapply(n.min:n.components, wrapper.2)
    suppressMessages(sfStop())
    for (n in 1:(n.components - n.min + 1)){
      X.des.mats[[n + n.min - 1]]      <- results[[n]][[1]]
      diag.mats[[n + n.min - 1]]       <- results[[n]][[2]]
      one.mats[[n + n.min - 1]]        <- results[[n]][[3]]
      periods[[2]][[n + n.min - 1]]    <- results[[n]][[4]]
      treatments[[2]][[n + n.min - 1]] <- results[[n]][[5]]
      clusters[[2]][[n + n.min - 1]]   <- results[[n]][[6]]
      means[[2]][[n + n.min - 1]]      <- results[[n]][[7]]
      Sigma.fact.mats[[n + n.min - 1]] <-
        Sigma.fact[1:((n + n.min - 1)*(Ti - t)), ]
      chol.Sigma.c.bar.mats[[n + n.min - 1]] <-
        chol.Sigma.c.bar[1:((n + n.min - 1)*(Ti - t)),
                         1:((n + n.min - 1)*(Ti - t))]
    }
    remove(results)
    sink()
    
  }
  Sigma.22 <- as(matrix(sigma.c^2, n.sw*t, n.sw*t),
                 "dspMatrix") + Diagonal(x = rep(sigma.e^2, n.sw*t))
  chol.Sigma.22 <- matrix(0, nrow = n.sw*t, ncol = n.sw*t)
  chol.Sigma.22[upper.tri(chol.Sigma.22, diag = TRUE)] <- chol(Sigma.22)@x
  chol.Sigma.22 <- as(chol.Sigma.22, "dtpMatrix")
  # Compute adjustment factors for the blinded method
  n        <- n.sw
  N.t      <- n*C*t
  if (t == 1){
    factor <- correction^2*sum(X[, 1:t]) - (correction^2/C)*sum(X[, 1:t])^2
  } else {
    factor <- correction^2*sum(X[, 1:t]) -
                (correction^2/C)*sum(colSums(X[, 1:t])^2)
  }
  factor   <- n*factor/(C*t - t)
  # Print summary of progress if desired
  if (summary == TRUE){
    print("Determining performance of SSRE procedure...")
  }
  # Simulate trials to assess performance of re-estimation procedure
  sink("NULL")
  suppressMessages(sfInit(parallel = parallel, cpus = cpus))
  suppressMessages(sfLibrary(lme4))
  suppressMessages(sfLibrary(stats))
  suppressMessages(sfLibrary(MASS))
  suppressMessages(sfLibrary(dplyr))
  suppressMessages(sfLibrary(Matrix))
  suppressMessages(sfLibrary(pbkrtest))
  sfExport("C", "Ti", "n", "alpha", "beta", "delta", "tau", "sigma.e", "seeds",
           "sigma.c", "pi", "n.min", "n.max", "t", "blinded", "correction",
           "assume.sigma.c.0", "n.components", "REML", "var.pi", "sigma.pi",
           "X.des.mats", "diag.mats", "one.mats", "periods", "treatments",
           "clusters", "kr", "means", "chol.Sigma.22", "factor", "reestnFinder",
           "singleTrial", "Sigma.fact.mats", "chol.Sigma.c.bar.mats", "X")
  results           <- sfLapply(1:replicates, wrapper)
  suppressMessages(sfStop())
  sink()
  results           <- matrix(unlist(results), nrow = replicates, ncol = 8,
                              byrow = TRUE)
  colnames(results) <- c("Reject", "N", "Under.P",
                         "Over.P", "est.sigma.c", "est.sigma.e",
                         "f.est.sigma.c", "f.est.sigma.e")
  # Compute summary statistics
  av.results        <- c(alpha, beta, delta, sigma.e.tilde, sigma.c.tilde, n,
                         tau, sigma.e, sigma.c, pi, n.min, n.max, t,
                         as.numeric(blinded), correction,
                         as.numeric(assume.sigma.c.0), as.numeric(REML),
                         as.numeric(kr), replicates, colMeans(results),
                         quantile(results[, 2], c(0.25, 0.5, 0.75)),
                         quantile(results[, 5], c(0.25, 0.5, 0.75)),
                         quantile(results[, 6], c(0.25, 0.5, 0.75)))
  names(av.results) <- c("alpha", "beta", "delta", "sigma.e.tilde",
                         "sigma.c.tilde", "n.sw", "tau", "sigma.e", "sigma.c",
                         paste("pi", 1:Ti, sep = ""), "n.min", "n.max", "t",
                         "blinded", "correction", "assume.sigma.c.0", "REML",
                         "kr", "replicates", "Mean Reject", "Mean N",
                         "Mean Under.P", "Mean Over.P", "Mean est.sigma.c",
                         "Mean est.sigma.e", "Mean f.est.sigma.c",
                         "Mean f.est.sigma.e", "25th N", "50th N", "75th N",
                         "25th est.sigma.c", "50th est.sigma.c",
                         "75th est.sigma.e", "25th est.sigma.e",
                         "50th est.sigma.e", "75th est.sigma.e")
  # Print summary of progress if desired
  if (summary == TRUE){
    print("Outputting...")
  }
  # Return result
  output            <- list(alpha = alpha, av.results = av.results, beta = beta,
                            blinded = blinded, correction = correction,
                            cpus = cpus, delta = delta, kr = kr,
                            n.components = n.components,
                            n.min = n.min, n.max = n.max, parallel = parallel,
                            pi = pi, REML = REML, replicates = replicates,
                            results = results, sigma.c = sigma.c,
                            sigma.c.tilde = sigma.c.tilde, sigma.e = sigma.e,
                            sigma.e.tilde = sigma.e.tilde, sigma.pi = sigma.pi,
                            t = t, tau = tau,
                            assume.sigma.c.0 = assume.sigma.c.0,
                            var.pi = var.pi, X = X)
  return(output)
}
