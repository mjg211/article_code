################################################################################
# File:        Figure_and_Table_Production.R                                   #
# Author:      Michael J. Grayling (mjg211@cam.ac.uk)                          #
# Last Edited: 21/02/2018                                                      #
################################################################################

library(ggplot2)

# Use optimal.gs.sw.R to find the optimal designs listed in Table 1
TDS.1      <- list()
TDS.1[[1]] <- optimal.gs.sw(w = c(1/3, 1/3, 1/3), set.T = c(3, 5), seed = 33300)
TDS.1[[2]] <- optimal.gs.sw(w = c(1/2,   0, 1/2), set.T = c(3, 5), seed = 20200)
TDS.1[[3]] <- optimal.gs.sw(w = c(  0, 1/2, 1/2), set.T = c(3, 5), seed = 2200)
TDS.2      <- list()
TDS.2[[1]] <- optimal.gs.sw(w = c(1/3, 1/3, 1/3), set.T = c(3, 6, 9), C = 20,
                            Ti = 9, sigma.c = sqrt(1/9), sigma.e = 1,
                            beta = 0.2, delta = 0.24, N.SW = 1260, n.max = 70,
                            seed = 333)
TDS.2[[2]] <- optimal.gs.sw(w = c(1/2,   0, 1/2), set.T = c(3, 6, 9), C = 20,
                            Ti = 9, sigma.c = sqrt(1/9), sigma.e = 1,
                            beta = 0.2, delta = 0.24, N.SW = 1260, n.max = 70,
                            seed = 202)
TDS.2[[3]] <- optimal.gs.sw(w = c(  0, 1/2, 1/2), set.T = c(3, 6, 9), C = 20,
                            Ti = 9, sigma.c = sqrt(1/9), sigma.e = 1,
                            beta = 0.2, delta = 0.24, N.SW = 1260, n.max = 70,
                            seed = 22)
# Table 1 is created using components of TDS.1 and TDS.2

# Use analyse.gs.sw.R to determine the performance of the naive and adjusted
# analysis procedures for Figures 1 and 2. First TDS1
TDS.1                <- list()
TDS.1[[1]]           <- list()
TDS.1[[1]]$X         <- matrix(0, nrow = 4, ncol = 5)
TDS.1[[1]]$X[1, 1:5] <- 1
TDS.1[[1]]$X[2, 2:5] <- 1
TDS.1[[1]]$X[3, 3:5] <- 1
TDS.1[[1]]$X[4, 5]   <- 1
TDS.1[[1]]$Ti        <- 5
TDS.1[[1]]$C         <- 4
TDS.1[[1]]$n         <- 69
TDS.1[[1]]$set.T     <- c(3, 5)
TDS.1[[1]]$e         <- c(2.274, 1.660)
TDS.1[[1]]$f         <- c(0.409, 1.660)
TDS.1[[1]]$sigma.c   <- sqrt(0.02)
TDS.1[[1]]$sigma.e   <- sqrt(0.51)
TDS.1[[1]]$alpha     <- 0.05
TDS.1[[2]]           <- list()
TDS.1[[2]]$X         <- matrix(0, nrow = 4, ncol = 5)
TDS.1[[2]]$X[1, 1:5] <- 1
TDS.1[[2]]$X[2, 2:5] <- 1
TDS.1[[2]]$X[3, 3:5] <- 1
TDS.1[[2]]$X[4, 5]   <- 1
TDS.1[[2]]$Ti        <- 5
TDS.1[[2]]$C         <- 4
TDS.1[[2]]$n         <- 70
TDS.1[[2]]$set.T     <- c(3, 5)
TDS.1[[2]]$e         <- c(2.951, 1.604)
TDS.1[[2]]$f         <- c(0.677, 1.604)
TDS.1[[2]]$sigma.c   <- sqrt(0.02)
TDS.1[[2]]$sigma.e   <- sqrt(0.51)
TDS.1[[2]]$alpha     <- 0.05
TDS.1[[3]]           <- list()
TDS.1[[3]]$X         <- matrix(0, nrow = 4, ncol = 5)
TDS.1[[3]]$X[1, 1:5] <- 1
TDS.1[[3]]$X[2, 2:5] <- 1
TDS.1[[3]]$X[3, 3:5] <- 1
TDS.1[[3]]$X[4, 5]   <- 1
TDS.1[[3]]$Ti        <- 5
TDS.1[[3]]$C         <- 4
TDS.1[[3]]$n         <- 69
TDS.1[[3]]$set.T     <- c(3, 5)
TDS.1[[3]]$e         <- c(2.124, 1.706)
TDS.1[[3]]$f         <- c(-5.046, 1.706)
TDS.1[[3]]$sigma.c   <- sqrt(0.02)
TDS.1[[3]]$sigma.e   <- sqrt(0.51)
TDS.1[[3]]$alpha     <- 0.05
tau                  <- seq(from = -0.3, to = 0.5, by = 0.05)
TDS.1.analysis       <- list()
TDS.1.analysis[[1]]  <- analyse.gs.sw(TDS.1[[1]], tau = tau, seed = 1,
                                      replicates = 100000)
TDS.1.analysis[[2]]  <- analyse.gs.sw(TDS.1[[2]], tau = tau, seed = 2,
                                      replicates = 100000)
TDS.1.analysis[[3]]  <- analyse.gs.sw(TDS.1[[3]], tau = tau, seed = 3,
                                      replicates = 100000)
# Create data.frame for use with ggplot
df.TDS.1 <- data.frame(Design = factor(rep(c("Design 1", "Design 2",
                                             "Design 3"), each = 2*length(tau)),
                                       levels = c("Design 1", "Design 2",
                                                  "Design 3")),
                       Method = factor(rep(rep(c("N", "SO"),
                                               each = length(tau)), 3),
                                       levels = c("N", "SO")),
                       tau = rep(tau, 6),
                       Bias = c(TDS.1.analysis[[1]]$bias.naive,
                                TDS.1.analysis[[1]]$bias.adj,
                                TDS.1.analysis[[2]]$bias.naive,
                                TDS.1.analysis[[2]]$bias.adj,
                                TDS.1.analysis[[3]]$bias.naive,
                                TDS.1.analysis[[3]]$bias.adj),
                       RMSE = c(TDS.1.analysis[[1]]$rmse.naive,
                                TDS.1.analysis[[1]]$rmse.adj,
                                TDS.1.analysis[[2]]$rmse.naive,
                                TDS.1.analysis[[2]]$rmse.adj,
                                TDS.1.analysis[[3]]$rmse.naive,
                                TDS.1.analysis[[3]]$rmse.adj),
                       Coverage = c(TDS.1.analysis[[1]]$coverage.naive,
                                    TDS.1.analysis[[1]]$coverage.adj,
                                    TDS.1.analysis[[2]]$coverage.naive,
                                    TDS.1.analysis[[2]]$coverage.adj,
                                    TDS.1.analysis[[3]]$coverage.naive,
                                    TDS.1.analysis[[3]]$coverage.adj))
df.TDS.1$Method.by.Design <- paste(df.TDS.1$Design, df.TDS.1$Method)
# Generate component plots of Figure 1
p1 <- ggplot(data = df.TDS.1, mapping = aes(x = tau, y = Bias,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
        geom_line() + geom_point() + xlab(expression(italic(tau))) +
        ylab(expression(paste(italic(b), "(", italic(tau), ")", sep = "")))

p2 <- ggplot(data = df.TDS.1, mapping = aes(x = tau, y = RMSE,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
        geom_line() + geom_point() + xlab(expression(italic(tau))) +
        ylab(expression(paste(italic(RMSE), "(", italic(tau), ")", sep = "")))

p3 <- ggplot(data = df.TDS.1, mapping = aes(x = tau, y = Coverage,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
        geom_line() + geom_point() + xlab(expression(italic(tau))) +
        ylab(expression(paste(italic(COV), "(", italic(tau), ")", sep = "")))
# Now repeat for TDS2
TDS.2                    <- list()
TDS.2[[1]]               <- list()
TDS.2[[1]]$X             <- matrix(0, nrow = 20, ncol = 9)
TDS.2[[1]]$X[1:3, 1:9]   <- 1
TDS.2[[1]]$X[4:6, 2:9]   <- 1
TDS.2[[1]]$X[7:8, 3:9]   <- 1
TDS.2[[1]]$X[9:10, 4:9]  <- 1
TDS.2[[1]]$X[11:12, 5:9] <- 1
TDS.2[[1]]$X[13:15, 6:9] <- 1
TDS.2[[1]]$X[16:18, 8:9] <- 1
TDS.2[[1]]$X[19, 9]      <- 1
TDS.2[[1]]$Ti            <- 9
TDS.2[[1]]$C             <- 20
TDS.2[[1]]$n             <- 7
TDS.2[[1]]$set.T         <- c(3, 6, 9)
TDS.2[[1]]$e             <- c(2.638, 2.140, 1.651)
TDS.2[[1]]$f             <- c(-0.071, 0.666, 1.651)
TDS.2[[1]]$sigma.c       <- sqrt(1/9)
TDS.2[[1]]$sigma.e       <- 1
TDS.2[[1]]$alpha         <- 0.05
TDS.2[[2]]               <- list()
TDS.2[[2]]$X             <- matrix(0, nrow = 20, ncol = 9)
TDS.2[[2]]$X[1:4, 1:9]   <- 1
TDS.2[[2]]$X[5:7, 2:9]   <- 1
TDS.2[[2]]$X[8:10, 3:9]  <- 1
TDS.2[[2]]$X[11, 4:9]    <- 1
TDS.2[[2]]$X[12:13, 5:9] <- 1
TDS.2[[2]]$X[14:15, 6:9] <- 1
TDS.2[[2]]$X[16:17, 7:9] <- 1
TDS.2[[2]]$X[18:19, 8:9] <- 1
TDS.2[[2]]$X[20, 9]      <- 1
TDS.2[[2]]$Ti            <- 9
TDS.2[[2]]$C             <- 20
TDS.2[[2]]$n             <- 7
TDS.2[[2]]$set.T         <- c(3, 6, 9)
TDS.2[[2]]$e             <- c(14.412, 12.929, 1.582)
TDS.2[[2]]$f             <- c(0.040, 0.769, 1.582)
TDS.2[[2]]$sigma.c       <- sqrt(1/9)
TDS.2[[2]]$sigma.e       <- 1
TDS.2[[2]]$alpha         <- 0.05
TDS.2[[3]]               <- list()
TDS.2[[3]]$X             <- matrix(0, nrow = 20, ncol = 9)
TDS.2[[3]]$X[1:4, 1:9]   <- 1
TDS.2[[3]]$X[5:7, 2:9]   <- 1
TDS.2[[3]]$X[8:10, 3:9]  <- 1
TDS.2[[3]]$X[11, 4:9]    <- 1
TDS.2[[3]]$X[12:13, 5:9] <- 1
TDS.2[[3]]$X[14:15, 6:9] <- 1
TDS.2[[3]]$X[16, 7:9]    <- 1
TDS.2[[3]]$X[17:18, 8:9] <- 1
TDS.2[[3]]$X[19:20, 9]   <- 1
TDS.2[[3]]$Ti            <- 9
TDS.2[[3]]$C             <- 20
TDS.2[[3]]$n             <- 7
TDS.2[[3]]$set.T         <- c(3, 6, 9)
TDS.2[[3]]$e             <- c(2.261, 2.048, 1.788)
TDS.2[[3]]$f             <- c(-5.549, -4.325, 1.788)
TDS.2[[3]]$sigma.c       <- sqrt(1/9)
TDS.2[[3]]$sigma.e       <- 1
TDS.2[[3]]$alpha         <- 0.05
TDS.2.analysis           <- list()
TDS.2.analysis[[1]]      <- analyse.gs.sw(TDS.2[[1]], tau = tau, seed = 1,
                                      replicates = 100000)
TDS.2.analysis[[2]]      <- analyse.gs.sw(TDS.2[[2]], tau = tau, seed = 2,
                                      replicates = 100000)
TDS.2.analysis[[3]]      <- analyse.gs.sw(TDS.2[[3]], tau = tau, seed = 3,
                                      replicates = 100000)
# Create data.frame for use with ggplot
df.TDS.2 <- data.frame(Design = factor(rep(c("Design 1", "Design 2",
                                             "Design 3"), each = 2*length(tau)),
                                       levels = c("Design 1", "Design 2",
                                                  "Design 3")),
                       Method = factor(rep(rep(c("N", "SO"),
                                               each = length(tau)), 3),
                                       levels = c("N", "SO")),
                       tau = rep(tau, 6),
                       Bias = c(TDS.2.analysis[[1]]$bias.naive,
                                TDS.2.analysis[[1]]$bias.adj,
                                TDS.2.analysis[[2]]$bias.naive,
                                TDS.2.analysis[[2]]$bias.adj,
                                TDS.2.analysis[[3]]$bias.naive,
                                TDS.2.analysis[[3]]$bias.adj),
                       RMSE = c(TDS.2.analysis[[1]]$rmse.naive,
                                TDS.2.analysis[[1]]$rmse.adj,
                                TDS.2.analysis[[2]]$rmse.naive,
                                TDS.2.analysis[[2]]$rmse.adj,
                                TDS.2.analysis[[3]]$rmse.naive,
                                TDS.2.analysis[[3]]$rmse.adj),
                       Coverage = c(TDS.2.analysis[[1]]$coverage.naive,
                                    TDS.2.analysis[[1]]$coverage.adj,
                                    TDS.2.analysis[[2]]$coverage.naive,
                                    TDS.2.analysis[[2]]$coverage.adj,
                                    TDS.2.analysis[[3]]$coverage.naive,
                                    TDS.2.analysis[[3]]$coverage.adj))
df.TDS.1$Method.by.Design <- paste(df.TDS.2$Design, df.TDS.2$Method)
# Generate component plots of Figure 2
p1 <- ggplot(data = df.TDS.2, mapping = aes(x = tau, y = Bias,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
  geom_line() + geom_point() + xlab(expression(italic(tau))) +
  ylab(expression(paste(italic(b), "(", italic(tau), ")", sep = "")))

p2 <- ggplot(data = df.TDS.2, mapping = aes(x = tau, y = RMSE,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
  geom_line() + geom_point() + xlab(expression(italic(tau))) +
  ylab(expression(paste(italic(RMSE), "(", italic(tau), ")", sep = "")))

p3 <- ggplot(data = df.TDS.2, mapping = aes(x = tau, y = Coverage,
                                            by = Method.by.Design,
                                            col = Design, linetype = Method)) +
  geom_line() + geom_point() + xlab(expression(italic(tau))) +
  ylab(expression(paste(italic(COV), "(", italic(tau), ")", sep = "")))
