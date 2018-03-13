##### E.g. determination of designs for Table 2 and Figure 1 ###################

# Load required packages
library(iterpc)
library(snowfall)
library(schoolmath)

# Set global parameter values as given in Section 2.4
K <- 2 ; alpha <- 0.15; beta <- 0.20; delta <- 0.15; pESS <- 0.70; r <- 1; 
n.max <- 70; n.init <- 40

# Create matrix of Fisher designs to determine as specified in Section 2.4
alpha.1      <- seq(from = 0.01, to = 0.14, by = 0.01)
beta.1       <- seq(from = 0.01, to = 0.19, by = 0.01)
poss.designs <- expand.grid(alpha.1, beta.1)

# The function two_stage_fisher() checks a single (alpha.1,beta.1) combination.
# Search over these in parallel using fisher_wrapper()
suppressMessages(sfInit(parallel = TRUE, cpus = 7))
sfLibrary(iterpc)
sfExport("K", "r", "alpha", "beta", "poss.designs", "delta", "n.max", "n.init",
         "pESS")
sfExport("two_stage_fisher")
results <- sfLapply(1:nrow(poss.designs), fisher_wrapper)
# For each (alpha.1,beta.1) combation, five .csv files are produced in the
# current wd containing information on the determined design.

# The function two_stage_binary() computed the exact binomial designs. 
# Parallelisation is performed internally, so we simply run:
two_stage_binary()
# For each considered value of n (between n.min and n.max) a .csv file is
# produced containing summary information on the performance of each considered
# design

# We can then consider the results.score.1.##.csv files to determine the optimal
# Fisher design, and do the same for the exact binomial designs using the
# exact.designs.##.csv files