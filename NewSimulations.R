#### Libraries ####

library(Rcpp)
library(data.table)
library(tidyr)
library(parallel)
library(ggplot2)
library(plotly)
library(pbapply)

sourceCpp("SimulationBackendFull.cpp")


# For reference and playing, here's my function to find the sample mode
dens.mode <- function(x, na.rm = FALSE) {
  dens <- density(x = x, na.rm = na.rm, n = 2^20)
  return(dens$x[which.max(dens$y)])
}


#### Run Simulations ####

# Define parameters for various simulations
TrapSpacing <- seq(from = 0.25, to = 6, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 6, by = 0.25)
Boarder <- seq(from = 3, to = 6, by = 6)
Density <- seq(from = 0.5, to = 5, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)

# Make array specifying each simulation to be done
SimRuns <- rep(1:nrow(Parameters), each = 1000)
SimRuns <- cbind(SimRuns, 1:length(SimRuns))

# Run simulations
# system.time(
# Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
#   s <- SimRuns[x, ]
#   paramset <- Parameters[Parameters$paramset == s[1], ]
#   res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 8, trueDensity = paramset$Density, nForages = 4)
#   return(as.data.frame(res))
# })
# )
# Simulations <- rbindlist(Simulations)
# names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


system.time(
  Simulations <- mclapply(seq(nrow(SimRuns)), function(x) {
    s <- SimRuns[x, ]
    paramset <- Parameters[Parameters$paramset == s[1], ]
    res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 8, trueDensity = paramset$Density, nForages = 4)
    return(as.data.frame(res))
  }, mc.cores = 2, mc.set.seed = TRUE)
)
Simulations <- rbindlist(Simulations)
names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


# save results
fwrite(Parameters, "data/ParametersNewSample3.csv")
fwrite(Simulations, "data/SimulationsNewSample3.csv")






