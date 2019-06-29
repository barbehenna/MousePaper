#### Libraries ####

library(Rcpp)
library(data.table)
library(parallel)
library(pbapply)

sourceCpp("cpp/SimulationBackendFull.cpp")


#### Run Simulations ####

# Define parameters for various simulations
TrapSpacing <- seq(from = 0.1, to = 4, by = .1)
CatchRadius <- seq(from = 0.1, to = 4, by = 0.1)
Boarder <- seq(from = 3, to = 6, by = 3)
Density <- seq(from = 0.5, to = 5, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)


# Random paramset, not grided
numsim = 4e+6
Parameters <- data.frame(Density = runif(numsim, min=0.1,max=5), 
                         TrapSpacing = runif(numsim, min=0.05, max=4),
                         CatchRadius = runif(numsim, min=0.05, max=4),
                         Boarder = rep(c(3,6), each=numsim/2),
                         paramset = 1:numsim)

# Make array specifying each simulation to be done
SimRuns <- rep(1:nrow(Parameters), each = 1)
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
  }, mc.cores = 10, mc.set.seed = TRUE)
)
Simulations <- rbindlist(Simulations)
names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


# save results
fwrite(Parameters, "data/ParametersNewSample-random.csv")
fwrite(Simulations, "data/SimulationsNewSample-random.csv")
 





