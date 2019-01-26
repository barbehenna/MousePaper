library(Rcpp)
library(data.table)
# library(pbapply)

sourceCpp("SimulationBackendFull.cpp")


#### Run Simulations ####
# 100 interations of (ts = 0.25->3.0, cr = 0.25->2.0, B = 3->6, D = 0.5->4.5) from takes about 23 minutes

# Define parameters for various simulations 
TrapSpacing <- seq(from = 0.25, to = 3, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 2, by = 0.25)
Boarder <- seq(from = 3, to = 6, by = 1)
Density <- seq(from = 0.5, to = 5, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)

# Make array specifying each simulation to be done
SimRuns <- rep(1:nrow(Parameters), each = 100)
SimRuns <- cbind(SimRuns, 1:length(SimRuns))

# Run simulations
set.seed(0) # I think this works
system.time(
Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
  s <- SimRuns[x, ]
  paramset <- Parameters[Parameters$paramset == s[1], ]
  res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 4, trueDensity = paramset$Density, nForages = 4)
  return(as.data.frame(res))
})
)
Simulations <- rbindlist(Simulations)
names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")




# library(microbenchmark)
# microbenchmark(Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
#   s <- SimRuns[x, ]
#   paramset <- Parameters[Parameters$paramset == s[1], ]
#   res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 4, trueDensity = paramset$Density, nForages = 4)
#   return(as.data.frame(res))
# }), times = 10)




















