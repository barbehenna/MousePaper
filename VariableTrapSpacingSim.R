# This script is an extension of Study Simulation.R. It is basically
# the same code but I'm looping through a few different trap spacings.
# I'm both storing the output in a list and saving the raw simulated statistics.

library(parallel)
library(ggplot2)
library(data.table)

## ----- Simulation Constants----
d <- 1 #mice/m^2
# ts <- 1.5 #m
fs <- 7*ts+7 #m
np <- as.integer(d*fs*fs) #mice
nv <- 4
delta <- 0.5
iter <- 100

if (delta > ts/2) {
  stop("Something weird may happen with this run becuase a mouse can be caught by two traps in the same forage")
}

rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)

TrapSpacing <- seq(from = -2, to = 5, by = 1)
TrapSpacing <- 2^TrapSpacing

print(Sys.time())
VariableSpacingSimulation <- lapply(TrapSpacing, function(ts) {
  print(ts)
  
  # Simulate the studies
  ncores <- detectCores()
  cl <- makeCluster(ncores-1, type = "FORK")
  Studies <- parLapply(cl, 1:iter, function(x) studySim())
  stopCluster(cl)
  print(Sys.time())
  
  # Parse the simulated studies
  Studies <- lapply(1:iter, function(x) {
    p1 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day <= 2] == y))
    p1 <- unlist(p1)
    p2 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day >= 3] == y))
    p2 <- unlist(p2)
    return(data.frame(simnum=rep(x,64), trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
  })
  print(Sys.time())
  
  aHat <- 1:4 #ring numbers
  aHat <- (ts*2*aHat)^2 #concentric ring areas
  aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4
  squares <- list(1, 1:2, 1:3, 1:4, 4)
  
  # Calculate statistics
  Stats <- lapply(1:iter, function(x) {
    p1 <- lapply(squares, function(y) sum(Studies[[x]]$period1[Studies[[x]]$ring %in% y]))
    p2 <- lapply(squares, function(y) sum(Studies[[x]]$period2[Studies[[x]]$ring %in% y]))
    nHat_dHat <- lapply(1:length(squares), function(x) {
      if (p1[[x]] == p2[[x]]) {
        nHat <- NA
        dHat <- NA
      } else {
        nHat <- (p1[[x]]^2)/(p1[[x]]-p2[[x]])
        dHat <- nHat/aHat[[x]]
      }
      return(c(nHat, dHat))
    })
    nHat_dHat <- do.call(rbind, nHat_dHat)
    pHat <- lapply(1:length(squares), function(x) {
      if (p1[[x]] == 0) {
        pHat <- NA
      } else {
        pHat <- 1 - sqrt(p2[[x]]/p1[[x]])
      }
      return(pHat)
    })
    pHat <- unlist(pHat)
    return(data.frame(simnum=rep(x,length(squares)),
                      dHat=nHat_dHat[,2],
                      nHat=nHat_dHat[,1],
                      pHat=pHat,
                      aHat=aHat))
  })
  print(Sys.time())
  
  # Save the data
  StatsDF <- as.data.frame(rbindlist(Stats))
  write.csv(StatsDF, paste0("/Users/Alton/Documents/Projects/MousePaper/Data/","delta_", delta, "ts_", ts, ".csv"))
  return(StatsDF)
})


# DOUBLE CHECK THAT AHAT IS BEING ASSIGNED CORRECTLY
# AND ADD A SQUARE NUMBER TO STATSDF BEFORE IT'S RETURNED.
















