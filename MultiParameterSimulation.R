# Copyright (c) 2017 Alton Barbehenn

# This script is an extension of Study Simulation.R. It is basically
# the same code but I'm looping through a few different trap spacings.
# I'm both storing the output in a list and saving the raw simulated statistics.

library(parallel)
library(data.table)
library(ggplot2)

## ----- Simulation Constants----
# d <- 1 #mice/m^2
# ts <- 1.5 #m
# fs <- 7*ts+7 #m
# np <- as.integer(d*fs*fs) #mice
# nv <- 4
# delta <- 0.5
# iter <- 100

# if (delta > ts/2) {
#   stop("Something weird may happen with this run becuase a mouse can be caught by two traps in the same forage")
# }

rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)

nv <- 4


# Generate a dataframe containing every parameter combination of interest
# Unit of Measure for distance is Sigma
TrapSpacing <- seq(from = 1, to = 7, by = 0.25)
CatchRadius <- seq(from = 0.1, to = 2, by = 0.1)
Border <- seq(from = 3, to = 6, by = 0.5)
Density <- seq(from = 0.25, to = 2, by = 0.25)
Parameters <- expand.grid(Density, Border, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Border", "CatchRadius", "TrapSpacing")
Parameters$FieldSize <- 7*Parameters$TrapSpacing + Parameters$Border
Parameters$NumMice <- as.integer(Parameters$Density*Parameters$FieldSize*Parameters$FieldSize)
# Parameters$Ventures <- 4

remove_rows <- which(Parameters$CatchRadius > Parameters$TrapSpacing/2)
Parameters <- Parameters[-remove_rows,]

Parameters_list <- split(Parameters, seq(nrow(Parameters)))

#For testing
Parameters <- Parameters[1:5,]
Parameters_list <- Parameters_list[1:5]


iter <- 1000
print(paste("start:", Sys.time()))
MultiParameterSimulation <- lapply(Parameters_list, function(param) {
  ts <- param$TrapSpacing
  fs <- param$FieldSize
  np <- param$NumMice
  delta <- param$CatchRadius
  d <- param$Density
  
  # Simulate the studies
  ncores <- detectCores()
  cl <- makeCluster(ncores-1, type = "FORK")
  Studies <- parLapply(cl, 1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
  # Studies <- lapply(1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
  stopCluster(cl)

  # Parse the simulated studies
  Studies <- lapply(1:iter, function(x) {
    p1 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day <= 2] == y))
    p1 <- unlist(p1)
    p2 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day >= 3] == y))
    p2 <- unlist(p2)
    return(data.frame(simnum=rep(x,64), trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
  })

  aHat <- 1:4 #ring numbers
  aHat <- (ts*2*aHat)^2 #concentric ring areas
  aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4
  squares <- list(1, 1:2, 1:3, 1:4, 4)
  
  # Calculate Calhoun-Zippen statistics
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
  
  # Save the data
  StatsDF <- as.data.frame(rbindlist(Stats))
  StatsDF$square <- factor(rep(1:5, times = iter))
  StatsDF$TrapSpacing <- ts
  StatsDF$FielSize <- fs
  StatsDF$CatchRadius <- delta
  StatsDF$NumVentures <- nv
  StatsDF$density <- d
  
  print(Sys.time())
  return(StatsDF)
})

StudyAggregate <- as.data.frame(rbindlist(MultiParameterSimulation))
write.csv(StudyAggregate, paste0("~/Documents/MousePaper/data/StudyAggregate_", Sys.time(), ".csv"), row.names = FALSE)


# lapply(MultiParameterSimulation[[sample(5, 1:length(MultiParameterSimulation))]], function(x) {
#   dHat_plt <- ggplot(x, aes(x = dHat, color = square)) + geom_density(na.rm = TRUE) + xlim(quantile(na.omit(x$dHat), 0.005)[[1]], quantile(na.omit(x$dHat), 0.995)[[1]])
#   return(dHat_plt)
# })













