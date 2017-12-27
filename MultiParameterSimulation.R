# Copyright (c) 2017 Alton Barbehenn

# This script has many of the same comonents as StudySimulation.R. Whereas in StudySimulation.R
# I am playing around with visualization and if and where we can detect the edge effect, in this
# script I'm simulating studies across a large parameter space. With each simulated study I'm
# I'm calculating the statistics used in a paper by Calhoun-Zippen in Barbehenn 1974. I'm then
# saving the data for later analyses. This way I can save CPU and understand the parameter space 
# better before making it finer.

library(parallel)
library(ggplot2)
library(data.table)


rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)


# Generate a dataframe containing every parameter combination of interest
# TrapSpacing <- seq(from = 1, to = 10, by = 0.25)
# CatchRadius <- seq(from = 0.1, to = 4, by = 0.1)
# Boarder <- seq(from = 3, to = 10, by = 1)
# Density <- seq(from = 0.25, to = 2, by = 0.25)
TrapSpacing <- seq(from = 1, to = 6, by = 1)
CatchRadius <- seq(from = 0.5, to = 4, by = 0.5)
Boarder <- seq(from = 3, to = 6, by = 1)
Density <- seq(from = 0.5, to = 2, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$FieldSize <- 7*Parameters$TrapSpacing + Parameters$Boarder
Parameters$NumMice <- as.integer(Parameters$Density*Parameters$FieldSize*Parameters$FieldSize)

remove_rows <- which(Parameters$CatchRadius > Parameters$TrapSpacing/2)
Parameters <- Parameters[-remove_rows,]

Parameters_list <- split(Parameters, seq(nrow(Parameters)))

# For testing
# Parameters <- Parameters[1:20,]
# Parameters_list <- Parameters_list[1:20]


iter <- 100 #number of simulations per study (set of parameters)
print(Sys.time())
VariableSpacingSimulation <- lapply(Parameters_list, function(param) {
  ts <- param$TrapSpacing
  fs <- param$FieldSize
  np <- param$NumMice
  delta <- param$CatchRadius
  nv <- 4
  d <- param$Density

  # Simulate the studies
  ncores <- detectCores()
  cl <- makeCluster(ncores-1, type = "FORK")
  Studies <- parLapply(cl, 1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
  # Studies <- lapply(1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
  stopCluster(cl)

  # Parse the simulated studies for the number of mice caught in each of the two periods, per trap
  Studies <- lapply(Studies, function(x) {
    p1 <- lapply(1:64, function(y) sum(x$trap[x$day <= 2] == y))
    p1 <- unlist(p1)
    p2 <- lapply(1:64, function(y) sum(x$trap[x$day >= 3] == y))
    p2 <- unlist(p2)
    # return(data.frame(simnum=rep(x,64), trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
    return(data.frame(trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
  })

  # Build a list of the rings in each square
  squares <- list(1, 1:2, 1:3, 1:4, 4)

  # Estimate the area in each of the squares
  aHat <- 1:4 #ring numbers
  aHat <- (ts*2*aHat)^2 #concentric ring areas
  aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4

  # Calculate Calhoun-Zippen and Barbehenn statistics
  Stats <- lapply(Studies, function(x) {
    p1 <- unlist(lapply(squares, function(y) sum(x$period1[x$ring %in% y])))
    p2 <- unlist(lapply(squares, function(y) sum(x$period2[x$ring %in% y])))

    nHat <- (p1^2)/(p1-p2)
    nHat[which(nHat == Inf)] <- NA #Keep as Inf? Store both or post process?
    dHat <- nHat/aHat

    pHat <- 1 - sqrt(p2/p1)
    pHat[which(pHat == Inf)] <- NA

    pHatDropNeg <- pHat
    pHatDropNeg[which(pHat < 0)] <- NA # NA<0 returns NA and the which only returns TRUE locations

    pHatZeroNeg <- pHat
    pHatDropNeg[which(pHat < 0)] <- 0 # NA<0 returns NA and the which only returns TRUE locations

    return(data.frame(#simnum=rep(x,length(squares)),
                      dHat=nHat,
                      nHat=dHat,
                      pHat=pHat,
                      pHatZeroNeg=pHatZeroNeg,
                      pHatDropNeg=pHatDropNeg,
                      aHat=aHat))
  })

  # Save the data
  StatsDF <- as.data.frame(rbindlist(Stats))
  StatsDF$square <- factor(rep(1:5, times = iter))
  StatsDF$TrapSpacing <- ts
  StatsDF$FieldSize <- fs
  StatsDF$CatchRadius <- delta
  StatsDF$NumVisits <- nv
  StatsDF$density <- d
  StatsDF$simnum <- rep(1:iter, each = 5)

  print(Sys.time())
  return(StatsDF)
})

# Save to the data directory for later use
# The working directory should be the mousepaper project directory
StudyAggregate <- as.data.frame(rbindlist(VariableSpacingSimulation))
StudyAggregate$paramnum <- rep(1:length(Parameters_list), each = iter*5)
write.csv(StudyAggregate, paste0("data/StudyAggregate_", format(Sys.time(), format = "%Y%m%d_%H%M%S"), ".csv"), row.names = FALSE)





# This is an old implementation and is included for reference and comparison

# VariableSpacingSimulation <- lapply(Parameters_list, function(param) {
#   ts <- param$TrapSpacing
#   fs <- param$FieldSize
#   np <- param$NumMice
#   delta <- param$CatchRadius
#   nv <- 4
#   d <- param$Density
# 
#   # Simulate the studies
#   ncores <- detectCores()
#   cl <- makeCluster(ncores-1, type = "FORK")
#   Studies <- parLapply(cl, 1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
#   # Studies <- lapply(1:iter, function(x) studySim(ts=ts, fs=fs, np=np, delta=delta, nv=nv, d=d))
#   stopCluster(cl)
# 
#   # Parse the simulated studies for the number of mice caught in each of the two periods, per trap
#   Studies <- lapply(Studies, function(x) {
#     p1 <- lapply(1:64, function(y) sum(x$trap[x$day <= 2] == y))
#     p1 <- unlist(p1)
#     p2 <- lapply(1:64, function(y) sum(x$trap[x$day >= 3] == y))
#     p2 <- unlist(p2)
#     # return(data.frame(simnum=rep(x,64), trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
#     return(data.frame(trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
#   })
# 
#   # Build a list of the rings in each square
#   squares <- list(1, 1:2, 1:3, 1:4, 4)
# 
#   # Estimate the area in each of the squares
#   aHat <- 1:4 #ring numbers
#   aHat <- (ts*2*aHat)^2 #concentric ring areas
#   aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4
# 
#   # Calculate Calhoun-Zippen and Barbehenn statistics
#   Stats <- lapply(Studies, function(x) {
#     p1 <- lapply(squares, function(y) sum(x$period1[x$ring %in% y]))
#     p2 <- lapply(squares, function(y) sum(x$period2[x$ring %in% y]))
#     nHat_dHat <- lapply(1:length(squares), function(x) {
#       if (p1[[x]] == p2[[x]]) {
#         nHat <- NA
#         dHat <- NA
#       } else {
#         nHat <- (p1[[x]]^2)/(p1[[x]]-p2[[x]])
#         dHat <- nHat/aHat[[x]]
#       }
#       return(c(nHat, dHat))
#     })
#     nHat_dHat <- do.call(rbind, nHat_dHat)
#     pHat <- lapply(1:length(squares), function(x) {
#       if (p1[[x]] == 0) {
#         pHat <- NA
#       } else {
#         pHat <- 1 - sqrt(p2[[x]]/p1[[x]])
#       }
#       return(pHat)
#     })
#     pHatDropNeg <- lapply(pHat, function(x) {
#       if (x < 0 && !is.na(x)) {
#         return(NA)
#       } else {
#         return(x)
#       }
#     })
#     pHatZeroNeg <- lapply(pHat, function(x) {
#       if (x < 0 && !is.na(x)) {
#         return(0)
#       } else {
#         return(x)
#       }
#     })
#     pHat <- unlist(pHat)
#     pHatDropNeg <- unlist(pHatDropNeg)
#     pHatZeroNeg <- unlist(pHatZeroNeg)
#     return(data.frame(#simnum=rep(x,length(squares)),
#       dHat=nHat_dHat[,2],
#       nHat=nHat_dHat[,1],
#       pHat=pHat,
#       pHatZeroNeg=pHatZeroNeg,
#       pHatDropNeg=pHatDropNeg,
#       aHat=aHat))
#   })
# 
#   # Save the data
#   StatsDF <- as.data.frame(rbindlist(Stats))
#   StatsDF$square <- factor(rep(1:5, times = iter))
#   StatsDF$TrapSpacing <- ts
#   StatsDF$FieldSize <- fs
#   StatsDF$CatchRadius <- delta
#   StatsDF$NumVisits <- nv
#   StatsDF$density <- d
# 
#   print(Sys.time())
#   return(StatsDF)
# })
# 
# 
# StudyAggregate <- as.data.frame(rbindlist(VariableSpacingSimulation))
# StudyAggregate$paramnum <- rep(1:length(Parameters_list), each = iter*5)
# write.csv(StudyAggregate, paste0("data/StudyAggregate_", format(Sys.time(), format = "%Y%m%d_%H%M%S"), ".csv"), row.names = FALSE)


