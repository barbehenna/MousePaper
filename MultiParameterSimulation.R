# Copyright (c) 2017 Alton Barbehenn

# This script has many of the same comonents as StudySimulation.R. Whereas in StudySimulation.R
# I am playing around with visualization and if and where we can detect the edge effect, in this
# script I'm simulating studies across a large parameter space. With each simulated study I'm
# I'm calculating the statistics used in a paper by Calhoun-Zippen in Barbehenn 1974. I'm then
# saving the data for later analyses. This way I can save CPU and understand the parameter space 
# better before making it finer.

# Libraries
library(parallel)
library(data.table)
library(pbapply) # not sure about rng reliability
library(pbmcapply) #pretty sure it handles rng well
library(Rcpp)
Rcpp::sourceCpp(paste0(getwd(), "/SimulationBackend.cpp"))


# Simulation constants
iterations <- 1000
nv <- 4
rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)
squares <- list(1, 1:2, 1:3, 1:4, 4) # List of the rings in each square


# List of viable parameters combinations 
# TrapSpacing <- seq(from = 1, to = 6, by = 1)
# CatchRadius <- seq(from = 0.5, to = 4, by = 0.5)
# Boarder <- seq(from = 3, to = 6, by = 1)
# Density <- seq(from = 0.5, to = 2, by = 0.5)
TrapSpacing <- seq(from = 0.25, to = 6, by = 0.25)
CatchRadius <- seq(from = 0.125, to = 3, by = 0.125)
Boarder <- seq(from = 3, to = 6, by = 1)
Density <- seq(from = 0.5, to = 2, by = 0.25)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$FieldSize <- 7*Parameters$TrapSpacing + 2*Parameters$Boarder
Parameters$NumMice <- as.integer(Parameters$Density*Parameters$FieldSize*Parameters$FieldSize)

remove_rows <- which(Parameters$CatchRadius > Parameters$TrapSpacing/2) # Not viable
Parameters <- Parameters[-remove_rows,]

Parameters$paramset <- 1:nrow(Parameters) # Label each set of parameters for later reference

# Parameters_list <- split(Parameters, seq(nrow(Parameters)))# Split dataframe by row -> lapply across list
# Parameters_list <- rep(Parameters_list, times = iterations)
# Parameters_list <- lapply(seq_along(Parameters_list), function(x) {
#   cbind(Parameters_list[[x]], "UniqueID" = x)
# })

Parameters_mat <- as.matrix(Parameters)
Parameters_mat <- replicate(iterations, Parameters_mat, simplify = FALSE)
Parameters_mat <- do.call(rbind, Parameters_mat)
Parameters_mat <- cbind(Parameters_mat, UniqueID=1:nrow(Parameters_mat))
Parameters_list <- split(Parameters_mat, Parameters_mat[,8]) # Split on UniqueID value

###############
# For testing
###############
# Parameters_list <- Parameters_list[c(1:4, 334:339)]


# Start cores for parallelization
ncores <- detectCores()
cl <- makeCluster(ncores-1, type = "FORK")


# Simulate trapping data
# Then look at the trapping data by trap in two periods (day={1,2} and day={3,4})
print("Simulating Trapping")
# TrapData <- parLapply(cl, Parameters_list, function(x) {
# TrapData <- lapply(Parameters_list, function(x) {
# TrapData <- pblapply(X = Parameters_list, cl = cl, FUN = function(x) {
TrapData <- pbmclapply(Parameters_list, mc.cores = ncores-1, mc.style = "ETA", FUN = function(x) {
  # sim <- trapSim1(ts=x$TrapSpacing, fs=x$FieldSize, np=x$NumMice, delta=x$CatchRadius, nv=nv)
  sim <- trapSim1(ts=x[4], fs=x[5], np=x[6], delta=x[3], nv=nv)
  sim <- as.data.frame(sim)
  names(sim) <- c("trap", "day")
  sim$trap <- sim$trap+1
  sim$day <- sim$day+1
  sim <- na.omit(sim)
  out <- data.frame(period1 = unlist(lapply(seq(64), function(y) sum(sim$trap[sim$day <= 2] == y))),
                    period2 = unlist(lapply(seq(64), function(y) sum(sim$trap[sim$day >= 3] == y))),
                    ring = rings)
  out$paramset <- x[7]
  out$UniqueID <- x[8]
  return(out)
})


# Analyze trapping data
print("Calculating Statistics")
# Stats <- parLapply(cl, TrapData, function(x) {
# Stats <- lapply(TrapData, function(x) {
Stats <- pblapply(TrapData, cl = cl, function(x) {
  ts <- Parameters$TrapSpacing[Parameters$paramset == x$paramset[1]]
  
  aHat <- 1:4 #ring numbers
  aHat <- (ts*2*aHat)^2 #concentric ring areas
  aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4
  
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
  
  # Save the data
  out <- data.frame(nHat, dHat, pHat, pHatDropNeg, pHatZeroNeg, aHat, square=seq(5))
  out$paramset <- x$paramset[1]
  out$UniqueID <- x$UniqueID[1]
  return(out)
})


print("Calculating Error in Density Estimates")
# error <- pblapply(unique(Sim$UniqueID), cl = cl, function(x) {
DensityError <- pblapply(unique(Sim$UniqueID), function(x) {
  avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == x][1]
  ts <- Sim$TrapSpacing[Sim$UniqueID == x][1]
  cr <- Sim$CatchRadius[Sim$UniqueID == x][1]
  tmp <- data.frame(den, avg-den, avg/den, ts, cr)
  return(tmp)
})
DensityError <- rbindlist(DensityError)
names(DensityError) <- c("den", "abs", "perc", "TrapSpacing", "CatchRadius")


Stats <- rbindlist(Stats)
TrapData <- rbindlist(TrapData)


# Stop clusters
stopCluster(cl)


# Save data
CompleteTime <- format(Sys.time(), format = "%Y%m%d_%H%M%S")

write.csv(Stats, paste0("data/", CompleteTime, "_Stats.csv"), row.names = FALSE)
write.csv(Parameters, paste0("data/", CompleteTime, "_Parameters.csv"), row.names = FALSE)
write.csv(DensityError, paste0("data/", CompleteTime, "_DensityError.csv"), row.names = FALSE)
write.csv(TrapData, paste0("data/", CompleteTime, "_TrapData.csv"), row.names = FALSE)



