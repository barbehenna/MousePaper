# Copyright (c) 2018 Alton Barbehenn

# This script takes heavily from MiltiParameterSimulation.R (and so from StudySimulation.R too). 
# The goal of this script is to simulate studies across a large parameter space. With each simulated 
# study I'm calculating the statistics used in a paper by Calhoun-Zippen in Barbehenn 1974. I'm then
# appending the generated statistics to a csv file for further analysis.

# For each of the relatively fine points in the parameter space, this script simulates 1000 studies
# and the corresponding statistics.

########################## Initialization ########################## 

###### Libraries ######
library(Rcpp)
# Rcpp::sourceCpp(paste0(getwd(), "/SimulationBackend.cpp"))

library(parallel)
library(pbapply)

###### Simulation Constants ######

# Simulation Start Time
SimulationTime <- Sys.time()
#print(paste("Starting Simulation at", SimulationTime))

# Simulation Constants
iterations <- 100
nv <- 4
rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)
rings <- c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
           8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,
           8,7,6,6,6,6,6,6,6,6,6,6,6,6,7,8,
           8,7,6,5,5,5,5,5,5,5,5,5,5,6,7,8,
           8,7,6,5,4,4,4,4,4,4,4,4,5,6,7,8,
           8,7,6,5,4,3,3,3,3,3,3,4,5,6,7,8,
           8,7,6,5,4,3,2,2,2,2,3,4,5,6,7,8,
           8,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,
           8,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,
           8,7,6,5,4,3,2,2,2,2,3,4,5,6,7,8,
           8,7,6,5,4,3,3,3,3,3,3,4,5,6,7,8,
           8,7,6,5,4,4,4,4,4,4,4,4,5,6,7,8,
           8,7,6,5,5,5,5,5,5,5,5,5,5,6,7,8,
           8,7,6,6,6,6,6,6,6,6,6,6,6,6,7,8,
           8,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,
           8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8)
# rings <- matrix(rings, nrow = 16, ncol = 16)
# image(rings, col = rainbow(8))
squares <- list(1, 1:2, 1:3, 1:4, 4) # List of the rings in each square
squares <- list(1, 1:2, 1:3, 1:4, 1:5, 1:6, 1:7, 1:8, 8)

###### Build Parameter Space ######

# List of allpossible combinations in the parameter space
TrapSpacing <- seq(from = 0.25, to = 2, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 1, by = 0.25)
Boarder <- seq(from = 3, to = 3, by = 1)
Density <- seq(from = 0.5, to = 5, by = 1)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$FieldSize <- 7*Parameters$TrapSpacing + 2*Parameters$Boarder
Parameters$NumMice <- round(Parameters$Density*Parameters$FieldSize*Parameters$FieldSize)

# Remove unviable rows
# remove_rows <- which(Parameters$CatchRadius > Parameters$TrapSpacing/2) # Not viable
# Parameters <- Parameters[-remove_rows,]

# Label each set of parameters for later reference
Parameters$paramset <- 1:nrow(Parameters)
N <- nrow(Parameters) # For ease later

# Save parameterspace used for this simulation
write.csv(Parameters, paste0("data/", SimulationTime, "_Parameters.csv"), row.names = FALSE)


########################## Simulation ########################## 
write.table(x = t(c("nHat", "dHat", "pHat", "pHatDropNeg", "pHatZeroNeg", 
                  "aHat", "square", "paramset", "UniqueID")),
            file = paste0("data/", SimulationTime, "_Stats.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)

# Make cluster and give it necessary parameters
cl <- makePSOCKcluster(names = rep("localhost", 3))
clusterCall(cl = cl, Rcpp::sourceCpp, file = paste0(getwd(),"/SimulationBackend.cpp"))
clusterExport(cl = cl, varlist = list("Parameters", "iterations", "nv", "rings", "squares", "SimulationTime"))

parLapply(cl = cl, X = sample(N), fun = function(param) {
  ###### Get parameters for simulation ######
  TrapSpacing <- Parameters$TrapSpacing[Parameters$paramset == param]
  CatchRadius <- Parameters$CatchRadius[Parameters$paramset == param]
  NumMice <- Parameters$NumMice[Parameters$paramset == param]
  FieldSize <- Parameters$FieldSize[Parameters$paramset == param]
  
  ###### Repeat simulation iteration times ######
  for (iter in seq(iterations)) {
    ###### Run simulation ######
    sim <- trapSim1(ts=TrapSpacing, fs=FieldSize, np=NumMice, delta=CatchRadius, nv=nv)
    
    # Correct simulation
    sim <- as.data.frame(sim)
    names(sim) <- c("trap", "day")
    sim$trap <- sim$trap+1
    sim$day <- sim$day+1
    sim <- na.omit(sim)
    
    # Count mice per trap in each period
    period1 <- unlist(lapply(seq(256), function(y) sum(sim$trap[sim$day <= 2] == y)))
    period2 <- unlist(lapply(seq(256), function(y) sum(sim$trap[sim$day >= 3] == y)))
    traps <- data.frame(period1 = period1,
                        period2 = period2,
                        ring = rings)
    
    ####### Calculate Statistics on the simulation #######
    # Area estimates
    aHat <- 1:8 #ring numbers
    aHat <- (TrapSpacing*2*aHat)^2 #concentric ring areas
    aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4
    
    # Catches per ring in each period
    p1 <- unlist(lapply(squares, function(y) sum(traps$period1[traps$ring %in% y])))
    p2 <- unlist(lapply(squares, function(y) sum(traps$period2[traps$ring %in% y])))
    
    # Density and abundance estamates
    nHat <- (p1^2)/(p1-p2)
    nHat[which(nHat == Inf)] <- NA #Keep as Inf? Store both or post process?
    dHat <- nHat/aHat
    
    # Catch rate estimates
    pHat <- 1 - sqrt(p2/p1)
    pHat[which(pHat == Inf)] <- NA
    pHatDropNeg <- pHat
    pHatDropNeg[which(pHat < 0)] <- NA # NA<0 returns NA and the which only returns TRUE locations
    pHatZeroNeg <- pHat
    pHatZeroNeg[which(pHat < 0)] <- 0 # NA<0 returns NA and the which only returns TRUE locations
    # There was a typo in the above line, pHatDropNeg instead of pHatZeroNeg
    
    # Save the data
    out <- data.frame(nHat, dHat, pHat, pHatDropNeg, pHatZeroNeg, aHat, 
                      square = seq(9), paramset = param, UniqueID = ((param-1)*iterations + iter))
    write.table(x = out, file = paste0("data/", SimulationTime, "_Stats.csv"), 
                sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
})

print(paste("Time Elapsed =", Sys.time()-SimulationTime))
stopCluster(cl)
