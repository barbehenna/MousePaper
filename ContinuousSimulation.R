# Copyright (c) 2018 Alton Barbehenn

# This script takes heavily from MiltiParameterSimulation.R (and so from StudySimulation.R too). 
# The goal of this script is to simulate studies across a large parameter space. With each simulated 
# study I'm calculating the statistics used in a paper by Calhoun-Zippen in Barbehenn 1974. I'm then
# appending the generated statistics to a csv file for further analysis.

# For a fixed, continuous range of parameters, this script samples uniformly across all the parameters
# and calculates the corresponding statistics.

########################## Initialization ########################## 

###### Libraries ######
library(Rcpp)
Rcpp::sourceCpp(paste0(getwd(), "/SimulationBackend.cpp"))


###### Define Parameter Space ######
TrapSpacingMin <- 0.1
TrapSpacingMax <- 6
CatchRadiusMin <- 0.0
# Maximum catch radius is half of the trap spacing
BoarderMin <- 3
BoarderMax <- 6
DensityMin <- 0.5
DensityMax <- 2

###### Simulation Constants ######

# Simulation Start Time
SimulationTime <- Sys.time()
print(paste("Starting Simulation at", SimulationTime))

# Simulation Constants
# ncores <- detectCores()
iterations <- 50
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


########################## Simulation ########################## 
write.table(x = t(c("nHat", "dHat", "pHat", "pHatDropNeg", "pHatZeroNeg", 
                    "aHat", "square", "UniqueID")),
            file = paste0("data/", SimulationTime, "_Stats.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
write.table(x = t(c("Density", "Boarder", "CatchRadius", "TrapSpacing", 
                    "FieldSize", "NumMice", "paramset")),
            file = paste0("data/", SimulationTime, "_Parameters.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE)
for (iter in seq(iterations)) {
  print(iter)
  ###### Get parameters for simulation ######
  TrapSpacing <- runif(n = 1, min = TrapSpacingMin, max = TrapSpacingMax)
  CatchRadius <- runif(n = 1, min = CatchRadiusMin, max = TrapSpacing/2)
  Boarder <- runif(n = 1, min = BoarderMin, max = BoarderMax)
  Density <- runif(n = 1, min = DensityMin, max = DensityMax)
  FieldSize <- 7*TrapSpacing + 2*Boarder
  NumMice <- FieldSize*FieldSize*Density

  ###### Run simulation ######
  sim <- trapSim1(ts=TrapSpacing, fs=FieldSize, np=NumMice, delta=CatchRadius, nv=nv)
  
  # Correct simulation
  sim <- as.data.frame(sim)
  names(sim) <- c("trap", "day")
  sim$trap <- sim$trap+1
  sim$day <- sim$day+1
  sim <- na.omit(sim)
  
  # Count mice per trap in each period
  period1 <- unlist(lapply(seq(64), function(y) sum(sim$trap[sim$day <= 2] == y)))
  period2 <- unlist(lapply(seq(64), function(y) sum(sim$trap[sim$day >= 3] == y)))
  traps <- data.frame(period1 = period1,
                      period2 = period2,
                      ring = rings)
  
  ####### Calculate Statistics on the simulation #######
  # Area estimates
  aHat <- 1:4 #ring numbers
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
  pHatDropNeg[which(pHat < 0)] <- 0 # NA<0 returns NA and the which only returns TRUE locations
  
  # Save the data
  out <- data.frame(nHat, dHat, pHat, pHatDropNeg, pHatZeroNeg, aHat,
                    square = seq(5), UniqueID = iter)
  write.table(x = out, file = paste0("data/", SimulationTime, "_Stats.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  paramset <- data.frame(Density, Boarder, CatchRadius, TrapSpacing, FieldSize, NumMice, UniqueID = iter)
  write.table(x = paramset, file = paste0("data/", SimulationTime, "_Parameters.csv"),
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}

