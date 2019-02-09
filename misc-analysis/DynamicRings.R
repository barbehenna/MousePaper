# Copyright (c) 2018 Alton Barbehenn

# I'm going to dynamically figure out which rings to consider. Then, use that
# to inform the analysis of trap spacing vs catch radius vs density analysis.

# To do this, I'll be using the discrete dataset. 

# For a given parameter set, I have 1000 simulations. Notice that the values of 
# ring 1 don't come from the same distributions as they're estimating different
# densities. However, if we normalize the density estimates, we can talk about
# the distribution of ring 1. 


# Discrete data
#Sim <- read.csv('./data/2018-02-18 02:01:29_Stats.csv')
#Sim <- read.csv('./data/discStats10000rows.csv') #smaller set for testing
Sim <- read.csv('./data/discStatsMiddlerows.csv') # bigger set, but still small
Param <- read.csv('./data/2018-02-18 02:01:29_Parameters.csv')

# Collect everything in one place 
Sim <- merge(Sim, Param, by="paramset")

# Free some memory
rm(Param)

# Normalize dHats
Sim <- na.omit(Sim) # remove unneccessary columns first, so I don't delete extra data
Sim$dHat <- (Sim$dHat - Sim$Density) / Sim$Density

# This seems to work well for well behaved parameter sets, however, for sets where
# the inside is totally over-run, like paramset=1, the inner most ring is different
# from the second ring, which would be fine with my intuition, but it also says the
# second and third rings are similar, so that confuses me a little. Maybe if I play
# with the app it'll help things make more sense. 

# Redo code in lapply
SimByParamset <- split(x = Sim, f = Sim$paramset)
# SimByParamset <- SimByParamset[1:10]

# TODO: test for normality, if normal, use t-test
squaresPerParamset <- pblapply(SimByParamset, function(SimSubset) {
  # Keep paramset for easy merging
  paramset <- SimSubset$paramset[1]
  
  # Rings/squares to compare
  rings <- c(2, 3, 4)
  
  out <- data.frame()
  for (ring in rings) {
    # Compare two distributions from consecutive consentric rings
    test <- wilcox.test(x = SimSubset$dHat[SimSubset$square <= (ring-1)], 
                        y = SimSubset$dHat[SimSubset$square <= ring],
                        alternative = "two.sided")
    
    # Store results
    out <- rbind(out, c((ring-1), test$p.value, paramset))
  }
  # Store results
  out <- rbind(out, c((ring), test$p.value, paramset))
  
  # Format output
  out <- as.data.frame(out)
  names(out) <- c("ring", "p.value", "paramset")

  return (out)
})

# Collect outputs outputs
squaresPerParamset <- as.data.frame(do.call("rbind", squaresPerParamset))

# Merge parameter set values for some understanding
# squaresPerParamset <- merge(x = squaresPerParamset, y = Param, by="paramset")

# Build reference table for which square to use for each parameter set
goodSquares <- aggregate(squaresPerParamset, by=list(squaresPerParamset$paramset), min)

# If the minimum p-value is not significant, each square is the same,
# so we set that ring to 4, for all traps.
goodSquares$ring[goodSquares$p.value > 0.05] = 4

# How many squares that are 4 have p-value less than 0.05?





