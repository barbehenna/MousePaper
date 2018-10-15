# Copyright (c) 2018 Alton Barbehenn

# I'm going to dynamically figure out which rings to consider. Then, use that
# to inform the analysis of trap spacing vs catch radius vs density analysis.

# To do this, I'll be using the discrete dataset. 

# For a given parameter set, I have 1000 simulations. Notice that the values of 
# ring 1 don't come from the same distributions as they're estimating different
# densities. However, if we normalize the density estimates, we can talk about
# the distribution of ring 1. 


# Load Libraries
library(nortest)
library(pbapply)
library(tidyr)
library(ggplot2)
library(plotly)
library(data.table)


# library(parallel)
# BROKEN

# rpicluster <- makeCluster(c("ubuntu@ubuntu1.local", "ubuntu@ubuntu1.local", "ubuntu@ubuntu1.local", "ubuntu@ubuntu1.local", 
#                             "ubuntu@ubuntu3.local", "ubuntu@ubuntu3.local", "ubuntu@ubuntu3.local", "ubuntu@ubuntu3.local"),
#                           homogeneous=TRUE, rscript="/usr/bin/Rscript")
# rpicluster <- makeCluster(rep("ubuntu@ubuntu3.local", 1))
# This will just hang indefinately unless you can log in without a password. 
# The easiest way to do this (when there's a password) is to copy your public 
# key (~/.ssh/id_rsa.pub) to the remote machine and then append it to the 
# authorized users file (cat id_dsa.pub >> ~/.ssh/authorized_keys). 
# I have only done this with ubuntu1 and ubuntu3 as of writing this.

# BROKEN

# Discrete data 
Sim <- read.csv('./data/2018-02-18 02:01:29_Stats.csv')
# Sim <- read.csv('./data/discStats10000rows.csv') #smaller set for testing 
# Sim <- read.csv('./data/discStatsMiddlerows.csv') # bigger set, but still small
Param <- read.csv('./data/2018-02-18 02:01:29_Parameters.csv')

# Sim <- Sim[, c("paramset", "dHat")]
Sim <- as.data.table(Sim)
Param <- as.data.table(Param)

# Free some memory
# rm(Param)
Param <- Param[, list(Density, CatchRadius, TrapSpacing, paramset)]
Sim <- Sim[, list(dHat, square, paramset)]

# Create Coordinate (Coord) field of simulation data to ignore any parameters outside
# of CatchRadius and TrapSpacing (mostly just ignore Boundary). 
# I CAN ALSO SLICE TO HOLD BOUDARY CONSTANT.
Param$Coord <- paste(Param$CatchRadius, Param$TrapSpacing, sep = " ")
# Param$Coord <- factor(Param$Coord)
# Param <- Param[, list(Density, TrapSpacing, CatchRadius, Coord, paramset)]
# Param <- Param[, c("Density", "TrapSpacing", "CatchRadius", "Coord", "paramset")]

# Collect everything in one place 
# Sim <- merge(Sim, Param, by="paramset")
Sim <- Sim[Param, on = "paramset"]

# Normalize dHats
Sim <- na.omit(Sim) # remove unneccessary rows first, so I don't delete extra data
Sim$dHat <- (Sim$dHat - Sim$Density) / Sim$Density


# More memory management
rm(Param)
gc()

# This seems to work well for well behaved parameter sets, however, for sets where
# the inside is totally over-run, like paramset=1, the inner most ring is different
# from the second ring, which would be fine with my intuition, but it also says the
# second and third rings are similar, so that confuses me a little. Maybe if I play
# with the app it'll help things make more sense. 

# This will need to be by trap spacing and catch radius pairs. Right now, I'm breaking
# by all parameter changes. An alternate approach is to do this for a fixed boundary. 
# I don't need to fix density because I'm normalizing it. 


# Break simulations up by parameters
# SimByParamset <- split(x = Sim, f = Sim$paramset)
SimByParamset <- split(x = Sim, f = Sim$Coord)
# SimByParamset <- SimByParamset[1:10]


# More memorty management
# rm(Sim)
Sim <- Sim[, list(dHat, TrapSpacing, CatchRadius, square)]
gc()


# For each set of parameters:
# # Test for normality of densities applying to each ring 
# # Starting from inside, going outside, test for difference in mean in ring densities
# # # All normal use t-test
# # # Some not normal, use Mann-Whitney test
# # # They should have small difference because N=1000, so their powers should similar
# # Record the inner most ring for which there is homogeneity
# squaresPerParamset <- pblapply(SimByParamset, cl = cl, function(SimSubset) {
squaresPerParamset <- pblapply(SimByParamset, function(SimSubset) {
  # Keep paramset for easy merging
  # paramset <- SimSubset$paramset[1]
  # paramset <- SimSubset$Coord[1] # Factor
  # paramset <- levels(paramset)[paramset] # Character
  
  # Rings/squares to compare
  rings <- c(2, 3, 4)
  
  # Get normality of the rings
  normality <- rep(0, times = 4)
  for (i in 1:length(normality)) {
    # Lillie Test is the Kolomorov-Smirnov test comparing the data to a normal distribution
    # of known mean and variance (substituting in the sample mean and variance). The NULL 
    # hypothesis is that the distributions are equal, i.e. the data is distributed normally
    # with it's sample mean and sample variance as its mean and variance parameters. Hence
    # for SMALL p-values we reject the null and conclude that the data is not normally 
    # distributed with it's sample mean and variance.
    test <- lillie.test(SimSubset$dHat[SimSubset$square == i])
    normality[i] <- test$p.value
  }
  
  out <- data.frame()
  if (max(normality) < 0.05) { # Chance data is non-normal
    for (ring in rings) {
      # Compare two distributions from consecutive consentric rings
      ############################################################################
      #
      # FIND MANN-WHITNEY TEST INSTEAD (for theoretical support) (in NSM3?)
      #
      ############################################################################
      test <- wilcox.test(x = SimSubset$dHat[SimSubset$square == (ring-1)], 
                          y = SimSubset$dHat[SimSubset$square == ring],
                          alternative = "two.sided")
      
      # Store results
      out <- rbind(out, c(ring, test$p.value))
    }
    # Store results
    # out <- rbind(out, c(ring, test$p.value))
  } else { # Consider the densities to be normally distributed
    for (ring in rings) {
      # Compare two distributions from consecutive consentric rings
      test <- t.test(x = SimSubset$dHat[SimSubset$square == (ring-1)], 
                     y = SimSubset$dHat[SimSubset$square == ring],
                     alternative = "two.sided")
      
      # Store results
      out <- rbind(out, c(ring, test$p.value))
    }
    # Store results
    # out <- rbind(out, c(ring, test$p.value))
  }
  
  # Format output
  out <- as.data.frame(out)
  names(out) <- c("ring", "p.value")
  # out$Coord <- levels(SimSubset$Coord)[SimSubset$Coord[1]]
  out$Coord <- SimSubset$Coord[1]
  # out <- cbind(out, "normal.pvalue" = normality)

  # return (list("diffenenceInRings" = out, "normalityOfRings" = normality))
  return(out)
})


# More memory management
# rm(SimByParamset)
# gc()

# Collect outputs outputs
# squaresPerParamset <- as.data.frame(do.call("rbind", squaresPerParamset))

# Merge parameter set values for some understanding
# squaresPerParamset <- merge(x = squaresPerParamset, y = Param, by="paramset")

# Build reference table for which square to use for each parameter set
# goodSquares <- aggregate(squaresPerParamset, by=list(squaresPerParamset$paramset), min)

# If the minimum p-value is not significant, each square is the same,
# so we set that ring to 4, for all traps.
# goodSquares$ring[goodSquares$p.value > 0.05] = 4



###############################################

# For each simulation, I have the p-value per ring comparison. 
# I want:
# # the largest ring for which every ring lower than it (and equal) has
# # a no difference in mean

###############################################

goodSquares <- pblapply(squaresPerParamset, function(squareData) {
  coord <- squareData$Coord[1]
  squareData <- squareData[squareData$p.value > 0.05, ]
  
  if (is.element(2, squareData$ring)) {
    # Is it possible that 1-2 and 3-4 have the same means but 2-3 do not?
    # This would mean that the max is 4, but three is not an option...
    # Would this be an issue?
    return(c(coord, max(squareData$ring)))
  }
  else {
    return(c(coord, 1))
  }
})

goodSquares <- do.call(rbind, goodSquares)
goodSquares <- as.data.frame(goodSquares)
names(goodSquares) <- c("Coord", "square")

#################
# 
# Therer are only three factors for square
# in goodSquares. Is this because 4 is always
# different? Or is it becuase I have an error?
# 
#################


# Return the coordinates back to trap spacing and catch radiuses
goodSquares <- tidyr::separate(goodSquares, col = "Coord", into = c("CatchRadius", "TrapSpacing"), sep = " ")
goodSquares$TrapSpacing <- as.numeric(goodSquares$TrapSpacing)
goodSquares$CatchRadius <- as.numeric(goodSquares$CatchRadius)
goodSquares$square <- as.numeric(levels(goodSquares$square)[goodSquares$square])
goodSquares <- as.data.table(goodSquares)


###########################################################
# 
# Now to do the analysis...
# 
###########################################################




# We already removed NA values from Sim
# We have also already transformed dHat from it's raw values to
# # (value - true)/true
# squares are concentric and inclusive in their estimates, so if we 
# use square 3, we're using the trap data from squares 1 and 2, too.

# Get the dHat values for every point, only using it's dynammically selected square
# goodSquares <- merge(x = Sim[, c("TrapSpacing", "CatchRadius", "square", "dHat")],
#                      y = goodSquares,
#                      by = c("TrapSpacing", "CatchRadius", "square"))
goodSquares <- goodSquares[Sim, on = list(TrapSpacing, CatchRadius, square)]


# More memory management
rm(Sim)
rm(squaresPerParamset)
gc()


# Aggregate dHat's per point by mean and standard deviation
# goodSquares <- aggregate(dHat ~ ., data = goodSquares, FUN = function(x) { c(mean = mean(x), sd = sd(x)) })
# goodSquares <- as.data.table(goodSquares)
goodSquares <- goodSquares[, list(dHat.mean = mean(dHat), 
                                  dHat.sd = sd(dHat)), 
                           by = list(TrapSpacing, CatchRadius, square)]

# Plot the outter most square. This is the largest concentric ring such that all
# the density estimations are homogenious. 

plot_squares <- goodSquares %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~square, color = ~square) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "Trap Spacing"),
                      yaxis = list(title = "Catch Radius"),
                      zaxis = list(title = "Square")))
plot_squares

plot_error <- goodSquares %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~-log10(abs(dHat.mean))) %>%
  add_markers()
plot_error

#####
# TODO:
# # Fix plots. At the very least there's a merge gone wrong 
# # # and all of the points are there for each square, not 
# # # just the good ones.
# # Improve data load speed
# # Check for correctness (any squares where 4 is good?) 
# # # (is it really just squares 1 for most points?)


