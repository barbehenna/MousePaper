# Copyright (c) 2017 Alton Barbehenn

# This script loads a course, set of simulation parameters and their resulting statistics 
# (From MultiParameterSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)

Sim <- read.csv("data/StudyAggregate_2017-11-28_18-44-35.csv", header = TRUE, row.names = NULL)


AggStats <- aggregate(cbind(dHat, nHat, pHat, pHatZeroNeg, pHatDropNeg, aHat) ~ square + TrapSpacing + FieldSize + CatchRadius + density, data = Sim, median)







# Copyright (c) 2017 Alton Barbehenn

# This script loads a course, set of simulation parameters and their resulting statistics 
# (From MultiParameterSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)

Sim <- read.csv("data/StudyAggregate_2017-11-28_18-44-35.csv", header = TRUE, row.names = NULL)


# Aggregate the statistics with the same parameters
# I chose median at this point because of it's robustness
AggStats <- aggregate(cbind(dHat, nHat, pHat, pHatZeroNeg, pHatDropNeg, aHat) ~ square + TrapSpacing + FieldSize + CatchRadius + density, data = Sim, median)

# Using the data.table structure might be faster
Sim <- data.table(Sim)
AggStats <- Sim[, .(med_dHat = median(dHat), med_nHat = median(nHat)), by = .(square, TrapSpacing, FieldSize, CatchRadius, density)]


























