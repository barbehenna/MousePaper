# Copyright (c) 2017 Alton Barbehenn

# This script loads a course, set of simulation parameters and their resulting statistics 
# (From MultiParameterSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)
library(ggplot2)

Sim <- read.csv("data/StudyAggregate_20171227_000337.csv", header = TRUE, row.names = NULL)


# Aggregate the statistics with the same parameters
# I chose median at this point because of it's robustness
# AggStats <- aggregate(cbind(dHat, nHat, pHat, pHatZeroNeg, pHatDropNeg, aHat) ~ square + TrapSpacing + FieldSize + CatchRadius + density, data = Sim, median)

# Using the data.table structure might be faster and more readable
# OMITING NA VALUES!!
Sim <- data.table(Sim)
AggStats <- Sim[, .(med_dHat = median(na.omit(dHat)),
                    avg_dHat = mean(na.omit(dHat)),
                    std_dHat = sd(na.omit(dHat)),
                    med_nHat = median(na.omit(nHat)),
                    avg_nHat = mean(na.omit(nHat)),
                    std_nHat = sd(na.omit(nHat)),
                    med_pHat = median(na.omit(pHat)),
                    avg_pHat = mean(na.omit(pHat)),
                    std_pHat = sd(na.omit(pHat)),
                    med_pHatZeroNeg = median(na.omit(pHatZeroNeg)),
                    avg_pHatZeroNeg = mean(na.omit(pHatZeroNeg)),
                    std_pHatZeroNeg = sd(na.omit(pHatZeroNeg)),
                    med_pHatDropNeg = median(na.omit(pHatDropNeg)),
                    avg_pHatDropNeg = mean(na.omit(pHatDropNeg)),
                    std_pHatDropNeg = sd(na.omit(pHatDropNeg)),
                    med_aHat = median(na.omit(aHat))), 
                by = .(square, TrapSpacing, FieldSize, CatchRadius, density)]

# Look at the correlation between density and dHat, colored by ring/square
p <- ggplot(data = AggStats) + geom_boxplot(mapping = aes(x = factor(density), y = med_dHat))
p
cor.test(x = AggStats$density, y = AggStats$med_dHat, method = "pearson")
cor.test(x = AggStats$density, y = AggStats$med_dHat, method = "kendall")

# So they appear to be correlated but look systematically biased
# We can see this by looking at the peaks in the density plots
p <- ggplot(data = AggStats) + geom_density(mapping = aes(x = med_dHat, colour = factor(square)))
p





















