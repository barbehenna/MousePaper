# Copyright (c) 2017 Alton Barbehenn

# This script loads a course, set of simulation parameters and their resulting statistics 
# (From MultiParameterSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

Sim <- read.csv("data/StudyAggregate_2017-11-28_18-44-35.csv", header = TRUE, row.names = NULL)