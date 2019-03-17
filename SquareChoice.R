# This script looks at different methods of evaluating the density estimator we use (dHat). We are using three
# metrics (Mean Squared Error as well as a median and mode of the squared error as we suspect the mean might not
# be as robust as we want) to choose which square to use to use for our measurement of dHat

#### Libraries ####

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)


# For reference and playing, here's my function to find the sample mode
dens.mode <- function(x, na.rm = FALSE, n = 2^20) {
  dens <- density(x = x, na.rm = na.rm, n = n)
  return(dens$x[which.max(dens$y)])
}


#### Load Data ####

Parameters <- fread("data/NewBackend-ParametersFull_B3.csv")
Simulations <- fread("data/NewBackend-SimulationsFull_B3.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Remove invalid simulations
# Approx 1.7% of simulations, also cleans up dHat
Simulations <- Simulations[nHat >= 0]
Simulations <- Simulations[is.finite(nHat)]


#### Process Data ####

Simulations <- as_tibble(Simulations)

# Compute summary statistics by catch radius, trap spacing, and square
Simulations %>% 
  mutate(L2loss = (dHat - Density)^2) %>%
  group_by(CatchRadius, TrapSpacing, square) %>%
  summarise(meanse = mean(L2loss), medianse = median(L2loss), bias = mean(dHat-Density)) -> dHatMetrics

# Check where MSE and median squared error agree
dHatMetrics %>% 
  group_by(CatchRadius, TrapSpacing) %>%
  summarise(equal = (which.min(meanse) == which.min(medianse))) %>%
  ggplot(aes(x = CatchRadius, y = TrapSpacing, colour = equal)) + geom_point()

# Count how many are equal
dHatMetrics %>% 
  group_by(CatchRadius, TrapSpacing) %>%
  summarise(equal = (which.min(meanse) == which.min(medianse))) %>%
  ungroup(CatchRadius, TrapSpacing) %>%
  summarise(propequal = sum(equal)/length(equal))

# Slice out square with smallest MSE
dHatMetrics %>%
  group_by(CatchRadius, TrapSpacing) %>%
  slice(which.min(meanse)) -> dHatMetrics

# Look at which squares minimize MSE
dHatMetrics %>% 
  ggplot(aes(x = CatchRadius, y = TrapSpacing)) +
  geom_point(aes(colour = factor(square)))



