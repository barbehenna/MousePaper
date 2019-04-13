#### Libraries ####

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

#### Load Data ####

Parameters <- fread("data/NewBackend-ParametersFull_B3.csv")
Simulations <- fread("data/NewBackend-SimulationsFull_B3.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Remove invalid simulations
# Approx 1.7% of simulations, also cleans up dHat
Simulations <- Simulations[nHat >= 0]
Simulations <- Simulations[is.finite(nHat)]

Simulations <- as_tibble(Simulations)
FullSimulations <- Simulations #backup for easy of testing

# prep data
# compute bias so that obs + bias = true
# only look at square 7 measurements
# Down sample Simulations for speed
Simulations %>% 
  # filter(square == 7) %>%
  mutate(bias = Density - dHat) %>%
  sample_frac(0.005) -> Simulations


#### Model bias ####

Simulations %>% 
  group_by(TrapSpacing, CatchRadius) %>%
  summarise(meanbias = mean(bias), varbias = var(bias)) %>%
  filter(varbias < 1) %>%
  ggplot(aes(x = TrapSpacing, y = CatchRadius)) +
    geom_point(aes(col = abs(meanbias)<0.05, size = log10(varbias))) +
    geom_abline(slope = 0.8, intercept = -0.5, colour = 'red')

Simulations %>%
  mutate(crts = 2*CatchRadius/TrapSpacing) %>%
  ggplot(aes(x = crts, y = bias)) + 
    geom_point(aes(colour = log(TrapSpacing/CatchRadius))) + 
    geom_smooth(method = "lm")


# What I really care about is predicting Density based on a study.

lmmodel <- lm(Density ~ I(CatchRadius/TrapSpacing) + square + dHat, Simulations)
rfmodel <- randomForest(Density ~ dHat + square + TrapSpacing, Simulations)



Simulations %>%
  ggplot(aes(x = TrapSpacing, y = dHat/Density)) +
    geom_point()

Simulations  %>%
  filter(CatchRadius <= 2) %>%
  group_by(CatchRadius, TrapSpacing) %>%
  summarise(meanacc = mean(dHat/Density), varacc = var(dHat/Density)) %>%
  ggplot(aes(x = TrapSpacing, y = meanacc)) +
    geom_line(aes(colour = factor(CatchRadius)))

Simulations %>%
  filter(CatchRadius == 0.5) %>%
  group_by(TrapSpacing, square) %>%
  mutate(meanacc = mean(dHat/Density)) %>%
  ggplot(aes(x = TrapSpacing, y = meanacc)) +
    geom_line(aes(colour = factor(square)))
  
Simulations %>% 
  filter(CatchRadius == 0.5) %>%
  sample_frac(0.25) %>%
  ggplot(aes(x = TrapSpacing, y = log(dHat/Density), colour = factor(square))) + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 5))

Simulations %>%
  filter(CatchRadius == 0.05) %>%
  mutate(acc = )
  