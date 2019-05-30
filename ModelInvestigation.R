#### Libraries ####

library(MASS)
library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)

#### Load Data ####

Parameters <- fread("data/NewBackend-ParametersFull_B3.csv")
Simulations <- fread("data/NewBackend-SimulationsFull_B3.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Remove invalid simulations
# Approx 1.7% of simulations, also cleans up dHat
Simulations <- Simulations[nHat >= 0]
Simulations <- Simulations[is.finite(nHat)]
Simulations <- Simulations[is.finite(pHat)]

Simulations <- as_tibble(Simulations)
FullSimulations <- Simulations #backup for easy of testing

# prep data
# compute bias so that obs + bias = true
# only look at square 7 measurements
# Down sample Simulations for speed
Simulations %>% 
  # filter(square == 7) %>%
  mutate(bias = Density - dHat, acc = dHat/Density) %>%
  sample_frac(0.5) -> Simulations


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
  sample_frac(0.1) %>%
  ggplot(aes(x = sqrt(TrapSpacing), y = log(dHat/Density), colour = factor(square))) + 
    geom_smooth()

Simulations %>%
  filter(CatchRadius == 0.5 & square == 6) %>%
  mutate(acc = dHat/Density) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE), data = SimModData)
rlmmod =  rlm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE), data = SimModData)


Simulations %>%
  filter(CatchRadius == 0.5 & square == 6) %>%
  mutate(acc = dHat/Density) %>%
  sample_frac(0.1)  %>%
  ggplot(aes(x = sqrt(TrapSpacing), y = log(acc))) +
    geom_point() +
    geom_smooth(aes(colour = "default")) +
    geom_smooth(aes(colour = "lm"), method = "lm", formula = y ~ poly(x, degree = 3, raw = TRUE)) + 
    geom_smooth(aes(colour = "rlm"), method = "rlm", formula = y ~ poly(x, degree = 3, raw = TRUE))
# notice that these models seem to do really well, but clearly the variance is not equal across the plots
  

Simulations %>%
  filter(CatchRadius == 0.5 & square == 6) %>%
  mutate(acc = dHat/Density) %>%
  sample_frac(0.005)  %>%
  ggplot(aes(x = sqrt(TrapSpacing), y = log(acc))) +
  geom_point(aes(colour = sqrt(pHat))) +
  geom_smooth()


Simulations %>%
  filter(CatchRadius == 0.5 & square == 7) %>%
  mutate(acc = dHat/Density) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE) + poly(pHat, degree = 2, raw = TRUE), SimModData)
summary(lmmod)


Simulations %>% 
  filter(CatchRadius == 0.5) %>%
  mutate(acc = dHat/Density) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE) + poly(pHat, degree = 2, raw = TRUE) + square, SimModData)
summary(lmmod)
# looks like this model works almost just as well when we factor in  square as a variable
# degree 3 pHat works well too...

Simulations %>% 
  filter(CatchRadius == 0.5) %>%
  mutate(acc = dHat/Density) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE) + poly(pHat, degree = 2, raw = TRUE), SimModData)
summary(lmmod)
# but... it also well when not using square as a predictor


Simulations %>% 
  filter(CatchRadius == 0.5) %>%
  mutate(acc = dHat/Density) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE) + poly(pHat, degree = 2, raw = TRUE) + square, SimModData)
summary(lmmod)
SimModData$resid =  residuals(lmmod)
ggplot(SimModData, aes(x = sqrt(TrapSpacing), y = resid)) +
  geom_point() + 
  geom_smooth()


Simulations %>% 
  filter(CatchRadius == 0.5) %>%
  mutate(acc = dHat/Density, 
         sqrtTrapSpacing = sqrt(TrapSpacing), 
         sqrtTrapSpacing3 = sqrt(TrapSpacing)^3,
         pHat2 = pHat^2,
         pHat3 = pHat^3) %>%
  select(-dHat, -Density, -paramset, -uuid, -pd1, -pd2) -> SimModData
lmmod = lm(log(acc) ~ poly(sqrt(TrapSpacing), degree = 3, raw = TRUE) + sqrt(pHat) + pHat, SimModData)
summary(lmmod)
ggfortify::autoplot(lmmod)  #takes a while


#### May 26 ####

Simulations %>% 
  filter(CatchRadius<=4, TrapSpacing<=4) %>%
  sample_frac(0.5) %>%
  group_by(TrapSpacing, CatchRadius, square) %>%
  summarise(mse = mean((Density - dHat)^2)) %>%
  ggplot(aes(x = TrapSpacing, y = CatchRadius, colour = cut(log10(mse),5))) +
    facet_wrap(~square) + 
    geom_point() + 
    geom_abline(slope = 0.5) + 
    geom_abline(slope = 1) +
    geom_abline(slope = 2) #also try slope = 1 and 2
# This seems to indicate to me that squares [4:7] are the best in general, maybe even 3 too

# Check the distributuion of one point for reference
Simulations %>%
  filter(CatchRadius == 0.5, TrapSpacing == 2) %>%
  mutate(l2error = (Density - dHat)^2) %>%
  ggplot(aes(x = log10(l2error))) + 
    geom_density()

# Lets model square 6 with CR/TS in [0.5, 2] (good square, good zone)
# First looking at the data, dHat works REALLY well
Simulations %>%
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5) %>%
  ggplot(aes(x = factor(Density), y = dHat)) +
    geom_violin()

# Make model, works well
mod = Simulations %>%
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5) %>%
  sample_n(10000) %>%
  lm(Density ~ dHat, .)
summary(mod)
autoplot(mod)
# Can be better? Looks like theirs a funky effect

autoplot(mod) + geom_point(aes(colour = factor(Density)))
# Seems like there's some effect that's being driven by the actual density

mod = Simulations %>%
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5) %>%
  sample_n(10000) %>%
  lm(rep(1,10000) ~ I(dHat/Density) - 1, .)
summary(mod)
autoplot(mod)
# Density is definately the factor startifiing the diagnostc plots

# Look at modeling the bias
Simulations %>% 
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5, CatchRadius<2) %>%
  sample_n(10000) %>%
  ggplot(aes(x = CatchRadius/TrapSpacing, y = Density-dHat)) +
    geom_point()

mod = Simulations %>%
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5, CatchRadius<2) %>%
  sample_n(10000) %>%
  lm(Density ~ dHat + I(pHat^2), .)
summary(mod)
autoplot(mod)


# Create new variable dHat2 = - pd2^2 / (pd1 - pd2), so that dHat1 + dHat2 = pd1+pd2 
# Don't throw out any information
Simulations %>% 
  mutate(dHat2 = -(pd2)^2 / (pd1 - pd2)) -> Simulations
mod = Simulations %>%
  filter(square == 6, CatchRadius/TrapSpacing < 2, CatchRadius/TrapSpacing > 0.5, CatchRadius<2) %>%
  sample_n(10000) %>%
  lm(Density ~ dHat + I(dHat^2) + dHat2 + I(dHat2^2) + dHat:dHat2, .)
summary(mod)
autoplot(mod)





