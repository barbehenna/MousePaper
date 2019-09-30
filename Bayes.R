########
#
# An attempt at Bayesian modeling of the trapping data
#
#######

library(mousesim)
library(tidyr)
library(ggplot2)


#### Play around with longer catch cycles ####


## We want to get to where we can estimate the true number of mice with high
## confidence. Before we get there, let's visualize how the catches fall off
## by day.


TRUE_DENSITY <- 1 # True density for simulation
TS <- 1.2 # Trap Spacing
CR <- 0.5 # Catch Radius
BORDER <- 3 # Field Border
NSQUARES <- 8 # Number of concentric squares
NUM_OBS <- 100 # Number of catching days/observations
NUM_REPS <- 20 # Number of simulations (below)

## Run NUM_REPS simulations with above constants (tweak others as needed)
## Collect the number of mice caught by day as a column of a matrix
reps <- replicate(NUM_REPS, {
  catches <- GenTrapData(
    trapSpacing = TS, catchRadius = CR, border = BORDER,
    nSquares = NSQUARES, trueDensity = TRUE_DENSITY,
    nForages = NUM_OBS
  )
  colSums(catches)
})

## Format matrix
repDF <- data.frame(reps)
names(repDF) <- paste0("rep", 1:NUM_REPS)
repDF$day <- 1:NUM_OBS
repDF <- gather(repDF, "rep", "numCaught", -day)

## Plot how each simulation's number caught falls off over the days
ggplot(repDF, aes(x = day, colour = rep)) +
  geom_path(aes(y = log(numCaught)))




## Notice a few things.
## 1) The fall off is very steep at first, then more normal
## 2) Because of the mice outside of the traps, there is a very heavy tail
## 3) We can check that each simulation is catching most of the mice as
##    number of observations -> inf with the below line
## 4) If you change to plot to log(numCaught) vs day, the lines still fall
##    off faster than linear in the begining, so numCaught falls off more
##    than exponentially

colSums(reps) / round(TRUE_DENSITY * ((2 * NSQUARES - 1) * TS + 2 * BORDER)^2)




#### Bayesian Estimate of p ####

## Model:
## X_i day's until mouse i is captured
## p common probabilty of capture on a given day
##
## X_i|p ~iid Geom(p)
## p ~ Unif(0,1)
## p|Xs ~ Beta(a,b)
##
## n = num observations
## sx = sum of X's
## a = n + 1
## b = sx - n + 1


## Get sample data
catches <- GenTrapData(trapSpacing = TS, catchRadius = CR, border = BORDER, nSquares = NSQUARES, trueDensity = TRUE_DENSITY, nForages = NUM_OBS)
x <- colSums(catches)

## Or generate some more pure data from model we're using
x <- as.numeric(table(rgeom(1000, 0.5)))

## Calculate estimates over time
n <- cumsum(x)
sx <- cumsum(x * (1:length(x)))

## Parameters
alphas <- n + 1
betas <- sx - n + 1

## Adjustment hack to correct to include the number of unobserved mice in the study
pHats <- alphas / (alphas + betas)
d <- 1:length(x)
betas <- sx - n + d * n * (1 / (1 - (1 - pHats)^d) - 1) / pHats + 1
# betas <- sx - n + d * n * (1 / (1 - (1 - pHats)^d) - 1) + 1


## Posterior statistics
postMean <- alphas / (alphas + betas)
postMode <- (alphas - 1) / (alphas + betas - 2)
postVar <- alphas * betas / ((alphas + betas)^2 * (alphas + betas + 1))



#### Alternative Methods to get pHat ####

## OLS Regression
ols.mod <- lm(-log(y) ~ x, df)
summary(ols.mod)

## Regression works quite well on the data from rgeom but it gets thrown off
## by the trapping simulation data, partly because the tail behavior (mice
## trickling in) is not geometric in distribution.


## MLE
## Basically the same as the Bayes estimator
pHatMLE <- n / sx
