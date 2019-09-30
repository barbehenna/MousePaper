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
x <- as.numeric(table(rgeom(1000, 0.25)))

## New model (unoptimized code)

## Priors: Beta(1,1) =d Unif(0,1)
## Can adjust if you want to shift prior mean of p
## prior mean of p = a0/(a0+b0)
a0 <- 1
b0 <- 1

## Initialize pHat vectors (expectation and variance)
pHatVar <- numeric(length(x))
pHat <- numeric(length(x))

## Start with prior expectation of p
pHatOld <- a0 / (a0 + b0)

## Iterativly approximate p, collecting expectation and variance
for (i in seq(x)) {
  ## Observations to date
  no <- sum(x[1:i])
  sxo <- sum(x[1:i] * (1:i))

  ## Estimate full n and sx based on previous best estimate of p
  n <- no * (1 + (1 - pHatOld)^i / (1 - (1 - pHatOld)^i))
  sx <- sxo + (n - no) * (i + 1 / pHatOld)

  ## Calculate posterior parameters
  a <- n + a0
  b <- sx - n + b0

  ## Estimates and variance
  pHat[i] <- a / (a + b)
  pHatVar[i] <- a * b / (a + b)^2 / (a + b + 1)

  ## Update pHatOld
  pHatOld <- pHat[i]
}





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
