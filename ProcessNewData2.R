# Continuing from ProcessNewData.R
# i) Get better estimate for a-hat
#   a) probably as a function of TS alone, maybe linear, maybe quadratic
#   b) eliminate possibility that it depends on CR as well
# 
# ii) Choose square
#   a) get variance by square across all cr and ts pairs -> naive guess
#   b) find most accurate square for each cr and ts pair (baseline for other methods’ performance)
#   c) find elbows using aggregated results on a per cr and ts basis
#     - Do they exist?
#     - Where are they (what value)
#     - How consistent are they per simulation
#     - How do they compare to the most accurate square (compare accuracy and just, which square)
# 
# iii) go back and understand the mean and variance of our measurements in an ideal setting 
#      (maybe too on independent test data) 
#   a) can we improve still?
# 
# NOTES:
# 1. To estimate both area and cr, we may need to take many different measurements with different 
#     trap spacings before being able to pin down both variables. WE do know the true cr, so maybe 
#     we can look at the average error to assume for the measurements.
# 2. Should we just slice the border? 6 sigma is effectively an open field, we’re currently averaging 
#     across borders of 3, 4, 5, and 6. 


#### Libraries ####
library(data.table)
library(ggplot2)
library(ggfortify)


#### Load Data ####
# Load simulations where traps were allowed to overlap
Parameters <- fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Parameters.csv")
Simulations <-  fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Stats.csv")

# Merge in "True" density
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset", all.x = TRUE)


#### Process Data ####
# Limit data to just the full squares
Simulations <- Simulations[square != 5]

# standardize density estimates. Perfect value is 0 (dHat/Density = 1, if perfect)
Simulations[, `:=`(dHat.acc = (dHat - Density)/Density)] # Standardized density measurement - 1, centered about 0

# Aggregate
# Get best square and it's corresponding accuracy
# reduce size of table
Simulations.Agg <- Simulations[, list(avg.acc.err = mean(dHat.acc, na.rm = TRUE),
                                      var.acc.err = var(dHat.acc, na.rm = TRUE), 
                                      avg.nHat = mean(nHat, na.rm = TRUE), 
                                      med.nHat = median(nHat, na.rm = TRUE)),
                               by = list(TrapSpacing, CatchRadius, square, aHat)]

Simulations.Agg[, `:=`(resp = avg.acc.err + 1)]



#### Improve a-hat ####

# Analysis
# Multiplicative area correction
# accuracy is centered around 0, so about half of the avgerage accuracies should be negative, for this reason,
# I'll be using resp (= 1 + accuracy), which should be positive a lot more. In the previous analysis, the most
# accurate square was always positive, but this is not true of each square in general (it seems) 

# holding n-hat fixed, set 1 = (nHat / (aHat * aHatCorrection))/Density => aHatCorrection = (nHat / aHat) / Density = dHat (standardized, resp)
# Now look at the data with this formualtion
ggplot(Simulations.Agg, aes(x = TrapSpacing, y = resp)) +
  geom_point() +
  geom_smooth(method = "loess")

# This like some sort of expontential decay, but it has negative values
# The negative values can be absorbed into an additive constant (in the linear model) or I can redo with reps.
# The model y ~ exp(-x) fell off too slowly compared to the loess, but this model, y ~ exp(-x^2) looks great
ggplot(Simulations.Agg, aes(x = TrapSpacing, y = resp)) +
  geom_point() +
  geom_smooth(method = "loess", aes(colour = 'loess')) +
  geom_smooth(method = "lm", formula = y ~ exp(-x^2), aes(colour = 'lm'))

model <- lm(resp ~ exp(- (TrapSpacing^2)), Simulations.Agg)
summary(model)
autoplot(model)
# holy cow! those are some non-normal errors. Unfortunately, because this is all squares, some have negative 
# resp values, so it's not obvious if I can use boxcox on full data

# Let's see if we need to consider a shift in our model
# Non-linear least squares (i.e. still assuming Gauss-Markov conditions on the errors)

# I'm really just fitting the mean of each TrapSpacing
nls(resp ~ addshift + multcoef * exp(-(TrapSpacing - meanshift)^2), 
    Simulations.Agg, 
    start = list(meanshift = 1, multcoef = 5, addshift = 1))

sum(residuals(model)^2) # So the non-linear least squares model has a lower RSS (which is a good thing) 

ggplot(Simulations.Agg, aes(x = TrapSpacing, y = resp)) +
  geom_point() +
  geom_smooth(method = "loess", aes(colour = 'loess')) +
  geom_smooth(method = "lm", formula = y ~ exp(-x^2), aes(colour = 'lm')) + 
  stat_function(fun = function(x) -0.3899 + 716.0860*exp(-(x - (-1.8638))^2), aes(colour = 'nls'))


# Additive area correction
# holding nHat constant, set 1 = (nHat / (aHat + aHatCorrection))/Density => aHatCorrection = (nHat - aHat*Density)/Density
Simulations2 <- Simulations
Simulations2[, `:=`(aHatCorrection = (nHat - aHat*Density)/Density)]
Simulations.Agg2 <- Simulations2[, .(avg.aHatCorrection = mean(aHatCorrection, na.rm = TRUE)), by = .(TrapSpacing, CatchRadius, square, aHat)]

ggplot(Simulations.Agg2, aes(x = TrapSpacing, y = avg.aHatCorrection)) + 
  geom_point()
# REALLY not much pattern here

summary(model <- lm(avg.aHatCorrection ~ TrapSpacing, Simulations.Agg2)) 
sum(residuals(model)^2) # HUGE RSS!

ggplot(Simulations.Agg2, aes(x = TrapSpacing, y = avg.aHatCorrection)) + 
  geom_point() +
  geom_smooth(method = "loess")


# Let's focus on using the mutliplicative error
# Using the fit nls model: 
# SHOULD I BOOTSTRAP OR CROSS VALIDATE THIS MODEL?
area.mult.adj <- function(ts) { # This function returns the multiplicative factor to adjust aHat to make resp = 1
  return(-0.3899 + 716.0860*exp(-(ts - (-1.8638))^2))
}

Simulations[, `:=`(new.aHat = aHat * area.mult.adj(TrapSpacing))]
Simulations[, `:=`(new.dHat = nHat / new.aHat)]
Simulations[, `:=`(new.dHat.acc = (new.dHat - Density)/Density)]
Simulations.Agg <- Simulations[, list(new.avg.acc.err = mean(new.dHat.acc, na.rm = TRUE),
                                      new.resp = mean(new.dHat.acc, na.rm = TRUE) + 1,
                                      var.acc.err = var(new.dHat.acc, na.rm = TRUE)),
                               by = list(TrapSpacing, CatchRadius, square, new.aHat)]

ggplot(Simulations.Agg, aes(x = log(CatchRadius/TrapSpacing), y = new.resp)) + geom_point()
# Crazy, but we haven't picked the best squares yet, so hopefully that clears things up...



# Find most accurate square 
# new.dHat.acc should be close to 0
Simulations.Agg[, `:=`(best.abs.avg.acc = min(abs(new.avg.acc.err))), by = list(TrapSpacing, CatchRadius)]
Simulations.Agg <- Simulations.Agg[best.abs.avg.acc == abs(new.avg.acc.err), ]
Simulations.Agg[, best.abs.avg.acc := NULL]


ggplot(Simulations.Agg, aes(x = log(CatchRadius/TrapSpacing), y = new.resp)) + geom_point()
# Ok... so this approach seems to have very much stabilized the regime where CR is high and TS is low but it
# probably has made most of the rest of the cr-ts pairs worse. I've overfit one set of outliers. 







