library(data.table)


# Load simulations where traps were allowed to overlap
Parameters <- fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Parameters.csv")
Simulations <-  fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Stats.csv")

# Merge in "True" density
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset", all.x = TRUE)

# Process Simulations
# remove columns that contain NA's (or are just unused and easily created)
Simulations <- Simulations[, c("pHatDropNeg", "pHatZeroNeg") := NULL]
# remove rows containing an NA
Simulations <- na.omit(Simulations)
# standardize density estimates. Perfect value is 0 (dHat/Density = 1, if perfect)
Simulations[, `:=`(dHat.acc = (dHat - Density)/Density)] # Standardized density measurement - 1, centered about 0

# Get best square and it's corresponding accuracy
# reduce size of table
Simulations.Agg <- Simulations[, list(avg.acc.err = mean(dHat.acc), var.acc.err = var(dHat.acc)), by = list(TrapSpacing, CatchRadius, square)] 
# create temporary column for comparisons
Simulations.Agg[, `:=`(abs.best.acc = min(abs(avg.acc.err))), by = list(TrapSpacing, CatchRadius)] 
# remove rows that weren't the best
Simulations.Agg <- Simulations.Agg[abs.best.acc == abs(avg.acc.err), ] 
# remove temporary column
Simulations.Agg[, abs.best.acc:=NULL] 



library(ggplot2)
library(ggfortify)
library(plotly)
library(MASS)

# plot 3-d relationsip between TS, CR, and best avg accuracy
accplt <- Simulations.Agg[TrapSpacing > 0.5, ] %>% 
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~avg.acc.err, color = ~avg.acc.err) %>% 
  add_markers()
accplt


# Create repsonse and predictor variables. 
# resp is the most accurate average of the standardized densities by ring
# cr.by.ts is the ratio of (CR)/(TS) to measure how far we are from the line (2 * CR = TS)
# (the 2 should drop out as a constant in the model, but can be inluded later if desired)
Simulations.Agg[, `:=`(resp = avg.acc.err+1, cr.by.ts = CatchRadius/TrapSpacing)] 
# We want resp to be 1 (so log(resp) = 0), the bulk of point seem good
# Since the good points are roughly off the line CR/(2 TS) I'm going to look at the ratio
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point() + 
  geom_smooth(formula = y ~ x, method = "loess")
# Hot damn, that looks really nice. but there's still the question:
#### In general, what trap spacings do you want to use as the practioner without any prior knowedge of the mice?
# It seems like you want to use a medium trap spacing (see below)
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point(aes(colour = TrapSpacing)) + 
  geom_smooth(formula = y ~ x, method = "loess")


# Let's test the seemingly cubic relationship more formally
model <- lm(log(resp) ~ I(log(cr.by.ts)) + I(log(cr.by.ts)^2) + I(log(cr.by.ts)^3), Simulations.Agg)
summary(model)
autoplot(model)
# It looks relatively heavy-tailed, so maybe in the end we should use a robust linear model
# there's also a little upwards trend in the variance with fitted values and some systematic error in residuals, but generally ok

# test reduced model
# Alternatively, we can just look at the p-value for the x variable in the first model (same thing)
# CONCLUSION: the 'x' term does not contriubte significantly to the model's explainitory ability
reduced.model <- lm(log(resp) ~ I(log(cr.by.ts)^2) + I(log(cr.by.ts)^3), Simulations.Agg)
anova(reduced.model, model)  
# a 5th degree polynomial looks like it may fit the data even better, but it starts to lose interpertablity and compactness


# Look at the model! It looks so good :)
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point(aes(colour = (CatchRadius/TrapSpacing))) + 
  # scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ I(x^2) + I(x^3), method = "lm") # reduced model


# I know the fit goes up a little bit if you increase the polynomial order, but let's compare it to a loess regression
# It looks like the CI overlap, which is reassuring
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point(aes(colour = factor(round(TrapSpacing)))) + 
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(mapping = aes(linetype = "solid"), formula = y ~ I(x^2) + I(x^3), method = "lm") + # reduced model
  geom_smooth(mapping = aes(linetype = "dashed"), formula = y~x, method = "loess")


# Compare the model on a per-TS basis, looks interesting, helps to parse out differerent, marginal behaviors
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp), colour = factor(round(TrapSpacing)))) + 
  geom_point() + 
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ I(x^2) + I(x^3), method = "lm", se = FALSE) # reduced model
  # geom_smooth(formula = y ~ x, method = "loess", se = FALSE) # loess model



# try to see if there are any TS that give good results more often than others 
# (maybe recommend using those if you have no prior knowledge)
Simulations.Agg[, `:=`(num.close = sum(abs(resp - 1) < 0.05)), by = .(TrapSpacing)]
ggplot(Simulations.Agg, aes(x = TrapSpacing, y = num.close)) + geom_point()
# It looks like any TS can give good reuslts, probably is just they all give good results if CR is close to line 2CR = TS.







#### Variance ####

# We were just fitting accuarcy of our method, i.e. how close we can get to measuring the true density
# Now we're interested in giving the practitioner an estimated variance for their measurement.

varplt <- Simulations.Agg[TrapSpacing > 0.5, ] %>% 
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~var.acc.err, color = ~var.acc.err) %>% 
  add_markers()
varplt


Simulations.Agg[, `:=`(cr.by.ts = CatchRadius/TrapSpacing)]


ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(var.acc.err))) + 
  geom_point()
# the log-log plot of variance also looks cubic, which is reasuring, I think?

# add trend-line
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(var.acc.err))) + 
  geom_point() +
  geom_smooth(formula = y ~ x, method = "loess")
# definately looks cubic-ish

model <- lm(log(var.acc.err) ~ log(cr.by.ts) + I(log(cr.by.ts)^2) + I(log(cr.by.ts)^3), Simulations.Agg)
summary(model)
autoplot(model)
# Interestingly, it looks like the first term is significant in this model

# Plot use the cubic trend-line
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(var.acc.err), colour = factor(round(CatchRadius)))) + 
  geom_point() +
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ x + I(x^2) + I(x^3), method = "lm")

# In the model for mean accuary we wanted the repsonce to be near one (so the log is near zero), but now that we're
# modeling the variance, we want it to be as close to zero as possible (so the log is low or negative). From these plots
# we see that the log(variance) is REALLY high, there are actually only a few points with low variance. Interestingly, 
# those points that have low variance seems to also have low CR/TS ratio, indicating that if the traps are very spread out
# the measurements are consistantly poor, but maybe we can correct for that?

# And, again, let's compare the two methods on one plot for reassuance
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(var.acc.err))) + 
  geom_point() +
  geom_smooth(mapping = aes(linetype = "solid"), formula = y ~ x + I(x^2) + I(x^3), method = "lm") + 
  geom_smooth(mapping = aes(linetype = "dashed"), formula = y ~ x, method = "loess")
# even better overlap than the mean accuracy model













#### Square ####

# It would also be great if we could tell the practicioner which square to use 


sqrplt <- Simulations.Agg[TrapSpacing > 0.5, ] %>% 
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~square, color = ~square) %>% 
  add_markers()
sqrplt

# Flatten
ggplot(Simulations.Agg) + 
  geom_point(aes(x = TrapSpacing, y = CatchRadius, colour = factor(square)))

# Well, shoot. It doesn't look like there's any relationship (by eye)
# let's try the ratio again

ggplot(Simulations.Agg) + geom_point(aes(x = log(cr.by.ts), y = square))
# yeah... doesn't look like there's much of a trend that I can predict

summary(lm(square ~ log(cr.by.ts), Simulations.Agg))
# A very pittiful R2 value




# I tried removing all square 5's becuase they aren't square (they're rings) and I don't want 
# to use them for density prediction. This resulted in TS=0.5, CR=4.5 using square=4 and having
# a negative resp value (i.e. it's avg.acc.err = -1.34, so resp = -0.34), so the log is not defined.
# I guess, I know that there are dHat estimates that under-estimate, but this is the only point, 
# who's most accuate measurement (average by square) is THAT negative. I don't think this is a bug,
# just an issue with my understanding of the problem. 



# Accross the whole parameter space (cr x ts) what is the variance in the accuracy by square
# we want to know the expected error in accuracy for each square so that we can recommed a square 
# in the case where no informatiuon is supplied or available. We should also compare this to 
# another method, such as elbow ("inflection point") to see if we can supply a better rule.
# If we're not seeing any elbow, you're at a special case, and just use the larger square
# We may also be interested in the mean error 

# By square, what factor do we need to adjust the area by to get perfect density measurements (1)
# this is a function of TS (maybe cr...)


