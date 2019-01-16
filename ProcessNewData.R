library(data.table)


# Load simulations where traps were allowed to overlap
Parameters <- fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Parameters.csv")
Simulations <-  fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Stats.csv")


# Set key of each table
# set the ON clause as keys of the tables
# setkey(Parameters, paramset)
# setkey(Simulations, paramset)

# Merge in "True" density
# Simulations <- merge(x = Simulations, y = Parameters, all.x = TRUE)
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset", all.x = TRUE)


# Process Simulations
Simulations <- Simulations[, c("pHatDropNeg", "pHatZeroNeg") := NULL]
Simulations <- na.omit(Simulations)
# Simulations$dHat.acc <- (Simulations$dHat - Simulations$Density) / Simulations$Density
Simulations[, `:=`(dHat.acc = (dHat - Density)/Density)]

# Get best square and it's corresponding accuracy
Simulations.Agg <- Simulations[, list(avg.acc.err = mean(dHat.acc)), by = list(TrapSpacing, CatchRadius, square)]
Simulations.Agg[, `:=`(abs.best.acc = min(abs(avg.acc.err))), by = list(TrapSpacing, CatchRadius)]
Simulations.Agg <- Simulations.Agg[abs.best.acc == abs(avg.acc.err), ]
Simulations.Agg[, abs.best.acc:=NULL]



library(ggplot2)
library(ggfortify)
library(plotly)
library(MASS)

accplt <- Simulations.Agg[TrapSpacing > 0.5, ] %>% 
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~avg.acc.err, color = ~avg.acc.err) %>% 
  add_markers()
accplt


Simulations.Agg[, `:=`(resp = avg.acc.err+1, cr.by.ts = CatchRadius/TrapSpacing)]
# We want resp to be 1 (so log(resp) = 0), the bulk of point seem good
# Since the good points are roughly off the line CR/(2 TS) I'm going to look at the ratio
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + geom_point() + geom_smooth(formula = y ~ x, method = "loess")
# Hot damn, that looks really nice. but there's still the question:
#### In general, what trap spacings do you want to use as the practioner without any prior knowedge of the mice
# It seems like you want to use a medium trap spacing (see below)
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + geom_point(aes(colour = TrapSpacing)) + geom_smooth(formula = y ~ x, method = "loess")


model <- lm(log(resp) ~ I(log(cr.by.ts)) + I(log(cr.by.ts)^2) + I(log(cr.by.ts)^3), Simulations.Agg)
summary(model)
autoplot(model)

# test reduced model
# Alternatively, we can just look at the p-value for the x variable in the first model (same thing)
# CONCLUSION: the 'x' term does not contriubte significantly to the model's explainitory ability
reduced.model <- lm(log(resp) ~ I(log(cr.by.ts)^1) + I(log(cr.by.ts)^2) + I(log(cr.by.ts)^3) + I(log(cr.by.ts)^4) + I(log(cr.by.ts)^5), Simulations.Agg)
anova(reduced.model, model)  
# a 5th degree polynomial looks like it may fit the data even better, but it starts to lose interpertablity and compactness


# Look at the model! It looks so good :)
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point(aes(colour = factor(round(TrapSpacing)))) + 
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ I(x^2) + I(x^3), method = "lm") # reduced model


# I know the fit goes up a little bit if you increase the polynomial order, but let's compare it to a loess regression
# It looks like the CI overlap, which is reassuring
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp))) + 
  geom_point(aes(colour = factor(round(TrapSpacing)))) + 
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ I(x^2) + I(x^3), method = "lm") + # reduced model
  geom_smooth(formula = y~x, method = "loess")


# Compare the model on a per-TS basis, looks interesting, helps to parse out differerent, marginal behaviors
ggplot(Simulations.Agg, aes(x = log(cr.by.ts), y = log(resp), colour = factor(round(TrapSpacing)))) + 
  geom_point() + 
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(formula = y ~ I(x^2) + I(x^3), method = "lm", se = FALSE) # reduced model
  # geom_smooth(formula = y ~ x, method = "loess", se = FALSE) # loess model



Simulations.Agg[, `:=`(num.close = sum(abs(resp - 1) < 0.05)), by = .(TrapSpacing)]
ggplot(Simulations.Agg, aes(x = TrapSpacing, y = num.close)) + geom_point()





























