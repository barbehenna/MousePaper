library(data.table)


# Load simulations where traps were allowed to overlap
Parameters <- fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Parameters.csv")
Simulations <-  fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Stats.csv")

# Merge in "True" density
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset", all.x = TRUE)

# Process Simulations
Simulations <- Simulations[, c("pHatDropNeg", "pHatZeroNeg") := NULL]
Simulations <- na.omit(Simulations)
Simulations[, `:=`(dHat.acc = (dHat - Density)/Density)] # Standardized density measurement - 1, centered about 0

# Get best square and it's corresponding accuracy
# reduce size of table
Simulations.Agg <- Simulations[, list(avg.acc.err = mean(dHat.acc)), by = list(TrapSpacing, CatchRadius, square)] 
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
  geom_point(aes(colour = factor(round(TrapSpacing)))) + 
  scale_color_brewer(palette = "Spectral") +
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




























