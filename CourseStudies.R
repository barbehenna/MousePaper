# Copyright (c) 2017 Alton Barbehenn

# This script loads a course, set of simulation parameters and their resulting statistics 
# (From MultiParameterSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)
library(ggplot2)
library(pbapply)
library(parallel)

Sim <- read.csv("data/20180115_185702_Stats.csv", header = TRUE, row.names = NULL)
Parameters <- read.csv("data/20180115_185702_Parameters.csv", header = TRUE, row.names = NULL)

Sim <- merge(Sim, Parameters, by = "paramset")

# Aggregate the statistics with the same parameters
# I chose median at this point because of it's robustness
# AggStats <- aggregate(cbind(dHat, nHat, pHat, pHatZeroNeg, pHatDropNeg, aHat) ~ square + TrapSpacing + FieldSize + CatchRadius + density, data = Sim, median)

# Using the data.table structure might be faster and more readable
# OMITING NA VALUES!!
Sim <- data.table(Sim)
AggStats <- Sim[, .(med_dHat = median(na.omit(dHat)),
                    avg_dHat = mean(na.omit(dHat)),
                    std_dHat = sd(na.omit(dHat)),
                    med_nHat = median(na.omit(nHat)),
                    avg_nHat = mean(na.omit(nHat)),
                    std_nHat = sd(na.omit(nHat)),
                    med_pHat = median(na.omit(pHat)),
                    avg_pHat = mean(na.omit(pHat)),
                    std_pHat = sd(na.omit(pHat)),
                    med_pHatZeroNeg = median(na.omit(pHatZeroNeg)),
                    avg_pHatZeroNeg = mean(na.omit(pHatZeroNeg)),
                    std_pHatZeroNeg = sd(na.omit(pHatZeroNeg)),
                    med_pHatDropNeg = median(na.omit(pHatDropNeg)),
                    avg_pHatDropNeg = mean(na.omit(pHatDropNeg)),
                    std_pHatDropNeg = sd(na.omit(pHatDropNeg)),
                    med_aHat = median(na.omit(aHat))),
                # by = .(square, TrapSpacing, FieldSize, CatchRadius, Density)]
                by = .(paramset, square)]

AggStats <- merge(AggStats, Parameters, by = "paramset")

# Want to understand what values affect the estimate of density
model <- lm(avg_dHat ~ square + TrapSpacing + FieldSize + CatchRadius + Density, data = AggStats)
summary(model)

# Looks like we can ignore the effects of Trap Spacing and Catch Radius on our density estimation
model <- lm(avg_dHat ~ square + CatchRadius + Density, data = AggStats)
summary(model)

# Look at the correlation between density and dHat, colored by ring/square
p <- ggplot(data = AggStats) + geom_boxplot(mapping = aes(x = factor(density), y = med_dHat))
p
cor.test(x = AggStats$Density, y = AggStats$med_dHat, method = "pearson")
cor.test(x = AggStats$Density, y = AggStats$med_dHat, method = "kendall")

# So they appear to be correlated but look systematically biased
# We can see this by looking at the peaks in the density plots
p <- ggplot(data = AggStats) + geom_density(mapping = aes(x = med_dHat, colour = factor(square)))
p





p <- ggplot(data = AggStats[AggStats$square < 2]) + geom_point(mapping = aes(x = density, y = med_dHat, color = factor(CatchRadius)))
p


p <- ggplot(data = AggStats[AggStats$square < 2]) + geom_density(mapping = aes(med_dHat, colour = factor(CatchRadius)))
p






##################################################################################
# Is the distribution of densitiy estimates normal for each ring?
##################################################################################

# If we're looking at our simulations by parameter set, the only possible thing they 
# very on (within a paramset) is the square. It looks like square 1 is a large reason 
# why the distribution of dHat not normal. 
# Regardless, our data does not look normally distributed.

for (i in sample(seq(336), 10)) {
  qqnorm(Sim$dHat[Sim$paramset == i], main = paste("paramset =", i))
  qqline(Sim$dHat[Sim$paramset == i])
  p <- ggplot(Sim[Sim$paramset == i]) + geom_point(stat = "qq", aes(sample = dHat, colour = factor(square)))
  plot(p)
}

# Check for parallelism in qqlines of squares 2 and 3 <==> do they have the same distribution?


##################################################################################
# Are the density estimates consistently wrong?
##################################################################################


median(Sim$dHat[Sim$paramset == 126 & Sim$square <= 3], na.rm = TRUE)
# Is this always true? No.
# error <- NULL
# # This foor loop can easily be done using vector arithmatic
# for (i in unique(Sim$UniqueID)) {
#   med <- median(Sim$dHat[Sim$UniqueID == i & Sim$square <= 3], na.rm = TRUE)
#   den <- Sim$Density[Sim$UniqueID == i][1] #should only be one value
#   error <- rbind(error, c(den, med/den, med-den))
# }
# error <- as.data.frame(error)
# names(error) <- c("den", "perc", "abs")

error <- lapply(unique(Sim$UniqueID), function(x){
  avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == i][1]
  tmp <- data.frame(den, med/den, med-den)
  return(tmp)
})
error <- rbindlist(error)
names(error) <- c("den", "perc", "abs")


ggplot(error[abs(error$abs) < 1,]) + geom_density(aes(x = abs, colour = factor(den)))
ggplot(error[abs(1 - error$perc) < 1,]) + geom_density(aes(x = perc, colour = factor(den)))


comb <- gtools::combinations(n = length(unique(error$den)), r = 2, v = unique(error$den))
pval <- NULL
for (i in seq(NROW(comb))) {
  test <- ks.test(x = error$perc[error$den == comb[i,1]], y = error$perc[error$den == comb[i,2]])
  pval <- c(pval, test$p.value)
}
print(min(pval)) 
#So the distributions of the percent error in the density estimates 
#seem to come from the same distribution regardless of the true density

mean(error$perc, na.rm = TRUE) #mean error of density estimates
median(error$perc, na.rm = TRUE) #different from mean implies a skewed distribution? Or just a small sample? (Or something else?)

#MAYBE?? THEY WERE THE SAME DISTRIBUTION BEFORE THE SAMPLE SIZE INCREASED, BUT NOW THEY'RE NOT
#THEY LOOK BASICALLY THE SAME. (MAYBE THE SAMPLE IS TOO LARGE?)


# Let's try bootstrapping the pvalues (I think I'm doing this correctly)
# use suppressWarnings(expr) if there are a lot of warnings related to ties
bootstrap_pval <- gtools::combinations(n = length(unique(error$den)), r = 2, v = unique(error$den))
dimnames(bootstrap_pval)[[2]] <- list("dens1", "dens2")
pval <- pblapply(seq(nrow(comb)), function (x) {
  unlist(pblapply(seq(100000), function(y) {
    ks.test(x = sample(error$perc[error$den == bootstrap_pval[x,1]], 100), 
            y = sample(error$perc[error$den == bootstrap_pval[x,2]], 100), 
            exact = TRUE)$p.value
  }))
})

# Calculate the lower bound for a 95% CI
low_95 <- lapply(pval, function(x) {
  sort(x)[floor(0.025 * length(x))]
})
bootstrap_pval <- cbind(bootstrap_pval, low_95=unlist(low_95))

# Calculate the upper bound for a 95% CI
high_95 <- lapply(pval, function(x) {
  sort(x)[floor(0.975 * length(x))]
})
bootstrap_pval <- cbind(bootstrap_pval, high_95=unlist(high_95))

# Estimate the value of the parameter
guess <- lapply(pval, function(x) {
  median(x, na.rm = FALSE)
})
bootstrap_pval <- cbind(bootstrap_pval, guess=unlist(guess))

# Probability of rejecting the NULL with the Kolmogorov-Smirnov Test (alpha = 0.05)
# the complement of this value is the probability that we conclude that the two percent 
# error distributions are coming from the same distributions, i.e. our method is density invarient 
# This allows us to figure out a constant to multiply our density estimate by to derive the true value.
p_reject <- lapply(pval, function(x) {
  sum(x < 0.05)/length(x)
})
bootstrap_pval <- cbind(bootstrap_pval, p_reject=unlist(p_reject))

# So it seems like we can't diffinatively say that the distributions of 
# the errors are different or not, but it seems like they're the same


##################################################################################
# Visualize metrics
##################################################################################

mode_alton <- function(x) {
  d <- density(x, na.rm=TRUE)
  return(d$x[which(d$y == max(d$y))])
}

plot(density(error$perc, na.rm = TRUE), xlim = c(-1, 3))
abline(v = median(error$perc, na.rm = TRUE), col = "red")
abline(v = mean(error$perc, na.rm = TRUE), col = "green")
abline(v = mode_alton(error$perc), col = "blue")


##################################################################################
# Is our measure of density "tight" and normal for a specific trap spacing?
##################################################################################

Sim <- Sim[TrapSpacing == 1.0] # Pretty sure this trap spacing is known to work well

for (i in unique(Sim$UniqueID)) {
  med <- median(Sim$dHat[Sim$UniqueID == i & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == i][1] #should only be one value
  error <- rbind(error, c(den, med/den, med-den))
}
error <- as.data.frame(error)
names(error) <- c("den", "perc", "abs")


ggplot(Sim) + geom_density(aes(x = dHat, colour = factor(Density)))


##################################################################################
# BY trap spacing AND Catch Radius (for a fixed catch radius), 
# what is the distribution of error in density estimates? 
# Which ones yield better results?
##################################################################################

ncores <- detectCores()
cl <- makeCluster(ncores-1, type = "FORK")

error <- pblapply(unique(Sim$UniqueID), cl = cl, function(x) {
  avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == x][1]
  ts <- Sim$TrapSpacing[Sim$UniqueID == x][1]
  cr <- Sim$CatchRadius[Sim$UniqueID == x][1]
  tmp <- data.frame(den, avg-den, avg/den, ts, cr)
  return(tmp)
})

stopCluster(cl)
error <- rbindlist(error)
names(error) <- c("den", "abs", "perc", "TrapSpacing", "CatchRadius")

ggplot(data = error[error$CatchRadius == 0.5]) + geom_density(aes(x = perc, colour = factor(TrapSpacing))) + xlim(-1,3)


# Attempt one at bi-variate analysis
error<-data.table(error)
error_proc <- error[, .(perc = mean(perc)), by = .(TrapSpacing, CatchRadius)]
ggplot(data = error_proc, mapping = aes(x = TrapSpacing, y = CatchRadius))+geom_point(aes(colour = perc))

# Looks like the predictions get better as the catch radius increases
tmp <- error[error$TrapSpacing == 6,]
ggplot(data = tmp) + geom_density(mapping = aes(x = perc, colour = factor(CatchRadius))) + xlim(-1,3)



