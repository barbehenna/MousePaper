#### Libraries ####

library(Rcpp)
library(data.table)
library(tidyr)
library(parallel)
library(ggplot2)
library(plotly)
library(pbapply)

sourceCpp("SimulationBackendFull.cpp")


# For reference and playing, here's my function to find the sample mode
dens.mode <- function(x, na.rm = FALSE) {
  dens <- density(x = x, na.rm = na.rm, n = 2^20)
  return(dens$x[which.max(dens$y)])
}


#### Run Simulations ####

# Define parameters for various simulations
TrapSpacing <- seq(from = 0.25, to = 6, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 6, by = 0.25)
Boarder <- seq(from = 3, to = 6, by = 6)
Density <- seq(from = 0.5, to = 5, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)

# Make array specifying each simulation to be done
SimRuns <- rep(1:nrow(Parameters), each = 1000)
SimRuns <- cbind(SimRuns, 1:length(SimRuns))

# Run simulations
# system.time(
# Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
#   s <- SimRuns[x, ]
#   paramset <- Parameters[Parameters$paramset == s[1], ]
#   res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 8, trueDensity = paramset$Density, nForages = 4)
#   return(as.data.frame(res))
# })
# )
# Simulations <- rbindlist(Simulations)
# names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


system.time(
  Simulations <- mclapply(seq(nrow(SimRuns)), function(x) {
    s <- SimRuns[x, ]
    paramset <- Parameters[Parameters$paramset == s[1], ]
    res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 8, trueDensity = paramset$Density, nForages = 4)
    return(as.data.frame(res))
  }, mc.cores = 2, mc.set.seed = TRUE)
)
Simulations <- rbindlist(Simulations)
names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


# save results
fwrite(Parameters, "data/ParametersNewSample3.csv")
fwrite(Simulations, "data/SimulationsNewSample3.csv")








#### Run analysis ####

# Load results
Parameters <- fread("data/NewBackend-ParametersFull_B3.csv")
Simulations <- fread("data/NewBackend-SimulationsFull_B3.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Computer true number of mice
Simulations[, `:=`(nSquares = max(square)), by = .(uuid)] # number of squares used in specific simulation run
Simulations[, `:=`(nMice = round(Density * (2 * TrapSpacing * nSquares + 2 * Boarder)^2))] # re-calculate number of mice used


square.results <- Simulations[!is.na(dHat) & is.finite(dHat) & dHat > 0, .(mean.dHat = mean(dHat / Density, na.rm = TRUE), median.dHat = median(dHat / Density, na.rm = TRUE), mode.dHat = dens.mode(dHat / Density, na.rm = TRUE), var.dHat = var(dHat / Density, na.rm = TRUE)), by = .(square)]

square.results <- square.results %>% 
  gather(key = statistic, value = value, -square) %>%
  arrange(square)

ggplot(square.results) + geom_point(aes(x = square, y = value, colour = statistic), size = 2) + geom_path(aes(x = square, y = value, colour = statistic))

# define a hopefully well-behaved subset of the data to explore
SimSubset <- na.omit(Simulations)
SimSubset <- SimSubset[square == 3, ]
SimSubset <- SimSubset[is.finite(nHat), ]
SimSubset <- SimSubset[nHat > 0, ]



# How good are we doing with our density estimates?
ggplot(SimSubset) + geom_density(aes(x = dHat, colour = factor(Density))) + xlim(0, 8)

# Damn that's pretty, but how does it look when we standardize our estimates?
# dHat/Density = 1 if we are accurate
ggplot(SimSubset) + geom_density(aes(x = dHat / Density, colour = factor(Density))) + xlim(0, 2)

# So it looks like regardless of how many mice we have (all things considered other than density)
# our measurements are very good (the tails on our standardized densities aren't symmetric)

# Why are the tails not symmetric? It could just be becuase they're bounded on the left
ggplot(SimSubset) + geom_density(aes(x = dHat / Density, colour = factor(CatchRadius))) + xlim(0, 2)
# So small CatchRadius is why we see heavy tails on the left, but they almost over correct for that issue/observation
ggplot(SimSubset) + geom_density(aes(x = dHat / Density, colour = factor(TrapSpacing))) + xlim(0, 2)
# And low TrapSpacing is why we see heavy tails on the right

# Look just at what should be a "good" combination
ggplot(SimSubset[TrapSpacing == 5 & CatchRadius == 4]) + geom_density(aes(x = dHat / Density))
# yeah... that looks pretty nice :)


# Ok, let's parse out this more in 3d
dHatSummary <- SimSubset[, .(dHat.mean = mean(dHat / Density), dHat.var = var(dHat / Density)), by = .(TrapSpacing, CatchRadius)]

# Using log10 dHat.mean for to center the colors around the correct value (dHat/Density = 1 => log(dHat/Density) = 0)
p <- dHatSummary %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~ log10(dHat.mean)) %>%
  add_markers()
p

# We have variance information, let's plot that too with size
p <- dHatSummary %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~ log10(dHat.mean), size = ~dHat.var, marker = list(symbol = "circle", sizemode = "diameter")) %>%
  add_markers()
p

# Looks good, but let's crop it a little so that it's easier to read
p <- dHatSummary[TrapSpacing > 0.25 & CatchRadius > 0.25] %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~ log10(dHat.mean), size = ~dHat.var, marker = list(symbol = "circle", sizemode = "diameter")) %>%
  add_markers()
p

# Just for fun, look at only those points whose mean dHat/Density value is within 5% of correct
p <- dHatSummary[abs(dHat.mean - 1) < 0.05] %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~ log10(dHat.mean), size = ~dHat.var, marker = list(symbol = "circle", sizemode = "diameter")) %>%
  add_markers()
p


# save any of these plots as interactive html's using
# htmlwidgets::saveWidget(widget = as_widget(p), file = "~/Desktop/3d-dHat_plot.html")




# So... dHat is a pretty good estimate of the true density for most trap spacings



# There is a true catch radius, how can we find that or how do we run an expirament around that?

# Repeat at a few differnt trap spacings (along one of the lines below)
ggplot(dHatSummary[is.element(CatchRadius, c(1, 2, 3, 4, 5))]) +
  geom_hline(yintercept = 0) +
  geom_path(aes(x = TrapSpacing, y = log(dHat.mean), colour = factor(CatchRadius), linetype = factor(CatchRadius))) +
  ggtitle("Samples of Cross-sections From 3d Plot", subtitle = "Goal is dHat.mean = 1, drawn in black")

# It looks like we expect: for each curve there's a region where it performs quite well and in the middle,
# but, in general, it looks something like a downward cubic function. To the left of the good region, the
# estimator is consistently way too high and to the right it is too low. Looking at how it changes by Catch Radius
# we see that the good region shifts to the right as Catch Radius increases. This makes sense because we expect
# the method to give the best results when CT ~= TS/2. A much larger set of Catch Radius and Trap Spacing
# pairs may be needed to define exactly where the method works best (by eye I'd guess 2CR <= TS <= 3CR)








