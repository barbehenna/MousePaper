library(Rcpp)
library(data.table)
library(ggplot2)
library(plotly)
library(pbapply)

sourceCpp("SimulationBackendFull.cpp")


#### Run Simulations ####
# 100 interations of (ts = 0.25->3.0, cr = 0.25->2.0, B = 3->6, D = 0.5->4.5) from takes about 20 minutes
# 1 interation of (ts = 0.25->360, cr = 0.25->5.0, B = 3->6, D = 0.5->4.5) from takes about 1.67 minutes

# Define parameters for various simulations 
TrapSpacing <- seq(from = 0.25, to = 3, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 2, by = 0.25)
Boarder <- seq(from = 3, to = 6, by = 1)
Density <- seq(from = 0.5, to = 5, by = 0.5)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)

# Make array specifying each simulation to be done
SimRuns <- rep(1:nrow(Parameters), each = 100)
SimRuns <- cbind(SimRuns, 1:length(SimRuns))

# Run simulations
set.seed(0) # I think this works
system.time(
Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
  s <- SimRuns[x, ]
  paramset <- Parameters[Parameters$paramset == s[1], ]
  res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 4, trueDensity = paramset$Density, nForages = 4)
  return(as.data.frame(res))
})
)
Simulations <- rbindlist(Simulations)
names(Simulations) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")


# save results
# fwrite(Parameters, "data/ParametersNewSample.csv")
# fwrite(Simulations, "data/SimulationsNewSample.csv")



# library(microbenchmark)
# microbenchmark(Simulations <- lapply(seq(nrow(SimRuns)), function(x) {
#   s <- SimRuns[x, ]
#   paramset <- Parameters[Parameters$paramset == s[1], ]
#   res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 4, trueDensity = paramset$Density, nForages = 4)
#   return(as.data.frame(res))
# }), times = 10)





#### Run analysis ####

# Load results
Parameters <- fread("data/ParametersNewSample.csv")
Simulations <- fread("data/SimulationsNewSample.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Computer true number of mice
Simulations[, `:=`(nSquares =  max(square)),  by = .(uuid)] #number of squares used in specific simulation run
Simulations[, `:=`(nMice = round(Density * (2 * TrapSpacing * nSquares + 2 * Boarder)^2))] # re-calculate number of mice used

# Standardize density estimate
Simulations[, `:=`(dHat.std = dHat/Density)]


# define a hopefully well-behaved subset of the data to explore
SimSubset <- na.omit(Simulations)
SimSubset <-  SimSubset[square == 3, ]
SimSubset <-  SimSubset[is.finite(nHat), ]
SimSubset <-  SimSubset[nHat > 0, ]



ggplot(SimSubset) + geom_density(aes(x = nHat/nMice))

# explore how close nHat is to the true number of mice
# nice and consistent modes
# pretty heavy tailed even the median isn't exactly at the mode
ggplot(SimSubset) +
  geom_density(aes(x = nHat/nMice,  colour = factor(Density))) +
  geom_vline(xintercept = median(SimSubset$nHat/SimSubset$nMice)) +
  geom_vline(xintercept = mean(SimSubset$nHat/SimSubset$nMice)) +
  xlim(0,2) +
  ggtitle("median and mean")


# estimate mode of numeric distribution
dist.mode <- function(vals, na.rm = FALSE) {
  dens <- density(vals, n = 2^20, na.rm = na.rm)
  return(dens$x[which.max(dens$y)])
}


# same as above, but drawring the estimated mode of nHat/nMice
# now that's pretty sexy... Consistent multiplicative error for any number of mice
ggplot(SimSubset) +
  geom_density(aes(x = nHat/nMice,  colour = factor(Density))) +
  geom_vline(xintercept = dist.mode(SimSubset$nHat/SimSubset$nMice)) +
  xlim(0,2) +
  ggtitle("mode")

dist.mode(SimSubset$nHat/SimSubset$nMice)
mean(SimSubset$nHat/SimSubset$nMice)
median(SimSubset$nHat/SimSubset$nMice)

SimSubset[, `:=`(nHat.adj.mode = nHat/dist.mode(nHat/nMice), 
                 nHat.adj.mean = nHat/mean(nHat/nMice),
                 nHat.adj.median = nHat/median(nHat/nMice))]

# Check out nHat accuracy by adjustment (mode = 1, mean = 1, or median = 1)
ggplot(SimSubset) + geom_density(aes(x = nHat.adj.mode/nMice, colour = factor(Density))) + xlim(0,2)
ggplot(SimSubset) + geom_density(aes(x = nHat.adj.mode/nMice)) + xlim(0,2)


# Now that we've adjusted nHat, how are we doing on density estimates?
ggplot(SimSubset) + geom_density(aes(x =  (nHat.adj.mode/aHat), colour = factor(Density))) + xlim(0,24)

# Scaled by true density
ggplot(SimSubset) + geom_density(aes(x =  nHat.adj.mode/(aHat*Density), colour = factor(Density))) + xlim(0,8)
# They're all off by a constant factor? 


SimSubset[, dist.mode(nHat.adj.mode/(aHat*Density))]
SimSubset[, mean(nHat.adj.mode/(aHat*Density))]
SimSubset[, median(nHat.adj.mode/(aHat*Density))]

ggplot(SimSubset) + 
  geom_density(aes(x =  nHat.adj.mode/(aHat*Density), colour = factor(Density))) +
  geom_vline(xintercept = SimSubset[, dist.mode(nHat.adj.mode/(aHat*Density))]) + 
  xlim(0,8)


# Adjust density estimate by mode?
SimSubset[, `:=`(dHat.adj = nHat.adj.mode / (aHat*dist.mode(nHat.adj.mode/(aHat*Density))))]

SimSubset[, dist.mode(nHat.adj.mode/(aHat*Density))] #factor scaling aHat

# plot
ggplot(SimSubset) + geom_density(aes(x = dHat.adj/Density, colour = factor(Density))) + xlim(0,2)




## From scratch now:

SimSubset <- na.omit(Simulations)
SimSubset <-  SimSubset[square == 3, ]
SimSubset <-  SimSubset[is.finite(nHat), ]
SimSubset <-  SimSubset[nHat > 0, ]

nHat.scale = 1/dist.mode(SimSubset$nHat/SimSubset$nMice)
aHat.scale = dist.mode((SimSubset$nHat*nHat.scale) / (SimSubset$aHat*SimSubset$Density))


# illustrate nHat.adj accuracy
ggplot(SimSubset) + geom_density(aes(x = nHat*nHat.scale/nMice)) + xlim(0,2)

SimSubset[, `:=`(nHat.new = nHat * nHat.scale, aHat.new = aHat * aHat.scale)]
ggplot(SimSubset) + geom_density(aes(x = nHat.new/aHat.new, colour = factor(Density))) + xlim(0,7)
ggplot(SimSubset) + geom_density(aes(x = nHat.new/(aHat.new * Density), colour = factor(Density))) + xlim(0,2)


# Mode -> highest probability of having near zero error in density prediction
# Mean -> average is closest to the true density if you keep repleating the expirament and the expiraments are identical
# Median -> median is closest to the true density if you keep repleating the expirament and the expiraments are identical (more robust, slower convergence?)


dHat.scale = nHat.scale/aHat.scale ## This is really close to 1!
# what does the actual  dHat look like? We haven't checked yet...
ggplot(SimSubset) + geom_density(aes(x = dHat/Density)) + xlim(0,2)

# hot damn!
# Let's try to tease appart that heavy tail

ggplot(SimSubset) + geom_density(aes(x = dHat/Density, colour = factor(Density))) + xlim(0,2) # nope
ggplot(SimSubset) + geom_density(aes(x = dHat/Density, colour = factor(Boarder))) + xlim(0,2)  # nope
#  So  it doesn't  depend on field size or number of mice, nice!

ggplot(SimSubset) + geom_density(aes(x = dHat/Density, colour = factor(TrapSpacing))) + xlim(0,5)
# That  looks like it fixes it. Notice how symmetric the distribution of the TS = 3.0 is 

ggplot(SimSubset) + geom_density(aes(x = dHat/Density, colour = factor(round(CatchRadius)))) + xlim(0,4)
# looks like that small bump on the left is due to small CatchRadii



## Normal dist?
ggplot(SimSubset, aes(sample = dHat/Density, colour = factor(TrapSpacing))) + stat_qq()
# not too much, very heavy tailed


# Look at the density accuracy by ts and cr
ggplot(SimSubset[, .(dHat.acc = mean(dHat/Density)), by = .(TrapSpacing, CatchRadius)], aes(x = TrapSpacing, y = CatchRadius, colour = factor(cut(dHat.acc, 15)), size = 4)) + 
  geom_point()

# sample plot but in 3d so I don't need to round means
plot_ly(data = SimSubset[, .(dHat.acc = mean(dHat/Density)), by = .(TrapSpacing, CatchRadius)], 
        x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.acc, color = ~dHat.acc) %>% 
  add_markers()

# same as above but with size proportional to the variance
plot_ly(data = SimSubset[, .(dHat.mean = mean(dHat/Density), dHat.var = var(dHat/Density)), by = .(TrapSpacing, CatchRadius)], 
        x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~dHat.mean, size = ~dHat.var,  marker = list(symbol = 'circle', sizemode = 'diameter')) %>% 
  add_markers()


# Look at a sinlge point
ggplot(SimSubset[TrapSpacing==2 & CatchRadius == 1.5], aes(x = dHat/Density)) + geom_density() # looks nice :)
ggplot(SimSubset[TrapSpacing==2 & CatchRadius == 1.5], aes(sample = dHat/Density)) + stat_qq() + stat_qq_line() # Heavy tails... t-dist?

# # compare dist with normal and t(df=1)
# # pretty close to normal
# Temp.dat <- SimSubset[TrapSpacing==2 & CatchRadius == 1.5]
# Temp.dat[, `:=`(samples = (dHat/Density - mean(dHat/Density))/sd(dHat/Density))]
# Temp.dat[, `:=`(rnorm.samples = rnorm(nrow(Temp.dat)))]
# Temp.dat[, `:=`(rt1.samples = rt(nrow(Temp.dat), df = 1))]
# ggplot(Temp.dat) + 
#   geom_density(aes(x = samples, col = "data")) + 
#   geom_density(aes(x = rnorm.samples, colour = "rnorm")) + 
#   geom_density(aes(x = rt1.samples, colour = "rt1")) + 
#   xlim(-6,6)
