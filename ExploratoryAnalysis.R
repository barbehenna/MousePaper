# This is a place for exploritory analyses

#### Libraries ####

library(data.table)
library(tidyr)
library(ggplot2)
library(plotly)
library(pbapply)


# For reference and playing, here's my function to find the sample mode
dens.mode <- function(x, na.rm = FALSE) {
  dens <- density(x = x, na.rm = na.rm, n = 2^20)
  return(dens$x[which.max(dens$y)])
}


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
  gather(key = statistic, value = value, -square)

ggplot(square.results, aes(x = square, y = value, colour = statistic)) + geom_point(size = 2) + geom_line()

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





#### Experimenting ####

Simulations <- Simulations[!is.na(dHat) & is.finite(dHat) & dHat>0, ]
Temp.Sim <- Simulations[, .(median.dHat = median(dHat/Density), var.dHat = var(dHat/Density)), by = .(square, TrapSpacing, CatchRadius)]


ggplot(Temp.Sim, aes(x = square, y = median.dHat, colour = factor(CatchRadius), shape = factor(TrapSpacing))) + geom_point() + geom_line() + geom_hline(yintercept = 1)
ggplot(Temp.Sim, aes(x = square, y = log(median.dHat), colour = factor(CatchRadius), shape = factor(TrapSpacing))) + geom_point() + geom_line() + geom_hline(yintercept = 0)
# too messy

ggplot(Temp.Sim, aes(x = square, y = median.dHat, colour = cut(log(CatchRadius/TrapSpacing), breaks = 6))) + geom_point() + geom_hline(yintercept = 1)
# So it looks like if  log(CR/TS) is within (-1.06, 1.06) ie CR/TS is within 0.35 and 2.88 (or it's probably 1/e to e) then our model
# is pretty reliable for most squares. When the ratio is outside that, the estimator doesn't work very well. 
# If you zoom in, you can see that there's still variation when the abs(log(cr/ts)) <= 2, it's pretty tight cutting off at +-one.

# There are too many points, let's just look at one good slice of TS
ts = 4
ggplot(Temp.Sim[TrapSpacing == ts, ], aes(x = square, y = median.dHat, colour = factor(CatchRadius))) + geom_point() + geom_line() + geom_hline(yintercept = 1) + ggtitle(paste("slice of plot by ts =", ts))
# so in some cases there is an observable elbow, but it doesn't look like it happens at the same places for every cr
# and it doesn't necessarily happen at all.



tmp <- Temp.Sim[TrapSpacing==4 & CatchRadius==2.5]
ggplot(tmp, aes(x = square, y = median.dHat)) + geom_point() + geom_smooth(method = "loess")
ggplot(tmp, aes(x = square, y = var.dHat)) + geom_point() + geom_smooth(method = "loess")




tmp <- tmp[sample(8)] # assume .DS is randomly ordered
tmp <- tmp[order(square)]

# method 1 of picking the square:
# sort by square, look for a change in concavity. As square increases, approximate the second order derivative
# and look to see if it changes. Becuase the second order derivative needs three points to define it, we can 
# use left recursive approx, right recursive approx, or the centered method (maybe average them all?).
# I'll  start with just the left sided method for ease of implementation.

# Since we're only taking steps of size one, the running difference approximates the slope (f(x+h)-f(x))/h, h=1
# Similarly, to approximate the second order deriavaticve, it suffices to look at the running difference of those
# differences (f'(x+h)-f'(x))/h, h=1, f' approximated above.

# We're essentially fitting a parabolla at every set of three consecutive points. Because we expect dHat to increase
# at the edge, maybe it makes more sense to be conservative and only use the right parabolla. Of course, if we do so,
# we will never select square 7 or 8. 

tmp <- tmp[sample(8)] # assume .DS is randomly ordered
tmp <- tmp[order(square)]

first.order <- diff(tmp$median.dHat) #approx first order derivative at squares 1-7
second.order <- diff(first.order) #approx second order derivative at squares 1-6
# I THINK if I sort the squares 8,7,...,1, this method will give me the second order derivatives at squares 3-8


Temp.Sim[, `:=`(crts = paste(CatchRadius, TrapSpacing, sep = "_"))]
sims.list <- split(x = Temp.Sim, f = Temp.Sim$crts)
sims.list <- lapply(sims.list, function(x) {
  x <- x[order(square)]
  second.order <- diff(diff(x[, median.dHat]))
  return(min(c(which(second.order > 0.05)[1], 7), na.rm = TRUE))
})
sims.list <- stack(sims.list)
names(sims.list) <- c("elbow.est","coord")

Temp.Sim <- merge(x = Temp.Sim, y = sims.list, by.x = "crts", by.y = "coord")

ggplot(Temp.Sim[TrapSpacing==3], aes(x = square, y = median.dHat, colour = factor(CatchRadius))) + geom_point() + geom_line()



Simulations <- Simulations[uuid %in% sample(x = max(uuid), size = 20)]
SimSplit <- split(Simulations, Simulations$uuid)

tmp <- SimSplit[[2]]

ggplot(tmp, aes(x = square, y = dHat)) + geom_point() + geom_line() 

tmp <- tmp[order(-square)]
sec.der <- rev(diff(diff(tmp[, dHat]))) # it was approx second order derivatives at squares 8,7,...,3, the rev makes it from 3,4,...,7,8
which(abs(sec.der) < 0.05)[1]
max(c(2, which(abs(sec.der) < 0.05)[1]), na.rm = TRUE)

lapply(SimSplit, function(x) {
  x <- x[order(-square)]
  sec.der <- rev(diff(diff(x[, dHat])))
  return(c("uuid" = x[1,uuid], "good.square" = max(c(2, which(abs(sec.der) < 0.05)[1]), na.rm = TRUE)))
})

