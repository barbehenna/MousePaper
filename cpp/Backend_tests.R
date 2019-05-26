# Test script for the full simulation backend


#### Libraries ####
library(Rcpp)
library(data.table)
library(ggplot2)


# Load full simulation backend
Rcpp::sourceCpp('SimulationBackendFull.cpp')




#### Run Simulations ####

# Test Trap generation
nrow(GenTraps(nSquares = 8, trapSpacing = 0.5)) # correct number of traps
nrow(GenTraps(nSquares = 4, trapSpacing = 0.5))
nrow(GenTraps(nSquares = 1, trapSpacing = 0.5))

plot(GenTraps(nSquares = 8, trapSpacing = 1))
points(GenTraps(nSquares = 6, trapSpacing = 1), col = 'green')
points(GenTraps(nSquares = 4, trapSpacing = 1), col = 'red') 
points(GenTraps(nSquares = 2, trapSpacing = 1), col = 'blue') # rings symmertric from middle

range(GenTraps(nSquares = 8, trapSpacing = 1)) # correct range



# Test mouse foraging
nsamples = 1000
plot(GenMouse(nforages = nsamples, fieldSize = 13.5)) # forages as a bi-variate normal
nortest::lillie.test(rnorm(nsamples)) # calabrate
nortest::lillie.test(runif(nsamples))
nortest::lillie.test(GenMouse(nforages = nsamples, fieldSize = 13.5)[,1]) # test normality of x coordinates
nortest::lillie.test(GenMouse(nforages = nsamples, fieldSize = 13.5)[,2]) # test normality of y coordinates




# Test mice home-ranges
nsamples = 1000
plot(t(sapply(seq(nsamples), function(x) GenMouse(nforages = 1, fieldSize = 13.5))))
abline(v = c(-6.75, 6.75), col = 'red') # add field boarder
abline(h = c(-6.75, 6.75), col = 'red')
# the reason some points are outside is because this is the first forage location, not the home,
# so they can be +- a few sigma from their home, most are inside square. If you modify the function
# to return just NumericVector of homes, it works perfectly.




# Test is mouse caught
mouse <- GenMouse(4, fieldSize = 13.5)
Traps <- GenTraps(nrings = 8, trapspacing = 0.5)

# day 1 (increase mouse index 1 for other days)
dxdy <- t(apply(Traps, MARGIN = 1, FUN = function(x) {return(abs(x - mouse[2,]))}))
which(apply(dxdy, MARGIN = 1, FUN = function(x) {max(x) < 0.3})) - 1 # which trap(s) catch a mouse (-1 to match c++ index)
Traps[apply(dxdy, MARGIN = 1, FUN = function(x) {max(x) < 0.3}), ]
isCaught(forages = mouse, Traps = Traps, catchRadius = 0.3)
# iterated enough times to check outputs that I'm happy the answers match for isCaught()




# test raw data simulation
test <- GenTrapData(trapSpacing = 0.5, catchRadius = 0.3, boarder = 3, nSquares = 8, trueDensity = 1, nForages = 4)
dim(test) # good
sum(test) # approx number of mice 
apply(test, MARGIN = 2, FUN = sum) # falls of as expected

range(replicate(1000, {
  test <- GenTrapData(trapSpacing = 0.5, catchRadius = 0.3, boarder = 3, nSquares = 8, trueDensity = 1, nForages = 4)
  p1 <- sum(test[,1:2]) 
  p2 <- sum(test[,3:4]) 
  return(p1^2 / (p1 - p2))
})) ## NEVER ESTIMATES NEGATIVE NHAT FOR SQUARE 8!!!

test <- t(replicate(1000, {
  test <- GenTrapData(trapSpacing = 0.5, catchRadius = 0.1, boarder = 3, nSquares = 8, trueDensity = 1, nForages = 4)
  p1 <- sum(test[c(120,121,136,137), 1:2]) 
  p2 <- sum(test[c(120,121,136,137), 3:4]) 
  return(c(p1, p2))
}))
nhat <- test[,1]^2 / (test[,1] - test[,2])
range(nhat[is.finite(nhat)]) # For this CR, we see no negative nHat's for square 1 either! You can see them if you allow 
# non-overlapping traps in the middle, i.e. if cr = 0.1 for instance, there will be negative nHat values (not too many)




# Test function calculating both periods by trap
test <- GenTrapData(trapSpacing = 0.5, catchRadius = 0.3, boarder = 3, nSquares = 8, trueDensity = 3, nForages = 4)
pds <- t(apply(test, MARGIN =  1, FUN = function(x) {return(c(x[1]+x[2], x[3]+x[4]))}))
all.equal(pds, calcPeriodsByTrap(test, 4))
all(pds[,1] == calcPeriodsByTrap(test, 4)[,1])
all(pds[,2] == calcPeriodsByTrap(test, 4)[,2])

# test different number of forages
test <- GenTrapData(trapSpacing = 0.5, catchRadius = 0.3, boarder = 3, nSquares = 8, trueDensity = 3, nForages = 6)
pds <- t(apply(test, MARGIN =  1, FUN = function(x) {return(c(x[1]+x[2]+x[3], x[4]+x[5]+x[6]))}))
all.equal(pds, calcPeriodsByTrap(test, 6))
all(pds[,1] == calcPeriodsByTrap(test, 6)[,1])
all(pds[,2] == calcPeriodsByTrap(test, 6)[,2])



# Test getRings (matrix form)
image(GenRingAssignmentMat(4), col = rainbow(4))
image(GenRingAssignmentMat(6), col = rainbow(6))
image(GenRingAssignmentMat(8), col = rainbow(8))

# Test getRings (vector form)
image(matrix(GenRingAssignmentVec(4), 8, 8), col = rainbow(4))
image(matrix(GenRingAssignmentVec(6), 12, 12), col = rainbow(6))
image(matrix(GenRingAssignmentVec(8), 16, 16), col = rainbow(8))

all(c(GenRingAssignmentMat(8)) == GenRingAssignmentVec(8))
all(c(GenRingAssignmentMat(5)) == GenRingAssignmentVec(5))




# To test ProcessResults, I repeated
ProcessResults(uuid = 12002, paramset = 37, trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 4)
# checking the results to see how calulations are handled with varius pd1 and pd2, looks good



# Test parameter checking
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 4)
checkParameters(trapSpacing = 0, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 4)
checkParameters(trapSpacing = 0.5, catchRadius = -0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 4)
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 0, trueDensity = 2, nForages = 4)
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 0.5, trueDensity = 2, nForages = 4)
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 0, nForages = 4)
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 5)
checkParameters(trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 0)

RunSimulation(uuid = 1, paramset = 2, trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 0)
RunSimulation(uuid = 1, paramset = 2, trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 1)
RunSimulation(uuid = 1.9, paramset = 2.4, trapSpacing = 0.5, catchRadius = 0.2, boarder = 3, nSquares = 8, trueDensity = 2, nForages = 2)



# Generate some tests at scale to look at the general results
library(data.table)
library(ggplot2)
Sims <- replicate(1000, RunSimulation(uuid = 1, paramset = 1, trapSpacing = 0.5, catchRadius = 0.5, boarder = 3, nSquares = 4, trueDensity = 2, nForages = 4), simplify = FALSE)
Sims <- lapply(Sims, as.data.frame)
Sims <- rbindlist(Sims)
names(Sims) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")

# what precent good values
length(Sims[is.finite(nHat) & nHat > 0, nHat]) / nrow(Sims) # most! :) 

ggplot(Sims) + geom_density(aes(x = nHat, colour = factor(square))) + geom_vline(xintercept = (2*4*0.5 + 3)^2 * 2)
ggplot(Sims) + geom_density(aes(x = dHat, colour = factor(square))) + geom_vline(xintercept = 2)

# generally looks pretty well behaved at this point and gives numbers closer to what we'd expect








# Full test run:
Rcpp::sourceCpp("SimulationBackendFull.cpp")
TrapSpacing <- seq(from = 0.25, to = 2, by = 0.25)
CatchRadius <- seq(from = 0.25, to = 1, by = 0.25)
Boarder <- seq(from = 3, to = 3, by = 1)
Density <- seq(from = 0.5, to = 5, by = 1)
Parameters <- expand.grid(Density, Boarder, CatchRadius, TrapSpacing)
names(Parameters) <- c("Density", "Boarder", "CatchRadius", "TrapSpacing")
Parameters$paramset <- 1:nrow(Parameters)


Simulations <- rep(1:nrow(Parameters), each = 1000)
Simulations <- cbind(Simulations, 1:length(Simulations))

system.time(
TestSims <- lapply(1:nrow(Simulations), function(x) {
  s <- Simulations[x, ]
  paramset <- Parameters[Parameters$paramset == s[1], ]
  res <- RunSimulation(uuid = s[2], paramset = s[1], trapSpacing = paramset$TrapSpacing, catchRadius = paramset$CatchRadius, boarder = paramset$Boarder, nSquares = 4, trueDensity = paramset$Density, nForages = 4)
})
)
TestSims <- lapply(TestSims, as.data.frame)
TestSims <- rbindlist(TestSims)
names(TestSims) <- c("uuid", "paramset", "square", "pd1", "pd2", "pHat", "nHat", "aHat", "dHat")



#### Test new approach ####

# Test for isCaught2

# generate a mouse
mouse = GenMouse(nForages = 20, fieldSize = 13)
colMeans(mouse)
# generate traps
Traps = GenTraps()

# For just the first day, look at the frequency of the traps that catch the mouse
table(replicate(1000, isCaught2(forages = mouse, Traps = Traps, catchRadius = 1)[1]))
# For just the tenth day, look at the frequency of the traps that catch the mouse
table(replicate(1000, isCaught2(forages = mouse, Traps = Traps, catchRadius = 1)[10]))



# Test for GenAllMice
GenAllMice(trapSpacing = 2, catchRadius = 0.5, boarder = 3, nSquares = 4, trueDensity = 1, nForages = 4)



