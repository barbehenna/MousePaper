library(Rcpp)
library(raster)

Rcpp::sourceCpp('Documents/Projects/MousePaper/RandomWalk.cpp')
Rcpp::sourceCpp('Documents/Projects/MousePaper/SimulationBackend.cpp')

trapspacing <- 2
catchradius <- 0.5
boarder <- 4
den <- 10
fieldsize <- 7*trapspacing + 2*boarder
nmice <- round(den*fieldsize*fieldsize)
nvisits <- 4

st <- Sys.time()

# library(microbenchmark)
# microbenchmark(trapSim1(ts=trapspacing, fs=fieldsize, np=10, delta=catchradius, nv=nvisits),
#                trap_rw(ts=trapspacing, fs=fieldsize, np=10, delta=catchradius, nv=nvisits, stepSize=0.25, maxiter=1000),
#                times = 1000)
# Looks like random walk is about 500-600 times longer than teleportation

sim_curr <- trapSim1(ts=trapspacing, fs=fieldsize, np=nmice, delta=catchradius, nv=nvisits)
sim_curr <- as.data.frame(sim_curr)
names(sim_curr) <- c("trap", "day")
sim_curr$trap <- sim_curr$trap+1
sim_curr$day <- sim_curr$day+1
sim_curr <- na.omit(sim_curr)

sim_rw <- trap_rw(ts=trapspacing, fs=fieldsize, np=nmice, delta=catchradius, nv=nvisits, stepSize=0.25, maxiter=10)
sim_rw <- as.data.frame(sim_rw)
names(sim_rw) <- c("trap", "day", "iter")
sim_rw$trap <- sim_rw$trap+1
sim_rw$day <- sim_rw$day+1
sim_rw$iter <- sim_rw$iter+1
sim_rw <- na.omit(sim_rw)


sim_curr_vec <- unlist(lapply(seq(64), function(x) sum(sim_curr$trap == x)))
sim_curr_mat <- matrix(sim_curr_vec, nrow = 8)
sim_curr_rast <- raster(nrows = 8, ncols = 8)
sim_curr_rast[] <- sim_curr_mat

sim_rw_vec <- unlist(lapply(seq(64), function(x) sum(sim_rw$trap == x)))
sim_rw_mat <- matrix(sim_rw_vec, nrow = 8)
sim_rw_rast <- raster(nrows = 8, ncols = 8)
sim_rw_rast[] <- sim_rw_mat

print(Sys.time() - st)

plot(sim_curr_rast, main = "teleport")
plot(sim_rw_rast, main = "random walk")


############### Find correct number of iterations ###############

# Maybe, instead of finding out how many steps of a specific step-size
# each mouse should take as a max, we just compair it to the expected
# distance a mouse teleports to and then calculate the number of steps
# from that (given a step-size).

# If a mouse has a 1/100 chance of traveling 3 sigma away...

# TODO: make sure the mouse is randomly walking in the same
# range as it teleports to. So if it teleports out to 2-3 sigma,
# the walks should go out to 2-3 sigma. ALSO, do the walks give
# the mouse a distribution that is normal about the origin when
# looking at a slice?















