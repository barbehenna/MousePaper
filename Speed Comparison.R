# This script just serves to compare the speed of the main two ways that 
# I've found to do parallel computing in R.

library(microbenchmark)
library(parallel)
library(foreach)
library(doParallel)

## ----- Alton's simulation ----
d <- 1 #mice/m^2
ts <- 1.5 #m
fs <- 7*ts+7 #m
np <- as.integer(d*fs*fs) #mice
nv <- 4
delta <- 0.5

if (delta > ts/2) {
  stop("Something weird may happen with this run becuase a mouse can be caught by two traps in the same forage")
}


# Common parameters for the parallel computing
iter = 100
ncores=detectCores()


#parLapply
parLapply_time <- Sys.time()
Catches <- NULL
cl <- makeCluster(ncores-1, type="FORK")
Catches <- parLapply(cl, 1:iter, function(x) matrixSim())
stopCluster(cl)
Catches <- Reduce('+', Catches)
parLapply_time <- Sys.time() - parLapply_time


#Foreach
foreach_time <- Sys.time()
Catches <- NULL
cl <- makeCluster(ncores-1, type = "FORK")
registerDoParallel(cl)
Catches <- rep(0, times = 64)
Catches <- matrix(Catches, nrow = 8, ncol = 8)
Catches <- foreach(i=1:iter, .combine = '+') %dopar% {
  tmp <- matrixSim()
  tmp
}
stopCluster(cl)
foreach_time <- Sys.time() - foreach_time


# Compare the times
parLapply_time
foreach_time


times <- data.frame(iteration = c(100, 1000, 10000), 
                    'parLapply (s)' = c(2.821888, 25.09182, 246.794), 
                    'foreach (s)' = c(4.072864, 25.33407, 209.6563))

# So it seems reasonable to use either method, but it may be worth using 
# foreach() for the really computations with a high number of iterations

# Let's be a little more rigorous and detailed.

ncores=detectCores()

f <- function(ncores) {
  Catches <- NULL
  cl <- makeCluster(ncores-1, type="FORK")
  Catches <- parLapply(cl, 1:iter, function(x) matrixSim())
  stopCluster(cl)
  Catches <- Reduce('+', Catches)
}

g <- function(ncores) {
  Catches <- NULL
  cl <- makeCluster(ncores-1, type = "FORK")
  registerDoParallel(cl)
  Catches <- rep(0, times = 64)
  Catches <- matrix(Catches, nrow = 8, ncol = 8)
  Catches <- foreach(i=1:iter, .combine = '+') %dopar% {
    tmp <- matrixSim()
    tmp
  }
  stopCluster(cl)
} 

# About 1.40 minutes on Alton's work laptop 7 cores used
iter = 100
st <- Sys.time()
microbenchmark(times = 10,
               a <- f(ncores),
               b <- g(ncores))
Sys.time() - st
# Unit: seconds
# expr      min      lq     mean   median       uq      max neval cld
# a <- f(ncores) 3.975289 3.98679 4.169670 4.082982 4.208818 4.912478    10   a
# b <- g(ncores) 3.975198 4.05166 4.164247 4.100875 4.281331 4.566719    10   a


# About 8.35 minutes on Alton's work laptop 7 cores used
iter = 1000
st <- Sys.time()
microbenchmark(times = 10,
               a <- f(ncores),
               b <- g(ncores))
Sys.time() - st
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# a <- f(ncores) 23.42245 23.67325 24.54033 24.01744 24.39299 28.59057    10   a
# b <- g(ncores) 23.39965 24.12750 25.54315 24.55611 27.63518 29.41864    10   a


# About ____~78____ minutes on Alton's work laptop 7 cores used
iter = 10000
st <- Sys.time()
microbenchmark(times = 10,
               a <- f(ncores),
               b <- g(ncores))
Sys.time() - st



# So we see that there isn't a noticable difference between the two methods 
# when we're doing low numbers of iterations


