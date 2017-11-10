# This script simulates a study (currently 2 periods made of two checks on each of two days is hard-coded)
# I first simualte a study, then calculate the important statistics. I store both the data and the statistics
# for later analyses. 


library(parallel)


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

rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)


#Simulate/generate the studies
Studies <-NULL
iter <- 1000
ncores <- detectCores()
genTime <- Sys.time()
cl <- makeCluster(ncores-1, type = "FORK")
Studies <- parLapply(1:iter, function(x) studySim())
stopCluster(cl)
genTime <- Sys.time() - genTime
genTime


Studies <- lapply(1:4, function(x) studySim())



# Use list of studies
# We're going to build a dataframe containing:
#     the study number,
#     p-value of anova test comparing the rings across the whole study, and
#     the statistics for each period (as per lisp file).

nperiods <- 2
periods <- list(c(1,2), c(3,4)) #the first two days are one period and the second two are the other

#Generate data frame with values: iteration, trap number, catches in period one, and catches in period two
#In post processing we can add ring labels and the total catches
Traps <- lapply(Studies, function(x) {
  
})
#Traps <- do.call("rbind", traps)
library(data.table)
Traps <- as.data.frame(rbindlist(Traps))


