# Copyright (c) 2017 Alton Barbehenn

# This script simulates a study (currently 2 periods made of two checks on each of two days is hard-coded)
# I first simualte a study, then calculate the important statistics. I store both the data and the statistics
# for later analyses. 


library(parallel)
library(ggplot2)
library(data.table)


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

print(Sys.time())
#Simulate/generate the studies
Studies <-NULL
iter <- 1000
ncores <- detectCores()
cl <- makeCluster(ncores-1, type = "FORK")
Studies <- parLapply(cl, 1:iter, function(x) studySim())
stopCluster(cl)
print(Sys.time())

# genTime <- Sys.time()
# Studies <- lapply(1:iter, function(x) studySim())
# genTime <- Sys.time() - genTime
# genTime

# Generate data frame with values: simulation number, trap number, catches in period one, 
# and catches in period two. In post processing we can add ring labels and the total catches.
# The first two days are one period and the second two are the other
Studies <- lapply(1:iter, function(x) {
  p1 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day <= 2] == y))
  p1 <- unlist(p1)
  p2 <- lapply(1:64, function(y) sum(Studies[[x]]$trap[Studies[[x]]$day >= 3] == y))
  p2 <- unlist(p2)
  return(data.frame(simnum=rep(x,64), trap=1:64, period1=p1, period2=p2, total=p1+p2, ring=rings))
})
# Traps <- do.call(rbind, Studies)
Traps <- as.data.frame(rbindlist(Studies)) #Faster than do.call?
print(Sys.time())

# Use list of studies
# We're going to build a dataframe containing:
#     the study number,
#     p-value of anova test comparing the rings across the whole study, and
#     the statistics for each period (as per lisp file).

# These are the statistics from Lisp file:
#     pHat: catch-rate estimator, dHat: density estimator
#     nHat: abundency estimator, aHat: Area estimator.
# In dad's simulation look at concentric rings ((1), (1,2), (1,2,3), (1,2,3,4), and (4)).

aHat <- 1:4 #ring numbers
aHat <- (ts*2*aHat)^2 #concentric ring areas
aHat <- c(aHat, aHat[4]-aHat[3]) #just ring 4

squares <- list(1, 1:2, 1:3, 1:4, 4)

# for each simulated study:
#    for each square (1-5 for convience, 5 being just ring 4)
#        calculate the study statistics
Stats <- lapply(1:iter, function(x) {
  p1 <- lapply(squares, function(y) sum(Studies[[x]]$period1[Studies[[x]]$ring %in% y]))
  p2 <- lapply(squares, function(y) sum(Studies[[x]]$period2[Studies[[x]]$ring %in% y]))
  nHat_dHat <- lapply(1:length(squares), function(x) {
    if (p1[[x]] == p2[[x]]) {
      nHat <- NA
      dHat <- NA
    } else {
      nHat <- (p1[[x]]^2)/(p1[[x]]-p2[[x]])
      dHat <- nHat/aHat[[x]]
    }
    return(c(nHat, dHat))
  })
  nHat_dHat <- do.call(rbind, nHat_dHat)
  pHat <- lapply(1:length(squares), function(x) {
    if (p1[[x]] == 0) {
      pHat <- NA
    } else {
      pHat <- 1 - sqrt(p2[[x]]/p1[[x]])
    }
    return(pHat)
  })
  pHat <- unlist(pHat)
  return(data.frame(simnum=rep(x,length(squares)),
                    dHat=nHat_dHat[,2],
                    nHat=nHat_dHat[,1],
                    pHat=pHat,
                    aHat=aHat))
})
print(Sys.time())


# Look at parameters by aHat/square not simulation number
StatsDF <- as.data.frame(rbindlist(Stats))
# write.csv(StatsDF, "/Users/Alton/Documents/Projects/MousePaper/Stats_temp.csv", row.names = FALSE)
# StatsDF <- read.csv("/Users/Alton/Documents/Projects/MousePaper/Stats_temp.csv", header = TRUE)
Stats_aHat <- split(StatsDF, aHat)
aHat <- sort(aHat)
dHat <- unlist(lapply(Stats_aHat, function(x) median(na.omit(x$dHat))))
nHat <- unlist(lapply(Stats_aHat, function(x) median(na.omit(x$nHat))))
pHat <- unlist(lapply(Stats_aHat, function(x) median(na.omit(x$pHat))))
Stats_aHat <- data.frame(dHat, nHat, pHat, aHat, row.names = NULL)



StatsDF$aHat <- factor(StatsDF$aHat)
dHat_plt <- ggplot(StatsDF, aes(x = dHat, color = aHat)) + geom_density(na.rm = TRUE) + xlim(0, 5)
dHat_plt


dHat_plt <- ggplot(StatsDF, aes(x = aHat, y = dHat)) + geom_boxplot(na.rm = TRUE)
dHat_plt










