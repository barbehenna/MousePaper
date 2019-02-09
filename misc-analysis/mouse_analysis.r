# Copyright (c) 2017 Alton Barbehenn

## ---- Extra Libraries ---------
library(doParallel)
library(foreach)
library(raster)
library(reshape2)
library(pgirmess)



## ------ Global Variables ------
np <- 50 #number of mice
nv <- 9 #number of trips
ts <- 1.5 ## unit = sigma
fs <- 7*ts+5 #field size
delta <- 0.5 #area where mouse is caught 


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

#Simulate by parallel computing
# iter = 10000 -> 6.55 minutes (on 3 cores of Alton's mac)
# iter = 100000 -> 1.160445 hours (on 3 cores of Alton's mac)
# iter = 1000
# cores=detectCores()
# start_time <- Sys.time()
# cl <- makeCluster(cores-1)
# registerDoParallel(cl)
# Catches <- rep(0, times = 64)
# Catches <- matrix(Catches, nrow = 8, ncol = 8)
# Catches <- foreach(i=1:iter, .combine = '+') %dopar% {
#   tmp <- oneRun(print = TRUE)
#   tmp
# }
# stopCluster(cl)
# Sys.time() - start_time

periods <- matrix(1:4, nrow = 2, ncol = 2) #period_i <- period[,i]
Exp <- oneExperiment()
#study[[i]] is period i (i=1,2), study[[3]] is the total study
study <- list()
for(i in 1:2) { #for each period
  per <- periods[,i]
  traps <- Exp$trap[is.element(Exp$day, per)]
  catches <- rep(0, 64)
  mice <- lapply(1:64, function(x) catches[x] <- catches[x] + sum(traps==x))
  mice <- unlist(mice)
  study[[i]] <- data.frame(count = mice, ring = factor(rings))
}
total <- rep(0, 64)
total <- unlist(lapply(1:64, function(x) total[x] <- total[x] + sum(Exp$trap==x)))
study[[i+1]] <- data.frame(count=total, ring = factor(rings))

#p-value
p <- summary(aov(formula = count~ring, data = study[[i+1]]))[[1]]$`Pr(>F)`[1]




#Alternate method of parallel computing to simulate
# library(parallel)
# st <- Sys.time()
# cl <- makeCluster(3, type = "FORK")
# iter <- rep(TRUE, 5000)
# li <- parLapply(cl, iter, function(x) oneRun(x))
# stopCluster(cl)
# Catches <- Reduce('+', li)
# Sys.time()-st


r <- raster(xmn = 0, xmx = 8, ymn = 0, ymx = 8, nrows = 8, ncols = 8)
r[] <- Catches
plot(r)
Catches


#rings <- c(rep(4,8), 4,rep(3,6),4, 4,3,rep(2,4),3,4, 4,3,2,1,1,2,3,4, 4,3,2,1,1,2,3,4, 4,3,rep(2,4),3,4, 4,rep(3,6),4, rep(4,8))
rings <- c(4,4,4,4,4,4,4,4, 
           4,3,3,3,3,3,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,1,1,2,3,4, 
           4,3,2,2,2,2,3,4, 
           4,3,3,3,3,3,3,4, 
           4,4,4,4,4,4,4,4)
mice <- melt(Catches)
names(mice) <- c("row", "col", "count")
mice$ring <- factor(rings)

#Test for differences in the rings
model <- aov(formula = count ~ ring, data = mice)
summary(model)
TukeyHSD(model)

#Nonparametic alternative
kruskal.test(count ~ ring, data = mice)
kruskalmc(resp = count ~ ring, data = mice, probs = 0.05)



#Estimate population density
#nmice = fieldsize * density
catchrate <- mice$count[is.element(mice$ring, 1:2)]
subfield <- iter*(4*ts)^2
den_approx <- sum(catchrate)/subfield
den_approx


# Change so that I can look at the data for each forage
    # aggrigate forages into two periods (sum of captures)
# Try to replicate the statistics used in the original paper
    # see how well they work
# See if we can do any better by not aggegating the data into days
    # if so, what happens if we sample for a different number of days