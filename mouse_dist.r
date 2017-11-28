# Copyright (c) 2017 Alton Barbehenn

#This script is designed to understand the distribution of counts 
#within each ring of the traps (and if that's the best way at all)

## ---- Extra Libraries ---------
library(foreach)
library(doParallel)
library(raster)
library(reshape2)
library(ggplot2)
library(nortest)


## ----- Alton's simulation ----
d <- 1 #mice/m^2
ts <- 1.5 #m
fs <- 7*ts+5 #m
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

groups <- c(7,6,5,5,5,5,6,7, 
            6,4,3,3,3,3,4,6, 
            5,3,2,2,2,2,3,5, 
            5,3,2,1,1,2,3,5, 
            5,3,2,1,1,2,3,5, 
            5,3,2,2,2,2,3,5, 
            6,4,3,3,3,3,4,6, 
            7,6,5,5,5,5,6,7)


#Simulate by parallel computing
#iter = 10000 -> 6.55 minutes (on 3 cores of Alton's mac)
#iter = 100000 -> 1.160445 hours (on 3 cores of Alton's mac)
iter = 10000000
cores=detectCores()
start_time <- Sys.time()
cl <- makeCluster(cores-1)
registerDoParallel(cl)
Catches <- rep(0, times = 64)
Catches <- matrix(Catches, nrow = 8, ncol = 8)
Catches <- foreach(i=1:iter, .combine = '+') %dopar% {
  tmp <- matrixSim()
  tmp
}
stopCluster(cl)
Sys.time() - start_time


#Save the data (1000000 runs took 4.664996 hours on 7 Alton work laptop cores)
# write.matrix(Catches, file = "", sep = ",")
# test = read.csv(file = "/Users/abarbehenn/Documents/MousePaper/data/CatchesMatrix_1000000.csv", header = FALSE)
# test = as.matrix(test)
# Catches <- test

r <- raster(xmn = 0, xmx = 8, ymn = 0, ymx = 8, nrows = 8, ncols = 8)
r[] <- Catches
plot(r)
Catches


mice <- melt(Catches)
names(mice) <- c("row", "col", "count")
mice$ring <- factor(rings)


plt <- ggplot(data = mice) + geom_density(mapping = aes(x=count, colour=ring))
plt


qqnorm(mice$count[mice$ring==1])
qqline(mice$count[mice$ring==1], col=2)
shapiro.test(mice$count[mice$ring==1]) #Null: dist is normal

qqnorm(mice$count[mice$ring==2])
qqline(mice$count[mice$ring==2], col=2)
shapiro.test(mice$count[mice$ring==2]) #Null: dist is normal
lillie.test(mice$count[mice$ring==2]) #Null: dist is normal (must be more than 4 samples)

qqnorm(mice$count[mice$ring==3])
qqline(mice$count[mice$ring==3], col=2)
shapiro.test(mice$count[mice$ring==3]) #Null: dist is normal
lillie.test(mice$count[mice$ring==3]) #Null: dist is normal

qqnorm(mice$count[mice$ring==4])
qqline(mice$count[mice$ring==4], col=2)
shapiro.test(mice$count[mice$ring==4]) #Null: dist is normal
lillie.test(mice$count[mice$ring==4]) #Null: dist is normal


#From the qqplots we can pretty clearly see that that rings 3 and 4 are definately multinodal.
#Here are their plots by group
mice$group <- factor(groups)

#ring 3 is groups 3 and 4
qqnorm(mice$count[mice$group==3])
qqline(mice$count[mice$group==3], col=2)
shapiro.test(mice$count[mice$group==3]) #Null: dist is normal
lillie.test(mice$count[mice$group==3]) #Null: dist is normal

qqnorm(mice$count[mice$group==4])
qqline(mice$count[mice$group==4], col=2)
shapiro.test(mice$count[mice$group==4]) #Null: dist is normal

#ring 4 is groups 5, 6, and 7
qqnorm(mice$count[mice$group==5])
qqline(mice$count[mice$group==5], col=2)
shapiro.test(mice$count[mice$group==5]) #Null: dist is normal
lillie.test(mice$count[mice$group==5]) #Null: dist is normal

qqnorm(mice$count[mice$group==6])
qqline(mice$count[mice$group==6], col=2)
shapiro.test(mice$count[mice$group==6]) #Null: dist is normal
lillie.test(mice$count[mice$group==6]) #Null: dist is normal

qqnorm(mice$count[mice$group==7])
qqline(mice$count[mice$group==7], col=2)
shapiro.test(mice$count[mice$group==7]) #Null: dist is normal


#So it seems that the only group that isn't normal is group 3.
#This is really interesting because it isn't super appearent from 
#the raster plot which ones are the diffent ones or if this is a
#statistical annomoly.

plot(density(mice$count[mice$group==3]))
#if it were the ones just in from the corners, it would be exactly half the group
#it could be that one group has a tighter cluster or that the groups are shifted.




## Seems close to normal
## Let's do both nonparametric gaussian based tests
model <- aov(formula = count~ring, data = mice)
pval <- summary(model)[[1]]$`Pr(>F)`[1]
if (pval) {
  model <- TukeyHSD(model)
  #comparisons labels(model$ring) 
}



# Some of the rings seem normal, others are definately multi-modal
for (i in as.numeric(levels(mice$ring))) {
  print(i)
  print(shapiro.test(mice$count[mice$ring==i])$p.value)
  if (length(mice$count[mice$ring==i])>4) {
    print(lillie.test(mice$count[mice$ring==i])$p.value)
  }
}


# The groups all seem normally distributed
for (i in as.numeric(levels(mice$group))) {
  print(i)
  print(shapiro.test(mice$count[mice$group==i])$p.value)
  if (length(mice$count[mice$group==i])>4) {
    print(lillie.test(mice$count[mice$group==i])$p.value)
  }
}






