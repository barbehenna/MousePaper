#Distribution of counts within each ring

## ---- Extra Libraries ---------
library(doParallel)
library(foreach)
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


#Simulate by parallel computing
#iter = 10000 -> 6.55 minutes (on 3 cores of Alton's mac)
#iter = 100000 -> 1.160445 hours (on 3 cores of Alton's mac)
iter = 100
cores=detectCores()
start_time <- Sys.time()
cl <- makeCluster(cores-1)
registerDoParallel(cl)
Catches <- rep(0, times = 64)
Catches <- matrix(Catches, nrow = 8, ncol = 8)
Catches <- foreach(i=1:iter, .combine = '+') %dopar% {
  tmp <- oneRun(print = TRUE)
  tmp
}
stopCluster(cl)
Sys.time() - start_time


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


## Seems close to normal
## Let's do both nonparametric gaussian based tests
model <- aov(formula = count~ring, data = mice)
pval <- summary(model)[[1]]$`Pr(>F)`[1]
if (pval) {
  model <- TukeyHSD(model)
  #comparisons labels(model$ring) 
}












