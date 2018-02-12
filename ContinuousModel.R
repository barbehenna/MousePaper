# Copyright (c) 2018 Alton Barbehenn

# This script loads a continuous set of simulation parameters and their resulting statistics 
# (From ContinuousSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(pbapply)
library(parallel)

Sim <- read.csv("data/2018-01-27 00:13:19_Stats.csv", header = TRUE, row.names = NULL)
Parameters <- read.csv("data/2018-01-27 00:13:19_Parameters.csv", header = TRUE, row.names = NULL)

Sim <- merge(Sim, Parameters, by = "paramset")
rm(Parameters)

# ncores <- detectCores()
# cl <- makeCluster(ncores-1, type = "FORK")

# error <- pblapply(unique(Sim$UniqueID), cl = cl, function(x) {
error <- pblapply(unique(Sim$UniqueID), function(x) {
  avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == x][1]
  ts <- Sim$TrapSpacing[Sim$UniqueID == x][1]
  cr <- Sim$CatchRadius[Sim$UniqueID == x][1]
  tmp <- data.frame(den, avg/den, ts, cr)
  return(tmp)
})

# stopCluster(cl)
error <- rbindlist(error)
names(error) <- c("den", "perc", "TrapSpacing", "CatchRadius")


# ggplot(data = error[error$CatchRadius == 0.5]) + geom_density(aes(x = perc, colour = factor(TrapSpacing))) + xlim(-1,3)


# Attempt one at bi-variate analysis
# error<-data.table(error)
# error_proc <- error[, .(perc = mean(perc)), by = .(TrapSpacing, CatchRadius)]
# ggplot(data = error_proc, mapping = aes(x = TrapSpacing, y = CatchRadius))+geom_point(aes(colour = perc, size = 4))

# Looks like the predictions get better as the catch radius increases
# tmp <- error[error$TrapSpacing == 6,]
# ggplot(data = tmp) + geom_density(mapping = aes(x = perc, colour = factor(CatchRadius))) + xlim(-1,3)


# This shouldn't be needed because I'm sampling sufficiantly percise real numbers
# There's a MUCH cleaner way to do this
# pts <- expand.grid(sort(unique(error$TrapSpacing)), sort(unique(error$CatchRadius)))
# names(pts) <- c("TrapSpacing", "CatchRadius")
# pts$mean_perc <- 0.0
# pts$size <- 0.0
# for (i in seq(nrow(pts))) {
#   ts <- pts$TrapSpacing[i]
#   cr <- pts$CatchRadius[i]
#   avg <- mean(error$perc[error$TrapSpacing == ts & error$CatchRadius == cr], na.rm = TRUE)
#   pts$size[i] <- length(error$perc[error$TrapSpacing == ts & error$CatchRadius == cr]) #just for reference
#   if (is.finite(avg)) {
#     pts$mean_perc[i] <- avg
#   }
# }


p <- plot_ly(x = error$TrapSpacing, y = error$CatchRadius, z = error$perc) %>% add_surface()
saveWidget(p, file = "perc_vs_ts_cr.html")





