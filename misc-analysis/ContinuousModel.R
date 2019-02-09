# Copyright (c) 2018 Alton Barbehenn

# This script loads a continuous set of simulation parameters and their resulting statistics 
# (From ContinuousSimulation.R). It then runs some of the interesting analyses to 
# determine whether or not the ranges are good and if we want to make the resolution finer 
# for more data points. 

library(data.table)
# library(ggplot2)
library(plotly)
library(htmlwidgets)
# library(pbapply)
# library(parallel)
 
# Sim <- read.csv("data/2018-01-27 00:13:19_Stats.csv", header = TRUE, row.names = NULL)
# Parameters <- read.csv("data/2018-01-27 00:13:19_Parameters.csv", header = TRUE, row.names = NULL)

# Sim <- merge(Sim, Parameters, by = "UniqueID")            
# rm(Parameters)

# Sim <- Sim[,c("UniqueID", "square", "Density", "CatchRadius", "TrapSpacing", "dHat")]
# write.csv(Sim, "data/SimAbv.csv", row.names = FALSE)

# Sim <- read.csv("data/SimAbv.csv", header = TRUE, row.names = NULL)
Sim <- read.csv("./data/SimAbv.csv", header = TRUE, row.names = NULL)

Sim <- Sim[Sim$square <= 3, ]
Sim <- data.table(Sim)
SimAgg <- Sim[, .(avg = mean(dHat, na.rm = TRUE)), by = .(UniqueID,Density,TrapSpacing,CatchRadius)]

# SimAgg <- read.csv("./data/SimAgg_Continuous.csv", header = TRUE)
write.csv(SimAgg, "./data/SimAgg_Continuous.csv", row.names = FALSE)

SimAgg$PercError <- SimAgg$avg / SimAgg$Density

# ncores <- detectCores()
# cl <- makeCluster(ncores-1, type = "FORK")

# error <- pblapply(unique(Sim$UniqueID), function(x) {
# error <- pblapply(unique(Sim$UniqueID), cl = cl, function(x) {
#   avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
#   den <- Sim$Density[Sim$UniqueID == x][1]
#   ts <- Sim$TrapSpacing[Sim$UniqueID == x][1]
#   cr <- Sim$CatchRadius[Sim$UniqueID == x][1]
#   tmp <- data.frame(den, avg/den, ts, cr)
#   return(tmp)
# })


i <- 0
error <- data.frame(den = as.numeric(), perc = as.numeric(), ts = as.numeric(), cr = as.numeric())
for (id in unique(Sim$UniqueID)) {
  avg <- mean(Sim$dHat[Sim$UniqueID == x & Sim$square <= 3], na.rm = TRUE)
  den <- Sim$Density[Sim$UniqueID == x][1]
  ts <- Sim$TrapSpacing[Sim$UniqueID == x][1]
  cr <- Sim$CatchRadius[Sim$UniqueID == x][1]
  error <- rbind(error, c(den, avg/den, ts, cr))
  if (i %% 100000 == 0) {
    print(paste(i, Sys.Time()))
  }
  i <- i + 1
}


# stopCluster(cl)
# error <- rbindlist(error)
names(error) <- c("den", "perc", "TrapSpacing", "CatchRadius")
write.csv(error, "ContinuousError.csv", row.names = FALSE)

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


Sim <- Sim[Sim$square <= 3, ]
Sim <- data.table(Sim)
SimAgg <- Sim[, .(avg = median(dHat, na.rm = TRUE)), by = .(UniqueID,Density,TrapSpacing,CatchRadius)]


p <- plot_ly(SimAgg, x = ~TrapSpacing, y = ~CatchRadius, z = ~PercError) %>% add_markers()

p <- plot_ly(x = error$TrapSpacing, y = error$CatchRadius, z = error$perc) %>% add_surface()
saveWidget(p, file = "perc_vs_ts_cr2.html")



SimAgg_Continuous <- read.csv("~/Documents/Projects/MousePaper/data/SimAgg_Continuous.csv")
SimAgg_Continuous$PercError <- SimAgg_Continuous$avg / SimAgg_Continuous$Density
SimAgg_Continuous <- na.omit(SimAgg_Continuous)

p <- plot_ly(x = SimAgg_Continuous$TrapSpacing, 
             y = SimAgg_Continuous$CatchRadius, 
             z = SimAgg_Continuous$PercError) %>% add_surface()
p



p <- plot_ly(tmp, x = ~TrapSpacing, y = ~CatchRadius, z = ~PercError,
             marker = list(color = ~PercError, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE)) %>% 
  add_markers()
p


SimAgg_Continuous$CR_TS_Ratio <- SimAgg_Continuous$CatchRadius / SimAgg_Continuous$TrapSpacing
# g <- ggplot(data = SimAgg_Continuous, mapping = aes(x = CR_TS_Ratio, y = PercError)) + geom_point()
# g








