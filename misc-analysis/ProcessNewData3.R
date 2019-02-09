library(data.table)
library(ggplot2)
library(tidyr)


dist.mode <- function(vals, na.rm = FALSE) {
  dens <- density(vals, n = 2^20, na.rm = na.rm)
  return(dens$x[which.max(dens$y)])
}


Parameters <- fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Parameters.csv")
Simulations <-  fread("~/Documents/Projects/MousePaper/data/2019-01-10 21_53_47_Stats.csv")

# Merge in "True" density
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset", all.x = TRUE)


# Limit data to just the full squares
Simulations <- Simulations[square != 5]
# table(sign(Simulations$nHat))
# sum(is.na(Simulations$nHat))
Simulations <- Simulations[nHat > 0]


Simulations[, `:=`(resp = dHat/Density)]

ggplot(Simulations[sample(x = nrow(Simulations), size = 5000)], aes(x = TrapSpacing, y = resp)) + geom_point()

Simulation.Agg <- Simulations[, .(resp.avg = mean(resp, na.rm = TRUE), 
                                  resp.med = median(resp, na.rm = TRUE),
                                  resp.mode = dist.mode(resp, na.rm = TRUE)), by = .(TrapSpacing)]
Simulation.Agg <- Simulation.Agg %>% 
  gather(-TrapSpacing, key = "method", value = "resp")

ggplot(Simulation.Agg, aes(x = TrapSpacing, y = resp, colour = method)) + geom_point() + ggtitle("No square 5, positive nHat only") + labs(subtitle = "One point per TrapSpacing")
# the last two points are ever so slightly bi-modal, the other mode falls in line with the other points
# ggplot(Simulations[TrapSpacing==5.75], aes(x = resp)) + geom_density() + xlim(c(0, 100))




# percent valid (nHat>0) points
Simulation.Agg <- Simulations[, .(total = .N, total.valid = sum(nHat > 0, na.rm = TRUE)), by = .(TrapSpacing, CatchRadius)]
Simulation.Agg[, `:=`(percent.valid = total.valid / total)]
ggplot(Simulation.Agg, aes(x = TrapSpacing, y = CatchRadius, colour = factor(cut(percent.valid, 6)))) + geom_point() + scale_color_brewer("Spectral")


ggplot(Simulations[TrapSpacing == 3 & CatchRadius == 1.5 & Density == 5, ], aes(x = nHat, colour = factor(Boarder))) + geom_density() + xlim(c(-1000, 500)) # Pretty similar
ggplot(Simulations[TrapSpacing == 3 & CatchRadius == 1.5 & Boarder == 5], aes(x = nHat, colour = factor(Density))) + geom_density() + xlim(c(-500, 500)) + ggtitle("TS=3, CR=1.5, B=5") # Pretty different


Simulations[square == 3, `:=`(nHat.acc = nHat/NumMice), ]
Simulation.Agg <- Simulations[, .(avg.nHat.acc = dist.mode(nHat.acc)), by = .(TrapSpacing, CatchRadius)]
ggplot(Simulation.Agg, aes(x = TrapSpacing, y = CatchRadius, colour = factor(cut(avg.nHat.acc, 5)))) + geom_point()
