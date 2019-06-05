#### Libraries ####

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)


#### Load Data ####

Parameters <- fread("data/NewBackend-ParametersFull_B3.csv")
Simulations <- fread("data/NewBackend-SimulationsFull_B3.csv")

# merge in parameters
Simulations <- merge(x = Simulations, y = Parameters, by = "paramset")

# Remove invalid simulations
# Approx 1.7% of simulations, also cleans up dHat
Simulations <- Simulations[nHat >= 0]
Simulations <- Simulations[is.finite(nHat)]
Simulations <- Simulations[is.finite(pHat)]


#### Prep Data ####

Simulations %>%
  mutate(square = paste0("square", square)) %>%
  select(uuid, square, dHat, Density, TrapSpacing, CatchRadius) -> Simulations

Simulations = dcast(Simulations, ... ~ square, fill = NA, value.var = "dHat")
Simulations = na.omit(Simulations)

head(Simulations)

fwrite(Simulations, "data/trainingdata.csv")
