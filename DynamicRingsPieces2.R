#### Libraries ####

library(data.table)
library(pbapply)

library(ggplot2)
library(plotly)



####  Process Simulations ####

setwd("~/Documents/Projects/MousePaper/data/DiscreteSimByTSCRPair/")
files <- list.files()
# files <- sample(files, size = 75)


# file <- files[60]
# file <- sample(files, 1)


processedSimulations <- pblapply(files, FUN = function(file) {
  # Load Sim point
  Sim <- fread(file)
  
  ## Prepare data
  # Remove rows with NA values
  Sim <- na.omit(Sim)
  
  # Remove unwanted squares
  # 1 is too small, 5 is just the outer ring
  # Sim <- Sim[!is.element(square, c(1,5))]
  
  # Standardize density measurements
  # Percent error / Accuaracy of method
  # This is the same as standardized dHat (dHat / Density - 1)  
  Sim$dHat <- (Sim$dHat - Sim$Density) / Sim$Density
  
  # ggplot(Sim) + geom_density(aes(x = dHat, colour = factor(square))) + xlim(c(-1,1)) + ggtitle(paste("CR = ", Sim$CatchRadius[1], "TS =", Sim$TrapSpacing[1]))
  
  # Mean of  accuarcy  by  square
  Sim <-  Sim[, list(dHat.mean = mean(dHat), dHat.var = var(dHat)), by=list(square, TrapSpacing, CatchRadius, aHat)]
  
  # Sim <-Sim[is.element(Sim$square, c(2,3,4)), ]
  
  return(Sim[which.min(abs(dHat.mean)), ])
})
processedSimulations <- rbindlist(processedSimulations)





#### Plot + Analyze #####


plot_squares <- processedSimulations %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~square, color = ~square) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = "Trap Spacing"),
    yaxis = list(title = "Catch Radius"),
    zaxis = list(title = "Square")
  ))
plot_squares



plot_accuracy <- processedSimulations %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean, color = ~square) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = "Trap Spacing"),
    yaxis = list(title = "Catch Radius"),
    zaxis = list(title = "Mean dHat")
  ))
plot_accuracy


# If a set of parameters has square 5 as its most accurate "square", then maybe the traps 
# are too far apart and the fence effect helps the accuracy by increasing the catchs. Just a though



model <- lm(dHat.mean ~ log(CatchRadius/TrapSpacing), processedSimulations)
summary(model)
autoplot(model)

processedSimulations$model.fitted <- predict(model, newdata = processedSimulations)
model.error <- processedSimulations %>% 
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~dHat.mean,  color = ~(dHat.mean - model.fitted)) %>%
  add_markers()
model.error

ggplot(processedSimulations,  aes(x = CatchRadius/TrapSpacing, y = log(dHat.mean))) + geom_point() + geom_smooth() + geom_smooth(method = "lm")




# Let's look at the amount of area needed to correct for the mean density error

processedSimulations$area.adjustment <- processedSimulations$dHat.mean * processedSimulations$aHat

plot_areaadj <- processedSimulations %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~area.adjustment, color = ~factor(square)) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = "Trap Spacing"),
    yaxis = list(title = "Catch Radius"),
    zaxis = list(title = "Area Adjustment")
  ))
plot_areaadj
# There is a lot more noise when you zoom in or crop or look at the info tab

ggplot(processedSimulations) + geom_density(aes(x = area.adjustment, colour = factor(square))) + xlim(c(-10, 10))
# Notice how large the range is for at least a few of the squares!

# Basic statistics on adjustment by square
processedSimulations[, list(mean = mean(area.adjustment), med = median(area.adjustment), 
                            var = var(area.adjustment), min = min(area.adjustment),  
                            max = max(area.adjustment)), 
                     by = square]

