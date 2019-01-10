library(data.table)
library(nortest)
library(pbapply)


setwd("~/Documents/Projects/MousePaper/data/DiscreteSimByTSCRPair/")
files <- list.files()

processedData <- pblapply(X = files, FUN = function(file) {
  # Constants
  rings <- c(2, 3, 4)

  # Load data
  Sim <- fread(file)

  # Prep data
  Sim <- na.omit(Sim)

  # Standardize density measurements
  # Percent error
  Sim$dHat <- (Sim$dHat - Sim$Density) / Sim$Density

  # fix subtle differences?
  # Sim <- Sim[Sim$Boarder == 6, ]

  # Test for normality of each ring independently
  normality <- rep(0, times = 4)
  for (i in 1:length(normality)) {
    # Lillie Test is the Kolomorov-Smirnov test comparing the data to a normal distribution
    # of known mean and variance (substituting in the sample mean and variance). The NULL
    # hypothesis is that the distributions are equal, i.e. the data is distributed normally
    # with it's sample mean and sample variance as its mean and variance parameters. Hence
    # for SMALL p-values we reject the null and conclude that the data is not normally
    # distributed with it's sample mean and variance.
    test <- lillie.test(Sim$dHat[Sim$square == i])
    normality[i] <- test$p.value
  }

  # Compare two distributions from consecutive consentric rings
  # If normally distributed use t-test, otherwise use Wilcoxon test (non-parametric alternative)
  #
  # Do I need to use Bon-Feroni Correction for this cut-off?
  # Yes. I want to know if all of the rings are normal, together, so
  # use alpha/n level test (0.05/4 = 0.0125). Conservative, but standard.
  #
  out <- data.frame()
  if (all(normality < 0.0125)) { # Chance data is non-normal
    for (ring in rings) {
      # Compare two distributions from consecutive consentric rings
      # FIND MANN-WHITNEY TEST INSTEAD (for theoretical support) (in NSM3?)
      test <- wilcox.test(
        x = Sim$dHat[Sim$square == (ring - 1)],
        y = Sim$dHat[Sim$square == ring],
        alternative = "two.sided"
      )
      # Store results
      out <- rbind(out, c(ring, test$p.value))
    }
  } else { # Consider the densities to be normally distributed
    for (ring in rings) {
      # Compare two distributions from consecutive consentric rings
      test <- t.test(
        x = Sim$dHat[Sim$square == (ring - 1)],
        y = Sim$dHat[Sim$square == ring],
        alternative = "two.sided"
      )

      # Store results
      out <- rbind(out, c(ring, test$p.value))
    }
  }
  names(out) <- c("ring", "p.value")

  # Get the largest homogeneous concentric ring
  # Is this the best rule?
  if (is.element(2, out$ring[out$p.value > 0.05])) {
    bestRing <- max(out$ring)
  } else {
    bestRing <- 1
  }

  # Process data
  Sim <- Sim[Sim$square == bestRing, ]

  return(list(
    "CatchRadius" = Sim$CatchRadius[1],
    "TrapSpacing" = Sim$TrapSpacing[1],
    "meanerror" = mean(Sim$dHat),
    "sderror" = sd(Sim$dHat),
    "squareUsed" = bestRing
  ))
})


processedData <- rbindlist(processedData)

# Save data
# fwrite(processedData, file = "~/Desktop/processedData.csv")

library(ggplot2)
library(plotly)


plot_squares <- processedData %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~squareUsed, color = ~squareUsed) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = "Trap Spacing"),
    yaxis = list(title = "Catch Radius"),
    zaxis = list(title = "Square")
  ))
plot_squares

plot_error <- processedData %>%
  plot_ly(x = ~TrapSpacing, y = ~CatchRadius, z = ~meanerror, color = ~ -log10(abs(sderror))) %>%
  add_markers()
plot_error
