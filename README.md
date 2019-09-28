# MousePaper

The goal of this project is to better understand one method of estimating population statistics of a group via catch-and-release. We do this by simulation. We have a simulation written in C++ to generate many random studies. Each of these studies are then aggregated and analyzed with an R script. There is also a Shiny app that can be updated to show what the raw trap  data looks like (on the indiviual and pooled scales).(To use the app you may need to roll back the SimulationBackend.R a few commits, before it was moved to  8 squares.)

## Package

For the ease of use, there is an R package that wraps all of the C++ code we use to run the simulations. You can either load the use Rcpp by sourcing the file called `mousesim/src/SimulationBackendFull.cpp` or you can build and install the package then use it like any other. There are many ways to build the package, the easiest may run the following commands from the main package's directory

```
R CMD build mousesim
R CMD INSTALL mousesim_1.0.tar.gz 
```

After that, you can load the package like any other: `library(mousesim)`.

