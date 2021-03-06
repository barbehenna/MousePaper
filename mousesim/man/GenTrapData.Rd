% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GenTrapData}
\alias{GenTrapData}
\title{Generate Raw Trapping Data}
\usage{
GenTrapData(trapSpacing, catchRadius, border, nSquares, trueDensity,
  nForages)
}
\arguments{
\item{trapSpacing}{(double >= 0) distance between traps}

\item{catchRadius}{(double >= 0) taxicab metric distance catch threshold}

\item{border}{(double >= 0) buffer outside trap grid}

\item{nSquares}{(int >= 1) num concentric rings of traps}

\item{trueDensity}{(double >= 0) true density of mice}

\item{nForages}{(int >= 1) number of days to simulate}
}
\value{
matrix of number of mice caught in each trap and each day
}
\description{
For a given set of simulation parameters, calculate the remaining parameters and then 
simulate the raw trap data. The return of this function is a NumericMatrix where each
row is a single trap and each day is a single day (ordered as usual). The value of each
entry is the number of mice caught in that trap on that day.
}
