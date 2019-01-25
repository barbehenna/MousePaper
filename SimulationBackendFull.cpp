#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;


// TODO: 
// function to calculate statistic by square
//    using: function to generate the square each trap corresponds to
// -- or --
// function to calculate statistic by square
//    using: function that calculates catches per square by day
// 
// function to save single simulation at a file name (optional for above function)



// For a given single mouse's forage coordinates, go through all traps and 
// determine which (if any) trap caught the mouse. Return -1 if the mouse is
// not caught and the traps index of the trap that caught the mouse if it was
// caught.
// The general algorithm (as implemented currently) does not implment recapture
// of the mouse, but could be extended to do so by returning an array of indicies
// (one per forage) of the traps the mouse was caught in.
// [[Rcpp::export]]
NumericVector isCaught(NumericMatrix forages, NumericMatrix Traps, double catchRadius) {
  std::vector<int> caughtin;
  NumericVector result(2); // first term is trap, second term is day
  double dx, dy;
  
  result(0) = -1;
  result(1) = -1;
  
  for (int day = 0; day < forages.nrow(); day++) { // for each forage
    for (int trap = 0; trap < Traps.nrow(); trap++) { // for each trap
      // calculate distance to trap
      dx = forages(day,0) - Traps(trap,0);
      dy = forages(day,1) - Traps(trap,1);
      
      // potentially catch mouse if close enough to trap
      // square catch area
      if ((std::abs(dx) < catchRadius) && (std::abs(dy) < catchRadius)) {
        // if mouse is close to trap, it could be caught
        caughtin.push_back(trap);
        // std::cout << "catch! ->  day " << day << " trap " << trap << ": " << Traps(trap,0) << ", " << Traps(trap,1) << std::endl;
      }
    }
    
    // If the mouse was caught in a trap on this forage
    if (!caughtin.empty()) {
      //return caughtin[(int) runif(1, 0, caughtin.size())(0)];
      result(0) = caughtin[(int) runif(1, 0, caughtin.size())(0)];
      result(1) = day;
      return result;
    }
  }
  
  return result;
}


// Generates a single mouse in a field of spcified size and returns all of the
// locations the mouse forages (number of forages given by nforages).
// In the returned matrix, each row is a forage coordinate (x,y).
// [[Rcpp::export]]
NumericMatrix GenMouse(int nForages, double fieldSize) {
  NumericMatrix mouse(nForages, 2); // forage locations for a single mouse
  
  NumericVector home = runif(2); // x and y home coordinates
  home = (home - 0.5) * fieldSize;
  
  mouse(_,0) = rnorm(nForages) + home[0]; // generate the forages
  mouse(_,1) = rnorm(nForages) + home[1];

  return mouse;
}


// This function returns a matrix of trap coordinates, as such, the returned
// matrix should be a nx2 matrix, (x,y) pairs
// [[Rcpp::export]]
NumericMatrix GenTraps(int nSquares = 8, double trapSpacing = 1.0) {
  const int trapsperside = 2*nSquares;
  double b[trapsperside];
  NumericMatrix trapCoordinates(trapsperside*trapsperside, 2);
  
  for (int i = 0; i < trapsperside; i++) {
    b[i] = (trapSpacing/2) * ((2*i) - (trapsperside - 1));
  }
  
  for (int i = 0; i < trapsperside; i++) {
    for (int j = 0; j < trapsperside; j++) {
      trapCoordinates((trapsperside*i)+j, 0) = b[j];
      trapCoordinates((trapsperside*i)+j, 1) = -b[i];
    }
  }
  
  return trapCoordinates;
}


// For a given set of simulation parameters, calculate the remaining parameters and then 
// simulate the raw trap data. The return of this function is a NumericMatrix where each
// row is a single trap and each day is a single day (ordered as ususal). The value of each
// entry is the number of mice caught in that trap on that day. 
// [[Rcpp::export]]
NumericMatrix GenTrapData(double trapSpacing, double catchRadius, double boarder, int nSquares, double trueDensity, int nForages) {
  double fieldSize = (((2 * nSquares) - 1) * trapSpacing) + (2 * boarder);
  int nmice = std::round(fieldSize * fieldSize * trueDensity);
  
  // std::cout << "fs = " << fieldSize << " nmice = " << nmice << std::endl; // check correct simulation values
  
  // Generate (place) traps 
  NumericMatrix Traps = GenTraps(nSquares, trapSpacing);
  
  // Data collection data-structure
  // One row per trap and one column for catch opprotunity
  NumericMatrix trapCountByDay(Traps.nrow(), nForages); // initialized to 0

  // for each of the nmice
  for (int mouse = 0; mouse < nmice; mouse++) {
    // generate the mouse
    NumericMatrix currMouse = GenMouse(nForages, fieldSize);
    
    // trap mouse (maybe)
    // recall that isCaught returns a length 2 NumericVector
    // first term is the trap and second term is day which it's caught
    // both terms are initialized to -1, to indicate not caught and both 
    // should be non-negative if the mouse is caught
    NumericVector mouseRes = isCaught(currMouse, Traps, catchRadius);
    
    if (mouseRes(0) >= 0 && mouseRes(1) >= 0) { // we caught the mouse
      trapCountByDay(mouseRes(0), mouseRes(1))++;
    }
  }
  
  return trapCountByDay;
}


// For a given number of square, return a NumericVector containing which square each trap is in
// [[Rcpp::export]]
NumericMatrix GenRingAssignment(int nSquares) {
  NumericMatrix rings(2*nSquares, 2*nSquares);
  
  // std::fill(rings.begin(), rings.end(), nSquares);
  
  for (int i = 0; i < rings.nrow(); i++) {
    for (int j = 0; j < rings.ncol(); j++) {
      //rings(i,j) = nSquares - std::min(i, j); // just top left quadrent
      //rings(i,j) = std::max(i, j) - nSquares + 1; // just bottom right quadrent
      rings(i,j) = std::max(nSquares - std::min(i, j), std::max(i, j) - nSquares + 1);  // All traps
    }
  }
  
  return rings;
}



// Calculate statistics by square
// [[Rcpp::export]]
void processResults(double trapSpacing, double catchRadius, double boarder, int nSquares, double trueDensity, int nForages) {
  NumericMatrix collectData = GenTrapData(trapSpacing, catchRadius, boarder, nSquares, trueDensity, nForages);
  NumericMatrix ringAssignment = GenRingAssignment(nSquares)
}






