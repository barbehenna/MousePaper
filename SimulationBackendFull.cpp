#include <Rcpp.h>
#include <vector>
#include <iostream>
using namespace Rcpp;


// For a given single mouse's forage coordinates, go through all traps and 
// determine which (if any) trap caught the mouse. Return -1 if the mouse is
// not caught and the traps index of the trap that caught the mouse if it was
// caught.
// The general algorithm (as implemented currently) does not implment recapture
// of the mouse, but could be extended to do so by returning an array of indicies
// (one per forage) of the traps the mouse was caught in.
// [[Rcpp::export]]
int isCaught(NumericMatrix forages, NumericMatrix Traps, double catchRadius) {
  std::vector<int> caughtin;
  double dx, dy;
  
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
        std::cout << "catch! ->  day " << day << " trap " << trap << ": " << Traps(trap,0) << ", " << Traps(trap,1) << std::endl;
      }
    }
    
    // If the mouse was caught in a trap on this forage
    if (!caughtin.empty()) {
      return caughtin[(int) runif(1, 0, caughtin.size())(0)];
    }
  }
  
  return -1;
}


// Generates a single mouse in a field of spcified size and returns all of the
// locations the mouse forages (number of forages given by nforages).
// In the returned matrix, each row is a forage coordinate (x,y).
// [[Rcpp::export]]
NumericMatrix GenMouse(int nforages, double fieldSize) {
  NumericMatrix mouse(nforages, 2); // forage locations for a single mouse
  
  NumericVector home = runif(2); // x and y home coordinates
  home = (home - 0.5) * fieldSize;
  
  mouse(_,0) = rnorm(nforages) + home[0]; // generate the forages
  mouse(_,1) = rnorm(nforages) + home[1];

  return mouse;
}


// This function returns a matrix of trap coordinates, as such, the returned
// matrix should be a nx2 matrix, (x,y) pairs
// [[Rcpp::export]]
NumericMatrix GenTraps(int nrings = 8, double trapspacing = 1.0) {
  const int trapsperside = 2*nrings;
  double b[trapsperside];
  NumericMatrix trapCoordinates(trapsperside*trapsperside, 2);
  
  for (int i = 0; i < trapsperside; i++) {
    b[i] = (trapspacing/2) * ((2*i) - (trapsperside - 1));
  }
  
  for (int i = 0; i < trapsperside; i++) {
    for (int j = 0; j < trapsperside; j++) {
      trapCoordinates((trapsperside*i)+j, 0) = b[j];
      trapCoordinates((trapsperside*i)+j, 1) = -b[i];
    }
  }
  
  return trapCoordinates;
}

