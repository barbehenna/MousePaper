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
// Matrix Output
// [[Rcpp::export]]
NumericMatrix GenRingAssignmentMat(int nSquares) {
  NumericMatrix rings(2*nSquares, 2*nSquares);
  
  for (int i = 0; i < rings.nrow(); i++) {
    for (int j = 0; j < rings.ncol(); j++) {
      //rings(i,j) = nSquares - std::min(i, j); // just top left quadrent
      //rings(i,j) = std::max(i, j) - nSquares + 1; // just bottom right quadrent
      rings(i,j) = std::max(nSquares - std::min(i, j), std::max(i, j) - nSquares + 1);  // All traps
    }
  }
  
  return rings;
}


// For a given number of square, return a NumericVector containing which square each trap is in
// Vector Output
// Don't need to worry too much about the order of indexing becasue the ring lables are invarient
// to flips and rotations.
// [[Rcpp::export]]
NumericVector GenRingAssignmentVec(int nSquares) {
  NumericVector rings(4*nSquares*nSquares);
  
  for (int i = 0; i < 2*nSquares; i++) { //2*nSquares = sqrt(length(rings))
    for (int j = 0; j < 2*nSquares; j++) {
      //rings(i,j) = nSquares - std::min(i, j); // just top left quadrent
      //rings(i,j) = std::max(i, j) - nSquares + 1; // just bottom right quadrent
      //rings(i,j) = std::max(nSquares - std::min(i, j), std::max(i, j) - nSquares + 1);  // All traps (matrix form)
      rings((2*nSquares*i) + j) = std::max(nSquares - std::min(i, j), std::max(i, j) - nSquares + 1);
    }
  }
  
  return rings;
}


// Sums the catchs at each trap into the sum of catches in the first half and the sum in the second half
// Sufficiently messy and maybe something we want to change up in the future.
// [[Rcpp::export]]
NumericMatrix calcPeriodsByTrap(NumericMatrix catchData, int nForages) {
  NumericMatrix periodSums(catchData.nrow(), 2);
  
  int cutoff = nForages/2;
  // std::cout << "cutoff: " << cutoff << " using nForages: " << nForages <<  std::endl;
  
  // Sum across first half of forages
  NumericMatrix::Column pd1 = periodSums(_,0);
  for (int day = 0; day < cutoff; day++) {
    pd1 = pd1 + catchData(_,day);
  }
  
  // Sum across first half of forages
  NumericMatrix::Column pd2 = periodSums(_,1);
  for (int day = cutoff; day < catchData.ncol(); day++) {
    pd2 = pd2 + catchData(_,day);
  }
  
  return periodSums;
}



// Calculate statistics by square
// [[Rcpp::export]]
NumericMatrix ProcessResults(int uuid, int paramset, double trapSpacing, NumericMatrix collectData) {
  const int nForages = collectData.ncol();
  const int nSquares = std::sqrt(collectData.nrow())/2;
  NumericVector ringAssignment = GenRingAssignmentVec(nSquares);
  
  // stats by column (in  order): uuid, paramset, square, pd1, pd2, pHat, nHat, aHat, dHat
  NumericMatrix Stats(nSquares, 9); 
  
  // Record the simulation's uuid
  // Reference the second column
  // Changes propagate to xx (same applies for Row)
  NumericMatrix::Column uuidCol = Stats(_,0);
  std::fill(uuidCol.begin(), uuidCol.end(), uuid);
  
  // Record the simulation's parameter set
  NumericMatrix::Column paramsetCol = Stats(_,1);
  std::fill(paramsetCol.begin(), paramsetCol.end(), paramset);
  
  // Record the square values
  for (int square = 0; square < nSquares; square++) {
    Stats(square,2) =  square+1;
  }
  
  // Calculate the number of mice caught in the first half and second half of forages
  NumericVector periods = calcPeriodsByTrap(collectData, nForages);
  
  // record sum of periods by square (pd1 and pd2)
  for (int trap = 0; trap < collectData.nrow(); trap++) { // for each trap
    for (int square = ringAssignment(trap)-1; square < nSquares; square++) { // for each square at least as big
      // rings are [1:n],  indexes are [0:n-1] 
      // ex: ring 1 in in every square, ring 4 is only in rings >=4
      // increment count of square and period by traps' value
      Stats(square,3) = Stats(square,3) + periods(trap,0);
      Stats(square,4) = Stats(square,4) + periods(trap,1);
    }
  }
  
  //
  // Now that the data is fully aggregated (for our purposes), we can compute our main statistics
  //
  
  // pHat
  NumericMatrix::Column pHat = Stats(_,5);
  pHat = 1 - sqrt(Stats(_,4) / Stats(_,3)); // 1-sqrt(pd2/pd1)
  // if pd1 = 0 and pd2 > 0, this implementation returns -Inf
  // if pd1 = pd2 = 0, this implementation returns NaN

  // nHat
  NumericMatrix::Column nHat = Stats(_,6);
  nHat = Stats(_,3)*Stats(_,3)/(Stats(_,3) - Stats(_,4)); // pd1^2 / (pd1 - pd2)
  // if pd1 = pd2 = 0, this implementation return NaN
  // if pd1 = pd2 > 0, this implementation returns Inf
  // if 0 < pd1 < pd2, this implementation returns a negative nHat
  
  // aHat
  NumericMatrix::Column aHat = Stats(_,7);
  aHat = 4*Stats(_,2)*Stats(_,2)*trapSpacing*trapSpacing; // (2 * square * ts)^2
  
  // dHat
  NumericMatrix::Column dHat = Stats(_,8);
  dHat = Stats(_,6) / Stats(_,7); // nHat/aHat
  // if pd1 = pd2 = 0, i.e. nHat = NaN, this implementation returns NaN
  // if pd1 = pd2 > 0, i.e. nHat = Inf, this implementation returns Inf
  // if 0 < pd1 < pd2, i.e. nHat < 0, this implementation returns a negative dHat
  
  return Stats;
}


// Check that the parameters input are valid and meaningful for the expirament
// type errors will be handled at runtime and cause an error
// [[Rcpp::export]]
bool checkParameters(double trapSpacing, double catchRadius, double boarder, int nSquares, double trueDensity, int nForages) {
  if (trapSpacing <= 0) {
    std::cout << "Invalid trapSpacing " << trapSpacing << " (<0)." << std::endl;
    return false;
  }
  
  if (catchRadius <= 0) {
    std::cout << "Invalid catchRadius " << catchRadius << " (<0)." << std::endl;
    return false;
  }
  
  if (boarder <= 0) {
    std::cout << "Invalid boarder " << boarder << " (<0)." << std::endl;
    return false;
  }
  
  if (nSquares <= 0) {
    std::cout << "Invalid nSquares " << nSquares << " (<0)." << std::endl;
    return false;
  }
  
  if (trueDensity <= 0) {
    std::cout << "Invalid trueDensity " << trueDensity << " (<0)." << std::endl;
    return false;
  }
  
  if (nForages <= 0) {
    std::cout << "Invalid nForages " << nForages << " (<0)." << std::endl;
    return false;
  }
  
  if ((nForages % 2) != 0) {
    std::cout << "Invalid nForages " << nForages << " (odd)." << std::endl;
    return false;
  }
  
  return true;
}


// Function to run simulation
// Separated to allow R to do both parts independently (to save or examine raw catch data)
//[[Rcpp::export]]
NumericMatrix RunSimulation(int uuid, int paramset, double trapSpacing, double catchRadius, double boarder, int nSquares, double trueDensity, int nForages) {
  // check if the inputs will cause an obvious error
  if (checkParameters(trapSpacing, catchRadius, boarder, nSquares, trueDensity, nForages)) {
    NumericMatrix collectData = GenTrapData(trapSpacing, catchRadius, boarder, nSquares, trueDensity, nForages);
    NumericMatrix stats = ProcessResults(uuid, paramset, trapSpacing, collectData);
    
    return stats;
  } else {
    return NumericMatrix(0,0);
  }
}


// For a given single mouse's forage coordinates, go through all traps and 
// determine which (if any) trap caught the mouse for a given forage. Return -1 
// if the mouse is not caught and the traps index of the trap that caught the 
// mouse if it was caught.
// The general algorithm (as implemented currently) does not implment recapture
// of the mouse, but could be extended to do so by returning an array of indicies
// (one per forage) of the traps the mouse was caught in.
// [[Rcpp::export]]
NumericVector isCaught2(NumericMatrix forages, NumericMatrix Traps, double catchRadius) {
  std::vector<int> caughtin;
  NumericVector result(forages.nrow()); // every entry is the trap that caught that mouse on that forage
  double dx, dy;
  
  // result = -1;
  
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
        // std::cout << "catch ->  day " << day << " trap " << trap << ": " << Traps(trap,0) << ", " << Traps(trap,1) << std::endl;
      }
    }
    
    // If the mouse was caught in a trap on this forage
    // record trap that caught the mouse
    if (!caughtin.empty()) {
      result(day) = caughtin[(int) runif(1, 0, caughtin.size())(0)];
    } else { // or record -1 for not caught
      result(day) = -1;
    }
    
    // reset caughtin for the next day
    caughtin.clear();
  }
  
  return result;
}


// This function generates all of the mice used in the simulation and returns a
// matrix (rows: mice, columns: forage) containing the trap the mouse is caught in.
// if the mouse is not caught, it's trap is -1. 
// [[Rcpp::export]]
NumericMatrix GenAllMice(double trapSpacing, double catchRadius, double boarder, int nSquares, double trueDensity, int nForages) {
  double fieldSize = (((2 * nSquares) - 1) * trapSpacing) + (2 * boarder);
  int nmice = std::round(fieldSize * fieldSize * trueDensity);
  
  // Initialize all mice martix
  // one row per mouse
  // columns are traps mouse is caught in
  NumericMatrix Mice(nmice, nForages);
  
  // Generate (place) traps
  NumericMatrix Traps = GenTraps(nSquares, trapSpacing);
  
  // generate mice
  for (int mouse = 0; mouse < nmice; mouse++) {
    // generate the mouse
    NumericMatrix currMouse = GenMouse(nForages, fieldSize);
    // record traps the mouse ends up in
    Mice(mouse,_) = isCaught2(currMouse, Traps, catchRadius);
  }
  
  return Mice;
}


// Taking in the mice matrix (rows: mouse, columns: trap caught in by forage),
// convert the data to a matrix (rows: traps, columns: number of mice caught on forage).
// The result of this function is a matrix where each entry is the number of mice first
// caught in a particular trap on a particular day, counting only the first catch. 
// [[Rcpp::export]]
NumericMatrix MiceDataToTrapData(NumericMatrix MiceData, int nSquares, int nForages) {
  int nTraps = 4 * nSquares * nSquares; // (2*nSquares)^2

  // Data collection data-structure
  // One row per trap and one column for catch opprotunity
  NumericMatrix trapCountByDay(nTraps, nForages); // initialized to 0
  
  // Convert to by trap and day counts
  for (int mouse = 0; mouse < MiceData.nrow(); mouse++) {
    for (int forage = 0; forage < nForages; forage++) {
      if (MiceData(mouse, forage) >= 0) { //caught the mouse
        trapCountByDay(MiceData(mouse, forage), forage)++;
        break;
      }
    }
  }
  
  return trapCountByDay;
}


// Taking in the traps each mouse is caught in on each of its forages,
// estimate the population size via catch and release. To aggregate the 
// observations, pass in the cut-off day (observation 1 < cut-off day, 
// observation 2 >= cut-off day). This function returns a single estimate
// of the total population size. 
