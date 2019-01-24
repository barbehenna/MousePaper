#include <Rcpp.h>
using namespace Rcpp;



// This function returns a matrix of trap coordinates, as such, the returned
// matrix should be a nx2 matrix, (x,y) pairs
// [[Rcpp::export]]
NumericMatrix GenerateTraps(int nrings = 8, double trapspacing = 1.0) {
  const int trapsperside = 2*nrings;
  double b[trapsperside];
  NumericMatrix trapCoordinates(std::pow(trapsperside, 2), 2);
  
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

