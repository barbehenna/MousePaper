// Copyright (c) 2017 Alton Barbehenn

// This is the same simulation as I have done in R, but converted to c++ so that I can compile
// it ahead of time and see large speed-ups, becuase c++ is faster than R.

#include <Rcpp.h>
using namespace Rcpp;

// returns a random double according to the uniform distribution on [0,1]
// double random_uniform() {
//   return((double)rand()/RAND_MAX);
// }

// returns a random double according to the standard normal distribution
// double random_norm() {
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::normal_distribution<float> dist(0.0, 1.0);
//   return(dist(gen));
// }

// Assume that no mouse during a visit can be caught by multiple traps
// [[Rcpp::export]]
NumericVector Traps(NumericMatrix trapCoords, NumericMatrix visits, int low, int nv, double delta){
  double dx, dy;
  bool caught = false;
  NumericVector mouse(2);
  
  for (int i = low; i < low + nv; i++) { //for a specific visit
    for (int j = 0; j < trapCoords.nrow(); j++) { // and each trap
      //calculate the distance to trap
      dx = sqrt(pow(trapCoords(j,0)-visits(i,0), 2.0));
      dy = sqrt(pow(trapCoords(j,1)-visits(i,1), 2.0));
      if (dx <= delta && dy <= delta) {
        mouse[0] = j;
        mouse[1] = i%nv;
        caught = true;
        break;
      }
    }
    if (caught) {
      break;
    }
  }
  
  if (!caught) {
    mouse[0] = NA_REAL;
    mouse[1] = NA_REAL;
  }
  
  //return(mouse);
  return mouse;
}

// [[Rcpp::export]]
NumericMatrix rngCpp(const int N) {
  NumericMatrix X(N, 4);
  X(_, 0) = runif(N);
  X(_, 1) = rnorm(N);
  X(_, 2) = rt(N, 5);
  X(_, 3) = rbeta(N, 1, 1);
  return X;
}

// [[Rcpp::export]]
NumericMatrix trapSim1(double ts, double fs, double np, double delta, int nv) {
  // ------ Generate Field Parameters ------
  // Constant time, maybe pass into function and only calculate once?
  double b[8];
  NumericMatrix trapCoords(64, 2); //64 rows, 2 columns (gx and gy: trap locations)
  
  for (int i = 0; i < 8; i++) {
    b[i] = (ts/2)*(2*i-7);
  }
  
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      trapCoords(8*i+j, 0) = b[j];
      trapCoords(8*i+j, 1) = b[i];
    }
  }
  
  // ------ Generate Mouse Visit Locations ------
  // Time scales linearly with number of mice (so quadratically with field size)
  double x, y;
  NumericMatrix visits(np*nv, 2); //where each mouse visits
  for (int mouse = 0; mouse < np; mouse++) {
    // x = fs * (0.5 - random_uniform());
    // y = fs * (0.5 - random_uniform());
    x = fs * (0.5 - runif(1)(0));
    y = fs * (0.5 - runif(1)(0));
    for (int trip = 0; trip < 4; trip++) {
      // visits(mouse*4+trip, 0) = x + random_norm();
      // visits(mouse*4+trip, 1) = y + random_norm();
      visits(mouse*4+trip, 0) = x + rnorm(1)(0);
      visits(mouse*4+trip, 1) = y + rnorm(1)(0);
    }
  }
  
  // ------ Figure Out if Each Mouse is Trapped and Where ------
  // Time scales linearly with number of mice (so quadratically with field size)
  // For now, traps are fixed
  int low;
  NumericVector TrapDay;
  NumericMatrix Catches(np, 2); // (trap, day) pairs
  for (int m = 0; m < np; m++) {
    low = m*nv;
    TrapDay = Traps(trapCoords, visits, low, nv, delta);
    Catches(m,0) = TrapDay[0];
    Catches(m,1) = TrapDay[1];
  }
  
  return(Catches);
}


