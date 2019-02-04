#include <Rcpp.h>
using namespace Rcpp;


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderCentral(NumericVector x, double h = 1) {
  NumericVector out(x.size());
  
  out[0] = NA_REAL;
  for (int i = 1; i < x.size()-1; i++) {
    out[i] = (x[i-1] - (2*x[i]) + x[i+1])/(h*h);
  }
  out[x.size()-1] = NA_REAL;
  
  return out;
}


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderForward(NumericVector x, double h = 1) {
  NumericVector out(x.size());
  
  for (int i = 0; i < x.size()-2; i++) {
    out[i] = (x[i+2] - 2*x[i+1] + x[i])/(h*h);
  }
  out[x.size()-2] = NA_REAL;
  out[x.size()-1] = NA_REAL;
  
  return out;
}


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderBackward(NumericVector x, double h = 1) {
  NumericVector out(x.size());
  
  out[0] = NA_REAL;
  out[1] = NA_REAL;
  for (int i = 2; i < x.size(); i++) {
    out[i] = (x[i] - 2*x[i-1] + x[i-2])/(h*h);
  }
  
  return out;
}

