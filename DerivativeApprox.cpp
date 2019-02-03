#include <Rcpp.h>
using namespace Rcpp;


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderCentral(NumericVector x, double h = 1) {
  NumericVector out(x.size()-2);
  
  for (int i = 1; i < x.size()-1; i++) {
    out[i-1] = (x[i+1] - 2*x[i] + x[i-1])/(h*h);
  }
  
  return out;
}


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderForward(NumericVector x, double h = 1) {
  NumericVector out(x.size()-2);
  
  for (int i = 0; i < x.size()-2; i++) {
    out[i-1] = (x[i+2] - 2*x[i+1] + x[i])/(h*h);
  }
  
  return out;
}


// vector must be of at least length 3 
// assuming constant step size
// [[Rcpp::export]]
NumericVector SecondOrderBackward(NumericVector x, double h = 1) {
  NumericVector out(x.size()-2);
  
  for (int i = 2; i < x.size(); i++) {
    out[i-1] = (x[i] - 2*x[i-1] + x[i-2])/(h*h);
  }
  
  return out;
}

