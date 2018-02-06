#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

// [[Rcpp::export]]
int sign(double x) {
  if(x > 0) {
    return 1;
  }
  else if (x < 0) {
    return -1;
  }
  else {
    return 0;
  }
}


// [[Rcpp::export]]
int nearTrap_Square(NumericVector curr, NumericMatrix trapCoords, double delta) {
  double dx, dy;
  
  for (int i = 0; i < trapCoords.nrow(); i++) {
    //calculate the distance to trap
    dx = sqrt(pow(curr(0) - trapCoords(i,0), 2.0));
    dy = sqrt(pow(curr(1) - trapCoords(i,1), 2.0));
    if (dx <= delta && dy <= delta) {
      return i;
    }
  }
  return -1;
}


// location is the mean (-pi, pi)
// concentration is similar to the 1/variance (0, inf)
// [[Rcpp::export]] 
double von_mises(double location, double concentration) {
  double tau, rho, r;
  double u1, u2, u3;
  double c, f, z;
  
  if (concentration == 0.0) { // concentration > 0
    return runif(1, -M_PI, M_PI)[0];
  }
  
  tau = 1.0 + sqrt(1.0 + 4.0*concentration*concentration);
  rho = (tau - sqrt(2.0*tau)) / (2.0*concentration);
  r = (1.0 + rho*rho) / (2.0*rho);
  
  while (1) {
    u1 = runif(1)[0];
    z = cos(M_PI*u1);
    f = (1.0 + r*z)/(r + z);
    c = concentration*(r - f);
    
    u2 = runif(1)[0];
    
    if (u2 <= c*(2.0 - c)) {
      break;
    }
    if (c <= log(c/u2) + 1.0) {
      break;
    }
  }
  
  u3 = runif(1)[0];
  
  return location + sign(u3-0.5)*acos(f);
}


// [[Rcpp::export]]
NumericVector step(NumericVector curr, NumericVector init, double walksize) {
  double angle_mean, angle_est, dist;
  double dx, dy;
  NumericVector out(2);
  //std::cout << "Curr = " << curr(0) << " " << curr(1) << std::endl;
  
  // dx = init(0) - curr(0);
  // dy = init(1) - curr(0);
  dx = curr(0) - init(0);
  dy = curr(1) - init(1);

  angle_mean = atan2(dy, dx) + M_PI;
  dist = sqrt(dx*dx + dy*dy);
  angle_est = von_mises(angle_mean, dist);

  out(0) = curr(0) + walksize*cos(angle_est);
  out(1) = curr(1) + walksize*sin(angle_est);
  //std::cout << "new curr = " << curr(0) << " " << curr(1) << std::endl;
  
  return out;
}

//If I want to export step() I think I need to make it so that it doesn't modify x and y in place, 
//rather I'll need to return an array/vector containing and x and a y component

// [[Rcpp::export]]
NumericMatrix trap_rw(double ts, double fs, double np, double delta, int nv, double stepSize = 0.25, int maxiter=1000) {
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
      trapCoords(8*i+j, 1) = -1*b[i];
    }
  }
  
  // ------ Generate Mouse Walks ------
  // Time scales linearly with number of mice (so quadratically with field size)
  NumericVector init(2); // mouse home (x0, y0)
  NumericVector curr(2); // forage locations (xc, yc)
  int trap;
  bool caught;
  NumericMatrix Catches(np,3); // (trap, day, iter) pairs
  for (int mouse = 0; mouse < np; mouse++) {
    init(0) = fs * (0.5 - runif(1)(0));
    init(1) = fs * (0.5 - runif(1)(0));
    caught = false;
    
    Catches(mouse, 0) = NA_REAL;
    Catches(mouse, 1) = NA_REAL;
    Catches(mouse, 2) = NA_REAL;
    
    for (int day = 0; day < 4; day++) { // I assume one forage per day
      if (caught) {
        break;
      }
      
      curr(0) = init(0);
      curr(1) = init(1);
      for (int iter = 0; iter < maxiter; iter ++) {
        curr = step(curr, init, stepSize);
        trap = nearTrap_Square(curr, trapCoords, delta);
        if(trap > -1) {
          Catches(mouse, 0) = trap;
          Catches(mouse, 1) = day;
          Catches(mouse, 2) = iter;
          caught = true;
          break;
        }
      }
    }
  }
  
  return Catches;
}


