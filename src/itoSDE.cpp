#include <Rcpp.h>
#include <random>
#include <chrono>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

std::vector<double> discretize(double a, double b, unsigned int n)
{
  double h = (b-a)/n;
  std::vector<double> y(n+1);
  for(unsigned int i = 0; i< n+1;i++)
  {
    y[i] = a+h*i;
  }
  return y;
}

double runiform(double a, double b)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator (seed);
  //static std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(generator);
}

double rgauss(double mean, double sd)
{
  // Box-Mueller transform
  double u1 = runiform(0, 1);
  double u2 = runiform(0, 1);
  double x = std::sqrt(-2 * std::log(u1));
  double y = std::cos(2 * std::atan(1.0)*4.0*u2);
  double z = x * y;
  return mean + sd * z;
}

//' Solve a constant-coefficient SDE for Ito diffusion
//'
//' @param x0 initial value
//' @param t the time period to solver over
//' @param drift the drift-rate
//' @param volat the volatility coefficient
//' @param n number of variates in time-discretization
//'
//' @description {Solve the SDE for a constant-coefficient Ito diffusion.}
//' @return numeric vector
// [[Rcpp::export]]
std::vector<double> em_ito(double x0, double t, double drift, double volat, unsigned int n)
{
  std::vector<double> x(n+1);
  double h = t/n;
  x[0] = x0;
  double z = 0;
  for(unsigned int i = 1; i < n+1; i++)
  {
    z = rgauss(0, 1);
    x[i] = x[i-1]+h*drift+volat*sqrt(h)*z;
  }
  return x;

}



