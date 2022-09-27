#include "potential.hpp"
#include <cmath>
#include <random>
// A Guassion distribution term
std::default_random_engine generator;
std::normal_distribution<double>distribution(0.0, 1.0);
extern const double tau_D, tau_m, tau;
//extern double h;
double x_update(double x, double xp, double y, double h){
    return 2.0*x - xp + tau_D/tau_m*grad_U(x)*h*h + pow(tau_D, 2.0)/(tau_m*tau)*(y-x)*h*h;
}
double y_update(double y, double x, double h){
   return y-tau_D*(y-x)*h/tau + sqrt(2.0*h)*distribution(generator);
}


