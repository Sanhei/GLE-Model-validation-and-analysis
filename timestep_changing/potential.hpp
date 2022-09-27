#ifndef POTENTIAL
#define POTENTIAL

#include <functional>
#include <math.h>
//extern double U_depth;
double grad_U(double x){
    double gradU = -U_depth*4*(pow(x, 3.0) - x);
    return gradU;
}
#endif
