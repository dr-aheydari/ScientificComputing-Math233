#include "XYVelocity.h"
#include <cmath>
#include <omp.h>

double velocity_X::operator()(double x, double y) const
{
    return -1*y ;
}

double velocity_Y::operator()(double x, double y) const
{
    return x;
}
