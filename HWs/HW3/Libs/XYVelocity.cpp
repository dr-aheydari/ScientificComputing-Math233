#include "XYVelocity.h"
#include <cmath>
#include <math.h>
#include <omp.h>

//double velocity_X::operator()(double x, double y) const
//{
//    return -1*y ;
//}

//double velocity_Y::operator()(double x, double y) const
//{
//    return x;
//}



double velocity_X::operator()(double x, double y) const
{
    return sin(y) ;
}

double velocity_Y::operator()(double x, double y) const
{
    return (rand() % 5 - rand() % 5) + 0.25 * cos(x);
}
