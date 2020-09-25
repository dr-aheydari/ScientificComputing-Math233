#include "libs.h"
#include <math.h>


double velocity_X::operator()(double x, double y) const
{
    return 10*x ;
}

double velocity_Y::operator()(double x, double y) const
{
    return 16*y;
}
