#include "libs.h"
#include <math.h>
#include <omp.h>

double velocity_X::operator()(double x, double y) const
{
    return 10*x ;
}

double velocity_Y::operator()(double x, double y) const
{
    return 16*y;
}

double xyFunction::operator ()(double x, double y) const
{
    return 12*x + 16*y;
}

void xyFunction::Assign(std::vector<double>& funct, double x_max, double y_max)
{
    // vector of values defined at the nodes of the grid
    funct.resize(x_max*y_max);
    // fill with 0.
    std::fill(funct.begin(), funct.end(), 0.0);
#pragma omp parallel for collapse(2)
    for (int i=0; i<x_max; i++)
        for (int j=0; j < y_max; j++)
        {
            double x  = double(i)/x_max;
            double y  = double(j)/y_max;
            funct[i+ j*y_max]=(*this)(x,y);
        }
}
