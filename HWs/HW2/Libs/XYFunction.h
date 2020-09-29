#ifndef XYFUNCTION_H
#define XYFUNCTION_H
#include <vector>
#include <cmath>
#include "XYVelocity.h"
#include "grid2d.h"

class xyFunction:CF_2
{
private:
    int N_x,N_y;
    double x_min,x_max,y_min,y_max;
    Grid2D my_grid;
   public:
    xyFunction(Grid2D grid);
    double operator()(double x, double y) const;
    std::vector<double> Assign();
    std::vector<double> InitialCondition();
    ~xyFunction(){};

};
#endif



