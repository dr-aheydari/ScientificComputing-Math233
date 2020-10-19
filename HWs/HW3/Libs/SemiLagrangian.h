#ifndef SEMILAGRANGIAN_H
#define SEMILAGRANGIAN_H

#include "grid2d.h"
#include "XYVelocity.h"
#include <vector>
#include <cmath>
#include <tuple>
class SemiLagrangian
{
private:

        Grid2D my_grid;
        std::vector<double> init_sols, fn_sols,Vx_n,Vy_n;
        velocity_X* field_x;
        velocity_Y* field_y;
        double dt;
public:
    SemiLagrangian();
    SemiLagrangian(const Grid2D grid, const std::vector<double>& initial_sols, velocity_X Vx, velocity_Y Vy, double dt_);
    std::tuple<double,double> trajectory_interpolation(int n);
    std::vector<double> Solve();
    std::vector<double> LS_Solve(std::vector<double>& fn);

    double Clip_xd(double x_d);
    double Clip_yd(double y_d);
    double Clip_x_star(double x_star);
    double Clip_y_star(double y_star);

};
#endif

