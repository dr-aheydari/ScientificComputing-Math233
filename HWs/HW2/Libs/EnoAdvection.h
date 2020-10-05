#ifndef ENOADVECTION_H
#define ENOADVECTION_H
#include <vector>
#include <cmath>
#include "grid2d.h"
#include "XYVelocity.h"

class Eno_Advection
{
private:
    Grid2D my_grid;
    std::vector<double> init_sols;
    std::vector<double> tn_sols;

    std::vector<double> x_sol_tn;
    std::vector<double> y_sol_tn;

    velocity_X* field_x;
    velocity_Y* field_y;
    double t_max;
    char flag;

public:
    Eno_Advection(Grid2D grid, std::vector<double>& initial_sols, velocity_X Vx, velocity_Y Vy);
//    void one_Step(double del_t);
    std::vector<double> one_Step(double del_t);
    void OutputVTK(char* name);
    double RK_4(double x0, double y0, double x, double h);
    double Condition(int n, char xOry);
    std::vector<double> Solve();
    double Solve_x(int n);
    double Solve_y(int n);
    ~Eno_Advection(){};
};
#endif



