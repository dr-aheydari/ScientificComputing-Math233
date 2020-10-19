#include "SemiLagrangian.h"
#include <math.h>
#include <cfloat>
#include <algorithm>
#include <omp.h>


template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

SemiLagrangian::SemiLagrangian(const Grid2D grid, const std::vector<double>& initial_sols, velocity_X Vx, velocity_Y Vy, double dt_)
{
    my_grid = grid;
    // explicitly setting = 0. to avoid any errors
    fn_sols.resize(my_grid.get_N() * my_grid.get_M(), 0.);
    Vx_n.resize(my_grid.get_N() * my_grid.get_M(), 0.);
    Vy_n.resize(my_grid.get_N() * my_grid.get_M(), 0.);
    fn_sols = initial_sols;
    field_x = &Vx;
    field_y = &Vy;
    int N_x = my_grid.get_N();
    int N_y = my_grid.get_M();
    dt = dt_;
#pragma omp parallel for
    for (int n = 0; n < N_x*N_y; n++)
    {
        double x = my_grid.x_from_n(n);
        double y = my_grid.y_from_n(n);
    }



}


///////// Returning the coordinates of deprature using RK2
std::tuple<double,double> SemiLagrangian::trajectory_interpolation(int n)
{
    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);

    double x_star = x - 0.5 * dt * field_x ->operator ()(x,y);
    double y_star = y - 0.5 * dt * field_y ->operator ()(x,y);

    x_star = Clip_x_star(x_star);
    y_star = Clip_y_star(y_star);

    double x_d = x - dt * field_x ->operator ()(x,y);
    double y_d = y - dt * field_y ->operator ()(x,y);

    // to check if our trajectory is not outside the domain
    x_d = Clip_xd(x_d);
    y_d = Clip_yd(y_d);

    return std::make_tuple(x_d,y_d);

}


// if we are just doing Semi-Lagrangian with no reinitialization
std::vector<double> SemiLagrangian::Solve()
{
    std::vector<double> fn_p1;
    fn_p1.resize(my_grid.get_N() * my_grid.get_M());

#pragma omp parallel
    omp_set_num_threads(8);

#pragma omp parallel for
    for (int i = 0; i <my_grid.get_N() * my_grid.get_M(); i++)
    {
        // we will actually get x_d, y_d
        double x,y;
        std::tie(x,y) = trajectory_interpolation(i);
        fn_p1[i] = my_grid.quadraticENO_interpolation(fn_sols, x, y,*field_x,*field_y);
    }

    fn_sols = fn_p1;

    return fn_sols;

}

// if we are doing advection with reinitialization
std::vector<double> SemiLagrangian::LS_Solve(std::vector<double>& fn)
{
    std::vector<double> fn_p1;
    fn_p1.resize(my_grid.get_N() * my_grid.get_M());

#pragma omp parallel
    omp_set_num_threads(8);

#pragma omp parallel for
    for (int i = 0; i <my_grid.get_N() * my_grid.get_M(); i++)
    {
        // we will actually get x_d, y_d
        double x,y;
        std::tie(x,y) = trajectory_interpolation(i);
        fn_p1[i] = my_grid.quadraticENO_interpolation(fn, x, y,*field_x,*field_y);
    }

    return fn_p1;

}



double SemiLagrangian::Clip_x_star(double x_star)
{
    // check that x and y star are in the domain, if not get the max or min with
    // min/max values of the domain
    double xmin = my_grid.get_xmin();
    double xmax = my_grid.get_xmax();

    x_star = std::max(x_star,xmin);
    x_star = std::min(x_star,xmax);

    return x_star;
}

double SemiLagrangian::Clip_y_star(double y_star)
{
    double ymin = my_grid.get_ymin();
    double ymax = my_grid.get_ymax();

    y_star = std::max(y_star,ymin);
    y_star = std::min(y_star,ymax);

    return y_star;
}



double SemiLagrangian::Clip_xd(double x_d)
{
    if (x_d < my_grid.get_xmin())
        x_d = my_grid.get_xmin();
    if (x_d > my_grid.get_xmax())
        x_d = my_grid.get_xmax();

    return x_d;
}


double SemiLagrangian::Clip_yd(double y_d)
{
    if(y_d < my_grid.get_ymin())
        y_d = my_grid.get_ymin();
    if (y_d > my_grid.get_ymax())
        y_d =  my_grid.get_ymax();

    return y_d;
}
