#include "XYFunction.h"
#include <cmath>
#include <omp.h>

xyFunction::xyFunction(Grid2D grid)
{
    my_grid = grid;
    N_x = grid.get_N();
    N_y = grid.get_M();
    std::tie(x_min,x_max) = grid.get_x_range();
    std::tie(y_min,y_max) = grid.get_y_range();

}


double xyFunction::operator ()(double x, double y) const
{
    // initial condition sqrt{(x-0.25)^2 + y^2} - 0.2

    return sqrt(std::pow(x - 0.25,2) + std::pow(y,2))- 0.2;
}

// assinging the function values to all the nodes in the mesh
std::vector<double> xyFunction::Assign()
{
    // vector of values defined at the nodes of the grid
    std::vector<double> funct;
    funct.resize(N_x*N_y);
    // fill with 0.
    std::fill(funct.begin(), funct.end(), 0.0);

// configure OMP stuff if we forget in the main
#pragma omp parallel
    omp_set_num_threads(8);

#pragma omp parallel for collapse(2)
    for (int i=0; i<N_x; i++)
        for (int j=0; j < N_y; j++)
        {
            int n = my_grid.n_from_ij(i, j);
            double x  = my_grid.x_from_n(n);
            double y  = my_grid.y_from_n(n);
            funct[i+ j*N_y]=(*this)(x,y);
        }
    return funct;
}


//for this problem we do not have the initial condition
std::vector<double> xyFunction::InitialCondition()
{
    // vector of values defined at the nodes of the grid
    std::vector<double> init_cond;
    init_cond.resize(N_x*N_y);
    // fill with 0.
    std::fill(init_cond.begin(), init_cond.end(), 0.0);
#pragma omp parallel for collapse(2)
    for (int i=0; i<N_x; i++)
        for (int j=0; j < N_y; j++)
        {
            int n = my_grid.n_from_ij(i, j);
            double x  = my_grid.x_from_n(n);
            double y  = my_grid.y_from_n(n);
            if ((*this)(x,y) <= 0)
                init_cond[i+ j*N_y]= 1;
            else
                init_cond[i+ j*N_y]= 0;
        }
    return init_cond;
}

