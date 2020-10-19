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
    double first_part = std::pow(x/7,2) * sqrt(abs(abs(abs(x)-3)/(abs(x) - 3))) + std::pow(y/3,2) * sqrt(abs(y+(3*sqrt(33)/7))/(y+3*sqrt(33)/7));
    double second_part = abs(x/2) - ((3*sqrt(33) - 7)/122)*std::pow(x,2) - 3 + sqrt(1-std::pow(abs(abs(x) - 2) - 1,2)) - y;

    double third_part = (9 * sqrt(    abs((abs(x) - 1) * (abs(x)-0.75))   /   ((1-abs(x))*(abs(x)-0.75))   ) - 8*abs(x) - y);

    double fourth_part = (3*abs(x) + 0.75*sqrt(    (abs((abs(x) - 0.75) * (abs(x) - 0.5)))
                                                   / ((0.75 - abs(x)) * (abs(x) - 0.5))
                                                   )

                                                   - y);

    double fifth_part = (2.25 * sqrt((abs((x - 0.5) * (x + 0.5))) / ((0.5 - x)*(0.5 + x))) - y);

    double sixth_part = (((6*sqrt(10))/7) + (1.5 - 0.5 * abs(x)) * sqrt(abs(abs(x) - 1)/(abs(x) - 1)) - ((6*sqrt(10))/14) * sqrt(4 - std::pow((abs(x) - 1),2)) - y);

//    if (first_part * second_part * third_part * fourth_part * fifth_part * sixth_part <= 0)
    if (third_part == 0)
        std::cout<<  third_part << std::endl;

    if (first_part * second_part <= 0 && first_part * second_part >= -3)
        return 1;
    else if ( sixth_part >= -2 && sixth_part <= 0)
        return 1;
    else if ( fourth_part > 2 && fourth_part <= 0)
        return 1;
//    else if ( third_part > -2 && third_part <= 0)
//        return 1;
    else if ( fifth_part > -2 && fifth_part <= 0)
        return 1;
    else
        return 0;
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

//#pragma omp parallel for collapse(2)
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
//#pragma omp parallel for collapse(2)
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

