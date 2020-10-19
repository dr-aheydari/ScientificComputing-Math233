#include "LevelSet.h"
#include <cmath>
#include <omp.h>
#include <algorithm>


template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}


LevelSet::LevelSet(const Grid2D& grid, const std::vector<double> Init_cond, double dt_)
{
    my_grid = grid;
    N_x = my_grid.get_N();
    N_y = my_grid.get_M();
    level_set_0.resize(N_x * N_y);
    level_set_0 = Init_cond;
    level_set_n.resize(N_x * N_y);
    dt = dt_;
}


double LevelSet::operator ()(double x, double y) const
{
    // level set function we want to do!
    // for now we keep it as the same one we did for advection
    return sqrt(std::pow(x - 0.25,2) + std::pow(y,2)) - 0.2;
}

void LevelSet::Assign_ls()
{
#pragma omp parallel for
    for (int n = 0; n < N_x * N_y; n++)
    {
        double x = my_grid.x_from_n(n);
        double y = my_grid.y_from_n(n);
        level_set_n[n] = (*this)(x,y);
    }

}


void LevelSet::Perturb_ls(double pertubation)
{
// to make it independent of main threads
#pragma omp parallel
    omp_set_num_threads(8);

#pragma omp parallel for
    for(int n=0; n<N_x * N_y; n++)
    {
        level_set_n[n] = pertubation * level_set_n[n];
    }

}

std::vector<double> LevelSet::get_level_set_n() const
{
    return level_set_n;
}


void LevelSet::set_level_set_n(std::vector<double> new_level_set)
{
    level_set_n = new_level_set;
}


std::vector<double> LevelSet::Reinitialize()
{
    std::vector<double> level_set_np1(N_x*N_y);
#pragma omp parallel
    omp_set_num_threads(8);

#pragma omp parallel for
    // Gudonov Scheme
    for(int n=0; n < (N_x * N_y)  ; n++)
    {
            // the functions will take care of the boundary stuff
            double dx_m = my_grid.dx_backward(level_set_n, n);
            double dx_p = my_grid.dx_forward(level_set_n, n);

            double dy_m = my_grid.dy_backward(level_set_n, n);
            double dy_p = my_grid.dy_forward(level_set_n, n);

            double phi_0 = level_set_0[n];
            // check for sign of Phi0
            double sgn = (phi_0 > 0) ? 1:-1;
            // check which derivative we should use
            double dx = Godunov_Check(sgn, dx_p, dx_m);
            double dy = Godunov_Check(sgn, dy_p, dy_m);

           level_set_np1[n] = level_set_n[n] + dt*sgn*(1 - sqrt(std::pow(dx,2)+std::pow(dy,2)));
        }

    return level_set_np1;

}


void LevelSet::reinitialize(int iteration = 20)
{
    // call the reinitialization scheme for the given iterations
    for(int it = 0; it < iteration; it++)
    {
        level_set_n = Reinitialize();
    }
}

/////////////////// TO WRITE OUT THE REINIT SOLUTIONS TO A VTK FILE

//void LevelSet::reinitialize(int iteration = 20)
//{
    //int counter = 0;
    //char name[250];

//    std::vector<double> PHI0(N_x,N_y);
//

//    int counter = 0;
//    char name[250];
//    for(int it = 0; it < iteration; it++)
//    {
//        level_set_n = Reinitialize();
//        sprintf(name,"/Users/aliheydari/Box/Course Material/Fall 2020/Scientific Computing/VTK_Sims/HW3/LevelSetTests/N=%i_reinitialized,t=%i.vtk",N_x, counter);
//        my_grid.initialize_VTK_file(name);
//        my_grid.print_VTK_Format(level_set_n, "LevelSet_at_nodes",name);
//        counter ++;
//    }
//}


double LevelSet::Godunov_Check(int sign, double d_p, double d_m) const
{
    if (sign*d_m <= 0 && sign*d_p <= 0)
        return d_p;
    else if (sign*d_m >= 0 && sign*d_p >= 0)
        return d_m;
    else if (sign*d_m <= 0 && sign*d_p >= 0)
        return 0;
    else if (sign*d_m >= 0 && sign*d_p <= 0)
    {
        if (abs(d_m) >= abs(d_p))
            return d_m;
        else if(abs(d_m) <= abs(d_p))
            return d_p;
    }

}

LevelSet::~LevelSet()
{
    std::cout<<"Destructor for Level Set"<<std::endl;
}





