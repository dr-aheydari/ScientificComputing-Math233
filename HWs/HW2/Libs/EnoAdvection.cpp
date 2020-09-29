#include "EnoAdvection.h"
#include <math.h>
#include <omp.h>
#include <cfloat>

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

 Eno_Advection::Eno_Advection(Grid2D grid, std::vector<double>& initial_sols, velocity_X Vx, velocity_Y Vy)
{
     my_grid = grid;
     tn_sols.resize(my_grid.get_N() * my_grid.get_M(), 0.);
     tn_sols = initial_sols;
     field_x = &Vx;
     field_y = &Vy;
     x_sol_tn.resize(my_grid.get_N() * my_grid.get_M(), 0.);
     y_sol_tn.resize(my_grid.get_N() * my_grid.get_M(), 0.);


}

void Eno_Advection::OutputVTK(char* name)
{
    //output to vtk file
    my_grid.initialize_VTK_file(name);
    my_grid.print_VTK_Format(tn_sols, "solution_at_nodes",name);
}


std::vector<double> Eno_Advection::Solve()
{
    double dx = my_grid.get_dx();
    double dt = 0.25 * dx;
//    #pragma omp for
        for(int q=0; q< my_grid.get_N() * my_grid.get_M(); q++)
        {
            Solve_x(q);
            Solve_y(q);
        }

        tn_sols = one_Step(dt);
    return tn_sols;
}
std::vector <double> Eno_Advection::one_Step(double del_t)
 {
        double dx = my_grid.get_dx();
        std::vector<double> tn_p1;
        tn_p1.resize(my_grid.get_M()*my_grid.get_N());
//#pragma omp for
        for(int n = 0; n < my_grid.get_M()*my_grid.get_N(); n++)
        {
            tn_p1[n] = tn_sols[n] + (del_t/(dx*2 + DBL_EPSILON)) * (Solve_x(n) + Solve_y(n));
//            std::cout<<"tn+1: "<<tn_p1[i]<<std::endl;
        }

        return tn_p1;
 }


double Eno_Advection::Solve_x(int n)
{
    double dx = my_grid.get_dx();
    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);
    double neededComponent_x = my_grid.Auto_dx(tn_sols,n,*field_x); //+ (field_x->operator()(x,y)) * my_grid.Auto_dxx(tn_sols,n,*field_x);

    return neededComponent_x;
}

double Eno_Advection::Solve_y(int n)
{
    double dy = my_grid.get_dy();
    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);
    double neededComponent_y = my_grid.Auto_dy(tn_sols,n,*field_y) ;//+ (field_x->operator()(x,y))* my_grid.Auto_dyy(tn_sols,n,*field_y);

    return neededComponent_y;

    // f_j+1/2
    // f_j-1/2
//    double f_plus_half = tn_sols[i+j*my_grid.get_N()] +(sgn(Condition(n,'y')))*my_grid.Auto_dy(tn_sols,n,field_y);
//    double f_minus_half = tn_sols[i - 1 +(j - 1)*my_grid.get_N()] +(sgn(Condition(n,'y')))*my_grid.Auto_dy(tn_sols,n,field_y);

//    double f_star = tn_sols[i+j*my_grid.get_N()] - (dt/dy) *field_y(x,y) * (f_plus_half - f_plus_half);

     // do something to get solution y
}


double Eno_Advection::Condition(int n, char xOry)
{
    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);
    if (xOry == 'x')
    {
        double dx = my_grid.get_dx();
        return (field_x->operator()(x,y) + field_x->operator()(x+dx,y));
    }
    else
    {
        double dy = my_grid.get_dy();
        return (field_x->operator()(x,y) + field_x->operator()(x,y+dy));
    }

}


