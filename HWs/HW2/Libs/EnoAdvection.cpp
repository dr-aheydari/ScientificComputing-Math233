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

#pragma omp parallel for
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

        std::vector<double> tn_p1;
        tn_p1.resize(my_grid.get_M()*my_grid.get_N());

        for(int n = 0; n < my_grid.get_M()*my_grid.get_N(); n++)
        {
             tn_p1[n] = tn_sols[n] - (del_t) * (Solve_x(n) + Solve_y(n));
        }

        return tn_p1;
 }


double Eno_Advection::Solve_x(int n)
{

    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);
    double neededComponent_x =  (field_x->operator()(x,y))* my_grid.Auto_dx(tn_sols,n,*field_x);// + (field_x->operator()(x,y))* my_grid.Auto_dxx(tn_sols,n,*field_x);

    return neededComponent_x;
}

double Eno_Advection::Solve_y(int n)
{

    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);
    double neededComponent_y = (field_y->operator()(x,y))* my_grid.Auto_dy(tn_sols,n,*field_y);// + (field_y->operator()(x,y)) * my_grid.Auto_dyy(tn_sols,n,*field_y);
    return neededComponent_y;

}


double Eno_Advection::Condition(int n, char xOry)
{
    double x = my_grid.x_from_n(n);
    double y = my_grid.y_from_n(n);

    if (xOry == 'x')
    {
        std::cout<<field_x->operator()(x,y)<<std::endl;
        return (field_x->operator()(x,y));
    }
    else
    {
        return (field_x->operator()(x,y));
    }

}


