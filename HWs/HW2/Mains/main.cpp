#include <iostream>
#include <Libs/XYVelocity.h>
#include <Libs/grid2d.h>
#include <Libs/EnoAdvection.h>
#include <Libs/XYFunction.h>
#include <omp.h>

int main()
{

    velocity_X V_x;
    velocity_Y V_y;
    int N_x = 200;
    int N_y = 200;
    double x_min =-1;
    double x_max = 1;
    double y_min = -1;
    double y_max = 1;

    //Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
    Grid2D my_grid(N_x,N_y,x_min,x_max,y_min,y_max);
    xyFunction F_xy(my_grid);
    std::cout<<"grid resolution in x :"<<my_grid.get_dx()<<std::endl;
    std::vector<double> IC = F_xy.InitialCondition();
    /// /// /// /// /// ///
    ///FOR NUMBER 13
//    std::vector<double> IC = F_xy.Assign();
    /// /// /// /// /// ///
//    my_grid.display(IC);
    std::cout<<"---------------------------------" << std::endl;


//    my_grid.display(IC);
    Eno_Advection sol(my_grid, IC, V_x, V_y);

//    // create a name for the vtk file, with the grid size in it

    double t = 0.;
    double t_max = 2*acos(0.0);
    double dt = my_grid.get_dx() * 0.25;
    // for testing
    t_max = 100*dt;

    std::vector<double> tn_sols;
    tn_sols.resize(N_x*N_y);
    while (t < 20*dt)
    {
//        std::cout<<dt<<std::endl << std::endl;
        char name[250];
        sprintf(name,"/Users/aliheydari/Documents/Scientific_Computing/VTK_Tests/t=%f.vtk",t);
        my_grid.initialize_VTK_file(name);
        if (t == 0)
            tn_sols = IC;
        else
            tn_sols = sol.Solve();

        my_grid.display(tn_sols);
        std::cout<<"---------------------------------" << std::endl;
        t += dt;
        my_grid.print_VTK_Format(tn_sols, "solution_at_nodes",name);
    }

    return 0;
    }

/*
// for level set section (Question 13)
std::cout<<IC[40]<<std::endl;
std::vector<double> funct = F_xy.Assign();

*/

