#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <Libs/XYVelocity.h>
#include <Libs/grid2d.h>
#include <Libs/EnoAdvection.h>
#include <Libs/XYFunction.h>
#include <Libs/VectorUtils.h>
#include <omp.h>

int main()
{
    
    velocity_X V_x;
    velocity_Y V_y;
    int N_x = 256;
    int N_y = 256;
    double x_min =-1;
    double x_max = 1;
    double y_min = -1;
    double y_max = 1;
    
    Grid2D my_grid(N_x,N_y,x_min,x_max,y_min,y_max);
    xyFunction F_xy(my_grid);
    std::cout<<"grid resolution in x :"<<my_grid.get_dx()<<std::endl;
    /// /// /// /// /// ///
    ///FOR NUMBER 11
    //    std::vector<double> IC = F_xy.InitialCondition();
    /// /// /// /// /// ///
    ///
    /// /// /// /// /// /// ///
    ///FOR NUMBER 13
    std::vector<double> IC = F_xy.Assign();
    /// /// /// /// /// ///
    my_grid.display(IC);
    std::cout<<"---------------------------------" << std::endl;
    
    
    Eno_Advection sol(my_grid, IC, V_x, V_y);
    
    double t = 0.;
    double t_max = 4*acos(0.0);
    /* for convergence analysis in space we need this dt
     //double dt = std::pow(my_grid.get_dx(),2);
     */
    double dt = 0.25 * my_grid.get_dx();
    std::cout<<"dt :"<< dt << std::endl;
    
    // for testing
    std::vector<double> tn_sols;
    tn_sols.resize(N_x*N_y);
    
    int counter = 0;
    
    while (t < t_max)
    {
        char name[250];
        sprintf(name,"/Users/aliheydari/Documents/Scientific_Computing/VTK_Tests/N=%i,t=%i.vtk",N_x,counter);
        my_grid.initialize_VTK_file(name);
        if (t == 0)
            tn_sols = IC;
        else
            tn_sols = sol.Solve();
        
        if (counter % 50 == 0)
            std::cout<<"->" << counter <<"<-" << std::endl;
        t += dt;
        my_grid.print_VTK_Format(tn_sols, "solution_at_nodes",name);
        counter ++;
    }
    
    /* TO GET ERROR
     std::cout<<"-> Calculating difference-<-" << std::endl;
     std::vector<double> sols_diff = Vector_Minus(IC,tn_sols);
     my_grid.display(sols_diff);
     std::cout<<"NORM:"<<std::endl;
     std::cout<<vectorNorm(sols_diff) << std::endl;
     */
    return 0;
}


