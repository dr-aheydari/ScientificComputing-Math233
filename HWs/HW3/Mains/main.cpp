#include <omp.h>
#include <vector>
#include <iostream>
#include <Libs/grid2d.h>
#include <Libs/LevelSet.h>
#include <Libs/XYFunction.h>
#include <Libs/VectorUtils.h>
#include <Libs/SemiLagrangian.h>


int main()
{
    // set number of threads
#pragma omp parallel
    omp_set_num_threads(8);
    int N_x = 128;
    int N_y = 128;

    // set x and y velocity fields
    velocity_X x_vfield;
    velocity_Y y_vfield;

    // create a square grid [-1,1]^2
    Grid2D my_grid(N_x,N_y,-7,7,-7,7);
    std::cout<<"grid resolution in x :"<<my_grid.get_dx()<<std::endl;
    std::cout<<"grid resolution in y :"<<my_grid.get_dy()<<std::endl;

    // get initial condition
    xyFunction init_cond(my_grid);
    std::vector<double> IC = init_cond.Assign();

    double t = 0.;
    double tmax = 8*acos(0.0);
    /* for convergence analysis in space we need this dt
    //double dt = std::pow(my_grid.get_dx(),2);
    */
    double dx = my_grid.get_dx();
    double dt = 0.5 * dx;
//    double dt = (1./10) * my_grid.get_dx();
    std::cout<<"dt :"<< dt << std::endl;
    int counter = 0;



////////////////////////////////////////////////////Semi-Lagrangian////////////////////////////////////////////////////////
//    SemiLagrangian advection(my_grid,IC, x_vfield, y_vfield, dt);
//    // for testing
//    std::vector<double> fn_sols;
//    fn_sols.resize(N_x*N_y);

//
//    while (t < tmax)
//    {
//        char name[250];
//        sprintf(name,"/Users/aliheydari/Box/Course Material/Fall 2020/Scientific Computing/VTK_Sims/HW3/Tests/N=%i,t=%i.vtk",N_x,counter);
//        my_grid.initialize_VTK_file(name);
//        if (t == 0)
//            fn_sols = IC;
//        else
//            fn_sols = advection.Solve();

//        if (counter % 10 == 0)
//            std::cout<<"->" << counter <<"<-" << std::endl;
//        t += dt;
//        std::cout<<"Current Time: " << t << std::endl;
//        my_grid.print_VTK_Format(fn_sols, "solution_at_nodes",name);
//        counter ++;
//    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Semi-Lagrangian



////////////////////////////////////////////////////Error-Analysis////////////////////////////////////////////////////////


//std::cout<<"-> Calculating difference-<-" << std::endl;
//std::vector<double> sols_diff = Vector_Minus(IC,fn_sols);
////my_grid.display(sols_diff);
//std::cout<<"NORM:"<<std::endl;
//std::cout<<vectorNorm(sols_diff) << std::endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////Error-Analysis






//////////////////////////////////////////////////Level-Set Reinitialize////////////////////////////////////////////////////////
//    char name[250];
//    sprintf(name,"/Users/aliheydari/Box/Course Material/Fall 2020/Scientific Computing/VTK_Sims/HW3/LevelSetTests/N=%i_actual_reinitialized.vtk",N_x);
//    my_grid.initialize_VTK_file(name);
//    my_grid.print_VTK_Format(IC, "LevelSet_at_nodes",name);

//    LevelSet ls_nodes(my_grid, dt);
//    ls_nodes.Assign_ls();

////    ls_nodes.Perturb_ls(10);
//    ls_nodes.Perturb_ls(10);
//    std::vector<double> ls_func = ls_nodes.get_level_set_n();

//    sprintf(name,"/Users/aliheydari/Box/Course Material/Fall 2020/Scientific Computing/VTK_Sims/HW3/LevelSetTests/N=%i_perturbed.vtk",N_x);
//    my_grid.initialize_VTK_file(name);
//    my_grid.print_VTK_Format(ls_func, "LevelSet_at_nodes",name);


//    ls_nodes.reinitialize(1);
//    ls_func = ls_nodes.get_level_set_n();
//    my_grid.display(ls_func);




////     OcTreeLevelSet ls_nodes(my_octree,level_set_n);
////     ls_nodes.reinitialize();
////     ls_nodes.perturb_Level_Function(1E-2*my_octree.dx_finest_resolution());



//////////////////////////////////////////////////Level-Set////////////////////////////////////////////////////////
       SemiLagrangian advection(my_grid,IC, x_vfield, y_vfield, dt);
       t= 0;
       counter = 0;
       LevelSet ls_nodes(my_grid,IC , (0.1*dx)/10);
//       ls_nodes.Assign_ls();
       std::vector<double> ls_func(N_x,N_y);

       while (t < tmax)
       {
           char name[250];
           sprintf(name,"/Users/aliheydari/Box/Course Material/Fall 2020/Scientific Computing/VTK_Sims/HW3/LevelSetTests/N=%i,t=%i.vtk",N_x,counter);
           my_grid.initialize_VTK_file(name);
           if (t == 0)
              ls_func = IC;
           else
           {
               // advect the solution
               ls_func = advection.LS_Solve(ls_func);
               // set this advection as the level_set function
               ls_nodes.set_level_set_n(ls_func);
               // reinitialize the solution
               ls_nodes.reinitialize(20);
               // get the reinitialized function as the solution
               ls_func = ls_nodes.get_level_set_n();
           }
//           my_grid.display(ls_func);
           if (counter % 10 == 0)
              std::cout<<"->" << counter <<"<-" << std::endl;
           my_grid.print_VTK_Format(ls_func, "LevelSet_at_nodes",name);
           t += dt;
           counter ++;
       }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Level-Set


 return 0;


}



/*
////////////////// Testing interpolation error
    std::vector<double> test_function;
    test_function.resize(N_x*N_y);
    xyFunction something(my_grid);

#pragma omp parallel for
    for (int i=0; i < N_x*N_y; i++)
    {
        double x = my_grid.x_from_n(i);
        double y = my_grid.y_from_n(i);
        test_function[i] = something(x,y);
    }

//    my_grid.display(test_function);

    double x_inter = 0.823;
    double y_inter = 0.543;

    double interp_value = my_grid.bilinear_interpolation(test_function,x_inter, y_inter);

    std::cout<<"interp :"<<interp_value<<std::endl;

    std::cout<<"My interpolation error: "<<abs(something(x_inter,y_inter) - interp_value) << std::endl;

//////////////////
*/
