#include <iostream>
#include <Libs/libs.h>
#include <Libs/grid2d.h>
#include <omp.h>

int main()
{
    
    velocity_X V_x;
    velocity_Y V_y;
    // a class for handling the function
    xyFunction F_xy;
    int N_x = 100;
    int N_y = 100;
    //Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
    Grid2D my_grid(N_x,N_y,0,1,0,1);
    
    std::cout<<"grid resolution in x :"<<my_grid.get_dx()<<std::endl;
    std::vector<double> funct;
    F_xy.Assign(funct, N_x, N_y);
    //    my_grid.display(funct);
    std::cout<<"Automatic x differentiation \n"<<my_grid.Auto_dx(funct, 0, V_x)<<std::endl;
    std::cout<<"Automatic y differentiation \n"<<my_grid.Auto_dy(funct, 0, V_y)<<std::endl;
    
    
    // create a name for the vtk file, with the grid size in it
    char name[250];
    sprintf(name,"/Users/aliheydari/Documents/Scientific_Computing/VTK_Tests/grid2D_N=%d_M=%d.vtk",N_x,N_y);
    
    //export to vtk file
    //create the vtk with all the grid2d information
    my_grid.initialize_VTK_file(name);
    
    // output the values at the nodes in the vtk file.
    //The quantity will show up in paraview as 'value_at_nodes'
    my_grid.print_VTK_Format(funct, "value_at_nodes",name);
    
    return 0;
}
