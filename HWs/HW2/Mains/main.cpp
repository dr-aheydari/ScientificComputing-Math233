#include <iostream>
#include <Libs/libs.h>
#include <Libs/grid2d.h>
#include <omp.h>

int main()
{
    
    velocity_X F_x;
    velocity_Y F_y;
    int N_x = 10;
    int N_y = 10;
    //Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
    Grid2D my_grid(N_x,N_y,0,1,0,1);
    
    std::cout<<"grid resolution in x :"<<my_grid.get_dx()<<std::endl;
    
    // vector of values defined at the nodes of the grid
    std::vector<double> funct;
    funct.resize(N_x*N_y);
    // fill with 0.
    std::fill(funct.begin(), funct.end(), 0.0);
#pragma omp parallel for collapse(2)
    for (int i=0; i<N_x; i++)
        for (int j=0; j < N_y; j++)
        {
            double x  = double(i)/N_x;
            double y  = double(j)/N_y;
            funct[i+ j*N_y]= F_y(x,y);
            //            funct[i+ j*N_y]= F_x(x,y);
        }
    
    my_grid.display(funct);
    //    std::cout<<"Automatic differentiation \n"<<my_grid.Auto_dx(funct, 0, F_x)<<std::endl;
    std::cout<<"Automatic differentiation \n"<<my_grid.Auto_dy(funct, 14, F_y)<<std::endl;
    
    
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
