#include <iostream>
#include <Libs/libs.h>
#include <Libs/grid2d.h>
#include <omp.h>

int main()
{

    velocity_X F_x;
    velocity_Y F_y;
    int N_x = 5;
    int N_y = 10;
    //Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
    Grid2D new_grid(N_x,N_y,0,10,0,5);

    std::cout<<"grid resolution in x :"<<new_grid.get_dx()<<std::endl;

    // vector of values defined at the nodes of the grid
    std::vector<std::vector<double>> funct(N_x, std::vector<double> (N_y, 0.));;
    #pragma omp parallel for collapse(2)
    for (int n=0; n<N_y; n++)
        for (int m=0; m < N_x; m++)
            funct[m][n] = F_x(m,n);

    std::cout<<"value at funct[0][3]: "<<funct[0][3]<<std::endl;
    return 0;
    }
