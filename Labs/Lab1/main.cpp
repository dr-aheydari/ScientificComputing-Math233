#include <iostream>
#include <grid2d.h>

using namespace std;

int main()
{

    cout << "Hello World!" << endl;

    int N = 2;
    int M = 3;
    double xmin = 0;
    double xmax = 2.2;
    double ymin = 0;
    double ymax = 2;


    // the command below creates an instance of the grid2d class
    Grid2D grid;

    // lets use the fancy constructor
    Grid2D new_grid(N,M,xmin,xmax,ymin,ymax);

    cout<<"this ios the grid resolution in x :"<<new_grid.get_dx()<<endl;


    return 0;
}
