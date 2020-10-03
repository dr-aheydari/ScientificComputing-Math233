#include "grid2d.h"
#include <vector>
#include <cmath>
#include <tuple>
#include <cfloat>

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}


////////////////automated differentiation functions
double Grid2D::Auto_dx(const std::vector<double> &funct, int n, velocity_X field)
{
    
    switch (CheckXBoundary(n))
    {
        case 'f':
            return 0.0;//dx_forward(funct,n);
        case 'b':
            return 0.0; //dx_backward(funct,n);
        case 'c':
            double x = x_from_n(n);
            double y = y_from_n(n);
            double condition = field(x,y);
            if (sgn(condition) == 1)
            {
                return dx_backward(funct,n);
            }
            else
            {
                return dx_forward(funct,n);
            }
    }
    std::cout<<"HELLO_x?"<<std::endl;
}

double Grid2D::Auto_dy(const std::vector<double> &funct, int n, velocity_Y field)
{
    
    switch (CheckYBoundary(n))
    {
        case 'f':
            return 0.0; //dy_forward(funct,n);
        case 'b':
            return 0.0; //dy_backward(funct,n);
            // check these conditions
        case 'c':
            double x = x_from_n(n);
            double y = y_from_n(n);
            double condition = field(x,y);
            if (sgn(condition) == 1)
            {
                return dy_backward(funct,n);
            }
            else
            {
                return dy_forward(funct,n);
            }
    }
    std::cout<<"HELLO_y?"<<std::endl;
}


char Grid2D::CheckXBoundary(int n)
{
    // check for x boundary
    if (x_from_n(n) == xmax)
    {
        return 'b';
    }
    else if (x_from_n(n) == xmin)
    {
        return 'f';
    }
    // all is good then!
    else
    {
        return 'c';
    }
    
}

char Grid2D::CheckYBoundary(int n)
{
    // check for y boundary
    if (y_from_n(n) == ymin)
    {
        return 'f';
    }
    else if(y_from_n(n) == ymax)
    {
        return 'b';
    }
    else
    {
        return 'c';
    }
}




double Grid2D::Auto_dxx(const std::vector<double> &funct, int n, velocity_X field)
{
    
    switch (CheckXBoundary(n))
    {
        case 'f':
            return dxx_forward(funct,n);
        case 'b':
            return dxx_backward(funct,n);
        case 'c':
            return MinMod(dxx_center(funct,n),dxx_center(funct,n-1));
    }
}


double Grid2D::Auto_dyy(const std::vector<double> &funct, int n, velocity_Y field)
{
    
    switch (CheckYBoundary(n))
    {
        case 'f':
            return dyy_forward(funct,n);
        case 'b':
            return dyy_backward(funct,n);
        case 'c':
            return MinMod(dyy_center(funct,n),dyy_center(funct,n-1));
    }
}

//##############automated differentiation functions

////////////////all the derivative function
double Grid2D::dx_forward (const std::vector<double>& funct,int n)const
{
    
    SafetyCheck(funct);
    return (funct[n + 1] - funct[n])/(dx + DBL_EPSILON);
}

double Grid2D::dx_backward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return (funct[n] - funct[n - 1])/(dx + DBL_EPSILON);
}

double Grid2D::dy_forward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return (funct[n + N] - funct[n])/(dy + DBL_EPSILON);
    
}
double Grid2D::dy_backward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return (funct[n] - funct[n - N])/(dy + DBL_EPSILON);
}

//##############first derivatives


/////////////// 2nd order derivatives
double Grid2D::dxx_center (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return (funct[n + 1] - 2* funct[n] + funct[n-1])/(std::pow(dx,2)+ DBL_EPSILON);
}

double Grid2D::dxx_backward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return 0.;
    //    return (funct[n] - 2*funct[n - 1]+funct[n - 2])/(std::pow(dx,2)+ DBL_EPSILON);
}

double Grid2D::dxx_forward(const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return 0.;
    //    return (funct[n+2] - 2*funct[n+1]+funct[n])/(std::pow(dx,2)+ DBL_EPSILON);
}


double Grid2D::dyy_center (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    // get the (i,j) coordinate
    return (funct[n + N] - 2* funct[n] + funct[n - N])/(std::pow(dy,2)+ DBL_EPSILON);
}

double Grid2D::dyy_forward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return 0.;
    //return (funct[n + 2*N] - 2* funct[n + N] + funct[n])/(std::pow(dy,2)+ DBL_EPSILON);
}

double Grid2D::dyy_backward (const std::vector<double>& funct,int n)const
{
    SafetyCheck(funct);
    return 0.;
    //return (funct[n] - 2* funct[n - N] + funct[n - 2*N])/(std::pow(dy,2)+ DBL_EPSILON);
}

//############# 2nd order derivatives




double Grid2D::get_dy() const
{
    return dy;
}

double Grid2D::get_N() const
{
    return N;
}

double  Grid2D::get_M() const
{
    return N;
}
std::tuple <double, double>  Grid2D::get_x_range() const
{
    return std::make_tuple(xmin,xmax);
}
std::tuple <double, double>  Grid2D::get_y_range() const
{
    return std::make_tuple(ymin,ymax);
}

void Grid2D::display(std::vector<double>& funct) const
{
    // can not make this parallel since it can mess up the visualization
    for (int i=0; i < N; i++)
    {
        std::cout<<"|";
        for (int j=0; j < M; j++)
        {
            std::cout<<funct[i+j*M] << " ";
        }
        std::cout<<"|"<<std::endl;
    }
}


double Grid2D::MinMod(double lval, double rval)
{
    if(lval*rval<=0)
        return 0;
    else
    {
        if((abs(lval))<(abs(rval)))
            return lval;
        else
            return rval;
    }
}


void Grid2D::SafetyCheck(const std::vector<double>& funct) const
{
    try
    {
        CheckDim(funct);
    } catch (const char* msg)
    
    {
        std::cout<<std::endl;
        std::cerr << "Length error: " << msg << std::endl;
        exit(-1);
    }
}

void Grid2D::CheckDim(const std::vector<double>& funct) const
{
    if (int(funct.size()) > M*N)
    {
        throw "Function has a larger size than the grid. Terminating...";
    }
}


/////////////////////////////////////
Grid2D::Grid2D()
{
}
Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{
    N = N_; // for x
    M = M_; // for y
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;
    dx = (xmax - xmin) / (double) (N - 1);
    dy = (ymax - ymin) / (double) (M - 1);
}
double Grid2D::get_dx() const
{
    return dx;
}

int Grid2D::i_from_n(int n) const
{
    // n = i+ j * N
    return  n % N;
}
int Grid2D::j_from_n(int n) const
{
    // n = i+ j * N
    return  n / N;
}

int Grid2D::n_from_ij(int i, int j)
{
    return i+ j * N;
}

double Grid2D::x_from_n(int n)
{
    return  xmin + dx*i_from_n(n);
}
double Grid2D::y_from_n(int n)
{
    return  ymin + dy*j_from_n(n);
}

// initialize the .vtk file at the specified address with all the grid information

void Grid2D::initialize_VTK_file(std::string file_name)
{
    int node_of_cell[4];
    
    FILE *outFile = fopen(file_name.c_str(),"w");
    
    fprintf(outFile,"# vtk DataFile Version 2.0 \n");
    fprintf(outFile,"Quadtree Mesh \n");
    fprintf(outFile,"ASCII \n");
    fprintf(outFile,"DATASET UNSTRUCTURED_GRID \n");
    
    
    //% first output the list of nodes
    fprintf(outFile,"POINTS %d double \n",N*M);
    for (int n=0; n<N*M; n++)
        fprintf(outFile,"%e %e %e\n",x_from_n(n), y_from_n(n), 0.0);
    
    
    // then output the list of cells. each cell is composed of four nodes
    fprintf(outFile,"CELLS %d %d \n",(N-1)*(M-1),5*(N-1)*(M-1));
    for (int i=0; i<N-1; i++)
        for (int j=0; j<M-1; j++)
        {
            node_of_cell[0] = n_from_ij(i  ,j  );
            node_of_cell[1] = n_from_ij(i+1,j  );
            node_of_cell[2] = n_from_ij(i+1,j+1);
            node_of_cell[3] = n_from_ij(i  ,j+1);
            
            fprintf(outFile,"%d %d %d %d %d\n",4,node_of_cell[0], node_of_cell[1], node_of_cell[2], node_of_cell[3]);
        }
    //  }
    fprintf(outFile,"CELL_TYPES %d \n",(N-1)*(M-1));
    for (int n=0; n<(N-1)*(M-1); n++)    fprintf(outFile,"%d \n",9);
    fprintf(outFile,"POINT_DATA %d \n",N*M);
    fclose (outFile);
}
// this function write the values of the vector F into the vtk file. before using it, the .vtk file must have been initialized with all the grid infos
void Grid2D::print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name )
{
    
    FILE *outFile;
    outFile = fopen(file_name.c_str(),"a");
    fprintf(outFile,"SCALARS %s double 1 \n",data_name.c_str());
    fprintf(outFile,"LOOKUP_TABLE default \n");
    for (int n=0; n<N*M; n++) fprintf(outFile,"%e \n",F[n]);
    fclose (outFile);
}
