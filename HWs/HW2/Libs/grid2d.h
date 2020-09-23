#ifndef GRID2D_H
#define GRID2D_H

#include <iostream>
#include <vector>

class Grid2D
{
private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    double dx,dy;
public:
    Grid2D();
//    ~Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double get_dx() const;
    double get_dy() const;
    int i_from_n(int n) const;
    int j_from_n(int n) const;
    double x_from_n(int n);
    double y_from_n(int n);
    int n_from_ij(int i, int j);
    void initialize_VTK_file(std::string file_name);
    void print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name );
    // New functions
    double dx_forward (const std::vector<std::vector<double>>& funct,int n)const;
    double dx_backward (const std::vector<std::vector<double>>& funct,int n)const;
    double dy_forward (const std::vector<std::vector<double>>& funct,int n)const;
    double dy_backward (const std::vector<std::vector<double>>& funct,int n)const;

    /*
     * with 1D function (I do not like this way that Maxime has asked us to do it
    double dx_forward (const std::vector<double>& funct,int n)const;

    */
};

#endif // GRID2D_H
