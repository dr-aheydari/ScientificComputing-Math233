#ifndef GRID2D_H
#define GRID2D_H

#include <iostream>
#include <vector>
#include <Libs/XYVelocity.h>
//#include <tuple>


class Grid2D
{
private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    double dx,dy;
public:
    // we need the default constructor for setting initial condition
    Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    //interpolation
    double bilinear_interpolation(const std::vector<double>& func, const double x, const double y) const;
    double quadraticENO_interpolation(const std::vector <double> & func, const double x, const double y, velocity_X field_x, velocity_Y field_y) const;
    // getters
    int get_N() const;
    int get_M() const;
    std::tuple <double, double> get_x_range() const;
    std::tuple <double, double> get_y_range() const;
    double get_dx() const;
    double get_dy() const;
    int get_xmin() const;
    int get_xmax() const;
    int get_ymin() const;
    int get_ymax() const;
    int i_from_n(int n) const;
    int j_from_n(int n) const;
    double x_from_n(int n) const;
    double y_from_n(int n) const;
    int n_from_ij(int i, int j) const;
    // VTK handling
    void initialize_VTK_file(std::string file_name);
    void print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name );
    // Automatic differentiation (different than in ML)
    double Auto_dx(const std::vector<double>& funct,int n, velocity_X field);
    double Auto_dy(const std::vector<double>& funct,int n,velocity_Y field);
    double Auto_dxx(const std::vector<double>& funct,int n, velocity_X field);
    double Auto_dyy(const std::vector<double>& funct,int n,velocity_Y field);    
    // Automatic MinMod
    double Auto_MinMod_x(const std::vector<double>& funct,int n,velocity_X field) const;
    double Auto_MinMod_y(const std::vector<double>& funct,int n,velocity_Y field) const;
    // MinMod
    double MinMod(double lval, double rval) const;
    // first derivatives
    double dx_forward (const std::vector<double>& funct,int n)const;
    double dx_backward (const std::vector<double>& funct,int n)const;
    double dy_forward (const std::vector<double>& funct,int n)const;
    double dy_backward (const std::vector<double>& funct,int n)const;
    // second derivatives
    double dxx_center (const std::vector<double>& funct,int n)const;
    double dyy_center (const std::vector<double>& funct,int n)const;
    double dxx_forward (const std::vector<double>& funct,int n)const;
    double dyy_forward (const std::vector<double>& funct,int n)const;
    double dxx_backward (const std::vector<double>& funct,int n)const;
    double dyy_backward (const std::vector<double>& funct,int n)const;
    // safety checks
    void SafetyCheck(const std::vector<double>& funct) const;
    void CheckDim(const std::vector<double>& funct) const;
    char CheckXBoundary(int n) const;
    char CheckYBoundary(int n) const;
    // display a 1D vector of N*M as a matrix of N*M
    void display(std::vector<double>& funct) const;
    // ~
    ~Grid2D(){};
};

#endif // GRID2D_H
