#include "grid2d.h"

Grid2D::Grid2D()
{

}
Grid2D::Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_)
{
    N = N_;
    M = M_;
    xmin = xmin_;
    xmax = xmax_;
    ymin = ymin_;
    ymax = ymax_;
    dx = (xmax - xmin) / (N - 1);
    dy = (ymax - ymin) / (M - 1);
}
double Grid2D::get_dx()
{
    return dx;
}
