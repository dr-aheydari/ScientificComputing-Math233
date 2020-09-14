#ifndef GRID2D_H
#define GRID2D_H


class Grid2D
{
private:
    int N,M;
    double xmin,xmax,ymin,ymax;
    double dx,dy;
public:
    Grid2D();
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    double get_dx();
};

#endif // GRID2D_H
