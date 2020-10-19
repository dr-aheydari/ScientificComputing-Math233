#ifndef LEVELSET_H
#define LEVELSET_H
#include "XYVelocity.h"
#include "grid2d.h"
#include <vector>

class LevelSet: public CF_2
{
    private:
        int N_x,N_y;
        Grid2D my_grid;
        velocity_X* field_x;
        velocity_Y* field_y;
        std::vector<double> level_set_0;
        std::vector<double> level_set_n;
        double dt;

    public:
        LevelSet();
        LevelSet(const Grid2D& grid, const std::vector<double> Init_cond, double dt_);
        double operator ()(double x, double y) const;

        void Assign_ls();
        void set_level_set_n(std::vector<double> new_level_set);
        std::vector<double> get_level_set_n() const;

        std::vector<double> Reinitialize();
        void reinitialize(int iteration);
        void Perturb_ls(double small_pertubation);
        double Godunov_Check(int sign, double d_p, double d_m) const;

        ~LevelSet();
};



#endif



