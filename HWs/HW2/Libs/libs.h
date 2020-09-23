#ifndef LIBS_H
#define LIBS_H
#include <vector>
#include <math.h>

class CF_2
{
public:
    // we need this to be pure virtual since we are going to overload them
    virtual double operator()(double x, double y) const=0;

    // destructor just for safety
    virtual ~CF_2() {}
};

class velocity_X:CF_2
{
   public:
    double operator()(double x, double y) const
    {
        return (x*x) + (y*y);
    }

};

class velocity_Y:CF_2
{
   public:

    double operator()(double x, double y) const
    {
        return (2*x) + (2*y);
    }
};

#endif



