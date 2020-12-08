// NaN finder for the values of the total masses
#include <lib/AliTools/isNan.h>

// isNan function def:

int IsNan(double mass_1, double mass_2)
{
    
    switch (isnan(mass_1))
    {
        case true :
            throw std::runtime_error("\n !!!! Psi is Nan!!!! \n ... terminating");
            exit (1);
        case false :
            //all is good
            break;
    }
    
    switch (isnan(mass_2))
    {
        case true :
            throw std::runtime_error("\n !!!!Zeta is Nan!!!! \n ... terminating ");
            exit (1);
        case false:
            //all is good
            break;
    }
    
    
    switch (isinf(mass_1))
    {
        case true :
            throw std::runtime_error("\n !!!! Psi is inf !!!! \n ... terminating ");
            exit (1);
        case false:
            //all is good
            break;
    }
    
    
    switch (isinf(mass_2))
    {
        case true :
            throw std::runtime_error("\n !!!! Zeta is inf!!!! \n ... terminating");
            exit (1);
        case false :
            //all is good
            break;
    }
    
}


