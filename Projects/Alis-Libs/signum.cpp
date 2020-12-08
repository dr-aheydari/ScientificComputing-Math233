// main file for signum.h

#include <lib/AliTools/signum.h>

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
