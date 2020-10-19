#include "VectorUtils.h"



std::vector<double> Vector_Minus(const std::vector<double>& vec1,const std::vector<double>& vec2)
{
    std::vector<double> diff;
    std::transform( std::begin(vec1), std::end(vec1),
                             std::begin(vec2),
                             std::back_inserter(diff),
                             []( double a, double b ) { return std::abs(a-b) ; } ) ;

    return diff;
}

double vectorNorm(const std::vector<double>& vec1)
{
  return sqrt(inner_product(std::begin(vec1), std::end(vec1), std::begin(vec1), 0.0L));
}
