#ifndef VECTORUTILS_H
#define VECTORUTILS_H
#include <numeric>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <cmath>

// to take the difference between two std::vector<double>
std::vector<double> Vector_Minus(const std::vector<double>& vec1,const std::vector<double>& vec2);
//template<typename Iter_T> double vectorNorm(const std::vector<double>& vec1);
double vectorNorm(const std::vector<double>& vec1);

#endif



