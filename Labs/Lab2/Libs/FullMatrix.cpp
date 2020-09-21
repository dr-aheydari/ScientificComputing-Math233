#include "FullMatrix.h"
#include <iostream>
#include <vector>

FullMatrix::FullMatrix()
{ }

FullMatrix::FullMatrix(int N)
{
    size = N;
    // resize to a squared matrix
    values.resize(N*N);
    // fill with 0.
    std::fill(values.begin(), values.end(), 0.0);
}

void FullMatrix::add_Element(int i, int j, double v)
{
    SafetyCheck(i);
    SafetyCheck(j);
    values[i+j*size] = v;
}

double FullMatrix::get_Values(int i, int j)
{
    SafetyCheck(i);
    SafetyCheck(j);
    return values[i+j*size];
}

// check the diplay function!!!
void FullMatrix::display()
{
    // can not make this parallel since it can mess up the visualization
    for (int i=0; i < size; i++)
    {
        std::cout<<"|";
        for (int j=0; j < size; j++)
        {
            std::cout<<values[i+j*size] << " ";
        }
        std::cout<<"|"<<std::endl;
    }
}

/*
 * Passing by reference is faster than passing by value, specially when it comes to
 * large vectors and matrices
 * So we pass in &x
*/

std::vector<double> FullMatrix::mat_Vec_Prod(std::vector<double> &x)
{

    std::vector<double> res;
    res.resize(size);

    #pragma omp parallel for collapse(2)
    for (int i=0;i<size;i++)
    {
            for (int j=0;j<size;j++)
            {
                res[i]+= values[i+j*size] * x[j];
            }
        }

    return res;
}

double FullMatrix::mat_Mat_Prod(FullMatrix &B)
{
    // add something
}

void FullMatrix::SafetyCheck(int n)
{
    try
    {
        SafeGaurd(n);
    } catch (const char* msg)
    {
        std::cerr << msg << std::endl;
        exit(-1);
     }
}

void FullMatrix::SafeGaurd(int n)
{
    if ( n < 0 || n > size-1 )
    {
        throw "Must enter an integer between 0 and dimension - 1 (inclusive) ";
    }
}

void FullMatrix::DimCheck(int n)
{
    //
}




