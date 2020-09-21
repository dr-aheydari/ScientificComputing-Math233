#ifndef FULLMATRIX_H
#define FULLMATRIX_H
#include <vector>

class FullMatrix
{
public:
    //constructors
    FullMatrix();
    FullMatrix(int N);
    // required functions
    void add_Element(int i, int j, double vec);
    double get_Values(int i, int j);
    void display();
    std::vector <double> mat_Vec_Prod(std::vector <double> &x);
    double mat_Mat_Prod(FullMatrix &B);
    // safety and sanity checks
    void SafetyCheck(int n);
    void SafeGaurd(int n);
    void DimCheck(int n);
private:
    int size;
    std::vector<double> values;


};
#endif

