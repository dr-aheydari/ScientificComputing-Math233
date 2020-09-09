#include <cmath>
#include <iostream>
// proto
double Legendre(double x, int n);
double* SampledLegendre(double a, double b, int N, int n);
double L1(double x);
double L2(double x);
double L3(double x);
double L4(double x);
double L5(double x);
double L6(double x);
void SafetyCheck(int n);
void SafeGaurd(int n);
void DisplayArray(double* p, int N);

int main()
{   // set precision for cout
    std::cout.precision(32);
    // initiailize
    double a,b = 0;
    int n,N = 0;
    // we assume to need dynamic arrays (to make it more interesting)
    std::cout<<"Legendre Polynomial Evaluator!" << std::endl;
    std::cout<<"Enter start point of interval (a): " << std::endl;
    std::cin >> a;
    std::cout<<"Enter end point of interval (b): " << std::endl;
    std::cin >> b;
    std::cout<<"Enter number of partitions (N): " << std::endl;
    std::cin >> N;
    std::cout<<"Enter n for Legendre Polynomials (0 <= n <= 6): " << std::endl;
    std::cin >> n;
    // check for correct degree
    SafetyCheck(n);
    double *point = SampledLegendre(a, b, N, n);
    // show results
    DisplayArray(point, N);
    // delete the pointer and point to null
    delete[] point;
    point = nullptr;
    return 0;
}
// defs
double Legendre(double x, int n)
{
    switch(n)
    {
    case 0:
        return 1.0;
    case 1:
        return L1(x);
    case 2:
        return L2(x);
    case 3:
        return L3(x);
    case 4:
        return L4(x);
    case 5:
        return L5(x);
    case 6:
        return L6(x);
    }
    // if there is something wrong (which should not since we have a safty function
    return -2;
}
double L1(double x) {return x;}
double L2(double x) {return 0.5*(3*std::pow(x,2) - 1);}
double L3(double x) {return 0.5*(5*std::pow(x,3) - 3*x);}
double L4(double x) {return 0.125*(35*std::pow(x,4)-30*x*x+3); }
double L5(double x) {return 0.125*(63*std::pow(x,5)-70*std::pow(x,3)+ 15*x); }
double L6(double x) {return 0.0625*(231*std::pow(x,6)-315*std::pow(x,4)
                                    + 105*std::pow(x,2) - 5); }

double* SampledLegendre(double a, double b, int N, int n)
{
    SafetyCheck(n);
    double*L = new double[N];
    double diff = b-a;
    double delta_t = diff/(N-1);
    double curr_point = a;
    for (int i=0; i < N; i++)
    {
        L[i] = Legendre(curr_point,n);
        curr_point = curr_point + delta_t;

    }
    std::cout<<"Done\n";
    return L;
}

void SafetyCheck(int n)
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

void SafeGaurd(int n)
{
    if ( n < 0 || n > 6 )
    {
        throw "Must enter an integer between 0 and 6 (inclusive)";
    }
}

void DisplayArray(double* p, int N)
{   std::cout<<"Diplaying Array Values: " << std::endl;
    for (int i = 0; i < N; i++)
    {
        std::cout<<"Array["<<i<<"] ="<< *(p + i)<< std::endl;
    }
}

