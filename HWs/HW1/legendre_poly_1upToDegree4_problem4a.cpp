#include <cmath>
#include <iostream>
// proto
double Legendre(double x, int n);
double L1(double x);
double L2(double x);
double L3(double x);
double L4(double x);
void SafetyCheck(int n);
void SafeGaurd(int n);

int main()
{   // set precision for cout
    std::cout.precision(32);
    // initiailize
    double x = 0;
    int n = 0;
    std::cout<<"Legendre Polynomial Evaluator!" << std::endl;
    std::cout<<"Enter a double x: " << std::endl;
    std::cin >> x;
    std::cout<<"Enter n for Legendre Polynomials (0 < n < 5): " << std::endl;
    std::cin>>n;
    SafetyCheck(n);
    std::cout<<"Legendre polynomial of degree " << n << " at "
            << x << " is = " << Legendre(x,n) << std::endl;
    return 0;
}

// defs
double Legendre(double x, int n)
{
    switch(n)
    {
    case 1:
        return L1(x);
    case 2:
        return L2(x);
    case 3:
        return L3(x);
    case 4:
        return L4(x);
    }
    // if there is something wrong (which should not since we have a safety function
    return -2;
}
double L1(double x) {return x;}
double L2(double x) {return 0.5*(3*std::pow(x,2) - 1);}
double L3(double x) {return 0.5*(5*std::pow(x,3) - 3*x);}
double L4(double x) {return 0.125*(35*std::pow(x,4)-30*x*x+3); }
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
    if ( n < 1 || n > 4 )
    {
        throw "Must enter an integer between 0 and 5 (non-inclusive)";
    }
}


