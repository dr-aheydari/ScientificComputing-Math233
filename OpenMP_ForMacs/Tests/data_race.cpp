// THIS AN EXAMPLE OF A CODE THAT IS WRONG!
#include <omp.h>
#include<math.h>
#include <limits>
#include <iostream>

using namespace std;

int main( int argc, char ** argv )
{

    long N=10000000;

    double *value = new double[N];

    double temp;

    omp_set_num_threads(40);
     cout << "OpenMP - Max number of threads : " << omp_get_max_threads() << endl;

// initialize the vector in parallel
   #pragma omp parallel for private(temp) 
   // #pragma omp parallel for public(temp)
    for (long i= 0; i<N; i++)
    {
        temp = exp(cos(i));
        value[i] =sqrt(temp);
	cout<<temp<<endl;
    }  

    // Check in serial
    for (long i= 0; i<N; i++)
    {
        double error = abs(value[i]-sqrt(exp(cos(i))));
        if (error>1E-10)
            cout<<"Wrong value at line "<<i<<"  error : "<<error<<endl;
    }


    /// Solution:
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    ///
    /// PRIVATE VARIABLES !!!!
    // #pragma omp parallel for private(temp)

    return 0;
}

