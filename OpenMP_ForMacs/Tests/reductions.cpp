
#include <omp.h>
#include<math.h>
#include <limits>
#include <iostream>



using namespace std;




int main( int argc, char ** argv )
{

    omp_set_num_threads(4);


    
    long N=100000000;
    double *value = new double[N];

    // set the number of thread and do not allow it to change ....
    omp_set_num_threads(12);
    omp_set_dynamic(false);


    // initialize vector
#pragma omp parallel for
    for (long i= 0; i<N; i++)
        value[i] = 1;


    int example =1;

    // lets compute the average of the value vector

    double correct_avg = 0;
    double avg = 0;


    for (long i= 0; i<N; i++)
        correct_avg+=value[i];

    correct_avg/=N;
    cout<<"the correct average is "<<correct_avg<<endl;

// Naive parallel implementation
#pragma omp parallel for
    for (long i= 0; i<N; i++)
        avg+=value[i];
    avg/=N;
    cout<<"average computed in parallel (naive approach): "<<avg<<endl;


    // doing it by hand

    avg = 0;

//    // create a avegs vector that will contain the avg (sum....) for each thread
    int nb_of_threads = omp_get_max_threads();
    double *avgs = new double[nb_of_threads];
    for (int i=0; i<nb_of_threads; i++)
        avgs[i] = 0;



#pragma omp parallel for
    for (long i= 0; i<N; i++)
        avgs[omp_get_thread_num()]+=value[i];

    for(int i=0; i<nb_of_threads; i++ )
        avg+=avgs[i];
    avg/=N;
    cout<<"average computed in parallel (by hand): "<<avg<<endl;


    // using reduction

    avg = 0;
    #pragma omp parallel for reduction(+:avg)
            for (long i= 0; i<N; i++)
                avg+=value[i];
            avg/=N;
            cout<<"average computed in parallel (reduction): "<<avg<<endl;





    return 0;
}
