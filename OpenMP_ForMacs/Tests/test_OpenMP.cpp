
#include <omp.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include <limits>


using namespace std;



int main( int argc, char ** argv )
{

    long N=10000000; // 1E8

    double *value = new double[N];

    int example=2;

    double ti;
    double tf;
    double duration;
    double duration_serial;
    omp_set_num_threads(10);

    char name[300];

    sprintf(name,"scaling.txt");
    std::ofstream fp_scaling;
    fp_scaling.open(name,ios_base::out);


    // parallel vector filing
    switch (example)
    {
    case 0 :

        ti = omp_get_wtime();
        for (long i= 0; i<N; i++)
            value[i] = exp(sqrt(exp(cos(i))));
        duration = omp_get_wtime() - ti;

        cout<<"Computational time in serial  "<<"  "<<duration<<" s"<<endl;

        break;
    case 1 :

        omp_set_num_threads(2);

        ti = omp_get_wtime();
#pragma omp parallel for
        for (long i= 0; i<N; i++)
            value[i] = exp(sqrt(exp(cos(i))));
        duration = omp_get_wtime() - ti;

        cout<<"Computational time using "<<omp_get_max_threads()<<" threads "<<duration<<" s"<<endl;

        break;
    case 2 :
        for (int nb_threads=1; nb_threads<81; nb_threads++)
        {
            omp_set_num_threads(nb_threads);
            long N=nb_threads*1000000; // 1E8

            double *value = new double[N];

            ti = omp_get_wtime();
#pragma omp parallel for
            for (long i= 0; i<N; i++)
            {
                value[i] = exp(sqrt(exp(cos(i))));
            }
            duration = omp_get_wtime() -  ti;
            if (nb_threads == 1)
                duration_serial = duration;
            cout<<duration_serial/duration<<endl;
            fp_scaling<<duration_serial/duration<<endl;
            //cout<<"Computational time using "<<omp_get_max_threads()<<" threads "<<duration<<" s   speedup : "<<duration_serial/duration<<endl;
        }
        break;


    case 3:  //

        omp_set_dynamic(true);

        ti = omp_get_wtime();

#pragma omp parallel for
        for (long i= 0; i<N; i++)
            value[i] = exp(sqrt(exp(cos(i))));
        duration = omp_get_wtime() - ti;
        cout<<"Computational time on "<<omp_get_max_threads()<<"  "<<duration<<" s"<<endl;
        break;



    case 4:        // don't specify anything

        ti = omp_get_wtime();
        for (long i= 0; i<N; i++)
            value[i] =0.;

        cout<<"Computational time on "<<omp_get_max_threads()<<"  "<<duration<<" s"<<endl;
        break;
    }

    delete[] value;
    return 0;
}
