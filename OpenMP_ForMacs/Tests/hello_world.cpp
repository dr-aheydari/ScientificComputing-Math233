#include <omp.h>
#include <iostream>

using namespace std;

int main( int argc, char ** argv )
{
    omp_set_num_threads(10);
    cout << "OpenMP - Max number of threads : " << omp_get_max_threads() << endl;
#pragma omp parallel
//   cout <<"Hello World"<<endl;
    cout<<"Hello World  my name is thread  "<<omp_get_thread_num()<< endl;

    return 0;
}
