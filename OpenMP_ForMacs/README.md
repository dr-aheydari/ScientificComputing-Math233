# OpenMP on OSX (High Sirerra and above )

Step 0: Make sure `brew` is installed and ready to go

Step 1: install `llvm` 
`$ brew install llvm`

Step 2: install `libomp`
`$ brew install libomp`

Step 3: Run a test code using the following explicit calls and flags
`$ /usr/local/opt/llvm/bin/clang++ -fopenmp <YOURCODE>.cpp -L/usr/local/opt/llvm/lib -I/usr/local/opt/llvm/include -o <execName> `

example code to run :
```` C++
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
````
if you save as `hello_world.cpp`, you can run :
`$ /usr/local/opt/llvm/bin/clang++ -fopenmp hello_world.cpp -L/usr/local/opt/llvm/lib -I/usr/local/opt/llvm/include -o hello`
***n.b. : make sure your output makes sense*** 
Step 4: Run the executable
`$ ./hello`
***n.b.: make sure your output makes sense*** 

Step 5 (optional): Make an alias for the compiler call (i.e. */usr/local/opt/llvm/bin/clang++*) and the include flags for easier compilation

Step 6 (optional and unnecessary) : Talk to Ali about how to get this working for your specific IDE. 

(Personally, I would recommend using command line for compiling and execution.) 

Let me know if you have any questions or concerns. 


