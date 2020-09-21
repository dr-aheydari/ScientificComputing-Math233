#include <iostream>
#include <Libs/FullMatrix.h>
#include <omp.h>

int main()
{
    int dim = 2;
    FullMatrix instance(dim);
    instance.add_Element(0,0,2);
    instance.add_Element(0,1,1);
    instance.add_Element(1,0,0);
    instance.add_Element(1,1,2);
    instance.display();

// checking for a simple multiplication case of x={1,2}
    std::vector<double> x = {1,2};
    std::vector<double> res = instance.mat_Vec_Prod(x);
    std::cout<<"[" << res[0] << ", " << res[1] << "]" << std::endl;

    return 0;
}


