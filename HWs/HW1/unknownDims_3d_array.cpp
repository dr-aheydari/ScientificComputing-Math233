#include <iostream>
// declaring a 3 dimensional array where the three dimensions are unknown!!
// strategy: we create a pointer to pointers of pointers through a function


int main()
{
    // we do not know the value of these at run time
    int N, M, L;
    std::cout<<"Dynamic N,M,L array-> enter dimensions:"<< std::endl;
    std::cin >> N >> M >> L;
    std::cout<<"N: "<< N << ", M: " << M <<", L: " << L << std::endl;
    // do a pointer to a pointer of pointers
    double ***array3d = new double**[N];
    // asign the values
    for (int i = 0; i < L; i++)
    {
        // make the x dims to be pointers to pointers
        array3d[i] = new double*[M];
        for(int j = 0; j < M; j++)
        {
            // the most inner set of pointers being initialized
            array3d[i][j] = new double[L];
        }

    }

    // asign a value to the last postion of the dynamic array for
    array3d[N-1][M-1][L-1] = 13.21;

    std::cout<<" value of the array at ("<< N << ", " << M<<", "<<L << "): " << array3d[N-1][M-1][L-1] << std::endl;
    
    // delete the memory allocation
    for (int i = 0; i < L; i++)
    {
        for(int j = 0; j < M; j++)
        {
            // the most inner set of pointers being initialized
            delete array3d[i][j];
        }
        delete[] array3d[i];
    }
    delete[] array3d;

   return 0;
}

