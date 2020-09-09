#include <iostream>

// declaring a 3 dimensional array where the three dimensions are known in advance

int main()
{
// when the dimensions are known in advanced, say 3,4,5, for an array of real numbers

double array[3][4][5];
// output the sizes to make sure all is good
int iter[] = {0,1,2};

for (int i=0; i < 3; i++)
    {
        std::cout<<"Dimension "<<iter[i] << ": ";
        if (i == 0)
            std::cout<<(sizeof array / sizeof array[0]);
        else if (i == 1)
            std::cout<<(sizeof array[0] / sizeof array[0][0]);
        else
            std::cout<<(sizeof array[0][0] / sizeof array[0][0][0]);
        std::cout << std::endl;
    }


   return 0;
}

