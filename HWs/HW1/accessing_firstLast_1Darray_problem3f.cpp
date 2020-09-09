#include <iostream>

int main()
{
   // a random vector
   int A[] = { 4, 5, 7, 13, 25, 65, 98, 15, 99, 199, 200};
   // get 1D array's size
   int N = sizeof(A) / sizeof(A[0]);
   // cout info
   std::cout<<"Size of the array: " << N<<std::endl;
   std::cout<<"first element: " << A[0]<<std::endl;
   std::cout<<"last element: " << A[N-1]<<std::endl;


   return 0;
}


