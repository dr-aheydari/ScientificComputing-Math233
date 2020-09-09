# HW 1

This homework is a very simple "reminder" of `C++`, and thus the code is not reflective of the coding level of the class.

## Notes on Memory Leak checks:

All homework codes which used pointers and/or references were checked with `valgrind` for leak detections. All of the tested codes showed no memory leaks or misallocation. To reproduce:
````
$ g++ -g <desired_code>.cpp
 
$ valgrind --tool=memcheck --track-origins=yes --leak-check=full ./a.out
```` 
