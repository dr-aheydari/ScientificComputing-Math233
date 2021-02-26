# Installing PETSc Library on macOS

Step 1: Download 
Download the tarball and extract

Step 2: Configure

Run `./Configure --with-clanguage=cxx --download-hypre=1 --download-f2cblaslapack=1 --download-mpich=1` with possibly these additional arguments:

for release:

````bash
--with-debugging=0
COPTFLAGS="-O2 -match=native"
CXXOPTFLAGS="-O2 -match=native"
FOPTFLATS="-O2 -match=native"
```` 

for debug:

````
--with-debugging=1
COPTFLAGS="-g"
CXXOPTFLAGS="-g"
FOPTFLATS="-g"
````

and to have it with OpenMP:
````
--with-openmp
```
