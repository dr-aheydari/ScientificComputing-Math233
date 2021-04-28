# Installing PETSc Library on macOS

Step 1: Download 
Download the tarball and extract

Step 2: Configure

Run: 
````
./configure prefix = <DIRECTORY/FOR/INSTALLATION> --with-clanguage=cxx --download-hypre=1 --download-f2cblaslapack=1 --download-mpich=1
````
Note: the installation directory (i.e. `prefix` arg) can not be the same as the downloaded directory!

Additional arguments:

For release:

````
--with-debugging=0
COPTFLAGS="-O2 -match=native"
CXXOPTFLAGS="-O2 -match=native"
FOPTFLATS="-O2 -match=native"
```` 

For debug:

````
--with-debugging=1
COPTFLAGS="-g"
CXXOPTFLAGS="-g"
FOPTFLATS="-g"
````

and to enable config with OpenMP:
````bash
--with-openmp
````

After a successful config, run `make install` as usual. Then, get the linked libs with the command
````bash 
make PETSC_DIR=path_to_petsc_install_dir PETSC_ARCH=arch_name getlinklibs
````

## Set up in QT Creator `.pro`:

We now need to add the paths to the `.pro` file (or makefile for a different IDE)

````
PETSC_DIR = /path/to/petsc
PETSC_INCLUDES = $${PETSC_DIR}/include $${PETSC_DIR}/arch_name/include
PETSC_LIBS = "getlinklibs results"
INCLUDEPATH += $${PETSC_INCLUDES}
LIBS += $${PETSC_LIBS}
````
and this should be it. Now depending on usage, the usual flags can be added in the makefile. 

Here is an actual example of these options in the `.pro` file:

````
PETSC_DIR =  /Users/aliheydari/Documents/Softwares/petsc-3.14.0
PETSC_ARCH = arch-darwin-cxx-debug
PETSC_INCLUDES = $${PETSC_DIR}/include $${PETSC_DIR}/$${PETSC_ARCH}/include
PETSC_LIBS = -L$${PETSC_DIR}/$${PETSC_ARCH}/lib -Wl,-rpath,$${PETSC_DIR}/$${PETSC_ARCH}/lib -lpetsc -lmpich -L /usr/local/lib --with-openmp # -L/usr/X11R6/lib  -lHYPRE -lpthread -lflapack -lfblas  -ldl -lmpichf90 -lpthread -lquadmath  -lm  -lstdc++  -lopa -lmpl -lpetsc
# or another set of arguments
#LIBS += -L/$${PETSC_DIR}/$${PETSC_ARCH}/lib -lpetsc -L/usr/X11R6/lib -lHYPRE -lpthread -lflapack -lfblas  -ldl -lmpichf90 -lpthread -lquadmath -l -lm -lstdc++  -lmpich -lopa -lmpl  -ldl

````


