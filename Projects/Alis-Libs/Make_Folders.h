// this will be the header file
// creates folders automatically + textfiles for mass conservations

#ifndef MAKE_FOLDER
#define MAKE_FOLDER

#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <stdio.h>
#include <sys/uio.h>

#include </Users/aliheydari/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include<unistd.h>

using namespace std;

void MakeFolder(char compart, char hole, char DiffRXN,char DiffOnly,double D_psi, double D_zeta, double gamma_AtoB, double gamma_BtoA, int max_level, int min_level, double tOrder,
                string* txtPathA,string* txtPathB,string* daughterPathA,string* daughterPathB,string* motherPathB,string* motherPathA,char* FolderPath,
                char initCond, string FullPath,string* Return_FolderPAth);


#endif //MAKE_FOLDER
