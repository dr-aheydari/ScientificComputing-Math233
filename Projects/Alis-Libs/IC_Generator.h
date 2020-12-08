// this is the header file
// creates initial conditions of random ('r'), Gaussian ('n'), or a mixture of the Gaussians ('m'), or exact solution for testing ('e')

#ifndef IC_GENERATOR
#define IC_GENERATOR

#include <iostream>
#include <lib/amr/OcTree.h>
#include <lib/arrays/ArrayV.h>
#include <lib/amr/OcTreeCellBased.h>
#include <lib/AliTools/Psi_Zeta_Init_Class.h>

using namespace std;
using namespace CASL;

struct IC_Generator {
    ArrayV<double> zeta, psi;
};


IC_Generator gen(char choice ,ArrayV<double> current_A,ArrayV<double> current_B,OcTreeCellBased my_octree, Initial_Solution psi_initial, zeta_Initial_Solution zeta_initial);
IC_Generator gen_exact(char choice ,ArrayV<double> current_A,OcTreeCellBased my_octree, Initial_Solution psi_initial);

#endif // IC_GENERATOR
