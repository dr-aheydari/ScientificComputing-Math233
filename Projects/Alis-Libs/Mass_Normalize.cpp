//// this will be the main file
//// creates initial conditions of random ('r'), Gaussian ('n'), or a mixture of the Gaussians ('m')
//
//
//// DO INCLUDES HERE
//
//
//// Structs for weight ditribution
//struct Mass_Normalize {
//    ArrayV<double> zeta, psi;
//
//    double total_mass;
//};
//
////typedef struct greaterSmaller SoftDist;
//
//Mass_Normalize rescale(double initial_totalMass,double current_total ,ArrayV<double> current_A,ArrayV<double> current_B,OcTreeCellBased my_octree)
//{
//
//
//
//    // for now we are not gonna mess with num_funcs since we assume its always 2
//
//    Mass_Normalize corrected;
//
//    corrected.zeta.resize_Without_Copy(my_octree.number_Of_Leaves());
//    corrected.psi.resize_Without_Copy(my_octree.number_Of_Leaves());
//
//
//    double ratio = initial_totalMass/current_total; //2 * current_B + ratio;
//
//    //corrected.total_mass = corrected.psi +  corrected.zeta  ;
//
//    for (int i = 0; i < current_A.size(); i++) //because we got adding array operator, just inbetween step
//    {
//        corrected.psi(i) = ABS(current_A(i) * ratio);
//        corrected.zeta(i) = ABS(current_B(i) * ratio);
//    }
//
//    return corrected;
//}
//
