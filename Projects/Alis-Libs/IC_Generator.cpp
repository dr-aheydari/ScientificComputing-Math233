# include <lib/AliTools/IC_Generator.h>

IC_Generator gen(char choice ,ArrayV<double> current_A,ArrayV<double> current_B,OcTreeCellBased my_octree, Initial_Solution psi_initial, zeta_Initial_Solution zeta_initial)
{

    double epsilon = EPSILON;

    // for now we are not gonna mess with num_funcs since we assume its always 2

    IC_Generator vec;

    vec.zeta.resize_Without_Copy(my_octree.number_Of_Leaves());
    vec.psi.resize_Without_Copy(my_octree.number_Of_Leaves());



    if (current_A.size() == current_B.size())
    {

        switch(choice)
        {
            case 'r':
                for (int i = 0; i < current_A.size(); i++) //because we got adding array operator, just inbetween step
                {
                    vec.psi(i) = rand() / (RAND_MAX + 1. + epsilon);

                    vec.zeta(i) = rand() / (RAND_MAX + 1. + epsilon);
                }
                break;

            case 'u':
                for (int i = 0; i < current_A.size(); i++) //because we got adding array operator, just inbetween step
                {
                    vec.psi(i) = 0.1;

                    vec.zeta(i) = 0.1;
                }
                break;


                // for normal distribution
            case 'n':
                #pragma omp parallel for
                for (CaslInt i = 0; i<my_octree.number_Of_Leaves(); i++) //where adaptive meshing is done
                {
                    double x = my_octree.x_fr_i(my_octree.get_Leaf(i).icenter());
                    double y = my_octree.y_fr_j(my_octree.get_Leaf(i).jcenter());
                    double z = my_octree.z_fr_k(my_octree.get_Leaf(i).kcenter());

                    // used to be psi_initial, but just to make srue they have the same initial condition!!

                    vec.psi(i) = psi_initial(x,y,z); //first element of psi_n is the initial condition

                    vec.zeta(i) = zeta_initial(x,y,z);
                }

                break;

                //for exact solution to of the convergence analysis
                // exact solution is phi(x,y,z;t) = exp(x+y+z+t)

                case 'e':
                #pragma omp parallel for
                    std::cout<<"MADE IT HERE" << std::endl;
                    for (CaslInt i = 0; i<my_octree.number_Of_Leaves(); i++) //where adaptive meshing is done
                        {
                            double x = my_octree.x_fr_i(my_octree.get_Leaf(i).icenter());
                            double y = my_octree.y_fr_j(my_octree.get_Leaf(i).jcenter());
                            double z = my_octree.z_fr_k(my_octree.get_Leaf(i).kcenter());

                            // used to be psi_initial, but kist to ,ake srue they have the same initial condition!!

                            vec.psi(i) = exp(x+y+z); //first element of psi_n is the initial condition
                            vec.zeta(i) = vec.psi(i);
                        }

                            break;

            default:
                cout << "it didn't go through any of the cases, check and run again!" << endl;
                throw invalid_argument("[IC_Generator Error]: your choice is not implemented yet");

        }


        // bracket for if statement
    }

    else if (current_A.size() != current_B.size())
    {
        throw length_error("[IC_Generator Error]: both arrays must have the same size : ) ");

    }

    return vec;
}




IC_Generator gen_exact(char choice ,ArrayV<double> current_A,OcTreeCellBased my_octree, Initial_Solution psi_initial)
{

    double epsilon = EPSILON;

    // for now we are not gonna mess with num_funcs since we assume its always 2

    IC_Generator vec;
    vec.psi.resize_Without_Copy(my_octree.number_Of_Leaves());


   switch(choice)
        {
            case 'r':
                for (int i = 0; i < current_A.size(); i++) //because we got adding array operator, just inbetween step
                {
                    vec.psi(i) = rand() / (RAND_MAX + 1. + epsilon);
                }
                break;

            case 'u':
                for (int i = 0; i < current_A.size(); i++) //because we got adding array operator, just inbetween step
                {
                    vec.psi(i) = 0.1;
                }
                break;


                // for normal distribution
            case 'n':
                #pragma omp parallel for
                for (CaslInt i = 0; i<my_octree.number_Of_Leaves(); i++) //where adaptive meshing is done
                {
                    double x = my_octree.x_fr_i(my_octree.get_Leaf(i).icenter());
                    double y = my_octree.y_fr_j(my_octree.get_Leaf(i).jcenter());
                    double z = my_octree.z_fr_k(my_octree.get_Leaf(i).kcenter());

                    // used to be psi_initial, but just to make srue they have the same initial condition!!

                    vec.psi(i) = psi_initial(x,y,z); //first element of psi_n is the initial condition
                }

                break;

                //for exact solution to of the convergence analysis
                // exact solution is phi(x,y,z;t) = exp(x+y+z+t)


                case 'e':
                #pragma omp parallel for
                    std::cout<<"MADE IT HERE" << std::endl;
                    for (CaslInt i = 0; i<my_octree.number_Of_Leaves(); i++) //where adaptive meshing is done
                        {
                            double x = my_octree.x_fr_i(my_octree.get_Leaf(i).icenter());
                            double y = my_octree.y_fr_j(my_octree.get_Leaf(i).jcenter());
                            double z = my_octree.z_fr_k(my_octree.get_Leaf(i).kcenter());

                            // used to be psi_initial, but kist to ,ake srue they have the same initial condition!!

                            vec.psi(i) = exp(x+y+z); //first element of psi_n is the initial condition
                        }

                            break;

            default:
                cout << "it didn't go through any of the cases, check and run again!" << endl;
                throw invalid_argument("[IC_Generator Error]: your choice is not implemented yet");

        }

    return vec;
}
