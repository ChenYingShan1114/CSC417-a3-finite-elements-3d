#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>

//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the neohookean energy
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {

    // std::cout << F.size() << " " << C << " " << D << " " << psi << std::endl;
    // std::cout << F << std::endl;

    // reference https://en.wikipedia.org/wiki/Neo-Hookean_solid

    // principal stretches
    // double lambda1 = F(0, 0), lambda2 = F(1, 1), lambda3 = F(2, 2);

    // left Cauchy-Green deformation tensor
    // Eigen::Matrix3d B;
    // B = F * F.transpose();

    // principal invariants
    // double I1 = pow(lambda1, 2) + pow(lambda2, 2) + pow(lambda3, 2);
    double I1 = F.squaredNorm();
    // double I2 = pow(lambda1 * lambda2, 2) + pow(lambda1 * lambda3, 2) + pow(lambda2 * lambda3, 2);
    // double I3 = pow(lambda1 * lambda2 * lambda3, 2);

    // Jacobian of F
    // double J = lambda1 * lambda2 * lambda3;
    double J = F.determinant();

    // compressible neo-Hookean strain energy density
    psi = C * (I1 * pow(J, -2.0 / 3.0) - 3) + D * pow((J - 1), 2);

    // incompressible neo-Hookean strain energy density
    // psi = C * (I1 - 3);

}