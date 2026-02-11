#include <dpsi_neo_hookean_dF.h>
#include <iostream>

//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the 9x1 gradient of the potential energy wrt to the deformation gradient
void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    dw.setZero();

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

    // F^{-T}
    Eigen::Matrix3d F_inv = F.inverse();
    Eigen::Matrix3d F_inv_T = F_inv.transpose();

    // compressible neo-Hookean strain energy density -> first Piola-Kirchhoff
    Eigen::Matrix3d dwdF = Eigen::Matrix3d::Zero();
    dwdF = 2 * C * pow(J, -2.0 / 3.0) * (F - 1.0 / 3.0 * I1 * F_inv_T) + 2 * D * J * (J - 1) * F_inv_T;

    // gradient of potnetial energy
    for (int i = 0; i < dwdF.rows(); i++) {
        for (int j = 0; j < dwdF.cols(); j++) {
            dw(3 * i + j) = dwdF(i, j);
        }
    }
}