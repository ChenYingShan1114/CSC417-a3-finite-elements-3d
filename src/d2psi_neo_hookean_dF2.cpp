#include <d2psi_neo_hookean_dq2.h>
#include <iostream>
//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  psi - the 9x9 Hessian of the potential energy wrt to the deformation gradient
void d2psi_neo_hookean_dF2(Eigen::Matrix99d &ddw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    ddw.setZero();

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

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            // dF^{-T}_{ij}
            Eigen::Matrix3d dF_inv_T = Eigen::Matrix3d::Zero();
            dF_inv_T = -F_inv_T.col(j) * F_inv_T.row(i);

            // delta
            Eigen::Matrix3d delta = Eigen::Matrix3d::Zero();
            delta(i, j) = 1;

            // compressible neo-Hookean strain energy density -> second Piola-Kirchhoff
            Eigen::Matrix3d dFF = Eigen::Matrix3d::Zero();
            dFF = 2 * C * pow(J, -2.0 / 3.0) * (-2.0 / 3.0 * F_inv_T * (F(i, j) - 1.0 / 3.0 * I1 * F_inv_T(i, j)) 
                            + (delta - 2.0 / 3.0 * F * F_inv_T(i, j) - 1.0 / 3.0 * I1 * dF_inv_T))
                + 2 * D * J * ((2 * J - 1) * F_inv_T * F_inv_T(i, j) + (J - 1) * dF_inv_T);

            // put in ddw
            for (int k1 = 0; k1 < dFF.rows(); k1++) {
                for (int k2 = 0; k2 < dFF.cols(); k2++) {
                    ddw(3 * i + j, 3 * k1 + k2) = dFF(k1, k2);
                }
            }
        }
    }

    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         // ex. ij
    //         std::cout << "ex. " << i << " " << j << std::endl;
    //         int i1, i2, j1, j2;
    //         i1 = (i + 1)%3;
    //         i2 = (i + 2)%3;
    //         j1 = (j + 1)%3;
    //         j2 = (j + 2)%3;
    //         // std::cout << i1 << " " << i2 << " " << j1 << " " << j2 << std::endl;

    //         // Derivatives of cofactors that contain F(i, j)
    //         Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
    //         tmp(i1, j1) = F(i2, j2);   // dC(i1,j1)/dF(i,j)
    //         tmp(i2, j2) = F(i1, j1);  // dC(i2,j2)/dF(i,j)
    //         tmp(i1, j2) = -F(i2, j1);  // dC(i1,j2)/dF(i,j)
    //         tmp(i2, j1) = -F(i1, j2);  // dC(i2,j1)/dF(i,j)
    //         // std::cout << tmp << std::endl << std::endl;

    //         Eigen::Matrix3d manual = Eigen::Matrix3d::Zero();
    //         manual = -F_inv_T * F_inv_T(i, j) + 1/J * tmp;
    //         // std::cout << F_inv_T << std::endl << std::endl;
    //         // std::cout << manual << std::endl << std::endl;

    //         // std::cout << F_inv_T.col(0) << " " << F_inv_T.row(2) << std::endl << std::endl;
    //         Eigen::Matrix3d identity = Eigen::Matrix3d::Zero();
    //         identity = -F_inv_T.col(j) * F_inv_T.row(i);
    //         // std::cout << -F_inv_T.col(j) * F_inv_T.row(i) << std::endl << std::endl;
    //         // std::cout << identity - manual << std::endl << std::endl;

    //         // std::cout << "check" << std::endl << std::endl;
    //     }
    // }
}