#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  v0 - the undeformed tetrahedron volume
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  dV - the 12x1 gradient of the potential energy for a single tetrahedron
void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

        // calculate dphi for deformation gradient F
        Eigen::Matrix43d dphi = Eigen::Matrix43d::Zero();
        dphi_linear_tetrahedron_dX(dphi, V, element, X);

        Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        Eigen::MatrixXd q_tet;
        q_tet.resize(3, 4);
        q_tet.setZero();
        for (int i = 0; i < q_tet.rows(); i++) {
            for (int j = 0; j < q_tet.cols(); j++) {
                q_tet(i, j) = q(3 * element(j) + i);
            }
        }
        F = q_tet * dphi;

        // calculate gradient of Neohookean strain energy
        // use F into the equation and put the vector dV
        Eigen::Vector9d dw_9 = Eigen::Vector9d::Zero(); // 9x1
        dpsi_neo_hookean_dF(dw_9, F, C, D);

        Eigen::MatrixXd B;
        B.resize(dw_9.size(), dV.size());
        B.setZero();
        for (int i = 0; i < dV.rows() / dphi.rows(); i++) {
            for (int j = 0; j < dV.size(); j++) {
                B(j % 3 + 3 * i, 3 * (j / 3) + i) = dphi(j / 3, j % 3);
            }
        }
        dV = B.transpose() * dw_9;

        // std::cout << "neohookean linear dV: " << dV << std::endl;
        
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    // std::cout << volume << " dV_linear_tetrahedron_dq dV: " << dV << std::endl;

}