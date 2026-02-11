#include <phi_linear_tetrahedron.h>
#include <iostream>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the undeformed space at which to compute the energy density
//Output:
//  phi - the 4x1 values the basis functions
void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < T.cols(); j++) {
            T(i, j) = V(element(j+1), i) - V(element(0), i);
        }
    }
    Eigen::Matrix3d T_inv = T.inverse();

    Eigen::Vector3d dX = Eigen::Vector3d::Zero();
    for (int i = 0; i < dX.size(); i++) {
        dX(i) = X(i) - V(element(0), i);
    }

    Eigen::Vector3d tmp_phi = Eigen::Vector3d::Zero();
    tmp_phi = T_inv * dX;
    phi << 1 - tmp_phi(0) - tmp_phi(1) - tmp_phi(2), tmp_phi(0), tmp_phi(1), tmp_phi(2);
}