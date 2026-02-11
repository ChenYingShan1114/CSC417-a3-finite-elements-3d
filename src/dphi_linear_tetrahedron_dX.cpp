#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x4 vertex indices for this tetrahedron
//  X - the position in the undeformed space at which to compute the energy density
//Output:
//  dphi - the 4x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < T.cols(); j++) {
            T(i, j) = V(element(j+1), i) - V(element(0), i);
        }
    }
    Eigen::Matrix3d T_inv = T.inverse();    

    Eigen::Vector3d one = Eigen::Vector3d::Zero();
    one << 1, 1, 1;

    Eigen::Vector3d tmp_phi0 = Eigen::Vector3d::Zero();
    tmp_phi0 = -one.transpose() * T_inv;

    dphi << tmp_phi0.transpose(), T_inv;
}