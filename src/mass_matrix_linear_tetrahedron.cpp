 
 #include <mass_matrix_linear_tetrahedron.h>
 #include <iostream>
//Input:
//  qdot - generalized velocity for the FEM system
//  element - the 1x4 vertex indices for this tetrahedron
//  density - density of material
//  volume - the undeformed tetrahedron volume
//Output:
//  M - dense 12x12 per-tetrahedron mass matrix
 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    const double mass = density * volume / 20;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            if (i == j) {
                M(i, j) = 2 * mass;
            } else if (((j % 3) - (i % 3)) == 0 ) {
                M(i, j) = mass;
            }
        }
    }
 }