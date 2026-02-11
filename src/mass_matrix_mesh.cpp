#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
//Input:
//  qdot - generalized velocity for the FEM system
//  T - the mx4 vertex indices for tet mesh
//  density - density of material
//  v0 - the undeformed tetrahedra volumes
//Output:
//  M - Sparse mass matrix for the whole mesh.
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    M.resize(qdot.size(), qdot.size());
    M.setZero();
    for ( int i = 0; i < T.rows(); i++ ) {

        Eigen::Matrix1212d M0 = Eigen::Matrix1212d::Zero();
        mass_matrix_linear_tetrahedron(M0, qdot, T.row(i), density, v0(i));

        for (int j = 0; j < 4; j++) {
            int index1 = T(i, j);
            for (int k = 0; k < 4; k++) {
                int index2 = T(i, k);
                M.coeffRef( 3 * index1    , 3 * index2    ) += M0(3 * j    , 3 * k    );
                M.coeffRef( 3 * index1    , 3 * index2 + 1) += M0(3 * j    , 3 * k + 1);
                M.coeffRef( 3 * index1    , 3 * index2 + 2) += M0(3 * j    , 3 * k + 2);
                M.coeffRef( 3 * index1 + 1, 3 * index2    ) += M0(3 * j + 1, 3 * k    );
                M.coeffRef( 3 * index1 + 1, 3 * index2 + 1) += M0(3 * j + 1, 3 * k + 1);
                M.coeffRef( 3 * index1 + 1, 3 * index2 + 2) += M0(3 * j + 1, 3 * k + 2);
                M.coeffRef( 3 * index1 + 2, 3 * index2    ) += M0(3 * j + 2, 3 * k    );
                M.coeffRef( 3 * index1 + 2, 3 * index2 + 1) += M0(3 * j + 2, 3 * k + 1);
                M.coeffRef( 3 * index1 + 2, 3 * index2 + 2) += M0(3 * j + 2, 3 * k + 2);      
            }
        }
    }
}