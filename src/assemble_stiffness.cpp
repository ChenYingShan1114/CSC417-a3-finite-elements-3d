#include <assemble_stiffness.h>
#include <iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  v0 - the mx1 vector of undeformed tetrahedron volumes
//  C,D - material parameters for the Neo-Hookean model
//Output:
//  K - the sparse, global stiffness matrix
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D) { 
        K.resize(q.size(), q.size());
        K.setZero();

        for (int i = 0; i < v0.size(); i++) {
            Eigen::Matrix1212d H = Eigen::Matrix1212d::Zero();
            d2V_linear_tetrahedron_dq2(H, q, V, T.row(i), v0(i), C, D);

            for (int j = 0; j < T.row(i).size(); j++) {
                for (int k = 0; k < T.row(i).size(); k++) {
                    int i1 = T(i, j), i2 = T(i, k);

                    K.coeffRef(3 * i1    , 3 * i2    ) -= H(3 * j    , 3 * k    );
                    K.coeffRef(3 * i1    , 3 * i2 + 1) -= H(3 * j    , 3 * k + 1);
                    K.coeffRef(3 * i1    , 3 * i2 + 2) -= H(3 * j    , 3 * k + 2);
                    K.coeffRef(3 * i1 + 1, 3 * i2    ) -= H(3 * j + 1, 3 * k    );
                    K.coeffRef(3 * i1 + 1, 3 * i2 + 1) -= H(3 * j + 1, 3 * k + 1);
                    K.coeffRef(3 * i1 + 1, 3 * i2 + 2) -= H(3 * j + 1, 3 * k + 2);
                    K.coeffRef(3 * i1 + 2, 3 * i2    ) -= H(3 * j + 2, 3 * k    );
                    K.coeffRef(3 * i1 + 2, 3 * i2 + 1) -= H(3 * j + 2, 3 * k + 1);
                    K.coeffRef(3 * i1 + 2, 3 * i2 + 2) -= H(3 * j + 2, 3 * k + 2);
                }
            }
        }    
    };
