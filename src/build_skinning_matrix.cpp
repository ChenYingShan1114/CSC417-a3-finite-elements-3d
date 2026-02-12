#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  T - the mx4 tetrahedron connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  V_skin - lx3 matrix of vertices of the display mesh
//Output:
//  N - the lxn sparse skinning matrix 
// n = V.rows(); m = T.rows(); l = V_skin.rows();

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    N.resize(V_skin.rows(), V.rows());
    N.setZero();

    for (int i = 0; i < V_skin.rows(); i++) {
        for (int j = 0; j < T.rows(); j++) {

            // use each skinning point to calculate the shade function phi of every tetrahedron 
            Eigen::Vector4d phi = Eigen::Vector4d::Zero();
            phi_linear_tetrahedron(phi, V, T.row(j), V_skin.row(i).transpose());

            // check whether this point i is inside this tetrahedron j
            if (phi(0) >= 0 and phi(1) >= 0 and phi(2) >= 0 and phi(3) >= 0) {
                
                // if yes, but the value in the sparse skinning matrix N
                for (int k = 0; k < T.row(j).size(); k++) {
                    N.coeffRef(i, T(j, k)) = phi(k);
                }
                break;
            }
        }
    }
}