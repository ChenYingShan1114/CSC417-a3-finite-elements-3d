#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

#include <iostream>
//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// volume - volume of tetrahedron
// C,D - material parameters
//Output:
//  energy - potential energy of tetrahedron
void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

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

        // calculate Neohookean strain energy
        // use F into the equation and put the energy value into e
        psi_neo_hookean(e, F, C, D);
        // std::cout << "auto neohooken linear tet energy: " << e << std::endl;
        
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    // std::cout << "V linear tetrahedron energy: " << energy << std::endl;
    
}