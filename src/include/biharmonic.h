#ifndef BIHARMONIC_DEFORMATION_H
#define BIHARMONIC_DEFORMATION_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

class BiharmonicDeformation
{
public:
    BiharmonicDeformation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

    void setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                        const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions);

    Eigen::MatrixXd computeDeformation();

private:
    Eigen::MatrixXd V, V_deformed;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> L, M; // Laplacian (L) and Mass (M) matrices
    Eigen::VectorXi anchor_indices, handle_indices;
    Eigen::MatrixXd anchor_positions, handle_positions;

    void computeLaplacian();
    void computeMassMatrix();
};

#endif // BIHARMONIC_DEFORMATION_H