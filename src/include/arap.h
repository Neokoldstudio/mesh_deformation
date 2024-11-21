#ifndef ARAP_H
#define ARAP_H

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

class ARAP
{
public:
    ARAP(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

    void setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                        const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions);

    Eigen::MatrixXd computeDeformation();

private:
    Eigen::MatrixXd V, V_deformed;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> L, M; // Laplacian (L) and Mass (M) matrices
    Eigen::VectorXi anchor_indices, handle_indices;
    Eigen::MatrixXd anchor_displacement, handle_displacement;

    void computeLaplacian();
    void computeMassMatrix();
};

#endif