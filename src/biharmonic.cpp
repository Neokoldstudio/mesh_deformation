#include "include/biharmonic.h"
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

BiharmonicDeformation::BiharmonicDeformation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    : V(V), F(F), V_deformed(V)
{
    computeLaplacian();
    computeMassMatrix();
}

void BiharmonicDeformation::setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                                           const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions)
{
    this->anchor_indices = anchor_indices;
    this->anchor_positions = anchor_positions;
    this->handle_indices = handle_indices;
    this->handle_positions = handle_positions;
}

void BiharmonicDeformation::computeLaplacian()
{
    igl::cotmatrix(V, F, L); // Compute the cotangent Laplacian matrix
}

void BiharmonicDeformation::computeMassMatrix()
{
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M); // Compute the mass matrix
}

Eigen::MatrixXd BiharmonicDeformation::computeDeformation()
{
    // Compute the right-hand side
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), V.cols());
    for (int i = 0; i < anchor_indices.size(); ++i)
    {
        B.row(anchor_indices(i)) = anchor_positions.row(i);
    }
    for (int i = 0; i < handle_indices.size(); ++i)
    {
        B.row(handle_indices(i)) = handle_positions.row(i);
    }

    // Construct the biharmonic operator
    Eigen::SparseMatrix<double> M_inv = M.cwiseInverse();
    Eigen::SparseMatrix<double> L2 = L.transpose() * M_inv * L;

    // Solve the linear system
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L2);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Decomposition failed" << std::endl;
        return V_deformed;
    }
    V_deformed = solver.solve(B);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Solving failed" << std::endl;
        return V_deformed;
    }

    return V_deformed;
}
