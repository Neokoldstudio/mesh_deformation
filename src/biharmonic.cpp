#include "include/biharmonic.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

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
    Eigen::SparseMatrix<double> K = L.transpose() * M.cwiseInverse() * L;


    for (int i = 0; i < anchor_indices.size(); i++)
    {
        Aeq.insert(i, anchor_indices(i)) = 1;
        Beq(i) = anchor_positions(i);
    }

    for (int i = 0; i < handle_indices.size(); i++)
    {
        Aeq.insert(i + anchor_indices.size(), handle_indices(i)) = 1;
        Beq(i + anchor_indices.size()) = handle_positions(i);
    }

    // Solve the system of equations
    Eigen::SparseMatrix<double> AeqT = Aeq.transpose();
    Eigen::SparseMatrix<double> AeqT_Aeq = AeqT * Aeq;
    Eigen::VectorXd AeqT_Beq = AeqT * Beq;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(AeqT_Aeq);
    xeq = solver.solve(AeqT_Beq);

    Eigen::VectorXd Beq_prime = Beq - Aeq * xeq;
    Eigen::VectorXd B_prime = B - A * xeq;

    solver.compute(A);
    x = solver.solve(B_prime);

    return V_deformed = x + xeq;
}
