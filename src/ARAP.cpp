#include "include/arap.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

ARAP::ARAP(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    : V(V), F(F), V_deformed(V)
{
    computeLaplacian();
    computeMassMatrix();
}

void ARAP::setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                          const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions)
{
    this->anchor_indices = anchor_indices;
    this->anchor_displacement = anchor_positions - V(anchor_indices, Eigen::all);
    this->handle_indices = handle_indices;
    this->handle_displacement = handle_positions - V(handle_indices, Eigen::all);
}

void ARAP::computeLaplacian()
{
    igl::cotmatrix(V, F, L); // Compute the cotangent Laplacian matrix
}

void ARAP::computeMassMatrix()
{
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M); // Compute the mass matrix
}

Eigen::MatrixXd ARAP::computeDeformation()
{
    // Step 1: Compute the bi-Laplacian K
    Eigen::SparseMatrix<double> K = L.transpose() * M.cwiseInverse() * L; // TODO fix this : we want to avoid computing the inverse of M explicitly

    // Step 2: Combine anchor and handle indices and positions
    Eigen::VectorXi b(anchor_indices.size() + handle_indices.size());
    b << anchor_indices, handle_indices;

    Eigen::MatrixXd bc(anchor_displacement.rows() + handle_displacement.rows(), anchor_displacement.cols());
    bc << anchor_displacement,
        handle_displacement;

    // Step 3: Compute the free indices
    Eigen::VectorXi all_indices = Eigen::VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1);
    Eigen::VectorXi free_indices;
    Eigen::VectorXi _;
    igl::setdiff(all_indices, b, free_indices, _);

    // Step 4: Partition K into blocks
    Eigen::SparseMatrix<double> K_ff, K_fc;
    igl::slice(K, free_indices, free_indices, K_ff); // Free-Free block
    igl::slice(K, free_indices, b, K_fc);            // Free-Constrained block

    // Step 5: Solve for free displacements
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_ff);
    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Decomposition failed!");
    }

    Eigen::MatrixXd d_f = solver.solve(-K_fc * bc);

    // Step 6: Assemble the full displacement vector
    Eigen::MatrixXd d(K.rows(), bc.cols());
    d.setZero();
    igl::slice_into(d_f, free_indices, 1, d);                   // Insert free displacements
    igl::slice_into(anchor_displacement, anchor_indices, 1, d); // Ensure anchors stay in place
    igl::slice_into(handle_displacement, handle_indices, 1, d); // Insert handle displacements

    // Step 7: Apply the displacement to the original vertices
    V_deformed = V + d;

    return V_deformed;
}
