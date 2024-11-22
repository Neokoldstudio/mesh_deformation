#include "include/biharmonic.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>

BiharmonicDeformation::BiharmonicDeformation(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    : V(V), F(F), V_deformed(V)
{
    computeLaplacian();
    computeMassMatrix();
    K = L.transpose() * M.cwiseInverse() * L; // precompute the bilaplacian operator
}

void BiharmonicDeformation::setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                                           const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions)
{
    this->anchor_indices = anchor_indices;
    this->anchor_displacement = anchor_positions - V(anchor_indices, Eigen::all);
    this->handle_indices = handle_indices;
    this->handle_displacement = handle_positions - V(handle_indices, Eigen::all);
}

void BiharmonicDeformation::computeLaplacian()
{
    igl::cotmatrix(V, F, L);
}

void BiharmonicDeformation::computeMassMatrix()
{
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
}

Eigen::MatrixXd BiharmonicDeformation::computeDeformation()
{
    // Combine anchor and handle indices and positions to form a constraint vector
    Eigen::VectorXi d_user_indices(anchor_indices.size() + handle_indices.size());
    d_user_indices << anchor_indices, handle_indices;

    // Combine anchor and handle displacements to form the vector d_user
    Eigen::MatrixXd d_user(anchor_displacement.rows() + handle_displacement.rows(), anchor_displacement.cols());
    d_user << anchor_displacement, handle_displacement;

    // Compute the free indices
    Eigen::VectorXi all_indices = Eigen::VectorXi::LinSpaced(V.rows(), 0, V.rows() - 1);
    Eigen::VectorXi free_indices;
    Eigen::VectorXi _;
    igl::setdiff(all_indices, d_user_indices, free_indices, _);

    // Partition K into blocks (as seen in the calculation part, we do not need more than K_ff and K_fc)
    Eigen::SparseMatrix<double> K_ff, K_fc;
    igl::slice(K, free_indices, free_indices, K_ff);
    igl::slice(K, free_indices, d_user_indices, K_fc);

    // Solve for free displacements
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_ff);
    if (solver.info() != Eigen::Success)
    {
        throw std::runtime_error("Decomposition failed!");
    }

    Eigen::MatrixXd d_f = solver.solve(-K_fc * d_user);

    // Assemble the full displacement vector
    Eigen::MatrixXd d(K.rows(), d_user.cols());
    d.setZero();
    igl::slice_into(d_f, free_indices, 1, d);                   // Insert free displacements
    igl::slice_into(anchor_displacement, anchor_indices, 1, d); // Ensure anchors stay in place
    igl::slice_into(handle_displacement, handle_indices, 1, d); // Insert handle displacements

    // Apply the displacement to the original vertices
    V_deformed = V + d;

    return V_deformed;
}
