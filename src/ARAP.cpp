#include "include/arap.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/setdiff.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <vector>

ARAP::ARAP(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    : V(V), F(F), V_deformed(V)
{
    computeLaplacian();
    computeMassMatrix();
    initialiseNeighbours(V, F);
}

void ARAP::setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                          const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions)
{
    this->anchor_indices = anchor_indices;
    this->anchor_positions = anchor_positions;
    this->handle_indices = handle_indices;
    this->handle_positions = handle_positions;
}

void ARAP::computeLaplacian()
{
    igl::cotmatrix(V, F, L); // Compute the cotangent Laplacian matrix
}

void ARAP::computeMassMatrix()
{
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M); // Compute the mass matrix
}

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <vector>
#include <iostream>

Eigen::MatrixXd ARAP::computeDeformation()
{
    const int max_iterations = 10;
    const double tolerance = 1e-5;
    Eigen::MatrixXd V_prev = V;

    // Combine anchor and handle indices
    Eigen::VectorXi constrained_indices(anchor_indices.size() + handle_indices.size());
    constrained_indices << anchor_indices, handle_indices;

    // Combine anchor and handle positionss
    Eigen::MatrixXd constrained_positions(anchor_positions.rows() + handle_positions.rows(), anchor_positions.cols());
    constrained_positions << anchor_positions, handle_positions;

    // Initialize the linear solver
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    for (int i = 0; i < max_iterations; i++)
    {
        std::vector<Eigen::Matrix3d> R(V.rows(), Eigen::Matrix3d::Identity());

        // Compute the B matrix
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(V.rows(), 3);
        for (int j = 0; j < V.rows(); j++)
        {
            for (int k : neighborsMatrix[j])
            {
                Eigen::Vector3d p = V_prev.row(j) - V_prev.row(k);
                double wij = L.coeff(j, k); // Use the cotangent laplacian to compute the weight
                B.row(j) += 0.5 * wij * ((R[j] + R[k]) * p).transpose();
            }
        }

        // Modify the B matrix to incorporate constraints
        igl::slice_into(-constrained_positions, constrained_indices, 1, B);

        // Create a modified Laplacian matrix
        Eigen::SparseMatrix<double> L_modified = L;
        for (int k = 0; k < constrained_indices.size(); k++)
        {
            int idx = constrained_indices[k];
            for (Eigen::SparseMatrix<double>::InnerIterator it(L_modified, idx); it; ++it)
            {
                it.valueRef() = 0.0;
            }
            L_modified.coeffRef(idx, idx) = 1.0;
        }

        // Recompute the decomposition with the modified Laplacian matrix
        solver.compute(L_modified);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Decomposition failed!" << std::endl;
            return V_deformed;
        }

        // Solve for the new vertex positions
        Eigen::MatrixXd P_prime = solver.solve(B);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Solving failed!" << std::endl;
            return V_deformed;
        }
        V_deformed = P_prime;

        // Update rotations
        for (int j = 0; j < V.rows(); j++)
        {
            Eigen::Matrix3d Cov = Eigen::Matrix3d::Zero();
            for (int k : neighborsMatrix[j])
            {
                Eigen::Vector3d p = V_prev.row(j) - V_prev.row(k);
                Eigen::Vector3d p_prime = V_deformed.row(j) - V_deformed.row(k);
                double wij = L.coeff(j, k); // Use the cotangent laplacian to compute the weight
                Cov += wij * p * p_prime.transpose();
            }

            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d U = svd.matrixU();
            Eigen::Matrix3d V = svd.matrixV();
            Eigen::Matrix3d signCorrection = Eigen::Matrix3d::Identity();

            double det = (V * U.transpose()).determinant();
            if (det < 0)
                signCorrection(2, 2) = -1;

            R[j] = U * signCorrection * V.transpose();
        }

        // Check for convergence
        if ((V_deformed - V_prev).norm() < tolerance)
        {
            break;
        }
        V_prev = V_deformed;
    }
    return V_deformed;
}

void ARAP::initialiseNeighbours(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    // Initialize the neighbor list as empty lists for each vertex
    neighborsMatrix.resize(V.rows());

    for (int i = 0; i < F.rows(); ++i)
    {
        for (int j = 0; j < F.cols(); ++j)
        {
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % F.cols());
            int v2 = F(i, (j + 2) % F.cols());

            // Add neighbors, ensuring no duplicates
            if (std::find(neighborsMatrix[v0].begin(), neighborsMatrix[v0].end(), v1) == neighborsMatrix[v0].end())
                neighborsMatrix[v0].push_back(v1);
            if (std::find(neighborsMatrix[v0].begin(), neighborsMatrix[v0].end(), v2) == neighborsMatrix[v0].end())
                neighborsMatrix[v0].push_back(v2);

            if (std::find(neighborsMatrix[v1].begin(), neighborsMatrix[v1].end(), v0) == neighborsMatrix[v1].end())
                neighborsMatrix[v1].push_back(v0);
            if (std::find(neighborsMatrix[v1].begin(), neighborsMatrix[v1].end(), v2) == neighborsMatrix[v1].end())
                neighborsMatrix[v1].push_back(v2);

            if (std::find(neighborsMatrix[v2].begin(), neighborsMatrix[v2].end(), v0) == neighborsMatrix[v2].end())
                neighborsMatrix[v2].push_back(v0);
            if (std::find(neighborsMatrix[v2].begin(), neighborsMatrix[v2].end(), v1) == neighborsMatrix[v2].end())
                neighborsMatrix[v2].push_back(v1);
        }
    }
}
