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

Eigen::VectorXd ARAP::SC_L1_Shrinkage(const Eigen::VectorXd &x, double lambda)
{
    // Ensure the input vector type is Eigen::VectorXd
    Eigen::VectorXd y = x;
    double norm_x = x.norm();
    if (norm_x > lambda)
    {
        y = (1.0 - lambda / norm_x) * x;
    }
    else
    {
        y.setZero();
    }
    return y;
}

Eigen::MatrixXd ARAP::computeDeformation()
{
    const int max_iterations = 4;
    const double tolerance = 1e-5;
    Eigen::MatrixXd V_prev = V;

    const double s = 1.0; // regularisation parameter

    // Combine anchor and handle indices
    Eigen::VectorXi constrained_indices(anchor_indices.size() + handle_indices.size());
    constrained_indices << anchor_indices, handle_indices;

    // Combine anchor and handle positions
    Eigen::MatrixXd constrained_positions(anchor_positions.rows() + handle_positions.rows(), anchor_positions.cols());
    constrained_positions << anchor_positions, handle_positions;

    // Create the mask of unconstrained indices
    std::vector<int> unconstrained_indices_vector;
    for (int i = 0; i < V.rows(); ++i)
    {
        if (std::find(constrained_indices.data(), constrained_indices.data() + constrained_indices.size(), i) == constrained_indices.data() + constrained_indices.size())
        {
            unconstrained_indices_vector.push_back(i);
        }
    }
    Eigen::VectorXi unconstrained_indices = Eigen::Map<Eigen::VectorXi>(unconstrained_indices_vector.data(), unconstrained_indices_vector.size());

    // Initialize the linear solver for the reduced system
    Eigen::SparseMatrix<double> L_reduced;
    igl::slice(L, unconstrained_indices, unconstrained_indices, L_reduced);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_reduced);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Decomposition failed!" << std::endl;
        return V_deformed;
    }

    // Initialize variables for ADMM
    Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(V.rows(), 3); // For SC-L1 shrinkage
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(V.rows(), 3); // For deformation gradient
    Eigen::MatrixXd V_update = V;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // Step 1: Update X (symmetric factor of deformation gradient)
        for (int j = 0; j < V.rows(); ++j)
        {
            Eigen::Matrix3d Cov = Eigen::Matrix3d::Zero();
            for (int k : neighborsMatrix[j])
            {
                Eigen::Vector3d p = V_prev.row(j) - V_prev.row(k);
                Eigen::Vector3d p_prime = V_deformed.row(j) - V_deformed.row(k);
                double wij = L.coeff(j, k);
                Cov += wij * p * p_prime.transpose();
            }

            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d U = svd.matrixU();
            Eigen::Matrix3d V = svd.matrixV();
            X.row(j) = (V * U.transpose()).row(0); // Taking the first row (the symmetric factor)
        }

        // Step 2: Update Z using SC-L1 shrinkage
        for (int j = 0; j < V.rows(); ++j)
        {
            Z.row(j) = SC_L1_Shrinkage(V_deformed.row(j) - V_prev.row(j), 0.1); // Update with SC-L1 shrinkage
        }

        // Step 3: Update vertex positions (V) with the linear system
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(unconstrained_indices.size(), 3);
        for (int i = 0; i < unconstrained_indices.size(); ++i)
        {
            int j = unconstrained_indices[i];
            for (int k : neighborsMatrix[j])
            {
                Eigen::Vector3d p = V_prev.row(j) - V_prev.row(k);
                double wij = L.coeff(j, k);
                B.row(i) += 0.5 * wij * ((X.row(j) + X.row(k)) * p);
            }
        }

        // Solve the reduced system
        Eigen::MatrixXd P_prime_unconstrained = solver.solve(B);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Solving failed!" << std::endl;
            return V_deformed;
        }

        // Reconstruct the full solution
        V_deformed = V_prev;
        igl::slice_into(P_prime_unconstrained, unconstrained_indices, 1, V_deformed);
        igl::slice_into(constrained_positions, constrained_indices, 1, V_deformed);

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
