#ifndef ARAP_H
#define ARAP_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

class ARAP
{
public:
    ARAP(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    void setConstraints(const Eigen::VectorXi &anchor_indices, const Eigen::MatrixXd &anchor_positions,
                        const Eigen::VectorXi &handle_indices, const Eigen::MatrixXd &handle_positions);
    Eigen::MatrixXd computeDeformation();

private:
    void computeLaplacian();
    void computeMassMatrix();
    void initialiseNeighbours(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_deformed;
    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> M;
    Eigen::VectorXi anchor_indices;
    Eigen::MatrixXd anchor_positions;
    Eigen::VectorXi handle_indices;
    Eigen::MatrixXd handle_positions;
    std::vector<std::vector<int>> neighborsMatrix;

    Eigen::VectorXd SC_L1_Shrinkage(const Eigen::VectorXd &, double);
};

#endif
