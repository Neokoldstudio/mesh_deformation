#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>
#include "include/viewer.h"

int main(int argc, char *argv[])
{
    {
        Eigen::MatrixXd V,U;
        Eigen::MatrixXi F;

        using namespace Eigen;
        using namespace std;
        using namespace igl;
        igl::read_triangle_mesh("../assets/bunny.obj", V, F);

        view(V, F);
    }
}