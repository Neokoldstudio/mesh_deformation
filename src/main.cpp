#include <igl/eigs.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <queue>
#include "include/viewer.h"

using namespace Eigen;
using namespace std;
using namespace igl;

int main(int argc, char *argv[])
{
    {
        MatrixXd V, U;
        MatrixXi F;

        read_triangle_mesh("../assets/Bunny.obj", V, F);

        view(V, F);
    }
}