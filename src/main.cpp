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

void usage()
{
    cout << "Usage: mesh_deformation [-i input_file] [-h | --help]" << endl;
    cout << "Options:" << endl;
    cout << "  -i input_file   Specify the input mesh file (default is bar.obj)" << endl;
    cout << "  -h, --help      Show this help message" << endl;
}

int main(int argc, char *argv[])
{
    string input_file = "../assets/bar.obj"; // default file

    for (int i = 1; i < argc; ++i)
    {
        if (string(argv[i]) == "-i" && i + 1 < argc)
        {
            input_file = argv[i + 1];
            break;
        }

        if (string(argv[i]) == "-h" || string(argv[i]) == "--help")
        {
            usage();
            return 0;
        }
    }

    MatrixXd V, U;
    MatrixXi F;

    if (!read_triangle_mesh(input_file, V, F))
    {
        cerr << "Failed to read mesh from " << input_file << endl;
        return 1;
    }

    view(V, F);

    return 0;
}