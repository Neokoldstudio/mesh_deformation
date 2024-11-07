#include <iostream>
#include "include/viewer.h"

using namespace Eigen;
using namespace std;
using namespace igl;

void view(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();
}