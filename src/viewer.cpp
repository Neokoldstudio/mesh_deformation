#include <iostream>
#include "include/viewer.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/readOBJ.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/arap.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>

using namespace Eigen;
using namespace std;
using namespace igl;

void view(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    Eigen::MatrixXd C = Eigen::MatrixXd::Constant(F.rows(), 3, 1); // color matrix initialized to white
    igl::opengl::glfw::Viewer viewer;

    // Attach a guizmo plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuizmoWidget guizmo;

    plugin.widgets.push_back(&guizmo);

    guizmo.T.block(0, 3, 3, 1) = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff()).transpose().cast<float>();
    // Update can be applied relative to this remembered initial transform
    const Eigen::Matrix4f T0 = guizmo.T;
    // Attach callback to apply imguizmo's transform to mesh

    std::vector<int> handle_indices;
    int anchor_index = 0;
    Eigen::MatrixXd V_orig = V;
    Eigen::MatrixXd V_handle = V;
    Eigen::VectorXi b;  // Handle vertices
    Eigen::VectorXi a;  // Anchor vertices
    Eigen::MatrixXd bc; // Handle target positions
    Eigen::MatrixXd ac; // Anchor positions

    guizmo.callback = [&](const Eigen::Matrix4f &T)
    {
        const Eigen::Matrix4d TT = (T * T0.inverse()).cast<double>().transpose();
        V_handle = (V_orig.rowwise().homogeneous() * TT).rowwise().hnormalized();
        viewer.data().set_vertices(V_handle);
        viewer.data().compute_normals();
    };

    // Maya-style keyboard shortcuts for operation
    viewer.callback_key_pressed = [&](decltype(viewer) &, unsigned int key, int mod)
    {
        switch (key)
        {
        case ' ':
            guizmo.visible = !guizmo.visible;
            return true;
        case 'W':
        case 'w':
            guizmo.operation = ImGuizmo::TRANSLATE;
            return true;
        case 'E':
        case 'e':
            guizmo.operation = ImGuizmo::ROTATE;
            return true;
        case 'R':
        case 'r':
            guizmo.operation = ImGuizmo::SCALE;
            return true;
        }
        return false;
    };

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    std::cout << R"(
    W,w  Switch to translate operation
    E,e  Switch to rotate operation
    R,r  Switch to scale operation
    )";

    static char filePath[128] = ""; // Buffer for file path input

    // Handle matrix to store selected vertex indices

    Eigen::Vector2f down_mouse_pos, up_mouse_pos;
    bool shift_pressed = false;

    viewer.callback_key_down = [&](decltype(viewer) &, unsigned int key, int mod) -> bool
    {
        if (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT)
        {
            shift_pressed = true;
            double x, y;
            glfwGetCursorPos(viewer.window, &x, &y);

            // Invert the y-coordinate
            int window_height;
            int window_width;
            glfwGetWindowSize(viewer.window, &window_width, &window_height);
            y = window_height - y;

            down_mouse_pos = Eigen::Vector2f(x, y);
            printf("%f %f\n", down_mouse_pos(0), down_mouse_pos(1));
            return true;
        }
        return false;
    };

    viewer.callback_key_up = [&](decltype(viewer) &, unsigned int key, int mod) -> bool
    {
        if (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT)
        {
            shift_pressed = false;
            double x, y;

            glfwGetCursorPos(viewer.window, &x, &y);

            // Invert the y-coordinate
            int window_height;
            int window_width;
            glfwGetWindowSize(viewer.window, &window_width, &window_height);
            y = window_height - y;

            up_mouse_pos = Eigen::Vector2f(x, y);
            printf("%f %f\n", up_mouse_pos(0), up_mouse_pos(1));
            // Define the rectangle
            Eigen::Vector2f min_pos = down_mouse_pos.cwiseMin(up_mouse_pos);
            Eigen::Vector2f max_pos = down_mouse_pos.cwiseMax(up_mouse_pos);

            // Project 3D points to screen space and check if they are inside the rectangle
            handle_indices.clear();
            for (int i = 0; i < V.rows(); ++i)
            {
                Eigen::Vector3f proj;
                igl::project(V.row(i).cast<float>(), viewer.core().view, viewer.core().proj, viewer.core().viewport, proj);
                if ((proj.head<2>().array() >= min_pos.array()).all() && (proj.head<2>().array() <= max_pos.array()).all())
                {
                    handle_indices.push_back(i);
                }
            }

            // Update anchor vertices
            a.resize(handle_indices.size());
            for (int i = 0; i < handle_indices.size(); ++i)
            {
                a(i) = handle_indices[i];
            }
            ac = V_orig(a, Eigen::all);

            return true;
        }
        return false;
    };
    viewer.callback_pre_draw = [&](decltype(viewer) &) -> bool
    {
        viewer.data().clear_points();
        for (const auto &idx : handle_indices)
        {
            viewer.data().add_points(V.row(idx), Eigen::RowVector3d(1, 0, 0));
        }
        return false;
    };
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.data().show_lines = false;
    viewer.launch();
}
