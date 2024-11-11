#include <iostream>
#include "include/viewer.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
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

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    static char filePath[128] = ""; // Buffer for file path input

    // Handle matrix to store selected vertex indices
    std::vector<int> handle_indices;
    std::vector<int> anchor_indices;
    Eigen::MatrixXd V_orig = V;
    Eigen::MatrixXd V_handle = V;
    Eigen::VectorXi b;  // Handle vertices
    Eigen::VectorXi a;  // Anchor vertices
    Eigen::MatrixXd bc; // Handle target positions
    Eigen::MatrixXd ac; // Anchor positions

    // ARAP data
    igl::ARAPData arap_data;
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
    arap_data.with_dynamics = true;

    bool arap_initialized = false;

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();
    };

    bool is_dragging = false;

    viewer.callback_mouse_down =
        [&V, &F, &C, &handle_indices, &anchor_indices, &V_handle, &V_orig, &b, &a, &bc, &ac, &arap_data, &is_dragging, &arap_initialized](igl::opengl::glfw::Viewer &viewer, int button, int modifier) -> bool
    {
        int fid;
        Eigen::Vector3f bc_f;
        // Cast a ray in the view direction starting from the mouse position
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                     viewer.core().proj, viewer.core().viewport, V, F, fid, bc_f))
        {
            // Get the vertex index from the barycentric coordinates
            Eigen::VectorXi vertices = F.row(fid);
            int selected_vertex;
            bc_f.maxCoeff(&selected_vertex);
            int vertex_index = vertices[selected_vertex];

            if (modifier == GLFW_MOD_SHIFT) // Check if Shift key is pressed for handles
            {
                // Add the selected vertex index to the handle matrix if not already present
                if (std::find(handle_indices.begin(), handle_indices.end(), vertex_index) == handle_indices.end())
                {
                    handle_indices.push_back(vertex_index);
                }

                // Update the list of handle vertices and their target positions
                b = Eigen::Map<const Eigen::VectorXi>(handle_indices.data(), handle_indices.size());
                bc.resize(handle_indices.size(), 3);
                for (int i = 0; i < handle_indices.size(); ++i)
                {
                    bc.row(i) = V_handle.row(handle_indices[i]);
                }

                // Reinitialize ARAP precomputation with updated handles and anchors
                Eigen::VectorXi combined_indices(b.size() + a.size());
                combined_indices << b, a;
                igl::arap_precomputation(V, F, 3, combined_indices, arap_data);
                arap_initialized = true;

                // Paint the selected face blue
                C.row(fid) << 0, 0, 1;
                viewer.data().set_colors(C);
                return true;
            }
            else if (modifier == GLFW_MOD_ALT) // Check if Alt key is pressed for anchors
            {
                // Add the selected vertex index to the anchor matrix if not already present
                if (std::find(anchor_indices.begin(), anchor_indices.end(), vertex_index) == anchor_indices.end())
                {
                    anchor_indices.push_back(vertex_index);
                }

                // Update the list of anchor vertices and their positions
                a = Eigen::Map<const Eigen::VectorXi>(anchor_indices.data(), anchor_indices.size());
                ac.resize(anchor_indices.size(), 3);
                for (int i = 0; i < anchor_indices.size(); ++i)
                {
                    ac.row(i) = V.row(anchor_indices[i]);
                }

                // Paint the selected face red
                C.row(fid) << 1, 0, 0;
                viewer.data().set_colors(C);
                return true;
            }
        }
        else if (modifier == GLFW_MOD_CONTROL) // Check if Control key is pressed for dragging
        {
            is_dragging = true; // Start dragging
            printf("dragging\n");
        }
        return false;
    };

    viewer.callback_mouse_up =
        [&is_dragging](igl::opengl::glfw::Viewer &, int, int) -> bool
    {
        is_dragging = false; // Stop dragging
        return false;
    };

    viewer.callback_mouse_move =
        [&V_handle, &b, &bc, &arap_data, &viewer, &is_dragging, &arap_initialized](igl::opengl::glfw::Viewer &, int mouse_x, int mouse_y) -> bool
    {
        if (!is_dragging || !arap_initialized || b.size() == 0)
            return false;

        // Update the position of all handle vertices based on mouse movement
        Eigen::Vector3f win(mouse_x, viewer.core().viewport(3) - mouse_y, 0.5);
        Eigen::Vector3f obj = igl::unproject(win, viewer.core().view, viewer.core().proj, viewer.core().viewport);

        // Calculate the translation vector
        Eigen::Vector3f lastVertexPosition = V_handle.row(b[0]).cast<float>();

        Eigen::Vector3f translation = obj - lastVertexPosition;

        // Apply the translation to all handle vertices
        for (int i = 0; i < b.size(); ++i)
        {
            bc.row(i) += translation.cast<double>();
        }

        // Solve ARAP with the updated handle positions
        igl::arap_solve(bc, arap_data, V_handle);

        // Update the mesh with the deformed vertices
        viewer.data().set_vertices(V_handle);
        viewer.data().compute_normals();
        return true;
    };

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.data().show_lines = false;
    viewer.launch();
}
