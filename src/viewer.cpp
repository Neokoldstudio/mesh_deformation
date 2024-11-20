#include <iostream>
#include "include/viewer.h"
#include "include/Constraints.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuizmoWidget.h>
#include <igl/readOBJ.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/arap.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>
#include "include/biharmonic.h"

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
    guizmo.operation = ImGuizmo::TRANSLATE;
    guizmo.visible = false;
    // Update can be applied relative to this remembered initial transform
    Eigen::Matrix4f T0 = guizmo.T;

    Constraints constraints;
    Eigen::MatrixXd V_orig = V;
    Eigen::MatrixXd V_handle = V;
    Eigen::VectorXi b;  // Handle vertices
    Eigen::VectorXi a;  // Anchor vertices
    Eigen::MatrixXd bc; // Handle target positions
    Eigen::MatrixXd ac; // Anchor positions

    // Store original positions of handle vertices
    Eigen::MatrixXd handle_positions;

    viewer.data().point_size = 10;

    // Biharmonic deformation object
    BiharmonicDeformation biharmonic(V, F); // Create the deformation object

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
        }
        return false;
    };

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    std::cout << R"(
    SPACE Toggle guizmo
    W,w  Switch to translate operation
    E,e  Switch to rotate operation
    SHIFT DOWN -> SHIFT UP Select anchored vertices by dragging a rectangle
    ALT + LEFT CLICK Select handle vertices
    )";

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
            // Define the rectangle
            Eigen::Vector2f min_pos = down_mouse_pos.cwiseMin(up_mouse_pos);
            Eigen::Vector2f max_pos = down_mouse_pos.cwiseMax(up_mouse_pos);

            // Project 3D points to screen space and check if they are inside the rectangle
            constraints.clearAnchorIndices();
            for (int i = 0; i < V_handle.rows(); ++i)
            {
                Eigen::Vector3f proj;
                igl::project(V_handle.row(i).cast<float>(), viewer.core().view, viewer.core().proj, viewer.core().viewport, proj);
                if ((proj.head<2>().array() >= min_pos.array()).all() && (proj.head<2>().array() <= max_pos.array()).all())
                {
                    constraints.addAnchorIndex(i);
                }
            }

            return true;
        }
        return false;
    };

    viewer.callback_mouse_down = [&](decltype(viewer) &, int button, int modifier) -> bool
    {
        if (button == GLFW_MOUSE_BUTTON_LEFT && (modifier & GLFW_MOD_ALT))
        {
            int fid;
            Eigen::Vector3f bc;
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                         viewer.core().proj, viewer.core().viewport,
                                         viewer.data().V, viewer.data().F, fid, bc))
            {
                // Find the closest vertex
                Eigen::MatrixXd::Index maxIndex;
                bc.maxCoeff(&maxIndex);
                int closest_vertex = viewer.data().F(fid, maxIndex);
                constraints.addHandleIndex(closest_vertex);
            }

            Eigen::MatrixXd V_constraints(constraints.getHandleSize(), V.cols());
            for (int i = 0; i < constraints.getHandleSize(); ++i)
            {
                V_constraints.row(i) = V.row(constraints.getHandleIndex(i));
            }

            // Store original positions of handle vertices
            handle_positions.resize(constraints.getHandleSize(), V.cols());
            for (int i = 0; i < constraints.getHandleSize(); ++i)
            {
                handle_positions.row(i) = V.row(constraints.getHandleIndex(i));
            }

            guizmo.T.block(0, 3, 3, 1) = 0.5 * (V_constraints.colwise().maxCoeff() + V_constraints.colwise().minCoeff()).transpose().cast<float>();
            guizmo.visible = true;

            return true;
        }

        return false;
    };

    guizmo.callback = [&](const Eigen::Matrix4f &T)
    {
        const Eigen::Matrix4d TT = (T * T0.inverse()).cast<double>().transpose();
        for (int i = 0; i < constraints.getHandleSize(); ++i)
        {
            int idx = constraints.getHandleIndex(i);
            V_handle.row(idx) = (handle_positions.row(i).homogeneous() * TT).hnormalized();
        }
        viewer.data().set_vertices(V_handle);
        viewer.data().compute_normals();

        // Update the transformation in the constraints object
        constraints.updateTransformation(TT.cast<float>());

        // After updating the handles, perform biharmonic deformation
        if (constraints.getHandleSize() > 0)
        {
            // Set the handle positions
            Eigen::MatrixXd handle_positions_updated(constraints.getHandleSize(), 3);
            for (int i = 0; i < constraints.getHandleSize(); ++i)
            {
                int idx = constraints.getHandleIndex(i);
                handle_positions_updated.row(i) = V_handle.row(idx);
            }

            // get anchor positions
            Eigen::MatrixXd ac(constraints.getAnchorSize(), 3);
            for (int i = 0; i < constraints.getAnchorSize(); ++i)
            {
                int idx = constraints.getAnchorIndex(i);
                ac.row(i) = V_orig.row(idx);
            }

            // Convert std::vector<int> to Eigen::VectorXi
            Eigen::VectorXi anchor_indices(constraints.getAnchorIndices().size());
            for (size_t i = 0; i < constraints.getAnchorIndices().size(); ++i)
            {
                anchor_indices(i) = constraints.getAnchorIndices()[i];
            }

            Eigen::VectorXi handle_indices(constraints.getHandleIndices().size());
            for (size_t i = 0; i < constraints.getHandleIndices().size(); ++i)
            {
                handle_indices(i) = constraints.getHandleIndices()[i];
            }

            biharmonic.setConstraints(anchor_indices, ac, handle_indices, handle_positions_updated);
            Eigen::MatrixXd V_deformed = biharmonic.computeDeformation();
            viewer.data().set_vertices(V_deformed);
            viewer.data().compute_normals();
        }
    };

    viewer.callback_pre_draw = [&](decltype(viewer) &) -> bool
    {
        viewer.data().clear_points();
        for (const auto &idx : constraints.getAnchorIndices())
        {
            viewer.data().add_points(V_handle.row(idx), Eigen::RowVector3d(1, 0, 0));
        }

        for (const auto &idx : constraints.getHandleIndices())
        {
            viewer.data().add_points(V_handle.row(idx), Eigen::RowVector3d(0, 1, 1));
        }
        return false;
    };

    menu.callback_draw_viewer_menu = [&constraints, &V_handle, &V_orig, &viewer, &guizmo, &V, &handle_positions, &biharmonic]()
    {
        // Add a text input field for loading constraints
        static char loadFilePath[256] = ""; // Separate buffer for load file path input
        ImGui::InputText("Load Constraints JSON File Path", loadFilePath, IM_ARRAYSIZE(loadFilePath));

        // Add a button for loading constraints
        if (ImGui::Button("Load Constraints"))
        {
            // Load the constraints when the button is pressed
            if (strlen(loadFilePath) > 0)
            {
                constraints.importConstraints(loadFilePath);

                const Eigen::Matrix4d TT = constraints.getTransformation().cast<double>();

                for (const auto &index : constraints.getHandleIndices())
                {
                    V_handle.row(index) = (V_orig.row(index).homogeneous() * TT).hnormalized();
                }

                Eigen::MatrixXd V_constraints(constraints.getHandleSize(), V.cols());
                for (int i = 0; i < constraints.getHandleSize(); ++i)
                {
                    V_constraints.row(i) = V.row(constraints.getHandleIndex(i));
                }

                // Store original positions of handle vertices
                handle_positions.resize(constraints.getHandleSize(), V.cols());
                for (int i = 0; i < constraints.getHandleSize(); ++i)
                {
                    handle_positions.row(i) = V.row(constraints.getHandleIndex(i));
                }

                guizmo.T.block(0, 3, 3, 1) = 0.5 * (V_constraints.colwise().maxCoeff() + V_constraints.colwise().minCoeff()).transpose().cast<float>();
                guizmo.visible = true;

                // Apply biharmonic deformation
                if (constraints.getHandleSize() > 0)
                {
                    // Set the handle positions
                    Eigen::MatrixXd handle_positions_updated(constraints.getHandleSize(), 3);
                    for (int i = 0; i < constraints.getHandleSize(); ++i)
                    {
                        int idx = constraints.getHandleIndex(i);
                        handle_positions_updated.row(i) = V_handle.row(idx);
                    }

                    // Get anchor positions
                    Eigen::MatrixXd ac(constraints.getAnchorSize(), 3);
                    for (int i = 0; i < constraints.getAnchorSize(); ++i)
                    {
                        int idx = constraints.getAnchorIndex(i);
                        ac.row(i) = V_orig.row(idx);
                    }

                    // Convert std::vector<int> to Eigen::VectorXi
                    Eigen::VectorXi anchor_indices(constraints.getAnchorIndices().size());
                    for (size_t i = 0; i < constraints.getAnchorIndices().size(); ++i)
                    {
                        anchor_indices(i) = constraints.getAnchorIndices()[i];
                    }

                    Eigen::VectorXi handle_indices(constraints.getHandleIndices().size());
                    for (size_t i = 0; i < constraints.getHandleIndices().size(); ++i)
                    {
                        handle_indices(i) = constraints.getHandleIndices()[i];
                    }

                    biharmonic.setConstraints(anchor_indices, ac, handle_indices, handle_positions_updated);
                    Eigen::MatrixXd V_deformed = biharmonic.computeDeformation();
                    viewer.data().set_vertices(V_deformed);
                    viewer.data().compute_normals();
                }
            }
        }

        // Add a space
        ImGui::Spacing();

        // Add a text input field for saving constraints
        static char saveFilePath[256] = ""; // Separate buffer for save file path input
        ImGui::InputText("Save Constraints JSON File Path", saveFilePath, IM_ARRAYSIZE(saveFilePath));

        // Add a button for saving constraints
        if (ImGui::Button("Save Constraints"))
        {
            // Save the constraints when the button is pressed
            if (strlen(saveFilePath) > 0)
            {
                constraints.exportConstraints(saveFilePath);
            }
        }
    };

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
    viewer.data().show_lines = false;
    viewer.launch();
}
