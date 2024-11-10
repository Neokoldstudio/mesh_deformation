#include <iostream>
#include "include/viewer.h"
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readOBJ.h>

using namespace Eigen;
using namespace std;
using namespace igl;

void view(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
    igl::opengl::glfw::Viewer viewer;

    // attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    double doubleVariable = 0.1f; // Shared between two menus

    static char filePath[128] = ""; // Buffer for file path input

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        ImGui::InputText("Mesh File Path", filePath, IM_ARRAYSIZE(filePath));
        if (ImGui::Button("Load Mesh"))
        {
            std::string filePathStr(filePath);
            if (igl::readOBJ(filePathStr, V, F)) // Change the read function based on the mesh format
            {
                viewer.data().clear();
                viewer.data().set_mesh(V, F);
            }
            else
            {
                std::cerr << "Failed to load mesh: " << filePathStr << std::endl;
            }
        }
    };

    viewer.data().set_mesh(V, F);
    viewer.launch();
}