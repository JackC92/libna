file(GLOB IMGUIFILEDIALOG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/*.h)
file(GLOB IMGUIFILEDIALOG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/*.cpp)

add_library(imgui_filedialog OBJECT
    ${IMGUIFILEDIALOG_HEADERS}
    ${IMGUIFILEDIALOG_SOURCES}
)
set_target_properties(imgui_filedialog PROPERTIES CXX_STANDARD 17)
target_link_libraries(imgui_filedialog PRIVATE glad glfw imgui)
target_compile_definitions(imgui_filedialog PRIVATE USE_STD_FILESYSTEM)
target_include_directories(imgui_filedialog INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog)
