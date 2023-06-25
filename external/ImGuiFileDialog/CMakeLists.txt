file(GLOB CUSTOM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/CustomImGuiFileDialogConfig.h)
file(GLOB RES_SOURCES 
		${CMAKE_CURRENT_SOURCE_DIR}/CustomFont.cpp 
		${CMAKE_CURRENT_SOURCE_DIR}/CustomFont.h)
file(GLOB IMGUIFILEDIALOG_SOURCES 
		${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/*.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/ImGuiFileDialog/*.h)
file(GLOB IMGUIIMPL_SOURCES 
		${CMAKE_SOURCE_DIR}/external/polyscope/deps/imgui/imgui/backends/imgui_impl_glfw.cpp
		${CMAKE_SOURCE_DIR}/external/polyscope/deps/imgui/imgui/backends/imgui_impl_glfw.h
		${CMAKE_SOURCE_DIR}/external/polyscope/deps/imgui/imgui/backends/imgui_impl_opengl3.cpp
		${CMAKE_SOURCE_DIR}/external/polyscope/deps/imgui/imgui/backends/imgui_impl_opengl3.h)
source_group(src\\Custom FILES ${CUSTOM_SOURCES})
source_group(src\\Res FILES ${RES_SOURCES})
source_group(src\\ImGuiFileDialog FILES ${IMGUIFILEDIALOG_SOURCES})
source_group(src\\ImGuiImpl FILES ${IMGUIIMPL_SOURCES})

add_library(imgui_filedialog STATIC
	${CUSTOM_SOURCES}
	${RES_SOURCES}
	${IMGUIFILEDIALOG_SOURCES}
	${IMGUIIMPL_SOURCES}
)
add_definitions(-DUSE_STD_FILESYSTEM)
set_target_properties(imgui_filedialog PROPERTIES CXX_STANDARD 17)
target_link_libraries(imgui_filedialog PRIVATE glad glfw imgui)
