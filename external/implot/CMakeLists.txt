file(GLOB IMPLOT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/implot/*.h)
file(GLOB IMPLOT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/implot/implot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot/implot_items.cpp)

add_library(implot OBJECT ${IMPLOT_HEADERS} ${IMPLOT_SOURCES})
set_target_properties(implot PROPERTIES CXX_STANDARD 17)
target_link_libraries(implot PRIVATE imgui)
target_compile_definitions(implot PRIVATE IMPLOT_DEBUG IMPLOT_DLL_EXPORT IMPLOT_BACKEND_ENABLE_OPENGL3 IMGUI_IMPL_OPENGL_LOADER_GLAD)
target_include_directories(implot INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/implot)
if(MSVC)
    target_compile_options(implot PRIVATE /W4 /WX /arch:AVX2 /fp:fast /permissive-) 
else()
    target_compile_options(implot PRIVATE -Wall -Wextra -pedantic -Werror -mavx2 -Ofast)
endif()
