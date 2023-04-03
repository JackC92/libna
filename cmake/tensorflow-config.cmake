# This minimal CMake Config Mode file for imporing TensorFlow should be copied to the root of the libtensorflow-(cpu|gpu) folder
get_filename_component(_TensorFlow_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)

add_library(TensorFlow::tensorflow STATIC IMPORTED)
set_target_properties(TensorFlow::tensorflow PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_TensorFlow_PREFIX}/include"
    IMPORTED_LOCATION "${_TensorFlow_PREFIX}/lib/tensorflow.lib")

set(TensorFlow_INCLUDES ${_TensorFlow_PREFIX}/include)
set(TensorFlow_LIBRARIES TensorFlow::tensorflow)

unset(_TensorFlow_PREFIX)
