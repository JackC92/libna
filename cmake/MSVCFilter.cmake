if (MSVC)
    foreach(header IN LISTS HEADER_FILES)
        get_filename_component(header_path "${header}" PATH)
        file(RELATIVE_PATH header_path_rel "${CMAKE_CURRENT_SOURCE_DIR}/include/na/" "${header_path}")
        string(REPLACE "/" "\\" header_path_rel_msvc "${header_path_rel}")
        source_group("Header Files\\${header_path_rel_msvc}" FILES "${header}")
    endforeach()
    foreach(source IN LISTS SRC_FILES)
        get_filename_component(source_path "${source}" PATH)
        file(RELATIVE_PATH source_path_rel "${CMAKE_CURRENT_SOURCE_DIR}/src/na/" "${source_path}")
        string(REPLACE "/" "\\" source_path_rel_msvc "${source_path_rel}")
        source_group("Source Files\\${source_path_rel_msvc}" FILES "${source}")
    endforeach()
    foreach(test IN LISTS TEST_FILES)
        if (${test} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp)
            continue()
        endif()
        message(STATUS ${test})
        get_filename_component(test_path "${test}" PATH)
        file(RELATIVE_PATH test_path_rel "${CMAKE_CURRENT_SOURCE_DIR}/tests/" "${test_path}")
        string(REPLACE "/" "\\" test_path_rel_msvc "${test_path_rel}")
        source_group("Source Files\\${test_path_rel_msvc}" FILES "${test}")
    endforeach()
endif()
