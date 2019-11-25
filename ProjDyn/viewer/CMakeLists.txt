add_library(viewer INTERFACE)
target_include_directories(viewer INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(viewer INTERFACE surface_mesh)
target_link_libraries(viewer INTERFACE nanogui ${NANOGUI_EXTRA_LIBS})