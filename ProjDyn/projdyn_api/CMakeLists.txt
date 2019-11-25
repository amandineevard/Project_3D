add_library(projdyn_api INTERFACE)
target_include_directories(projdyn_api INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(projdyn_api INTERFACE projdyn)
target_link_libraries(projdyn_api INTERFACE surface_mesh)
target_link_libraries(projdyn_api INTERFACE tetgen)
target_link_libraries(projdyn_api INTERFACE viewer)
target_link_libraries(projdyn_api INTERFACE nanogui ${NANOGUI_EXTRA_LIBS})