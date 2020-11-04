IF(LIBIGL_WITH_OPENGL AND LIBIGL_WITH_OPENGL_GLFW AND LIBIGL_WITH_OPENGL_GLFW_IMGUI)
    # Viewer APP
    add_executable(falcon_igl_view_field view_field.cpp )
    target_link_libraries(falcon_igl_view_field igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
    INSTALL( TARGETS  falcon_igl_view_field)
ENDIF()

IF(HAVE_ZLIB)
    add_definitions(-DHAVE_ZLIB)
ENDIF()

add_executable(falcon_igl_field_eig field_eig.cpp )
target_link_libraries(falcon_igl_field_eig igl::core )

add_executable(falcon_igl_field_depth field_depth.cpp depth_potential.cpp)
target_link_libraries(falcon_igl_field_depth igl::core )

add_executable(falcon_igl_field_resample field_resample.cpp )
target_link_libraries(falcon_igl_field_resample igl::core )

add_executable(falcon_igl_mesh_split mesh_split.cpp )
target_link_libraries(falcon_igl_mesh_split igl::core )

add_executable(falcon_igl_multilabel_overlap multilabel_overlap.cpp )
target_link_libraries(falcon_igl_multilabel_overlap igl::core )


#
INSTALL( TARGETS 
        falcon_igl_mesh_split
        falcon_igl_field_resample
        falcon_igl_field_depth
        falcon_igl_field_eig
       )
