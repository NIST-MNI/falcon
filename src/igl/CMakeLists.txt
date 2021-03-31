IF( LIBIGL_WITH_OPENGL_GLFW_IMGUI) # LIBIGL_WITH_OPENGL AND LIBIGL_WITH_OPENGL_GLFW AND
    message("Building interactive viewer")
    # Viewer APP
    add_executable(falcon_igl_view_field view_field.cpp )
    target_link_libraries(falcon_igl_view_field igl::csv igl::core igl::opengl_glfw igl::opengl_glfw_imgui)

    INSTALL( TARGETS  falcon_igl_view_field)
ENDIF()

IF(HAVE_ZLIB)
    add_definitions(-DHAVE_ZLIB)
ENDIF()


if(NOT TARGET igl::png)
    add_definitions(-DSTB_PNG)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/png ${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/stb/)
    add_library(stb_image STATIC 
        ../3rd_party/stb/stb_image.h 
        ../3rd_party/stb/stb_image_write.h
        png/igl_stb_image.h 
        png/igl_stb_image.cpp 
    )
    if(LIBIGL_USE_STATIC_LIBRARY)
        set(SOURCES_IGL_PNG
           png/readPNG.cpp  
           png/readPNG.h  
           png/writePNG.cpp  
           png/writePNG.h
        )
        add_library(igl_png STATIC ${SOURCES_IGL_PNG})
        target_link_libraries(igl_png ${IGL_SCOPE} igl_common)
        target_link_libraries(igl_png ${IGL_SCOPE} igl)
        target_link_libraries(igl_png ${IGL_SCOPE} stb_image)
    else()
        add_library(igl_png INTERFACE)
        target_link_libraries(igl_png INTERFACE stb_image)
    endif()
    target_include_directories(igl_png INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/png ${CMAKE_CURRENT_SOURCE_DIR}/../3rd_party/stb)
    add_library(igl::png ALIAS igl_png)
    set_property(TARGET igl_png PROPERTY EXPORT_NAME igl::png)
endif()

add_executable(falcon_igl_field_eig field_eig.cpp )
target_link_libraries(falcon_igl_field_eig igl::csv igl::core )

add_executable(falcon_igl_field_depth field_depth.cpp depth_potential.cpp)
target_link_libraries(falcon_igl_field_depth igl::csv igl::core )

add_executable(falcon_igl_field_resample field_resample.cpp )
target_link_libraries(falcon_igl_field_resample igl::csv igl::core )

add_executable(falcon_igl_mesh_split mesh_split.cpp )
target_link_libraries(falcon_igl_mesh_split igl::core )

add_executable(falcon_igl_multilabel_overlap multilabel_overlap.cpp )
target_link_libraries(falcon_igl_multilabel_overlap igl::csv igl::core )

add_executable(falcon_igl_surface_demons surface_demons.cpp depth_potential.cpp)
target_link_libraries(falcon_igl_surface_demons igl::csv igl::core )

add_executable(falcon_igl_mesh_decimate mesh_decimate.cpp )
target_link_libraries(falcon_igl_mesh_decimate igl::core )

add_executable(falcon_igl_mesh_smooth mesh_smooth.cpp )
target_link_libraries(falcon_igl_mesh_smooth igl::core )

add_executable(falcon_igl_mesh_debug mesh_debug.cpp )
target_link_libraries(falcon_igl_mesh_debug igl::core )

add_executable(falcon_igl_mesh_render mesh_render.cpp )
target_link_libraries(falcon_igl_mesh_render igl::png igl::csv igl::core  igl::embree )


add_executable(falcon_igl_marching_cubes marching_cubes.cpp)
target_link_libraries(falcon_igl_marching_cubes igl::core minc2-simple )

#
INSTALL( TARGETS 
        falcon_igl_mesh_split
        falcon_igl_field_resample
        falcon_igl_field_depth
        falcon_igl_field_eig
        falcon_igl_multilabel_overlap
        falcon_igl_surface_demons
        falcon_igl_mesh_render
        falcon_igl_marching_cubes
        falcon_igl_mesh_smooth
        falcon_igl_mesh_decimate
       )
