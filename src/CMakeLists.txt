
#adding 3rd party subprojects
# Eigen 
#ADD_SUBDIRECTORY(3rd_party/eigen)
#add_library(igl_eigen INTERFACE)
#target_include_directories(igl_eigen SYSTEM INTERFACE
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/eigen>
#  $<INSTALL_INTERFACE:include>
#)
#set_property(TARGET igl_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
#add_library(Eigen3::Eigen ALIAS igl_eigen)

# libigl
ADD_SUBDIRECTORY(3rd_party/libigl)

# STB
#ADD_SUBDIRECTORY(3rd_party/stb)

INCLUDE( ${LIBMINC_USE_FILE} )

IF(HAVE_BICPL)
  INCLUDE_DIRECTORIES(${BICPL_INCLUDE_DIRS})
  LINK_DIRECTORIES(${BICPL_LIBRARY_DIRS})
ENDIF(HAVE_BICPL)

LINK_LIBRARIES(${LIBMINC_LIBRARIES})

# IF(HAVE_LIBTIFF)
#   INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
#   LINK_LIBRARIES(${TIFF_LIBRARY})
# ENDIF(HAVE_LIBTIFF)

IF(HAVE_ZLIB)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
  add_definitions(-DHAVE_ZLIB=1)
ENDIF(HAVE_ZLIB)

# default verbosity level, change for debugging?
SET(FALCON_VERBOSE 0)

SET(FALCON_DATA_DEFAULT ${CMAKE_INSTALL_PREFIX} )
SET(FALCON_DATA_SOURCE  ${CMAKE_SOURCE_DIR})
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/falcon_config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/falcon_config.h @ONLY)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
   3rd_party/rply-1.1.4
   3rd_party/libcsv
   3rd_party/NrrdIO
   3rd_party/
   3rd_party/zstr/src
   3rd_party/spectra/include
   3rd_party/stb
   ${CMAKE_CURRENT_BINARY_DIR}/3rd_party/NrrdIO
)

SET(RPLY_LIB_SRC
    3rd_party/rply-1.1.4/rply.c
    3rd_party/rply-1.1.4/rplyfile.h
    3rd_party/rply-1.1.4/rply.h
)

SET(CSV_LIB_SRC
    3rd_party/libcsv/csv.h  3rd_party/libcsv/libcsv.c
)

# NrrdIO
ADD_SUBDIRECTORY(3rd_party/NrrdIO)

# RapidCSV
ADD_SUBDIRECTORY(3rd_party/csv)


# minc2-simple
#ADD_SUBDIRECTORY(3rd_party/minc2-simple/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/minc2-simple/src)

add_library(minc2-simple STATIC 
             3rd_party/minc2-simple/src/minc2-simple.c
             3rd_party/minc2-simple/src/minc2-matrix-ops.c
             3rd_party/minc2-simple/src/minc2-simple.h
             3rd_party/minc2-simple/src/minc2-simple-int.h
)

target_include_directories(minc2-simple INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/minc2-simple/src)
set_property(TARGET minc2-simple PROPERTY EXPORT_NAME minc2-simple)

# LibIGL based tools
ADD_SUBDIRECTORY(igl)



SET(FALCON_LIB_HDR
    nifti1.h
    nifti1_io.h
    falcon.h
    falcon_seg.h
    falcon_surfaces.h
    falcon_cortex.h
    falcon_inline.h
    falcon_point_inline.h
    falcon_sph_inline.h
    falcon_tritri_inline.h
    ${CMAKE_CURRENT_BINARY_DIR}/falcon_config.h
)

SET(FALCON_LIB_SRC
    falcon_stb.c
    falcon_config.c
    znzlib.c 
    nifti1_io.c 
    falcon.c 
    falcon_morph.c 
    falcon_nelder_mead.c 
    falcon_bla.c 
    falcon_point.c 
    falcon_interpolate.c 
    falcon_affine.c 
    falcon_aregister.c 
    falcon_aregister2.c 
    falcon_bounds.c 
    falcon_bseg.c 
    falcon_threshold.c 
    falcon_off.c 
    falcon_off_cuberille.c 
    falcon_off_remesh.c 
    falcon_off_bbox.c
    falcon_off_bbox_sph.c 
    falcon_off_obj2img.c 
    falcon_off_shrinkwrap.c 
    falcon_sobel.c 
    falcon_laplace_map.c 
    falcon_spmat.c 
    falcon_median.c 
    falcon_nbcsr.c 
    falcon_nregister.c 
    falcon_nregister_bspline.c 
    falcon_nregister_demons.c 
    falcon_off_fsasc2off.c 
    falcon_venseg.c 
    falcon_histogram.c 
    falcon_cortex.c 
    falcon_distance_map.c 
    falcon_inorm.c 
    falcon_dbc.c 
    falcon_lesion_fill.c 
    falcon_feature_map.c 
    falcon_bimodal_fit.c 
    falcon_segment.c 
    falcon_PABIC.c 
    falcon_cortex_edge.c 
    falcon_jacobian_map.c 
    falcon_siena.c 
    falcon_curve_fitting.c 
    falcon_off_ocv.c 
    falcon_pseudoT2.c 
    falcon_nregister_intrasubject_template.c 
    falcon_interpacket.c 
    falcon_mnc2nii.c 
    falcon_nii2mnc.c 
    falcon_morph_gray.c 
    falcon_rand.c
    falcon_nlseg.c
    mincio.c
    falcon_off_curvature.c
    falcon_cortex_debug.c
    falcon_ply.c
    falcon_slice.c
    falcon_text_io.c
    ${RPLY_LIB_SRC}
    ${CSV_LIB_SRC}
)

ADD_LIBRARY(falcon_lib ${FALCON_LIB_HDR}  ${FALCON_LIB_SRC})

LINK_LIBRARIES(falcon_lib)

ADD_EXECUTABLE(falcon_calc    falcon_calc.c)

ADD_EXECUTABLE(falcon_mean    niik_fmean.c)

ADD_EXECUTABLE(falcon_median  niik_fmedian.c)

ADD_EXECUTABLE(falcon_stdv    niik_fstdv.c)

ADD_EXECUTABLE(falcon_math    niikmath.c )

ADD_EXECUTABLE(falcon_bpf     niik_bpf.c )

ADD_EXECUTABLE(falcon_stx_registration  niik_stx_reg.c )

ADD_EXECUTABLE(falcon_xfm2mat   falcon_xfm2niikmat.c )

ADD_EXECUTABLE(falcon_mat2xfm   falcon_niikmat2xfm.c )

ADD_EXECUTABLE(falcon_inpaint   niik_inpaint.c )

ADD_EXECUTABLE(falcon_calc_pT2  niik_calc_pT2.c )

IF(BUILD_NVK AND GTK2_FOUND)
  ADD_EXECUTABLE(NVK falcon_viewer.c )
  add_definitions(${GTK2_DEFINITIONS})
  INCLUDE_DIRECTORIES(${GTK2_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(NVK ${GTK2_LIBRARIES} )
  INSTALL( TARGETS NVK  DESTINATION bin)
ENDIF(BUILD_NVK AND GTK2_FOUND)

ADD_EXECUTABLE(falcon_simu_focalpath_test     niik_simu_focalpath_test.c)

ADD_EXECUTABLE(falcon_test1                   niik_classifyEM.c)

ADD_EXECUTABLE(falcon_skull_segment           niik_skull_segment.c)

ADD_EXECUTABLE(falcon_nls_skull               niik_nls_skull.c)

ADD_EXECUTABLE(falcon_simu_cortical_atrophy   niik_simu_cortical_atrophy.c)

ADD_EXECUTABLE(falcon_cortex_refine           niikcortex_refine.c falcon_cortex_deformation.c )

ADD_EXECUTABLE(falcon_surface_refine          niikcortex_surface_refine.c falcon_cortex_deformation.c )

ADD_EXECUTABLE(falcon_cortex_smooth           niikcortex_smooth.c  falcon_cortex_smoothing.c falcon_cortex_deformation.c )

ADD_EXECUTABLE(falcon_cortex_remesh           niikcortex_remesh.c  falcon_cortex_deformation.c )

ADD_EXECUTABLE(falcon_cortex_initics          niikcortex_initics.c falcon_cortex_initics.c )

ADD_EXECUTABLE(falcon_surface_check           niikcortex_check_surface.c )

ADD_EXECUTABLE(falcon_cortex_initocs          niikcortex_initocs.c falcon_cortex_initocs.c )

ADD_EXECUTABLE(falcon_cortex_initocs_simple   niikcortex_initocs_simple.c falcon_cortex_initocs.c )

ADD_EXECUTABLE(falcon_cortex_shrinkwrap       niikcortex_shrinkwrap.c  )

ADD_EXECUTABLE(falcon_cortex_mod              niikcortex_mod.c  )

ADD_EXECUTABLE(falcon_off2asc                 niikcortex_off2asc.c)

ADD_EXECUTABLE(falcon_surface_split           niikcortex_surface_split.c)

ADD_EXECUTABLE(falcon_surface_smooth          niikcortex_surface_smooth.c)

ADD_EXECUTABLE(falcon_transform_surface       niikcortex_transform_surface.c)

ADD_EXECUTABLE(falcon_transform_surface_ants  niikcortex_transform_surface_ants.c)
TARGET_LINK_LIBRARIES(falcon_transform_surface_ants NrrdIO)

ADD_EXECUTABLE(falcon_avg                     niik_avg.c)

ADD_EXECUTABLE(falcon_cortex_calc_thickness   niikcortex_calc_thickness.c) # TODO: depricate in favor of falcon_cortex_measure

ADD_EXECUTABLE(falcon_field_smooth            falcon_field_smooth.c)

ADD_EXECUTABLE(falcon_cortex_measure          niikcortex_measure.c)

ADD_EXECUTABLE(falcon_cortex_calc_curvature   niikcortex_calc_curvature.c  )

ADD_EXECUTABLE(falcon_resample_field_sph      niik_resample_field_sph.c)

ADD_EXECUTABLE(falcon_surface_signals         falcon_surface_signals.c
                                              volume_tools.c  
                                              volume_tools.h
               )

ADD_EXECUTABLE(falcon_surface_colour          falcon_surface_colour.c)

ADD_EXECUTABLE(falcon_midsag                  falcon_midsag.c)

ADD_EXECUTABLE(falcon_slice_extract           falcon_slice_extract.c)

ADD_EXECUTABLE(falcon_csv_stats               falcon_csv_stats.c)

#test prog
ADD_EXECUTABLE(falcon_test_sphere             falcon_test_sphere.c )


IF(BICPL_FOUND)
  ADD_EXECUTABLE(falcon_obj2off              falcon_mniobj2off.c)
  ADD_EXECUTABLE(falcon_off2obj              falcon_off2mniobj.c)
  
  TARGET_LINK_LIBRARIES(falcon_off2obj      ${BICPL_LIBRARIES})
  TARGET_LINK_LIBRARIES(falcon_obj2off      ${BICPL_LIBRARIES})

  INSTALL( TARGETS
    falcon_off2obj
    falcon_obj2off
    DESTINATION bin
  )  
ENDIF(BICPL_FOUND)

ADD_EXECUTABLE(falcon_off2ply        falcon_off2ply.c)
ADD_EXECUTABLE(falcon_ply2off        falcon_ply2off.c)
ADD_EXECUTABLE(falcon_off2pov        falcon_off2pov.c)
ADD_EXECUTABLE(falcon_ply_info       falcon_ply_info.c)

INSTALL( TARGETS 
#  falcon_calc
#  falcon_mean
#  falcon_median
#  falcon_stdv
#  falcon_avg
#  falcon_bpf
#  falcon_calc_pT2
  falcon_field_smooth
  falcon_csv_stats
  falcon_cortex_calc_thickness
  falcon_cortex_measure
  falcon_cortex_refine
  falcon_surface_refine
  falcon_cortex_initics
  falcon_cortex_initocs
  falcon_cortex_shrinkwrap
  falcon_surface_colour
  falcon_slice_extract
  falcon_resample_field_sph
  falcon_cortex_smooth
  falcon_surface_smooth
  falcon_surface_check
  falcon_cortex_mod
  falcon_surface_split
  falcon_transform_surface
  falcon_cortex_calc_curvature
  falcon_resample_field_sph
  falcon_surface_signals
  falcon_midsag
  falcon_math
  falcon_off2pov
  falcon_off2ply
  falcon_ply2off
#  falcon_inpaint
#  falcon_mat2xfm
#  falcon_nls_skull
#  falcon_simu_cortical_atrophy
#  falcon_simu_focalpath_test
#  falcon_skull_segment
#  falcon_stx_registration
#  falcon_test1
#  falcon_xfm2mat
  DESTINATION bin)


#ADD_SUBDIRECTORY(surface_tools)
