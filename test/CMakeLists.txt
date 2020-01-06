SET(TESTDIR        ${CMAKE_SOURCE_DIR}/test)
SET(FALCON_BIN     ${CMAKE_BINARY_DIR}/src)
SET(FALCON_SCRIPTS ${CMAKE_SOURCE_DIR}/bin)
SET(FALCON_DATA    ${CMAKE_BINARY_DIR})

SET(MODELDIR  ${CMAKE_SOURCE_DIR}/data/icbm152_model_09c)
SET(OUTDIR    ${CMAKE_CURRENT_BINARY_DIR})

# for running proper NIIK command
SET(FALCON_TEST_ENVIRONMENT "PATH=${FALCON_BIN}:${FALCON_SCRIPTS}:${MINC_TOOLKIT_DIR}/bin:$ENV{PATH}")

IF(NOT APPLE)
    LIST(APPEND FALCON_TEST_ENVIRONMENT "LD_LIBRARY_PATH=${FALCON_BIN}:${MINC_TOOLKIT_DIR}/lib:$ENV{LD_LIBRARY_PATH}")
ELSE(NOT APPLE)
    LIST(APPEND FALCON_TEST_ENVIRONMENT "DYLD_LIBRARY_PATH=${FALCON_BIN}:${MINC_TOOLKIT_DIR}/lib:$ENV{DYLD_LIBRARY_PATH}")
ENDIF(NOT APPLE)

LIST(APPEND FALCON_TEST_ENVIRONMENT "FALCON_HOME=${CMAKE_SOURCE_DIR}")
LIST(APPEND FALCON_TEST_ENVIRONMENT "FALCON_DATA=${CMAKE_SOURCE_DIR}")
LIST(APPEND FALCON_TEST_ENVIRONMENT "FALCON_BIN=${FALCON_BIN}")
LIST(APPEND FALCON_TEST_ENVIRONMENT "FALCON_SCRIPTS=${FALCON_SCRIPTS}")

SET(FALCON_TEST_ENVIRONMENT_PRL ${FALCON_TEST_ENVIRONMENT})
LIST(APPEND FALCON_TEST_ENVIRONMENT "OMP_NUM_THREADS=1")
LIST(APPEND FALCON_TEST_ENVIRONMENT_PRL "OMP_NUM_THREADS=4")

message("FALCON_TEST_ENVIRONMENT:${FALCON_TEST_ENVIRONMENT}")


INCLUDE( ${LIBMINC_USE_FILE} )
# simple tests of basic functionality
INCLUDE_DIRECTORIES(../src  ${CMAKE_CURRENT_BINARY_DIR}/../src)

LINK_LIBRARIES(falcon_lib ${LIBMINC_LIBRARIES})

IF(HAVE_LIBTIFF)
    INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIRS})
    LINK_LIBRARIES(${TIFF_LIBRARY})
ENDIF(HAVE_LIBTIFF)

IF(HAVE_ZLIB)
    INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
    LINK_LIBRARIES(${ZLIB_LIBRARIES})
ENDIF(HAVE_ZLIB)


ADD_EXECUTABLE(test_text_io  test_text_io.c )
ADD_TEST(NAME text_io COMMAND test_text_io ${CMAKE_CURRENT_SOURCE_DIR}/table.csv ${CMAKE_CURRENT_BINARY_DIR}/table_test.csv)

ADD_EXECUTABLE(test_off  test_off.c )
ADD_TEST(NAME test_off COMMAND test_off ${CMAKE_CURRENT_SOURCE_DIR}/test.ply ${CMAKE_CURRENT_BINARY_DIR}/test.ply)


MACRO(ADD_FALCON_TEST name cmd)
   ADD_TEST( NAME ${name}
       COMMAND ${cmd} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
        ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
        ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
   )
   set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${FALCON_TEST_ENVIRONMENT}")
   
ENDMACRO(ADD_FALCON_TEST)

MACRO(ADD_FALCON_TEST_PRL name cmd)
   ADD_TEST( NAME ${name}
       COMMAND ${cmd} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
        ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
        ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
   )
   set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${FALCON_TEST_ENVIRONMENT_PRL}")
   
ENDMACRO(ADD_FALCON_TEST_PRL)


# # create zero 
# ADD_FALCON_TEST(CREATE_ZERO minccalc -express 0 -byte -clobber ${SAMPLE_BRAIN} ${OUTDIR}/zero.mnc)
# 
# # create one
# ADD_FALCON_TEST(CREATE_ONE minccalc -express 1 -byte -clobber ${SAMPLE_BRAIN} ${OUTDIR}/one.mnc)
# 
# # extract GW interface
# ADD_FALCON_TEST(EXTRACT_GW ${FALCON_BIN}/niikmath thresh  -in=${SAMPLE_BRAIN} -thresh=250 -out=${OUTDIR}/GW.mnc )
# 
# # extract 'brain' 
# ADD_FALCON_TEST(EXTRACT_BRAIN ${FALCON_BIN}/niikmath thresh  -in=${SAMPLE_BRAIN} -thresh=20 -out=${OUTDIR}/BRAIN.mnc )
# 
# # close
# ADD_FALCON_TEST(CLOSE ${FALCON_BIN}/niikmath closebrain -in=${OUTDIR}/GW.mnc -out=${OUTDIR}/GW_close.mnc -radius=1 )
# 
# #add_dependencies(CLOSE EXTRACT_GW)
# 
# # dilate WM
# ADD_FALCON_TEST(DILATE ${FALCON_BIN}/niikmath dilate -in=${OUTDIR}/GW_close.mnc -out=${OUTDIR}/GW_dilate.mnc -radius=1)
# 
# # dilate brain
# # dilate
# ADD_FALCON_TEST(DILATE2 ${FALCON_BIN}/niikmath dilate -in=${OUTDIR}/BRAIN.mnc -out=${OUTDIR}/BRAIN_dilate.mnc -radius=1)
# 
# 
# # laplace
# ADD_FALCON_TEST(LAPLACE ${FALCON_BIN}/niikmath laplacemap -imglist=${OUTDIR}/GW.mnc,${OUTDIR}/GW_dilate.mnc -out=${OUTDIR}/GW_laplace.mnc -xyz=0.5,0.5,0.5 )
# 
# # shrink-wrap
# ADD_FALCON_TEST(SHRINK_WRAP ${FALCON_BIN}/niikmath shrinkwrap -mask=${OUTDIR}/GW_dilate.mnc -out=${OUTDIR}/GW_dilate_obj.off -iter=100 -radius=35 -val=10 -elen=1 )
# 
# # initics
# #ADD_FALCON_TEST(INITICS ${FALCON_BIN}/niikcortex_initics ${SAMPLE_BRAIN} ${OUTDIR}/GW.mnc ${OUTDIR}/GW_laplace.mnc ${OUTDIR}/GWI_init_ics.off ${fn}_GWI_mask_init_ics_${hemi}.off
# 
# # extract pial surface
# ADD_FALCON_TEST(EXTRACT_PIAL ${FALCON_BIN}/niikmath niikcortex-initocs -in=${SAMPLE_BRAIN} -imglist=${OUTDIR}/BRAIN_dilate.mnc,${OUTDIR}/zero.mnc,${OUTDIR}/GW.mnc 
#   -obj=${OUTDIR}/GW_dilate_obj.off -out=${OUTDIR}/GW_init_ocs.off -iter=100 )
# 
# # refine (deform) both WM and pial surfaces
# ADD_FALCON_TEST(REFINE ${FALCON_BIN}/niikcortex_deform ${SAMPLE_BRAIN} ${OUTDIR}/BRAIN_dilate.mnc ${OUTDIR}/zero.mnc ${OUTDIR}/GW.mnc ${OUTDIR}/zero.mnc ${OUTDIR}/zero.mnc 
#   ${OUTDIR}/GW_dilate_obj.off ${OUTDIR}/GW_init_ocs.off ${OUTDIR}/ics.off ${OUTDIR}/ocs.off -nonctx-mask ${OUTDIR}/zero.mnc)
# 
# # calculate thickness 
# ADD_FALCON_TEST(THICKNESS ${FALCON_BIN}/niikcortex_calc_thickness ${OUTDIR}/ics.off ${OUTDIR}/ocs.off ${OUTDIR}/thickness.txt )

#ADD_FALCON_TEST(SULCUS_PHANTOM ${CMAKE_CURRENT_SOURCE_DIR}/simple_cortex_thickness.sh ${FALCON_BIN} )

ADD_FALCON_TEST(BINARY_MORPHOLOGY ${TESTDIR}/test_morphology.sh ${FALCON_BIN} ${OUTDIR}/morphology_temp ) 

ADD_FALCON_TEST(SULCUS_PHANTOM    ${TESTDIR}/simple_cortex_thickness.sh ${FALCON_BIN} ${TESTDIR}/sulcus_2layers.mnc sulcus_phantom 250 90 4.0 ${OUTDIR}/sulcus_phantom_temp ) 

ADD_FALCON_TEST(SULCUS_PHANTOM_NOISY ${TESTDIR}/simple_cortex_thickness.sh ${FALCON_BIN} ${TESTDIR}/sulcus_2layers_noisy.mnc sulcus_phantom_noisy 250 90 4.0 ${OUTDIR}/sulcus_phantom_temp_noisy ) 

set_tests_properties(SULCUS_PHANTOM PROPERTIES TIMEOUT 1200)

set_tests_properties(SULCUS_PHANTOM_NOISY PROPERTIES TIMEOUT 1200)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/complex_phantom.mnc
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/make_complex_phantom.sh ${MINC_TOOLKIT_DIR} ${CMAKE_CURRENT_BINARY_DIR}/complex_phantom.mnc
  COMMENT "Generate complex phantom"
  VERBATIM
)

add_custom_command(
  OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/fg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bin_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/dist_phantom.mnc
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/make_simple_phantom.sh ${MINC_TOOLKIT_DIR} ${CMAKE_CURRENT_BINARY_DIR}/fg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bin_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/dist_phantom.mnc
  COMMENT "Generate simple phantom"
  VERBATIM
)

add_custom_target(make_phantoms ALL 
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/complex_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/fg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bin_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/dist_phantom.mnc
                  )

ADD_FALCON_TEST(COMPLEX_PHANTOM ${TESTDIR}/simple_cortex_thickness.sh 
    ${FALCON_BIN} ${CMAKE_CURRENT_BINARY_DIR}/complex_phantom.mnc complex_phantom 120 50 3.2 ${OUTDIR}/complex_phantom_temp 0.1 -complex ) 

set_tests_properties(COMPLEX_PHANTOM PROPERTIES TIMEOUT 1500)


ADD_FALCON_TEST(SURFACE_REFINE ${TESTDIR}/run_surface_refine.sh 
    ${MINC_TOOLKIT_DIR} ${FALCON_BIN} ${CMAKE_CURRENT_BINARY_DIR}/fg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bg_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/bin_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/dist_phantom.mnc ${CMAKE_CURRENT_BINARY_DIR}/surface_refine   )

#set_tests_properties(SURFACE_REFINE PROPERTIES TIMEOUT 1500)



if(HAVE_POVRAY)
  ADD_FALCON_TEST(COMPLEX_PHANTOM_QC 
                ${TESTDIR}/thickness_qc_test.sh ${FALCON_BIN} ${FALCON_SCRIPTS} ${CMAKE_CURRENT_BINARY_DIR}/complex_phantom ${OUTDIR}/complex_phantom_qc.png ) 

  SET_TESTS_PROPERTIES(COMPLEX_PHANTOM_QC PROPERTIES DEPENDS COMPLEX_PHANTOM)

endif(HAVE_POVRAY)

# 2. work on "real" brain

FUNCTION(DOWNLOAD_AND_UNPACK_ARCHIVE url md5 archive_file output )
  get_filename_component(archive_name "${archive_file}" NAME_WE)
  
  file(DOWNLOAD "${url}" "${archive_file}"  EXPECTED_MD5  "${md5}" SHOW_PROGRESS )

  add_custom_command( 
    OUTPUT "${output}"
    COMMAND cmake -E tar zxf ${archive_file} 
    DEPENDS ${archive_file}
    COMMENT "Unpacking ${archive_name}" VERBATIM
  )
  add_custom_target(unpack_${archive_name} ALL 
                  DEPENDS "${output}"
                  )
    
ENDFUNCTION(DOWNLOAD_AND_UNPACK_ARCHIVE)

OPTION(RUN_LONG_TESTS "Run long tests (~ 2 hrs)" OFF)

IF(RUN_LONG_TESTS)
message("Configuring long test")
#SET(CTEST_TEST_TIMEOUT

# need to download data from the web
#SET(REMOTE "http://www.bic.mni.mcgill.ca/~vfonov/temp/brainweb_sample.tar.gz")
#SET(LOCAL "file:///Users/vfonov/src/cache/brainweb_sample.tar.gz")
DOWNLOAD_AND_UNPACK_ARCHIVE( "http://www.bic.mni.mcgill.ca/~vfonov/temp/brainweb_sample.tar.gz"
                            "9458bd3158fa22b7a7b06b70090b2cc4"
                            "${CMAKE_CURRENT_BINARY_DIR}/brainweb_sample.tar.gz"
                            "${CMAKE_CURRENT_BINARY_DIR}/brainweb"
                            )

# run complete cortex thickness estimation
ADD_FALCON_TEST_PRL( BRAINWEB_THICKNESS_TEST ${TESTDIR}/run_falcon.sh 
                ${FALCON_BIN} ${MODELDIR} 
                ${CMAKE_CURRENT_BINARY_DIR}/brainweb/brainweb_t1_icbm152_2009c_normal_1mm_pn0_rf0 
                ${CMAKE_CURRENT_BINARY_DIR}/brainweb_measure 
                ${CMAKE_CURRENT_BINARY_DIR}/brainweb_work # work dir
                )
                
set_tests_properties(BRAINWEB_THICKNESS_TEST PROPERTIES TIMEOUT 7200)


                            

ENDIF(RUN_LONG_TESTS)

