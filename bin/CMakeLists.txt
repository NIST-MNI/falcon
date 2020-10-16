SET(FALCON_HOME_DEFAULT  "${CMAKE_INSTALL_PREFIX}")

configure_file(falcon_run.sh
  ${CMAKE_CURRENT_BINARY_DIR}/falcon_run.sh @ONLY )

# alias for compatibility
configure_file(falcon_run_v2.sh
  ${CMAKE_CURRENT_BINARY_DIR}/falcon_run_v2.sh @ONLY )

INSTALL(PROGRAMS
          ${CMAKE_CURRENT_BINARY_DIR}/falcon_run.sh
          ${CMAKE_CURRENT_BINARY_DIR}/falcon_run_v2.sh
          falcon_slice_qc.sh
          falcon_off_qc.sh
          falcon_off_qc_2.sh
          falcon_antsAtroposN4.sh
          falcon_antsBrainExtraction.sh
          falcon_face_render.sh
          falcon_make_face.sh
        DESTINATION
          bin )
