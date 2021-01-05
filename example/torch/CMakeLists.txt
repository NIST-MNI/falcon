find_package(Torch QUIET)


IF(TORCH_FOUND)
message("Building TORCH example")

add_executable(minc_torch_segment minc_torch_segment.cpp)
set_property(TARGET minc_torch_segment PROPERTY CXX_STANDARD 11)
message("LIBMINC_LIBRARIES=${LIBMINC_LIBRARIES}")

target_link_libraries(minc_torch_segment ${TORCH_LIBRARIES} minc2-simple ${LIBMINC_LIBRARIES})


ENDIF(TORCH_FOUND)