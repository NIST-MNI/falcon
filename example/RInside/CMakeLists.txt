
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}")

FIND_PACKAGE(RInside)

IF(RINSIDE_FOUND)

message("Building RInside example")
message("RINSIDE_LIBS=${RINSIDE_LIBS}")
message("RINSIDE_LIB_DIR=${RINSIDE_LIB_DIR}")
message("RINSIDE_INCLUDE_DIR=${RINSIDE_INCLUDE_DIR}")


link_directories(${RINSIDE_LIB_DIR})
include_directories(${RINSIDE_INCLUDE_DIR})

add_executable(example_minc_rinside example_minc_rinside.cpp)
target_link_libraries(example_minc_rinside ${RINSIDE_LIBS} minc2-simple ${LIBMINC_LIBRARIES})

ENDIF(RINSIDE_FOUND)
