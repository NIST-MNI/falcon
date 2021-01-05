# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.8)
PROJECT(minc2-simple C)

FIND_PATH( MINC_TOOLKIT_DIR "lib${LIB_SUFFIX}/LIBMINCConfig.cmake" HINTS "$ENV{MINC_TOOLKIT}" "/opt/minc" )
FIND_PACKAGE( LIBMINC REQUIRED HINTS ${MINC_TOOLKIT_DIR}/lib${LIB_SUFFIX})

INCLUDE( ${LIBMINC_USE_FILE} )

IF(LUAROCKS_PREFIX)
    MESSAGE(STATUS "Installing Torch through Luarocks")
    STRING(REGEX REPLACE "(.*)lib/luarocks/rocks.*" "\\1" CMAKE_INSTALL_PREFIX  "${LUAROCKS_PREFIX}")
    MESSAGE(STATUS "Prefix inferred from Luarocks: ${CMAKE_INSTALL_PREFIX}")
ENDIF()
FIND_PACKAGE(Torch REQUIRED)


MESSAGE(STATUS "Prefix : ${CMAKE_INSTALL_PREFIX}")

INCLUDE_DIRECTORIES(../src)

#FILE(GLOB luasrc *.lua)
SET(SRC
 ../src/minc2-simple.c
 ../src/minc2-matrix-ops.c
)
#SET(CMAKE_C_FLAGS "-fPIC ${CMAKE_C_FLAGS}")

ADD_LIBRARY(minc2-simple SHARED ${SRC})
target_link_libraries(minc2-simple ${LIBMINC_LIBRARIES})
install(TARGETS minc2-simple DESTINATION lib)

SET(LUASRC 
 init.lua
)

ADD_TORCH_PACKAGE(minc2_simple "" "${LUASRC}" "minc2-simple Package")
#kate: indent-width 2; replace-tabs on; 
