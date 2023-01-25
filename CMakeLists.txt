cmake_minimum_required(VERSION 3.0)


if(POLICY CMP0069)
cmake_policy(SET CMP0069 NEW)
endif()

PROJECT(FALCON VERSION 0.9.18)

enable_language(C CXX)

SET(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

OPTION(USE_OPENMP  "Use OpenMP multi-threading" OFF)
OPTION(BUILD_NVK   "Build viewer, requires GTK2" OFF)
OPTION(HAVE_POVRAY "Is POVRAY available?" OFF)

#LIBIGL options
option(LIBIGL_USE_STATIC_LIBRARY      "Use libigl as static library" OFF)
option(LIBIGL_COPYLEFT_COMISO         "Use CoMiso"                   OFF)
option(LIBIGL_COPYLEFT_TETGEN         "Use Tetgen"                   OFF)
option(LIBIGL_COPYLEFT_CGAL           "Use CGAL"                     OFF)

option(LIBIGL_RESTRICTED_TRIANGLE     "Use Triangle"                 OFF)
option(LIBIGL_RESTRICTED_MOSEK        "Use MOSEK"                    OFF)
option(LIBIGL_RESTRICTED_MATLAB       "Use MATLAB"                   OFF)

option(LIBIGL_OPENGL                  "Use OpenGL"                   OFF)
option(LIBIGL_EMBREE                  "Use Embree"                   OFF)
option(LIBIGL_PNG                     "Use PNG"                      OFF)
option(LIBIGL_XML                     "Use XML"                      OFF)
option(LIBIGL_PREDICATES              "Use exact predicates"         OFF)
option(LIBIGL_GLFW                    "Use GLFW"                     OFF)
option(LIBIGL_IMGUI                   "Use ImGui"                    OFF)

if (USE_OPENMP)
    FIND_PACKAGE(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        set(HAVE_OPENMP  ON)
    else()
        message(WARNING "OpenMP NOT found")
    endif()
endif()

IF(${CMAKE_VERSION} VERSION_GREATER "3.9.4") 
  OPTION(USE_LTO "Use link time optimization" OFF)

  IF(USE_LTO)
    cmake_policy(SET CMP0069 NEW)
    # enable interprocedural and link-time optimization
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT error)

    IF(IPO_SUPPORTED)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    ENDIF(IPO_SUPPORTED)
  ENDIF(USE_LTO)
ENDIF()

# configuring for different systems
INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)

FIND_PATH(MINC_TOOLKIT_DIR names bin/mincresample HINTS $ENV{MINC_TOOLKIT} /opt/minc /opt/minc-itk4 )

FIND_PACKAGE(LIBMINC REQUIRED HINTS ${MINC_TOOLKIT_DIR}/lib )
FIND_PACKAGE(BICPL   REQUIRED HINTS ${MINC_TOOLKIT_DIR}/lib )
#FIND_PACKAGE(TIFF    QUIET )
FIND_PACKAGE(ZLIB    REQUIRED )

IF(BUILD_NVK)
    FIND_PACKAGE(GTK2 QUIET )
ENDIF(BUILD_NVK)

IF( NOT CMAKE_INSTALL_PREFIX )
   SET( CMAKE_INSTALL_PREFIX ${MINC_TOOLKIT_DIR} )
ENDIF()

#SET(HAVE_LIBTIFF OFF)
SET(HAVE_GTK2    OFF)
SET(HAVE_BICPL   OFF)

IF(HAVE_MINC1)
    ADD_DEFINITIONS(-DHAVE_MINC -DHAVE_MINC1)
ENDIF(HAVE_MINC1)

#IF(TIFF_FOUND)
#    SET(HAVE_LIBTIFF ON)
#ENDIF(TIFF_FOUND)

IF(ZLIB_FOUND)
    SET(HAVE_ZLIB ON)
ENDIF(ZLIB_FOUND)

IF(BUILD_NVK AND GTK2_FOUND)
    SET(HAVE_GTK2 ON)
ENDIF(BUILD_NVK AND GTK2_FOUND)

IF(BICPL_FOUND)
  SET(HAVE_BICPL ON)
  ADD_DEFINITIONS(-DHAVE_BICPL)
ENDIF(BICPL_FOUND)

SET(FALCON_DATA_SOURCE   ${CMAKE_SOURCE_DIR})
SET(FALCON_DATA_DEFAULT "${CMAKE_INSTALL_PREFIX}/share/falcon")

SET(NIIK_MAJOR_VERSION "0")
SET(NIIK_MINOR_VERSION "9")
SET(NIIK_MICRO_VERSION "16")

ENABLE_TESTING()
INCLUDE(CTest)

ADD_DEFINITIONS(-DHAVE_CONFIG_H)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(bin)

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(test)
ENDIF(BUILD_TESTING)

# install falcon data 
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/data
                    DESTINATION share/falcon )

# kate: space-indent on; hl cmake;indent-width 4; tab-width 4
