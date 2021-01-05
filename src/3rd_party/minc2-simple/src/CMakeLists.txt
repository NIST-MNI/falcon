ADD_LIBRARY(minc2-simple SHARED
            minc2-simple.c
            minc2-matrix-ops.c
            minc2-simple.h
            minc2-simple-int.h
            )

ADD_LIBRARY(minc2-simple-static STATIC
             minc2-simple.c
             minc2-matrix-ops.c
             minc2-simple.h
             minc2-simple-int.h
           )
            
TARGET_LINK_LIBRARIES(minc2-simple ${LIBMINC_LIBRARIES})
TARGET_LINK_LIBRARIES(minc2-simple-static ${LIBMINC_LIBRARIES})

set_target_properties(minc2-simple minc2-simple-static PROPERTIES MACOSX_RPATH ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
set_target_properties(minc2-simple minc2-simple-static PROPERTIES PUBLIC_HEADER "minc2-simple.h;minc2-simple-int.h")

INSTALL(TARGETS minc2-simple minc2-simple-static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
