if(LIBIGL_USE_STATIC_LIBRARY)
 set(SOURCES_IGL_CSV rapidcsv.h  readCSV.cpp  readCSV.h  writeCSV.cpp  writeCSV.h)
 add_library(igl_csv STATIC ${SOURCES_IGL_CSV})
else()
 add_library(igl_csv INTERFACE)   
endif()


target_link_libraries(igl_csv ${IGL_SCOPE} igl_common)
target_link_libraries(igl_csv ${IGL_SCOPE} igl)

add_library(igl::csv ALIAS igl_csv)
set_property(TARGET igl_csv PROPERTY EXPORT_NAME igl::csv)
