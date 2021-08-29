#pragma once
#ifndef __READ_MNIOBJ_CPP__
#define __READ_MNIOBJ_CPP__

#include <iostream>
#include <fstream>
#include <bicpl.h>

namespace igl {

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedC
  >
inline bool readMNIObj(
  const std::string & fname,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedN> & N,
  Eigen::PlainObjectBase<DerivedC> & C
  )
{

  VIO_File_formats         format;
  object_struct        **object_list;
  int n_objects=0;
  int n_points=0;
  int n_triangles=0;

  polygons_struct      *polygons;

  if( input_graphics_file( fname.c_str(), &format, &n_objects,
                           &object_list ) != VIO_OK ) {
    std::cerr<<"Can't open file"<<std::endl;
    return false;
  }

  if( n_objects != 1) {
    std::cerr<< "File "<< fname.c_str() << 
        " contains "<< n_objects << " objects, only one is supported presently"<<std::endl;
    return false;
  }

  if(get_object_type( object_list[0] ) != POLYGONS ) {
    std::cerr<< "File "<<fname.c_str() << 
        " contains unsupported object, only polygons are supported presently"<<std::endl;
    return false;
  }

  polygons = get_polygons_ptr(object_list[0]);

  n_points = polygons->n_points;

  V.resize(n_points,3);

  for(int i=0; i<n_points; i++) {
    V(i,0)=polygons->points[i].coords[0];
    V(i,1)=polygons->points[i].coords[1];
    V(i,2)=polygons->points[i].coords[2];
  }

  N.resize(n_points,3);

  for(int i=0; i<n_points; i++) {
    N(i,0)=polygons->normals[i].coords[0];
    N(i,1)=polygons->normals[i].coords[1];
    N(i,2)=polygons->normals[i].coords[2];
  }

  n_triangles = polygons->n_items;

  if(polygons->colour_flag != ONE_COLOUR) {
    C.resize(n_triangles,3);

    for(int i=0; i<n_triangles; i++) {

      C(i,0)=get_Colour_r_0_1(polygons->colours[i]);
      C(i,1)=get_Colour_g_0_1(polygons->colours[i]);
      C(i,2)=get_Colour_b_0_1(polygons->colours[i]);
    }
  }

  F.resize(n_triangles,3);
  /* read/write triangles */
  for(int i=0,idx=0; i<n_triangles; i++) {
    /*convert at most 3 points per face.... hack?*/
    for(int j=idx,k=0; j<polygons->end_indices[i] && k<3; j++,k++) {
      F(i,k)=polygons->indices[j];
    }
    idx = polygons->end_indices[i];
  }

  delete_object_list( n_objects, object_list );

  return true; 
} 

};//igl

#endif 
