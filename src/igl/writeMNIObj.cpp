#pragma once
#ifndef __WRITE_MNIOBJ_CPP__
#define __WRITE_MNIOBJ_CPP__

#include <iostream>
#include <fstream>
#include  <bicpl.h>

namespace igl {

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedC
>
bool writeMNIObj(
  const std::string & fname,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedN> & N,
  const Eigen::MatrixBase<DerivedC> & C,
  FileEncoding encoding
   )
{
  object_struct   *object_list;
  polygons_struct      *polygons;

  object_list=create_object(POLYGONS);
  polygons = get_polygons_ptr(object_list);

  VIO_Colour single_colour=make_rgba_Colour(200,100,100,255); /*reddish?*/
  if(C.rows()==1)
    single_colour=make_rgba_Colour(C(0,0)*255,C(0,1)*255,C(0,2)*255,255);

  initialize_polygons_with_size(polygons,
                                single_colour,  
                                NULL,
                                V.rows(),F.rows(),3);

  for(int i=0;i<V.rows();++i) {
    /*assign XYZ coords*/
    fill_Point(polygons->points[i], V(i,0),V(i,1),V(i,2));
  }
  /*assign normals?*/
  if(N.rows()==V.rows()) {
    for(int i=0;i<N.rows();++i) {
      fill_Point(polygons->normals[i],N(i,0),N(i,1),N(i,2));
    }
  } //TODO: display warning?

  /*assign indeces*/
  for(int i=0;i<F.rows();++i) {
    polygons->indices[i*3]  =F(i,0);
    polygons->indices[i*3+1]=F(i,1);
    polygons->indices[i*3+2]=F(i,2);
  }

  /*assume that end_indeces are properly set already*/

  if(C.rows()==F.rows()) {
    /*assign colours*/
    set_polygon_per_item_colours(polygons);
    for(int i=0;F.rows();++i)
      polygons->colours[i]=make_rgba_Colour(F(i,0)*255,F(i,1)*255,F(i,2)*255,255);
  }

  bool success= ( output_graphics_file( (char*) fname.c_str(),encoding==FileEncoding::Binary?BINARY_FORMAT:ASCII_FORMAT, 1, &object_list ) == VIO_OK);
  delete_object(object_list);

  return success;
}

}; //igl

#endif
