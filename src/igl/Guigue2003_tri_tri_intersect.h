#ifndef __TRI_TRI_INTERSECT_H__
#define __TRI_TRI_INTERSECT_H__

#include <igl/igl_inline.h>

int tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3], 
          double p2[3], double q2[3], double r2[3]);


// Three-dimensional Triangle-Triangle Overlap Test
// additionaly computes the segment of intersection of the two triangles if it exists. 
// coplanar returns whether the triangles are coplanar, 
// source and target are the endpoints of the line segment of intersection 
int tri_tri_intersection_test_3d(double p1[3], double q1[3], double r1[3], 
								 double p2[3], double q2[3], double r2[3],
								 int * coplanar, 
								 double source[3],double target[3]);


int coplanar_tri_tri3d(double  p1[3], double  q1[3], double  r1[3],
           double  p2[3], double  q2[3], double  r2[3],
           double  N1[3]);


// Two dimensional Triangle-Triangle Overlap Test
int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
          double p2[2], double q2[2], double r2[2]);


#ifndef IGL_STATIC_LIBRARY
#  include "Guigue2003_tri_tri_intersect.cpp"
#endif

#endif //__TRI_TRI_INTERSECT_H__