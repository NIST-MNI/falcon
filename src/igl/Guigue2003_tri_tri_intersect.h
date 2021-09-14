#ifndef __TRI_TRI_INTERSECT_H__
#define __TRI_TRI_INTERSECT_H__

#include <igl/igl_inline.h>
#include <Eigen/Core>

namespace igl {

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
IGL_INLINE bool tri_tri_overlap_test_3d(
  const Eigen::MatrixBase<DerivedP1> &  p1, 
  const Eigen::MatrixBase<DerivedQ1> &  q1, 
  const Eigen::MatrixBase<DerivedR1> &  r1, 
  const Eigen::MatrixBase<DerivedP2> &  p2, 
  const Eigen::MatrixBase<DerivedQ2> &  q2, 
  const Eigen::MatrixBase<DerivedR2> &  r2);


// Three-dimensional Triangle-Triangle Overlap Test
// additionaly computes the segment of intersection of the two triangles if it exists. 
// coplanar returns whether the triangles are coplanar, 
// source and target are the endpoints of the line segment of intersection 
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedS,typename DerivedT>
IGL_INLINE bool tri_tri_intersection_test_3d(
    const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
    const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
    bool & coplanar, 
    Eigen::MatrixBase<DerivedS> & source, 
    Eigen::MatrixBase<DerivedT> & target );


template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1>
IGL_INLINE bool coplanar_tri_tri3d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &normal_1);


// Two dimensional Triangle-Triangle Overlap Test
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
IGL_INLINE bool tri_tri_overlap_test_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2);


};
#ifndef IGL_STATIC_LIBRARY
#  include "Guigue2003_tri_tri_intersect.cpp"
#endif

#endif //__TRI_TRI_INTERSECT_H__