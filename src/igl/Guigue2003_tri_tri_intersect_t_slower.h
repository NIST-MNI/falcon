
/*
*         
*  Triangle-Triangle Overlap Test Routines        
*  July, 2002                                                          
*  Updated December 2003                                                
*                                                                       
*  This file contains C implementation of algorithms for                
*  performing two and three-dimensional triangle-triangle intersection test 
*  The algorithms and underlying theory are described in                    
*                                                                           
* "Fast and Robust Triangle-Triangle Overlap Test 
*  Using Orientation Predicates"  P. Guigue - O. Devillers
*                                                 
*  Journal of Graphics Tools, 8(1), 2003                                    
*                                                                           
*  Several geometric predicates are defined.  Their parameters are all      
*  points.  Each point is an array of two or three double precision         
*  floating point numbers. The geometric predicates implemented in          
*  this file are:                                                            
*                                                                           
*    int _tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)                         
*    int _tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)                         
*                                                                           
*    int _tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
*                                     coplanar,source,target)               
*                                                                           
*       is a version that computes the segment of intersection when            
*       the triangles overlap (and are not coplanar)                        
*                                                                           
*    each function returns 1 if the triangles (including their              
*    boundary) intersect, otherwise 0                                       
*                                                                           
*                                                                           
*  Other information are available from the Web page                        
*  http://www.acm.org/jgt/papers/GuigueDevillers03/                         
*                                                                           
*/

#ifndef IGL_GUIGUE2003_TRI_TRI_INTERSECT_T_H
#define IGL_GUIGUE2003_TRI_TRI_INTERSECT_T_H

#include <Eigen/Core>

/* function prototype */

// Three-dimensional Triangle-Triangle Overlap Test
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
inline bool _tri_tri_overlap_test_3d(
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
inline bool _tri_tri_intersection_test_3d(
    const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
    const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
    bool & coplanar, 
    Eigen::MatrixBase<DerivedS> & source, 
    Eigen::MatrixBase<DerivedT> & target );


template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1>
inline bool _coplanar_tri_tri3d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &normal_1);


// Two dimensional Triangle-Triangle Overlap Test
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
inline bool _tri_tri_overlap_test_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2);




/* some 3D macros */

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
inline bool igl_check_min_max(
  const Eigen::MatrixBase<DerivedP1> &  p1, 
  const Eigen::MatrixBase<DerivedQ1> &  q1, 
  const Eigen::MatrixBase<DerivedR1> &  r1, 
  const Eigen::MatrixBase<DerivedP2> &  p2, 
  const Eigen::MatrixBase<DerivedQ2> &  q2, 
  const Eigen::MatrixBase<DerivedR2> &  r2)
{ 
  using Scalar=typename DerivedP1::Scalar;

  auto v1=p2-q1; 
  /*IGL_SUB(v1,p2,q1)*/ 
  auto v2=p1-q1;
  /*IGL_SUB(v2,p1,q1)*/ 
  auto N1=v1.cross(v2); 
  /*IGL_CROSS(N1,v1,v2)*/ 
  v1=q2-q1; 
  /*IGL_SUB(v1,q2,q1)*/
  if(v1.dot(N1) > Scalar(0.0)) return false; 
  /*if (IGL_DOT(v1,N1) > 0.0f) return false;*/
  v1=p2-p1; 
  /*IGL_SUB(v1,p2,p1)*/
  v2=r1-p1; 
  /*IGL_SUB(v2,r1,p1)*/
  auto N2=v1.cross(v2); 
  /*IGL_CROSS(N1,v1,v2)*/ 
  v1=r2-p1; 
  /*IGL_SUB(v1,r2,p1)*/ 
  if(v1.dot(N2)>Scalar(0.0)) return false; 
  /*if (IGL_DOT(v1,N1) > 0.0f) return false;*/ 
  else return true; }



/* Permutation in a canonical form of T2's vertices */

#define IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
  if (dp2 > 0.0f) { \
     if (dq2 > 0.0f) return igl_check_min_max(p1,r1,q1,r2,p2,q2); \
     else if (dr2 > 0.0f) return igl_check_min_max(p1,r1,q1,q2,r2,p2);\
     else return igl_check_min_max(p1,q1,r1,p2,q2,r2); }\
  else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) return igl_check_min_max(p1,q1,r1,r2,p2,q2);\
    else if (dr2 < 0.0f) return igl_check_min_max(p1,q1,r1,q2,r2,p2);\
    else return igl_check_min_max(p1,r1,q1,p2,q2,r2);\
  } else { \
    if (dq2 < 0.0f) { \
      if (dr2 >= 0.0f)  return igl_check_min_max(p1,r1,q1,q2,r2,p2);\
      else return igl_check_min_max(p1,q1,r1,p2,q2,r2);\
    } \
    else if (dq2 > 0.0f) { \
      if (dr2 > 0.0f) return igl_check_min_max(p1,r1,q1,p2,q2,r2);\
      else  return igl_check_min_max(p1,q1,r1,q2,r2,p2);\
    } \
    else  { \
      if (dr2 > 0.0f) return igl_check_min_max(p1,q1,r1,r2,p2,q2);\
      else if (dr2 < 0.0f) return igl_check_min_max(p1,r1,q1,r2,p2,q2);\
      else return _coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);\
     }}}
  


/*
*
*  Three-dimensional Triangle-Triangle Overlap Test
*
*/

template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2> 
inline bool _tri_tri_overlap_test_3d(
  const Eigen::MatrixBase<DerivedP1> &  p1, 
  const Eigen::MatrixBase<DerivedQ1> &  q1, 
  const Eigen::MatrixBase<DerivedR1> &  r1, 
  const Eigen::MatrixBase<DerivedP2> &  p2, 
  const Eigen::MatrixBase<DerivedQ2> &  q2, 
  const Eigen::MatrixBase<DerivedR2> &  r2)
{
  using Scalar    = typename DerivedP1::Scalar;
  using RowVector = typename Eigen::Matrix<Scalar,1,3>;

  Scalar dp1, dq1, dr1, dp2, dq2, dr2;
  RowVector v1, v2;
  RowVector N1, N2; 
  
  /* Compute distance signs  of p1, q1 and r1 to the plane of
     triangle(p2,q2,r2) */

  v1=p2-r2;
  //IGL_SUB(v1,p2,r2)
  v2=q2-r2;
  //IGL_SUB(v2,q2,r2)

  N2=v1.cross(v2);
  //IGL_CROSS(N2,v1,v2)
  
  v1=p1-r2;
  //IGL_SUB(v1,p1,r2)
  dp1 = v1.dot(N2);
  //IGL_DOT(v1,N2);
  v1=q1-r2;
  //IGL_SUB(v1,q1,r2)
  dq1=v1.dot(N2);
  //dq1 = IGL_DOT(v1,N2);
  v1=r1-r2;
  //IGL_SUB(v1,r1,r2)
  dr1=v1.dot(N2);
  //dr1 = IGL_DOT(v1,N2);
  
  if (((dp1 * dq1) > Scalar(0.0) ) && ((dp1 * dr1) > Scalar(0.0)))  return false; 

  /* Compute distance signs  of p2, q2 and r2 to the plane of
     triangle(p1,q1,r1) */
  v1=q1-p1;
  //IGL_SUB(v1,q1,p1)
  v2=r1-p1;
  //IGL_SUB(v2,r1,p1)
  N1=v1.cross(v2);
  //IGL_CROSS(N1,v1,v2)
  v1=p2-r1;
  //IGL_SUB(v1,p2,r1)
  dp2 = v1.dot(N1);
  //dp2 = IGL_DOT(v1,N1);
  v1=q2-r1;
  //IGL_SUB(v1,q2,r1)
  dq2=v1.dot(N1);
  //dq2 = IGL_DOT(v1,N1);
  v1=r2-r1;
  //IGL_SUB(v1,r2,r1)
  dr2=v1.dot(N1);
  //dr2 = IGL_DOT(v1,N1);
  
  if (((dp2 * dq2) > Scalar(0.0)) && ((dp2 * dr2) > Scalar(0.0))) return false;

  /* Permutation in a canonical form of T1's vertices */


  if (dp1 > 0.0f) {
    if (dq1 > 0.0f) IGL_TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
    else if (dr1 > 0.0f) IGL_TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)  
    else IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
  } else if (dp1 < 0.0f) {
    if (dq1 < 0.0f) IGL_TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
    else if (dr1 < 0.0f) IGL_TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    else IGL_TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
  } else {
    if (dq1 < 0.0f) {
      if (dr1 >= 0.0f) IGL_TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
      else IGL_TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
    }
    else if (dq1 > 0.0f) {
      if (dr1 > 0.0f) IGL_TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
      else IGL_TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
    }
    else  {
      if (dr1 > 0.0f) IGL_TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
      else if (dr1 < 0.0f) IGL_TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
      else return _coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);
    }
  }
};


template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1>
inline bool _coplanar_tri_tri3d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &normal_1)
{

  using Scalar= typename DerivedP1::Scalar;
  using RowVector2D = typename Eigen::Matrix<Scalar,1,2>;

  RowVector2D P1,Q1,R1;
  RowVector2D P2,Q2,R2;

  Scalar n_x, n_y, n_z;

  n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
  n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
  n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


  /* Projection of the triangles in 3D onto 2D such that the area of
     the projection is maximized. */


  if (( n_x > n_z ) && ( n_x >= n_y )) {
    // Project onto plane YZ

      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1]; 
    
      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1]; 

  } else if (( n_y > n_z ) && ( n_y >= n_x )) {
    // Project onto plane XZ

    P1[0] = q1[0]; P1[1] = q1[2];
    Q1[0] = p1[0]; Q1[1] = p1[2];
    R1[0] = r1[0]; R1[1] = r1[2]; 
 
    P2[0] = q2[0]; P2[1] = q2[2];
    Q2[0] = p2[0]; Q2[1] = p2[2];
    R2[0] = r2[0]; R2[1] = r2[2]; 
    
  } else {
    // Project onto plane XY

    P1[0] = p1[0]; P1[1] = p1[1]; 
    Q1[0] = q1[0]; Q1[1] = q1[1]; 
    R1[0] = r1[0]; R1[1] = r1[1]; 
    
    P2[0] = p2[0]; P2[1] = p2[1]; 
    Q2[0] = q2[0]; Q2[1] = q2[1]; 
    R2[0] = r2[0]; R2[1] = r2[1]; 
  }

  return _tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
    
};



/*
*                                                                
*  Three-dimensional Triangle-Triangle Intersection              
*
*/

/*
   This macro is called when the triangles surely intersect
   It constructs the segment of intersection of the two triangles
   if they are not coplanar.
*/

// NOTE: a faster, but possibly less precise, method of computing
// point B is described here: https://github.com/erich666/jgt-code/issues/5
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1,typename DerivedN2,
typename DervideS,typename DerivedT>
inline bool igl_construct_intersection(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2,
  const Eigen::MatrixBase<DerivedN1> &N1, const Eigen::MatrixBase<DerivedN2> &N2,
        Eigen::MatrixBase<DervideS> & source,  Eigen::MatrixBase<DerivedT> & target
  )
{ 
  using Scalar= typename DerivedP1::Scalar;
  using RowVector3D = typename Eigen::Matrix<Scalar,1,3>;

  Scalar alpha;

  RowVector3D v1=q1-p1; 
  RowVector3D v2=r2-p1; 
  RowVector3D N=v1.cross(v2); 
  RowVector3D v=p2-p1;

  if (v.dot(N) > 0.0f) {
    v1=r1-p1; 
    N=v1.cross(v2); 
    if (v.dot(N) <= 0.0f) { 
      v2=q2-p1; 
      N=v1.cross(v2); 
      if (v.dot(N) > 0.0f) { 
        v1=p1-p2; 
        v2=p1-r1; 
        alpha = v1.dot(N2) / v2.dot(N2); 
        v1=alpha*v2; 
        source=p1-v1; 
        v1=p2-p1; 
        v2=p2-r2; 
        alpha = v1.dot(N1) / v2.dot(N1); 
        v1=alpha*v2; 
        target=p2-v1; 
        return true; 
      } else { 
        v1=p2-p1; 
        v2=p2-q2; 
        alpha = v1.dot(N1) / v2.dot(N1); 
        v1=alpha*v2; 
        source=p2-v1; 
        v1=p2-p1; 
        v2=p2-r2; 
        alpha = v1.dot(N1) / v2.dot(N1); 
        v1=alpha*v2; 
        target=p2-v1; 
        return true; 
      } 
    } else { 
      return false; 
    } 
  } else { 
    v2=q2-p1; 
    N=v1.cross(v2); 
    if (v.dot(N) < 0.0f) { 
      return false; 
    } else { 
      v1=r1-p1; 
      N=v1.cross(v2); 
      if (v.dot(N) >= 0.0f) { 
        v1=p1-p2; 
        v2=p1-r1; 
        alpha = v1.dot(N2) / v2.dot(N2); 
        v1=alpha*v2; 
        source=p1-v1; 
        v1=p1-p2; 
        v2=p1-q1; 
        alpha = v1.dot(N2) / v2.dot(N2); 
        v1=alpha*v2; 
        target=p1-v1; 
        return true; 
      } else { 
        v1=p2-p1; 
        v2=p2-q2; 
        alpha = v1.dot(N1) / v2.dot(N1); 
        v1=alpha*v2; 
        source=p2-v1; 
        v1=p1-p2; 
        v2=p1-q1; 
        alpha = v1.dot(N2) / v2.dot(N2); 
        v1=alpha*v2; 
        target=p1-v1; 
        return true; 
      }
    }
  }
} 

                
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedN1,typename DerivedN2,
typename DerivedS,typename DerivedT>
inline bool igl_tri_tri_inter_3d(
  const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
  const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
  const Eigen::MatrixBase<DerivedN1> & N1, const Eigen::MatrixBase<DerivedN2> & N2,
  typename  DerivedP1::Scalar dp2, typename  DerivedP1::Scalar dq2, typename  DerivedP1::Scalar dr2,
  Eigen::MatrixBase<DerivedS> &source, Eigen::MatrixBase<DerivedT> &target,
  bool &coplanar) {
  coplanar=false;
  if (dp2 > 0.0f) { 
     if (dq2 > 0.0f) return igl_construct_intersection(p1,r1,q1,r2,p2,q2,N1,N2,source,target); 
     else if (dr2 > 0.0f) return igl_construct_intersection(p1,r1,q1,q2,r2,p2,N1,N2,source,target);
     else return igl_construct_intersection(p1,q1,r1,p2,q2,r2,N1,N2,source,target); }
  else if (dp2 < 0.0f) { 
    if (dq2 < 0.0f) return igl_construct_intersection(p1,q1,r1,r2,p2,q2,N1,N2,source,target);
    else if (dr2 < 0.0f) return igl_construct_intersection(p1,q1,r1,q2,r2,p2,N1,N2,source,target);
    else return igl_construct_intersection(p1,r1,q1,p2,q2,r2,N1,N2,source,target);
  } else { 
    if (dq2 < 0.0f) { 
      if (dr2 >= 0.0f)  return igl_construct_intersection(p1,r1,q1,q2,r2,p2,N1,N2,source,target);
      else return igl_construct_intersection(p1,q1,r1,p2,q2,r2,N1,N2,source,target);
    } 
    else if (dq2 > 0.0f) { 
      if (dr2 > 0.0f) return igl_construct_intersection(p1,r1,q1,p2,q2,r2,N1,N2,source,target);
      else  return igl_construct_intersection(p1,q1,r1,q2,r2,p2,N1,N2,source,target);
    } 
    else  { 
      if (dr2 > 0.0f) return igl_construct_intersection(p1,q1,r1,r2,p2,q2,N1,N2,source,target);
      else if (dr2 < 0.0f) return igl_construct_intersection(p1,r1,q1,r2,p2,q2,N1,N2,source,target);
      else { 
        coplanar = true; 
        return _coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);
     } 
  }} }
  

/*
   The following version computes the segment of intersection of the
   two triangles if it exists. 
   coplanar returns whether the triangles are coplanar
   source and target are the endpoints of the line segment of intersection 
*/
template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2,
typename DerivedS,typename DerivedT>
inline bool _tri_tri_intersection_test_3d(
    const Eigen::MatrixBase<DerivedP1> & p1, const Eigen::MatrixBase<DerivedQ1> & q1, const Eigen::MatrixBase<DerivedR1> & r1, 
    const Eigen::MatrixBase<DerivedP2> & p2, const Eigen::MatrixBase<DerivedQ2> & q2, const Eigen::MatrixBase<DerivedR2> & r2,
    bool & coplanar, 
    Eigen::MatrixBase<DerivedS> & source, 
    Eigen::MatrixBase<DerivedT> & target )        
{
  using Scalar= typename DerivedP1::Scalar;
  using RowVector3D = typename Eigen::Matrix<Scalar,1,3>;

  Scalar dp1, dq1, dr1, dp2, dq2, dr2;
  RowVector3D v1, v2, v;
  RowVector3D N1, N2, N;
  Scalar alpha;

  // Compute distance signs  of p1, q1 and r1 
  // to the plane of triangle(p2,q2,r2)


  v1=p2-r2;
  v2=q2-r2;
  N2=v1.cross(v2);

  v1=p1-r2;
  dp1 = v1.dot(N2);
  v1=q1-r2;
  dq1 = v1.dot(N2);
  v1=r1-r2;
  dr1 = v1.dot(N2);
  
  if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > Scalar(0.0)))  return false; 

  // Compute distance signs  of p2, q2 and r2 
  // to the plane of triangle(p1,q1,r1)

  
  v1=q1-p1;
  v2=r1-p1;
  N1=v1.cross(v2);

  v1=p2-r1;
  dp2 = v1.dot(N1);
  v1=q2-r1;
  dq2 = v1.dot(N1);
  v1=r2-r1;
  dr2 = v1.dot(N1);
  
  if (((dp2 * dq2) > Scalar(0.0)) && ((dp2 * dr2) > Scalar(0.0))) return false;

  // Permutation in a canonical form of T1's vertices

  if (dp1 > Scalar(0.0)) {
    if (dq1 > Scalar(0.0)) return igl_tri_tri_inter_3d(r1,p1,q1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
    else if (dr1 > Scalar(0.0)) return igl_tri_tri_inter_3d(q1,r1,p1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
  
    else return igl_tri_tri_inter_3d(p1,q1,r1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
  } else if (dp1 < Scalar(0.0)) {
    if (dq1 < Scalar(0.0)) return igl_tri_tri_inter_3d(r1,p1,q1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
    else if (dr1 < Scalar(0.0)) return igl_tri_tri_inter_3d(q1,r1,p1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
    else return igl_tri_tri_inter_3d(p1,q1,r1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
  } else {
    if (dq1 < Scalar(0.0)) {
      if (dr1 >= Scalar(0.0)) return igl_tri_tri_inter_3d(q1,r1,p1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
      else return igl_tri_tri_inter_3d(p1,q1,r1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
    }
    else if (dq1 > Scalar(0.0)) {
      if (dr1 > Scalar(0.0)) return igl_tri_tri_inter_3d(p1,q1,r1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
      else return igl_tri_tri_inter_3d(q1,r1,p1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
    }
    else  {
      if (dr1 > Scalar(0.0)) return igl_tri_tri_inter_3d(r1,p1,q1,p2,q2,r2,N1,N2,dp2,dq2,dr2,source,target,coplanar);
      else if (dr1 < Scalar(0.0)) return igl_tri_tri_inter_3d(r1,p1,q1,p2,r2,q2,N1,N2,dp2,dr2,dq2,source,target,coplanar);
      else {
        // triangles are co-planar
        coplanar = true;
        return _coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1);
      }
    }
  }
};





/*
*
*  Two dimensional Triangle-Triangle Overlap Test    
*
*/


/* some 2D macros */

#define IGL_ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


#define IGL_INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (IGL_ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (IGL_ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (IGL_ORIENT_2D(P1,P2,Q1) > 0.0f) {\
  if (IGL_ORIENT_2D(P1,Q2,Q1) <= 0.0f) return true; \
  else return false;} else {\
  if (IGL_ORIENT_2D(P1,P2,R1) >= 0.0f)\
    if (IGL_ORIENT_2D(Q1,R1,P2) >= 0.0f) return true; \
    else return false;\
  else return false;}\
    else \
      if (IGL_ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
  if (IGL_ORIENT_2D(R2,Q2,R1) <= 0.0f)\
    if (IGL_ORIENT_2D(Q1,R1,Q2) >= 0.0f) return true; \
    else return false;\
  else return false;\
      else return false;\
  else\
    if (IGL_ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (IGL_ORIENT_2D(Q1,R1,R2) >= 0.0f)\
  if (IGL_ORIENT_2D(P1,P2,R1) >= 0.0f) return true;\
  else return false;\
      else \
  if (IGL_ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
    if (IGL_ORIENT_2D(R2,R1,Q2) >= 0.0f) return true; \
    else return false; }\
  else return false; \
    else  return false; \
 };



#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (IGL_ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (IGL_ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (IGL_ORIENT_2D(P1,Q1,R2) >= 0.0f) return true; \
        else return false;} else { \
      if (IGL_ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
  if (IGL_ORIENT_2D(R1,P1,P2) >= 0.0f) return true; else return false;} \
      else return false; } \
  } else {\
    if (IGL_ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (IGL_ORIENT_2D(P1,P2,R1) >= 0.0f) {\
  if (IGL_ORIENT_2D(P1,R1,R2) >= 0.0f) return true;  \
  else {\
    if (IGL_ORIENT_2D(Q1,R1,R2) >= 0.0f) return true; else return false;}}\
      else  return false; }\
    else return false; }}


template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
inline bool ccw_tri_tri_intersection_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2) 
{
  using Scalar= typename DerivedP1::Scalar;
  using RowVector3D = typename Eigen::Matrix<Scalar,1,2>;

  if ( IGL_ORIENT_2D(p2,q2,p1) >= 0.0 ) {
    if ( IGL_ORIENT_2D(q2,r2,p1) >= 0.0 ) {
      if ( IGL_ORIENT_2D(r2,p2,p1) >= 0.0 ) return true;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( IGL_ORIENT_2D(r2,p2,p1) >= 0.0 ) 
        INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else IGL_INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( IGL_ORIENT_2D(q2,r2,p1) >= 0.0 ) {
      if ( IGL_ORIENT_2D(r2,p2,p1) >= 0.0 ) 
        INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  IGL_INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else IGL_INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
};



template <typename DerivedP1,typename DerivedQ1,typename DerivedR1,
typename DerivedP2,typename DerivedQ2,typename DerivedR2>
inline bool _tri_tri_overlap_test_2d(
  const Eigen::MatrixBase<DerivedP1> &p1, const Eigen::MatrixBase<DerivedQ1> &q1, const Eigen::MatrixBase<DerivedR1> &r1,
  const Eigen::MatrixBase<DerivedP2> &p2, const Eigen::MatrixBase<DerivedQ2> &q2, const Eigen::MatrixBase<DerivedR2> &r2) 
{
  if ( IGL_ORIENT_2D(p1,q1,r1) < 0.0 )
    if ( IGL_ORIENT_2D(p2,q2,r2) < 0.0 )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( IGL_ORIENT_2D(p2,q2,r2) < 0.0 )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

};
#endif 
// IGL_GUIGUE2003_TRI_TRI_INTERSECT_T_H