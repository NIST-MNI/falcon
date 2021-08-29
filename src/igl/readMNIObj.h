#ifndef __READ_MNI_OBJ_H__
#define __READ_MNI_OBJ_H__

#include <igl/igl_inline.h>
#include <igl/FileEncoding.h>

#include <string>
#include <Eigen/Core>
namespace igl
{

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
  );
};


#ifndef IGL_STATIC_LIBRARY
#  include "readMNIObj.cpp"
#endif

#endif //__READ_MNI_OBJ_H__