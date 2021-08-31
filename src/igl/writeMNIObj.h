#ifndef __WRITE_MNI_OBJ_H__
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
bool writeMNIObj(
  const std::string & fname,
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedN> & N,
  const Eigen::MatrixBase<DerivedC> & C,
  FileEncoding encoding
   );

};

#ifndef IGL_STATIC_LIBRARY
#  include "writeMNIObj.cpp"
#endif

#endif //