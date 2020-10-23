#ifndef IGL_READ_CSV2_H
#define IGL_READ_CSV2_H

//#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace igl 
{
  // read a matrix from a csv file into a Eigen matrix
  // Templates:
  //   Scalar  type for the matrix
  // Inputs:
  //   str  path to .csv file
  // Outputs:
  //   M  eigen matrix 
template <typename Derived>
 bool readCSV(
  const std::string & csv_file, 
  Eigen::PlainObjectBase<Derived> & M,
  std::vector<std::string> &header,
  bool skip_header=false
  );
}

#ifndef IGL_STATIC_LIBRARY
#  include "readCSV.cpp"
#endif

#endif
