#ifndef IGL_WRITE_CSV2_H
#define IGL_WRITE_CSV2_H

//#include "igl/igl_inline.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace igl 
{
  // write a matrix to a csv file 
  // Templates:
  //   Scalar  type for the matrix
  // Inputs:
  //   str  path to .csv file
  //   M  eigen matrix     
  //   header  column names
template <typename Derived>
 bool writeCSV(
  const std::string &csv_file, 
  const Eigen::PlainObjectBase<Derived> & M,
  const std::vector<std::string> &header={""}
  );

  bool check_ext(const std::string &path, const std::string &ext);

}

bool inline igl::check_ext(const std::string &path, const std::string &ext)
{
  auto idx=path.rfind('.');
  if(idx != std::string::npos)
  {
      std::string _extension = path.substr(idx+1);
      return _extension==ext;
  }
  return false;
}


#ifndef IGL_STATIC_LIBRARY
#  include "writeCSV.cpp" 
#endif

#endif
