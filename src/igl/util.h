#ifndef __UTIL_H__
#define __UTIL_H__


#include <vector>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

inline bool extract_psi_the(const std::vector<std::string> &header, 
                     const Eigen::MatrixXd& D, 
                           Eigen::MatrixXd& out_psi_the )
{
  auto idx_psi=std::find(header.begin(),header.end(),"psi");
  auto idx_the=std::find(header.begin(),header.end(),"the");
  if(idx_psi!=header.end() && idx_the!=header.end())
  {
    size_t _idx_psi=idx_psi-header.begin();
    size_t _idx_the=idx_the-header.begin();

    out_psi_the.resize(D.rows(),2);

    out_psi_the<<D.col(_idx_psi),D.col(_idx_the);
    return true;
  } else {
    return false;
  }
}

inline void sph_to_xyz(const Eigen::MatrixXd& sph, 
                      Eigen::MatrixXd& xyz )
{

  xyz.resize(sph.rows(), 3);

  xyz << sph.col(1).array().sin() * sph.col(0).array().cos(), 
         sph.col(1).array().sin() * sph.col(0).array().sin(),
         sph.col(1).array().cos() ;

}

inline void xyz_to_sph(const Eigen::MatrixXd& xyz, 
                             Eigen::MatrixXd& sph)
{
  sph.resize(sph.rows(), 2);
  sph << 
    Eigen::VectorXd::NullaryExpr(xyz.rows(),
      [&xyz](auto r) { double the=acos(xyz(r,2)); return the>M_PI_2?the-M_PI_2:the; }),
    Eigen::VectorXd::NullaryExpr(xyz.rows(),
      [&xyz](auto r) { double psi=atan2(xyz(r,1), xyz(r,0)); return psi>M_PI_2?psi-M_PI_2:psi; });
}

#endif //__UTIL_H__