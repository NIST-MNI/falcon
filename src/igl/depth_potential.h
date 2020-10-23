#ifndef DEPTH_POTENTIAL_H
#define DEPTH_POTENTIAL_H

#include <Eigen/Core>
#include <Eigen/Sparse>


bool depth_potential(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,double alpha,Eigen::VectorXd &dp_out);
bool depth_potential_benchmark(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,double alpha,Eigen::VectorXd &dp_out);

#endif 
