#include <iostream>
#include <unistd.h>
#include <Eigen/Sparse>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "cxxopts.hpp"

#include "util.h"

template <typename DerivedF, typename DerivedV, typename T>
IGL_INLINE void length_matrix(
  const Eigen::MatrixBase<DerivedF> & F,
  const Eigen::MatrixBase<DerivedV> & V,
  Eigen::SparseMatrix<T>& A)
{
  using namespace std;
  using namespace Eigen;
  typedef typename DerivedF::Scalar Index;
  typedef typename DerivedV::Scalar Scalar;

  typedef Triplet<T> IJV;
  vector<IJV > ijv;
  ijv.reserve(F.size()*2);
  // Loop over **simplex** (i.e., **not quad**)
  for(int i = 0;i<F.rows();i++)
  {
    // Loop over this **simplex**
    for(int j = 0;j<F.cols();j++)
    for(int k = j+1;k<F.cols();k++)
    {
      // Get indices of edge: s --> d
      Index s = F(i,j);
      Index d = F(i,k);
      Scalar l = (V.row(s)-V.row(d)).norm();

      ijv.push_back(IJV(s,d,l));
      ijv.push_back(IJV(d,s,l));
    }
  }

  const Index n = V.rows();
  A.resize(n,n);
  switch(F.cols())
  {
    case 3:
      A.reserve(6*(F.maxCoeff()+1));
      break;
    case 4:
      A.reserve(26*(F.maxCoeff()+1));
      break;
  }
  A.setFromTriplets(ijv.begin(),ijv.end());
}


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Debug mesh");

  options
      .positional_help("<input> <output>")
      .show_positional_help();
  
  options.add_options()
    ("i,input",   "Input mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output mesh ",   cxxopts::value<std::string>())
    ("help", "Print help") ;
  
  options.parse_positional({"input"});
  auto par = options.parse(argc, argv);

  if( par.count("input") )
  {

    Eigen::MatrixXd V,N,UV,D;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header;

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, D, header))
    {
      std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
      std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
      std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
      std::cout << "header:";

      for(auto h:header)
        std::cout<<h<<"\t";
      std::cout << std::endl;

      // calculate edge lengths
      Eigen::SparseMatrix<double> xyz_lengths, sph_lengths; 
      
      Eigen::MatrixXd  PT;
      Eigen::MatrixXd  sph;

      // convert spherical coordinates to x,y,z on sphere 
      if( extract_psi_the(header, D, PT))
      {
        sph_to_xyz(PT, sph);
      } else {
        std::cerr<<"Can't get spherical coordinates!"<<std::endl;
        return 1;
      }
      length_matrix(F,V,xyz_lengths);
      xyz_lengths.makeCompressed();
      Eigen::VectorXd _xyz_lengths=Eigen::Map<Eigen::VectorXd>(xyz_lengths.valuePtr(),xyz_lengths.nonZeros());

      std::cout << "XYZ length min:" << _xyz_lengths.minCoeff() << " max:" << _xyz_lengths.maxCoeff() << std::endl;

      length_matrix(F,sph,sph_lengths);
      //sph_lengths.makeCompressed();
      Eigen::VectorXd _sph_lengths=Eigen::Map<Eigen::VectorXd>(sph_lengths.valuePtr(),sph_lengths.nonZeros());

      std::cout << "SPH length min:" << _sph_lengths.minCoeff() << " max:" << _sph_lengths.maxCoeff() << std::endl;

      if( par.count("output"))
      {
        // find zero length
        // Iterate over outside
        for(int k=0; k<sph_lengths.outerSize(); ++k)
        {
          // Iterate over inside
          for(typename Eigen::SparseMatrix<double>::InnerIterator it(sph_lengths,k); it; ++it)
          {
            //assert(it.value() != 0);
            if(sph_lengths.coeffRef(it.row(),it.col()) < 1e-8 )
               sph_lengths.coeffRef(it.row(),it.col()) = 1.0;
            else
               sph_lengths.coeffRef(it.row(),it.col()) = 0.0;
          }
        }
        // zero edge count
        Eigen::VectorXd zero_length_count = sph_lengths * Eigen::VectorXd::Ones(sph_lengths.cols());
        header.push_back("zc");
        Eigen::MatrixXd D1(D.rows(), D.cols()+1);
        D1<<D,zero_length_count;

        igl::writePLY(par["output"].as<std::string>(), V, F, E, N, UV, D1, header );
      }
             
    } else {
      std::cerr<<"Error reding ply:"<<argv[1]<<std::endl;
      return 1;
    }
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
