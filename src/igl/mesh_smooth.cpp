#include <iostream>
#include <unistd.h>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/adjacency_matrix.h"
#include "igl/sum.h"

#include "cxxopts.hpp"


int uniform_laplacian(const Eigen::MatrixXi &F,Eigen::SparseMatrix<double> &U)
{
  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);
  Eigen::SparseVector<double> Asum;
  igl::sum(A, 1, Asum);
  Eigen::SparseMatrix<double> Adiag;
  igl::diag(Asum, Adiag);
  U = A - Adiag ;
  return 0;
}


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Smooth mesh using Taubin algorithm");

  options
      .positional_help("<input> <output>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("i,input",   "Input mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output mesh ",   cxxopts::value<std::string>())
    ("iter",      "Number of iterations", cxxopts::value<int>()->default_value("100"))
    ("l,lambda",  "Lambda (>0)", cxxopts::value<double>()->default_value("0.1"))
    ("m,mju",     "Mju (<0)", cxxopts::value<double>()->default_value("-0.11"))
    ("help", "Print help") ;
  
  options.parse_positional({"input", "output"});
  auto par = options.parse(argc, argv);

  if( par.count("input") && 
      par.count("output") )
  {

    double mju=par["mju"].as<double>();
    double lambda=par["lambda"].as<double>();
    int iter=par["iter"].as<int>();

    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> headerV,headerF,headerE,comments;

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, 
                    D, headerV,
                    FD,headerF, 
                    ED,headerE, 
                    comments))
    {
      if(par["verbose"].as<bool>()) {
        std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
        std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
        std::cout<<"header:";

        for(auto h:headerV)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

      // Compute Laplace-Beltrami operator: #V by #V
      Eigen::SparseMatrix<double> L;
      igl::cotmatrix(V, F, L);
      //uniform_laplacian(F, L);
      // mass matrix
      //Eigen::SparseMatrix<double> M;
      //igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
      Eigen::SparseMatrix<double> I = Eigen::SparseMatrix<double>(Eigen::VectorXd::Ones(V.rows()).asDiagonal());
      //
      
      // apply smoothing
      auto S1 = I - lambda*L;
      auto S2 = I - mju*L;
      

      for(int i=0;i<iter;++i)
      {
        auto V1 = S1*V;
        auto V2 = S2*V1;
        //V = V2.eval();
        V = V2.eval();
      }

      if(!igl::writePLY(par["output"].as<std::string>(), V, F, E, N, UV, 
                  D,  headerV,
                  FD, headerF, 
                  ED, headerE, comments, igl::FileEncoding::Binary ))
      {
        std::cerr<<"Error writing:" << par["output"].as<std::string>() << std::endl;
        return 1;
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