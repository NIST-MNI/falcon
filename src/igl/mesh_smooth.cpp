#include <iostream>
#include <unistd.h>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/adjacency_matrix.h"
#include "igl/sum.h"

#include "cxxopts.hpp"


int uniform_laplacian(const Eigen::MatrixXi &F, Eigen::SparseMatrix<double> &U)
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
    ("i,input",   "Input mesh ",    cxxopts::value<std::string>())
    ("o,output",  "Output mesh ",   cxxopts::value<std::string>())

    ("iter",      "Number of iterations", cxxopts::value<int>()->default_value("100"))
    ("f,fast",    "Fast mode, using adjacency matrix instead of cotmatrix, less stable", cxxopts::value<bool>()->default_value("false"))
    ("implicit",  "Implicit filtering, using lambda only", cxxopts::value<bool>()->default_value("false"))

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
      
      if(par["fast"].as<bool>())
      {
        Eigen::SparseMatrix<double> L;
        uniform_laplacian(F,L);
        // uniform laplacian does not depend on coordinates, no need to recompute
        for(int i=0;i<iter;++i)
        {
          V -= lambda*L*V;
          V -= mju*L*V;
        }
      } else if(par["implicit"].as<bool>()) {
        Eigen::SparseMatrix<double> L;
        igl::cotmatrix(V, F, L);

        for(int i=0;i<iter;++i) {
          Eigen::SparseMatrix<double> M;
          igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC,M);        
          
          // Solve (M-delta*L) U = M*U
          const auto & S = (M - lambda*L);

          Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
          assert(solver.info() == Eigen::Success);
          V = solver.solve(M*V).eval();
        }
        
        //TODO: adjust for centroid?        
      } else {
        for(int i=0;i<iter;++i)
        {
          Eigen::SparseMatrix<double> L;
          igl::cotmatrix(V, F, L);
          V -= lambda*L*V;

          igl::cotmatrix(V, F, L);
          V -= mju*L*V;
        }
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
