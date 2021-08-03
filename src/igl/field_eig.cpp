#include <iostream>

// IGL
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <igl/eigs.h>

#include <Eigen/Sparse>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "csv/writeCSV.h"

// Eigen solver
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

// option parser
#include "cxxopts.hpp"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Calculate mesh eigen vectors");
  options
      .positional_help("[optional args]")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    // laplacian type
    //("rw", "Random walker normalization", cxxopts::value<bool>()->default_value("false"))

    // laplacian normalizations
    ("rw", "Random walker normalization", cxxopts::value<bool>()->default_value("false"))
    ("sym", "Symmetric normalization", cxxopts::value<bool>()->default_value("false"))
    // 
    ("k", "Number of eigen vectors", cxxopts::value<int>()->default_value("5"))
    ("eps", "Epsilon", cxxopts::value<float>()->default_value("1e-10"))
    // IO
    ("i,input", "Input file name (PLY)", cxxopts::value<std::string>())
    ("o,output", "Output file name (PLY|CSV)", cxxopts::value<std::string>())
    ("help", "Print help")
  ;

  options.parse_positional({"input", "output"});

  try {
    auto par = options.parse(argc, argv);

    if (par.count("help") || !par.count("input") || !par.count("output"))
    {
      std::cout << options.help({"", "Group"}) << std::endl;
      return 1;
    }

    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header,eig_header;
    std::vector<std::string> headerF, headerE, comments;

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, D, header,
              FD, headerF, ED, headerE, comments))
    {
      size_t k = par["k"].as<int>();

      if(par.count("verbose"))
      {
        std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Vert Data: " << D.rows() << "x"<< D.cols() << std::endl;

        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Face Data:    " << FD.rows() << "x"<< FD.cols() << std::endl;

        std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
        std::cout << "Edge Data:    " << ED.rows() << "x"<< ED.cols() << std::endl;

        std::cout << "Normals:    " << N.rows() << "x"<< N.cols() << std::endl;

        std::cout << "UV:    " << UV.rows() << "x"<< UV.cols() << std::endl;
      }

      Eigen::MatrixXd U;
      Eigen::VectorXd S;

      Eigen::SparseMatrix<double> Cot;
      Eigen::VectorXd Cot_Sum;

      igl::cotmatrix(V,F,Cot);

      Eigen::SparseMatrix<double> Laplacian;
      
      // if norm == 'sym': # make symmetric normalized laplacian
      //     Disqrt = sparse.dia_matrix((np.power(s, -0.5), 0), shape=(N, N)) # TODO am i doing it right?
      //     dia = sparse.dia_matrix((np.repeat(1.0,N), 0), shape=(N, N))
      //     L = dia - Disqrt @ weights @ Disqrt
      //     L = sparse.lil_matrix( L )
      // elif norm == 'rw': 
      //     Dinv = sparse.dia_matrix((np.reciprocal(s), 0), shape=(N, N)) # TODO am i doing it right?
      //     dia  = sparse.dia_matrix((np.repeat(1.0,N), 0.0), shape=(N, N))
      //     L = dia - Dinv @ weights
      //     L = sparse.lil_matrix( L )
      // else: # default no normalization
      //     dia = sparse.dia_matrix((s, 0), shape=(N, N))
      //     L = dia - weights
      //     L = sparse.lil_matrix( L )

      if(par.count("sym"))  {
        // make symmetric normalized laplacian
        // HACK : set diagonal to zero
        Cot.diagonal().setZero();
        igl::sum(Cot, 0, Cot_Sum);

        Eigen::SparseMatrix<double> norm = Eigen::SparseMatrix<double>( Cot_Sum.cwiseSqrt().cwiseInverse().asDiagonal() );
        Laplacian =  Eigen::SparseMatrix<double>(Eigen::VectorXd::Ones(norm.rows()).asDiagonal()) - norm.transpose() * Cot * norm;
      } else if(par.count("rw")) {
        // HACK : set diagonal to zero
        Cot.diagonal().setZero();
        igl::sum(Cot, 0, Cot_Sum);
        // randow walker normalization
        Eigen::SparseMatrix<double> norm = Eigen::SparseMatrix<double>( Cot_Sum.cwiseInverse().asDiagonal() );
        Laplacian =  Eigen::SparseMatrix<double>(Eigen::VectorXd::Ones(norm.rows()).asDiagonal()) - norm.transpose() * Cot;
      }
      else {
        // default normalization : use original Cot matrix
        Laplacian =  Cot;
      }
     
      try {
        Laplacian.makeCompressed();
        Spectra::SparseSymShiftSolve<double> Laplacian_op(Laplacian);

        Spectra::SymEigsShiftSolver< Spectra::SparseSymShiftSolve<double> >
            eigs(Laplacian_op, k+1, (k+1)*4, par["eps"].as<float>());
        
        eigs.init();
        eigs.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestMagn);

        if(eigs.info() == Spectra::CompInfo::Successful)
        {
            S = eigs.eigenvalues().bottomRows(k);
            std::cout << "Eigenvalues found:\n" << S << std::endl;
            U = eigs.eigenvectors().rightCols(k);
            // 0th eigen value should be 0 and the vector is mostly rounding errors
            for(auto i=0;i<k;++i)
            {
              char tmp[100];
              sprintf(tmp,"eig_%03d",i+1);
              eig_header.push_back(tmp);
            }
            std::string out=par["output"].as<std::string>();

            if(igl::check_ext(out, "ply") ) 
            {
              //concatenate with existing data
              Eigen::MatrixXd ND(std::max(D.rows(),U.rows()), D.cols()+U.cols());
              ND << D,U;

              for(auto const &h:eig_header)
                header.push_back(h);

              igl::writePLY(out, V, F, E, N, UV, ND, header,
                  FD, headerF, ED, headerE, comments, igl::FileEncoding::Binary );
            } else {
              igl::writeCSV(out, U, eig_header ); 
            }
        }
      } catch (const std::invalid_argument& e) {
        std::cerr << "invalid_argument: " << e.what() << std::endl;
        return 1;
      }
    } else {
      std::cerr<<"Error reding ply:" << par["input"].as<std::string>() << std::endl;
      return 1;
    }
  } catch (const cxxopts::OptionException& e)
  {
    std::cerr << "error parsing options: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
