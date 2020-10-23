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

    Eigen::MatrixXd V,N,UV,D;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header,eig_header;
 
    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, D, header))
    {
      size_t k = par["k"].as<int>();

      if(par.count("verbose"))
      {
        std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
      }

      Eigen::MatrixXd U;
      Eigen::VectorXd S;

      Eigen::SparseMatrix<double> Cot,M;
      Eigen::VectorXd Cot_Sum;

      igl::cotmatrix(V,F,Cot);
      //igl::massmatrix(V,F, igl::MASSMATRIX_TYPE_DEFAULT, M);
      //Cot *= -0.5; // convert to notation from 
      igl::sum(Cot, 1, Cot_Sum);

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
        Eigen::SparseMatrix<double> norm = Eigen::SparseMatrix<double>( Cot_Sum.cwiseSqrt().cwiseInverse().asDiagonal() );
        Laplacian =  Eigen::SparseMatrix<double>(Eigen::VectorXd::Ones(norm.rows()).asDiagonal()) - norm.transpose() * Cot * norm;
      } else if(par.count("rw")) {
        // randow walker normalization
        Eigen::SparseMatrix<double> norm = Eigen::SparseMatrix<double>( Cot_Sum.cwiseInverse().asDiagonal() );
        Laplacian =  Eigen::SparseMatrix<double>(Eigen::VectorXd::Ones(norm.rows()).asDiagonal()) - norm.transpose() * Cot;
      }
      else {
        // default laplacian, no normalization
        Laplacian =  Eigen::SparseMatrix<double>(Cot_Sum.asDiagonal()) - Cot;
      }
     
      try {
        Laplacian.makeCompressed();
        Spectra::SparseSymShiftSolve<double> Laplacian_op(Laplacian);
        Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> >
            eigs(&Laplacian_op, k+1, (k+1)*4, 0.0);
        
        eigs.init();
        eigs.compute(1000, 1e-10, Spectra::SMALLEST_MAGN);

        if(eigs.info() == Spectra::SUCCESSFUL)
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

              igl::writePLY(out, V, F, E, N, UV, ND, header );
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
