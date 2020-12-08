#include <igl/read_triangle_mesh.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <iostream>

#include <Eigen/Sparse>

#include <igl/readPLY.h>
#include <igl/writePLY.h>

#include "depth_potential.h"
#include "cxxopts.hpp"
#include "csv/writeCSV.h"

// See Boucher et al 2009 "Depth potential function for folding pattern representation, registration and analysis"
// See http://dx.doi.org/10.1016/j.media.2008.09.001

bool inline check_ext(const std::string &path, const std::string &ext)
{
  auto idx=path.rfind('.');
  if(idx != std::string::npos)
  {
      std::string _extension = path.substr(idx+1);
      return _extension==ext;
  }
  return false;
}


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Calculate depth potential for mesh");
  options
      .positional_help("[optional args]")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("a,alpha", "Alpha parameter", cxxopts::value<double>()->default_value("0.01"))
    ("i,input", "Input file name", cxxopts::value<std::string>())
    ("o,output", "Output file name", cxxopts::value<std::string>())
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
    std::vector<std::string> header;
    std::vector<std::string> headerF,headerE,comments;

    double alpha=par["alpha"].as<double>();

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, D, header, 
       FD,headerF, ED,headerE, comments))
    {

      if(par.count("verbose"))
      {
        std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;

        if(!header.empty())
        {
          std::cout<<"header:";
          for(auto h:header)
            std::cout<<h<<"\t";

          std::cout<<std::endl;
        }
      }
      Eigen::VectorXd dp;
      
      if(!depth_potential(V,F,alpha,dp)) {
          std::cerr<<"Solving failed"<<std::endl;
      } else {
        
        std::string out=par["output"].as<std::string>();

        // if this is supposed to be ply
        if(check_ext(out, "ply") ) {
          Eigen::MatrixXd ND(V.rows(), D.cols()+1);
          ND << D,dp;
          header.push_back("dp");

          igl::writePLY(out, V, F, E, N, UV, ND, header, FD, headerF, ED, headerE, comments, igl::FileEncoding::Binary );
        } else {
          igl::writeCSV(out, dp, {"dp"}); 
        }
      }
    } else {
      std::cerr<<"Error reding ply:"<<argv[1]<<std::endl;
      return 1;
    }
  } catch (const cxxopts::OptionException& e)
  {
    std::cerr << "error parsing options: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
