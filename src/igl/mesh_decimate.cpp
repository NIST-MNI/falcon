#include <iostream>
#include <unistd.h>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "igl/decimate.h"

#include "cxxopts.hpp"

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Decimate mesh, default : half");

  options
      .positional_help("<input> <output>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("i,input",   "Input mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output mesh ",   cxxopts::value<std::string>())
    ("n",         "Desired number of new vertices", cxxopts::value<int>())
    ("f",         "Desired fraction of new vertices", cxxopts::value<double>())
    ("help",      "Print help") ;
  
  options.parse_positional({"input", "output"});
  auto par = options.parse(argc, argv);

  if( par.count("input") && 
      par.count("output")  )
  {

    Eigen::MatrixXd V,N,UV,VD,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> headerV;
    std::vector<std::string> headerF,headerE,comments;

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, VD, headerV,
                    FD,headerF, ED,headerE, comments))
    {
      if(par["verbose"].as<bool>()) {
        std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << VD.rows() << "x"<< VD.cols() << std::endl;
        std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
        std::cout<<"header:";

        for(auto h:headerV)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

      Eigen::MatrixXd V1,N1,UV1,FD1,ED1;
      Eigen::MatrixXi F1,E1;
      Eigen::VectorXi J,I;
      std::vector<std::string> header1F,header1E;

      int n=V.rows()/2;

      if(par.count("n")) n=par["n"].as<int>();
      else if(par.count("f")) n=par["f"].as<double>()*V.rows();

      //decimating
      if(!igl::decimate(V, F, par["n"].as<int>(), V1, F1, J, I))
        std::cerr<<"Failed achieving desired number of target vertices"<<std::endl;

      // map data to the new vertices
      Eigen::MatrixXd D1 = 
        Eigen::MatrixXd::NullaryExpr( V1.rows(), VD.cols(), 
          [ &I, &VD]( auto r, auto c) {
            return VD(I(r),c);
          });

      //TODO: populate E1!
      if(!igl::writePLY(par["output"].as<std::string>(), V1, F1, E1, N1, UV1, D1, headerV,
                  FD1, header1F, ED1, header1E, comments, igl::FileEncoding::Binary ))
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
