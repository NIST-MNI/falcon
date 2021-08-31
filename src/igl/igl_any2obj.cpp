#include <iostream>

// IGL
//#include "readMNIObj.h"
#include "writeMNIObj.h"
#include "cxxopts.hpp"
#include <igl/read_triangle_mesh.h>
#include <igl/readPLY.h>
//#include <igl/edges.h>
#include <igl/per_vertex_normals.h>

int main(int argc,char **argv)
{
  cxxopts::Options options(argv[0], "Convert any format to MNI obj file");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",         cxxopts::value<bool>()->default_value("false"))
    ("i,input",  "Input mesh  ",           cxxopts::value<std::string>())
    ("o,output", "Output mesh  ",           cxxopts::value<std::string>())

    ("clobber", "Clobber output file ",      cxxopts::value<bool>()->default_value("false"))
    //("falcon", "Assume we are dealing with .ply file from falcon ",      cxxopts::value<bool>()->default_value("false"))

    ("help", "Print help") ;
  
  options.parse_positional({"input","output"});

  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();
  //bool for_falcon=par["falcon"].as<bool>();

  Eigen::MatrixXd V,N,C;
  Eigen::MatrixXi F;

  
  if(igl::read_triangle_mesh(par["input"].as<std::string>(),V,F))
  {
    if(verbose)
      std::cout<<"Vertices:"<<V.rows()<<" Faces:"<<F.rows()<<std::endl;

    // generate normals
    igl::per_vertex_normals(V,F,N);

    // if(for_falcon)
    // {
    //   if(verbose)
    //     std::cout<<"Generating list of edges, to be compatible with falcon, not saving normals"<<std::endl;
      
    //   Eigen::MatrixXi E;
    //   igl::edges(F,E);

    //   if(verbose)
    //     std::cout<<"Edges:"<<E.rows()<<std::endl;
      
    //   //TODO: save color if present
    //   if(!igl::writePLY(par["output"].as<std::string>(), V,F,E))
    //     std::cerr<<"Error writing:"<<par["output"].as<std::string>()<<std::endl;

    // } else {
      if(!igl::writeMNIObj(par["output"].as<std::string>(),V,F,N, Eigen::MatrixXd::Ones(1,3), igl::FileEncoding::Binary))
        std::cerr<<"Error writing:"<<par["output"].as<std::string>()<<std::endl;
//    }
  } else {
    std::cerr<<"Error reading:"<<par["input"].as<std::string>()<<std::endl;
  }

  return 0;
}