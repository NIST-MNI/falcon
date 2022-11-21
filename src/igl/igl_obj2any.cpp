#include <iostream>
#include <unistd.h>
// IGL
#include "readMNIObj.h"
//#include "writeMNIObj.h"
#include "cxxopts.hpp"
#include <igl/write_triangle_mesh.h>
#include <igl/writePLY.h>
#include <igl/edges.h>

int main(int argc,char **argv)
{
  cxxopts::Options options(argv[0], "Convert MNI obj files to other formats");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",         cxxopts::value<bool>()->default_value("false"))
    ("i,input",  "Input mesh  ",           cxxopts::value<std::string>())
    ("o,output", "Output mesh  ",           cxxopts::value<std::string>())

    ("clobber", "Clobber output file ",      cxxopts::value<bool>()->default_value("false"))
    ("falcon", "Try to generate .ply file compatible with falcon ", cxxopts::value<bool>()->default_value("false"))

    ("help", "Print help") ;
  
  options.parse_positional({"input","output"});

  auto par = options.parse(argc, argv);

  if( par.count("input") && 
      par.count("output")  )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }

    bool verbose=par["verbose"].as<bool>();
    bool for_falcon=par["falcon"].as<bool>();

    Eigen::MatrixXd V,N,C;
    Eigen::MatrixXi F;

    if(igl::readMNIObj(par["input"].as<std::string>(),V,F,N,C))
    {
      if(verbose)
      {
        std::cout<<"Vertices:"<<V.rows()<<" Faces:"<<F.rows()<<" Normals:"<<N.rows()<<" Colors:"<<C.rows()<<std::endl;
        if(C.rows()==1)
          std::cout<<"Object colour:"<<C<<std::endl;
      }

      if(for_falcon)
      {
        if(verbose)
          std::cout<<"Generating list of edges, to be compatible with falcon, not saving normals"<<std::endl;
        
        Eigen::MatrixXi E;
        igl::edges(F,E);

        if(verbose)
          std::cout<<"Edges:"<<E.rows()<<std::endl;
        
        //TODO: save color if present
        if(!igl::writePLY(par["output"].as<std::string>(), V, F, E, igl::FileEncoding::Binary))
          std::cerr<<"Error writing:"<<par["output"].as<std::string>()<<std::endl;

      } else {
        if(!igl::write_triangle_mesh(par["output"].as<std::string>(), V, F, igl::FileEncoding::Binary))
          std::cerr<<"Error writing:"<<par["output"].as<std::string>()<<std::endl;
      }
    } else {
      std::cerr<<"Error reading:"<<par["input"].as<std::string>()<<std::endl;
    }
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }

  return 0;
}