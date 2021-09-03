#include <iostream>
#include <unistd.h>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "cxxopts.hpp"

#include "minc2-simple.h"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Apply .xfm transform to a mesh");

  options
      .positional_help("<input> <output>")
      .show_positional_help();
  
  options.add_options()
    ("verbose", "Verbose",         cxxopts::value<bool>()->default_value("false"))
    ("i,input",   "Input mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output mesh ",   cxxopts::value<std::string>())
    ("t,transform","XFM transform",   cxxopts::value<std::string>())
    ("invert_transform", "Invert transform",         cxxopts::value<bool>()->default_value("false"))
    ("help", "Print help") ;
  
  options.parse_positional({"input","output"});
  auto par = options.parse(argc, argv);

  if( par.count("input") && par.count("output"))
  {

    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header,headerF,headerE;
    std::vector<std::string> comments;

    if(igl::readPLY(par["input"].as<std::string>(), V, F, E, N, UV, D, header, FD, headerF, ED, headerE, comments))
    {

      if(par["verbose"].as<bool>()) {
        std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
        std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
        std::cout << "header:";

        for(auto h:header)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

      if(par.count("transform"))
      {
        minc2_xfm_file_handle h=minc2_xfm_allocate0();
        if(!minc2_xfm_open(h,par["transform"].as<std::string>().c_str()))
        {
          bool invert_xfm=par["invert_transform"].as<bool>();

          for(int i=0;i<V.rows();++i)
          {
            Eigen::RowVector3d in=V.row(i);
            Eigen::RowVector3d out;
            if(invert_xfm)
              minc2_xfm_inverse_transform_point(h,in.data(),out.data());
            else
              minc2_xfm_transform_point(h,in.data(),out.data());
            V.row(i)=out;
          }
        } else {
          std::cerr<<"Error reading xfm:"<<par["transform"].as<std::string>()<<std::endl;
          return 1;
        }
        minc2_xfm_free(h);
      }

      comments.push_back("mesh_transform"); //TODO: make a proper stamp

      igl::writePLY(par["output"].as<std::string>(), V, F, E, N, UV, D, header, FD, headerF, ED, headerE, comments, igl::FileEncoding::Binary  );
             
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
