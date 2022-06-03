#include <iostream>
#include <unistd.h>
#include <numeric>
#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "igl/sum.h"

#include "cxxopts.hpp"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Average meshes (by vertex) preserve fields only from 1st mesh");

  options
      .positional_help("<inputs>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("i,input",   "Input meshes ",  cxxopts::value<std::vector<std::string>>())
    ("o,output",  "Output mesh",    cxxopts::value<std::string>())

    ("help", "Print help") ;
  
  options.parse_positional({"input"});
  auto par = options.parse(argc, argv);

  if( par.count("input") && 
      par.count("output") )
  {

    auto& v = par["input"].as<std::vector<std::string>>();
    std::cout << "Number of inputs :"<<size(v) <<std::endl;

    std::vector< Eigen::MatrixXd > Vs;
    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> headerV,headerF,headerE,comments;

    if(igl::readPLY(v[0], 
            V, F, E, N, UV, 
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
    } else {
      std::cerr<<"Error reading ply:"<<v[0]<<std::endl;
      return 1;
    }

    Vs.push_back(V);

    for(int i=1;i<size(v);++i)
    {
      Eigen::MatrixXd V_,N_,UV_,D_,FD_,ED_;
      Eigen::MatrixXi F_,E_;
      std::vector<std::string> headerV_,headerF_,headerE_,comments_;

      if(igl::readPLY(v[i], 
              V_, F_, E_, N_, UV_, 
              D_, headerV_,
              FD_,headerF_, 
              ED_,headerE_, 
              comments_))
      {
        // sanity check
        if(V_.rows()!=V.rows() || V_.cols()!=V.cols())
        {
          std::cerr<<"Inconsistent number of vertices:"<<v[i]<<" Expected:"<<V.rows()<<" Got:"<<V_.rows()<<std::endl;
          return 1;
        } else if(F_.rows()!=F.rows() || F_.cols()!=F.cols()) {
          std::cerr<<"Inconsistent number of faces:"<<v[i]<<" Expected:"<<F.rows()<<" Got:"<<F_.rows()<<std::endl;
          return 1;
        } else if(E_.rows()!=E.rows() || E_.cols()!=E.cols()) {
          std::cerr<<"Inconsistent number of edges:"<<v[i]<<" Expected:"<<E.rows()<<" Got:"<<E_.rows()<<std::endl;
          return 1;
        }
        Vs.push_back(V_);
      } else {
        std::cerr<<"Error reading ply:"<<v[i]<<std::endl;
        return 1;
      }
    }

    //Eigen::MatrixXd Vn=std::accumulate(Vs.begin(),Vs.end(),Eigen::MatrixXd::Zero(Vs[0].rows(),Vs[0].cols()));
    for(size_t i=1;i<size(v);++i)
      Vs[0]+=Vs[i];

    Vs[0]/=size(v);
    //TODO: store std ?

    if(!igl::writePLY(par["output"].as<std::string>(), Vs[0], F, E, N, UV, 
                D,  headerV,
                FD, headerF, 
                ED, headerE, comments, igl::FileEncoding::Binary ))
    {
      std::cerr<<"Error writing:" << par["output"].as<std::string>() << std::endl;
      return 1;
    }
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
