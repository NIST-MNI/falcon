#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
//#include <Eigen/Sparse>

//#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/sum.h>
// #include <igl/AABB.h>
// #include <igl/signed_distance.h>

#include <igl/read_triangle_mesh.h>
#include <igl/writePLY.h>
#include <igl/readPLY.h>
#include <unistd.h>

#include "csv/writeCSV.h"
#include "cxxopts.hpp"

#include "minc_volume.h"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Field sampler, extract values from a volume correspoding to mesh(es) nodes");

  options
      .positional_help("<input mesh> <source volume> <output measurement>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",          cxxopts::value<bool>()->default_value("false"))
    ("i,input",   "Input mesh ",             cxxopts::value<std::string>())
    ("s,source",  "Source volume ",          cxxopts::value<std::string>())
    ("o,output",  "Output measurement",      cxxopts::value<std::string>())

    ("second",    "Second input mesh ",      cxxopts::value<std::string>())
    ("distance",  "Sampling distance along normal or along link (voxels)",      cxxopts::value<double>()->default_value("0.0"))
    ("default",   "Default value ",      cxxopts::value<double>()->default_value("0.0"))

    ("clobber", "Clobber output file ",         cxxopts::value<bool>()->default_value("false"))
    ("labels",  "Treat source as label volume", cxxopts::value<bool>()->default_value("false"))
    ("header",  "Output header",                cxxopts::value<std::string>()->default_value("sample"))

    ("help", "Print help") ;
  
  options.parse_positional({"input", "source", "output"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();
  double  default_val = par["default"].as<double>();
  double  distance = par["distance"].as<double>();

  Eigen::MatrixXd V,N,UV,D,FD,ED;
  Eigen::MatrixXi F,E;
  std::vector<std::string> headerV,headerF,headerE,comments;

  Eigen::VectorXd M;

  if( par.count("input") &&  par.count("source") && par.count("output") )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }

    if(!igl::readPLY(par["input"].as<std::string>().c_str(), 
            V, F, E, N, UV, 
            D, headerV,
            FD,headerF, 
            ED,headerE, 
            comments) )
    {
        std::cerr<<"Error loading:"<<par["input"].as<std::string>().c_str()<<std::endl;
        return 1;
    }

    minc_volume vol;

    if( !load_volume(par["source"].as<std::string>().c_str(),vol, false ) )
    {
        std::cerr<<"Error loading:"<<par["like"].as<std::string>().c_str()<<std::endl;
        return 1;
    }
    Eigen::MatrixXd VV;
    VV.resize(V.rows(),3);

    // convert to voxel coordinates
    for(int i=0;i<V.rows();++i)
      VV.row(i)=vol.world_to_voxel(V.row(i));

    // calculate normals in voxel space, if needed
    if(distance>0.0) // need to search along normal both in and out
    {
      Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
      Eigen::RowVector3d hi(vol.dims[0]-1,vol.dims[1]-1,vol.dims[2]-1);

      igl::per_vertex_normals(V,F,N);
      if(par["labels"].as<bool>()) 
      {
        M=sample_values_nn<Eigen::VectorXd>(vol, VV, default_val);
      } else {
        M=sample_values<Eigen::VectorXd>(vol, VV, default_val);
      }

      // assume that zero is undefined value and should be replaced with something else
      // TODO: provide a way to configure it
      int bad_values=0;
      int fixed_values=0;
      for(int i=0;i<V.rows();++i)
      {
        if(M(i)==0.0)
        {
          bad_values++;

          double S=M(i);
          for(double d=0.0;d<=distance;d+=1.0)
          {
            Eigen::RowVector3d C1=VV.row(i)+N.row(i)*d; //out 
            Eigen::RowVector3d C2=VV.row(i)-N.row(i)*d; //in
            
            //std::cout<<"V="<<VV.row(i)<<" C1="<<C1<<" C2="<<C2<<std::endl;
            if((C1.array()>=lo.array() && C1.array()<hi.array()).all() ) //inside array
            {
              S=vol.sample_nn(C1);
              if(S!=0.0) break;
            }
            if(S==0.0 && (C2.array()>=lo.array() && C2.array()<hi.array()).all())
            {
              S=vol.sample_nn(C2);
              if(S!=0.0) break;
            }
          }
          if(S!=0.0)
            fixed_values++;
          M(i) = S; //hopefully better value (?)
        }
      }
      std::cout<<"Bad values:"<<bad_values<<" fixed:"<<fixed_values<<std::endl;
    } else {
      if(par["labels"].as<bool>()) 
      {
        M=sample_values_nn<Eigen::VectorXd>(vol, VV, default_val);
      } else {
        M=sample_values<Eigen::VectorXd>(vol, VV, default_val);
      }
    }

    if(par.count("output"))
    {
      if(par["labels"].as<bool>())
      {
        Eigen::VectorXi M_=M.cast<int>();
        if(! igl::writeCSV(par["output"].as<std::string>().c_str(),M_,{par["header"].as<std::string>()}))
        {
          std::cerr<<"Error writing:"<<par["output"].as<std::string>().c_str()<<std::endl;
          return 1;
        }
      }
      else
        if(! igl::writeCSV(par["output"].as<std::string>().c_str(),M,{par["header"].as<std::string>()}))
        {
          std::cerr<<"Error writing:"<<par["output"].as<std::string>().c_str()<<std::endl;
          return 1;
        }
    }

  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
