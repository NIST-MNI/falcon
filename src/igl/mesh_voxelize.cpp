#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/sum.h>
#include <igl/AABB.h>
#include <igl/signed_distance.h>

#include <igl/read_triangle_mesh.h>
#include <igl/writePLY.h>
#include <unistd.h>

#include "cxxopts.hpp"

#include "minc_volume.h"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Mesh voxelizer, calculate signed distance to closest mesh face per voxel");

  options
      .positional_help("<input mesh> <reference volume> <output volume>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",         cxxopts::value<bool>()->default_value("false"))
    ("i,input", "Source mesh ",           cxxopts::value<std::string>())

    ("o,output", "Output volume ",             cxxopts::value<std::string>())
    ("l,like",   "Reference volume ",          cxxopts::value<std::string>())
    ("clobber", "Clobber output file ",      cxxopts::value<bool>()->default_value("false"))
    ("help", "Print help") ;
  
  options.parse_positional({"input", "like", "output"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();

  Eigen::MatrixXd V;
  Eigen::MatrixXd W;
  Eigen::MatrixXi F;

  if( par.count("input") &&  par.count("like") )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }

    if(!igl::read_triangle_mesh(par["input"].as<std::string>(),W,F))
    {
        std::cerr<<"Error loading:"<<par["input"].as<std::string>().c_str()<<std::endl;
        return 1;
    }

    minc_volume vol;
    if( !load_volume(par["like"].as<std::string>().c_str(),vol, true ) )
    {
        std::cerr<<"Error loading:"<<par["like"].as<std::string>().c_str()<<std::endl;
        return 1;
    }
    
    V.resize(W.rows(),3);
    // convert to voxel coordinates
    for(int i=0;i<V.rows();++i)
      V.row(i)=vol.world_to_voxel(W.row(i));

    Eigen::MatrixXd GV;
    GV.resize( vol.dims.prod() ,3);

    for(int zi=0;zi<vol.dims[2];++zi)
        for(int yi=0;yi<vol.dims[1];++yi)
            for(int xi=0;xi<vol.dims[0];++xi)
                GV.row(xi+yi*vol.dims[0]+zi*vol.dims[0]*vol.dims[1]) = Eigen::RowVector3d(xi,yi,zi);

    /*       
    igl::AABB<Eigen::MatrixXd,3> tree;
    igl::FastWindingNumberBVH fwn_bvh;

    tree.init(V,F);*/
    Eigen::VectorXi CF;
    Eigen::MatrixXd CV,CN;

    igl::signed_distance(GV,V,F,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,vol.volume,CF,CV,CN);

    if(par.count("output"))
      save_volume(par["output"].as<std::string>().c_str(),
        par["like"].as<std::string>().c_str(),vol,false);
        
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
