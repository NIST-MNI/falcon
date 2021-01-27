#include <iostream>
#include <Eigen/Core>
#include <igl/marching_cubes.h>
#include <igl/writePLY.h>
#include <unistd.h>

#include "cxxopts.hpp"

#include "minc2-simple.h"

bool load_volume(const char*minc_file,Eigen::MatrixXd &GV, Eigen::VectorXd &V,int &x,int &y,int &z )
{
    struct minc2_dimension * _dims;
    int ndim, nelement;
    // load minc file
    minc2_file_handle h = minc2_allocate0();

    if(minc2_open(h, minc_file)!=MINC2_SUCCESS)
    {
        std::cerr << "Can't open " << minc_file  << " for reading" << std::endl;
        minc2_free(h);
        return false;
    }
    minc2_setup_standard_order(h);
    minc2_ndim(h, &ndim);
    minc2_nelement(h, &nelement);
    minc2_get_representation_dimensions(h, &_dims);

    //std::cout<<"load_float_tensor:" << _dims[0].length << "," << _dims[1].length << "," << _dims[2].length << std::endl;
    //torch::Tensor output  = torch::empty( { _dims[2].length, _dims[1].length, _dims[0].length }, torch::kFloat32);
    x=_dims[0].length;
    y=_dims[1].length;
    z=_dims[2].length;
    //Eigen::Array<float, Eigen::Dynamic, 1> _V(_dims[2].length* _dims[1].length* _dims[0].length);
    V.resize(_dims[2].length* _dims[1].length* _dims[0].length);

    if( minc2_load_complete_volume(h, V.data(), MINC2_DOUBLE) != MINC2_SUCCESS )
    {
        std::cerr << "Error reading data from minc file" << std::endl;
        minc2_close(h);
        minc2_free(h);
        return false;
    }
    minc2_close(h);
    minc2_free(h);

    GV.resize( _dims[2].length * _dims[1].length*_dims[0].length ,3);

    Eigen::RowVector3d start(_dims[0].start, _dims[1].start, _dims[2].start);
    Eigen::RowVector3d step(_dims[0].step, _dims[1].step, _dims[2].step);


    // initialize coordinates
    // GV = Eigen::MatrixXd::NullaryExpr( _dims[2].length * _dims[1].length*_dims[0].length ,3,
    //     [&_dims](int r,int c) { 
    //         switch(c) 
    //         {
    //             case 2:return   r/(_dims[1].length*_dims[0].length);
    //             case 1:return  (r%(_dims[1].length*_dims[0].length))/_dims[0].length;
    //             case 0:
    //             default:return (r%(_dims[1].length*_dims[0].length))%_dims[0].length;
    //         };
    //      });
    for(int zi=0;zi<z;++zi)
        for(int yi=0;yi<y;++yi)
            for(int xi=0;xi<x;++xi)
                GV.row(xi+yi*x+zi*x*y)=(Eigen::RowVector3d(xi,yi,zi).array()*step.array()).matrix()+start;
    // TODO: convert to world coordinates
    return true;
}


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Resample field");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",     cxxopts::value<bool>()->default_value("false"))
    ("i,input", "Source volume ",        cxxopts::value<std::string>())
    ("o,output", "Output mesh ",        cxxopts::value<std::string>())
    ("clobber", "Clobber output file ", cxxopts::value<bool>()->default_value("false"))
    ("t,threshold",    "Threshold ", cxxopts::value<double>()->default_value("0.5"))

    ("help", "Print help") ;
  
  options.parse_positional({"input", "output" });
  auto par = options.parse(argc, argv);

  if( par.count("input") && 
      par.count("output")  )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }
    Eigen::MatrixXd GV;
    Eigen::VectorXd V;
    int x,y,z;
    if( !load_volume(par["input"].as<std::string>().c_str(),GV,V,x,y,z) )
        return 1;

    std::cout<<"Grid:"<<GV.rows()<<"x"<<GV.cols()<<std::endl;
    std::cout<<" Vol:"<< V.rows()<<"x"<< V.cols()<<std::endl;

    Eigen::MatrixXd SV;
    Eigen::MatrixXi SF;
    igl::marching_cubes(V,GV, x,y,z, par["threshold"].as<double>(), SV,SF);
    std::cout<<"Extracted surface V:"<<SV.rows()<<" F:"<<SF.rows()<<std::endl;

    igl::writePLY(par["output"].as<std::string>(), SV, SF );
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
