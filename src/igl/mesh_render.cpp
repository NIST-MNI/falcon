#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <igl/jet.h>
#include <igl/parula.h>
#include <igl/png/writePNG.h>

#include <Eigen/Geometry> 

// embree
#include <igl/embree/EmbreeRenderer.h>


#include "igl/readPLY.h"
#include "csv/readCSV.h"

#include "cxxopts.hpp"

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Render mesh fields");

  options
      .positional_help("<mesh> <output>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose",     "Verbose output",  cxxopts::value<bool>()->default_value("false"))
    ("m,mesh",        "Input mesh ",      cxxopts::value<std::string>())
    ("o,output",      "Output png file ", cxxopts::value<std::string>())
    ("c,csv",         "Input csv ",     cxxopts::value<std::string>())
    ("f,field",       "Field name ",    cxxopts::value<std::string>())
    ("n",             "Field number ",  cxxopts::value<int>()->default_value("0"))
    ("w,width",       "Frame width ",  cxxopts::value<int>()->default_value("1280"))
    ("h,height",      "Frame height ", cxxopts::value<int>()->default_value("960"))
    ("z,zoom",        "Zoom ", cxxopts::value<double>()->default_value("1.0"))

    ("rx",        "Rotation around X ", cxxopts::value<double>()->default_value("0.0"))
    ("ry",        "Rotation around Y ", cxxopts::value<double>()->default_value("0.0"))
    ("rz",        "Rotation around Z ", cxxopts::value<double>()->default_value("0.0"))

    ("view",      "preset view: front back top bottom left right",     cxxopts::value<std::string>())

    ("flat",    "Flat color per face (average of corners)",  cxxopts::value<bool>()->default_value("false"))
    ("ortho",   "Use orthographic projection",  cxxopts::value<bool>()->default_value("false"))

    // colormaps
    ("cmap", "Color map: inferno,jet,magma,parula,plasma,viridis,turbo ", cxxopts::value<std::string>()->default_value("jet"))
    ("min",  "Min value ", cxxopts::value<double>())
    ("max",  "Max value ", cxxopts::value<double>())
    ("relevel", "Relevel discrete labels befor coloring ", cxxopts::value<bool>()->default_value("false"))


    ("help", "Print help")
  ;
  
  options.parse_positional({"mesh","output"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();
  if(par.count("mesh")&&par.count("output"))
  {
    Eigen::MatrixXd V,N,UV,D;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header;
 
    if(igl::readPLY(par["mesh"].as<std::string>(), V, F, E, N, UV, D, header))
    {
      const size_t k = 5;
      int idx_field = par["n"].as<int>();

      if(verbose) {
        std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
      }

      if(par.count("csv"))
      {
        Eigen::MatrixXd Dcsv;
        std::vector<std::string> csv_header;
        if(verbose) 
          std::cout<<"Reading from :"<<par["csv"].as<std::string>()<<std::endl;
        
        igl::readCSV(par["csv"].as<std::string>(), Dcsv, csv_header);
        if(Dcsv.cols()>0)
        {
          Eigen::MatrixXd D_(V.rows(), D.cols()+Dcsv.cols());
          D_ << D,Dcsv;
          for(auto const &h:csv_header)
            header.push_back(h);
          D = D_;
        }
      }

      if(verbose && !header.empty())
      {
        std::cout<<"Data:";
        for(auto h:header)
          std::cout<<h<<"\t";

        std::cout<<std::endl;
      }

      if(par.count("field"))
      {
        int i;
        for(i=0;i<header.size();++i)
          if(par["field"].as<std::string>()==header[i])
          {
            idx_field=i;
            break;
          }
        if(i==header.size())
        {
          std::cerr<<"Can't find field:\""<<par["field"].as<std::string>()<<"\""<<std::endl;
          return 0;
        }
      }

      bool flat_color  =par["flat"].as<bool>();
      bool orthographic=par["ortho"].as<bool>();

      // embree object
      igl::embree::EmbreeRenderer er;
      er.set_mesh(V,F,true);

      if(idx_field < D.cols() && idx_field>=0)
      {
        Eigen::VectorXd field=D.col(idx_field);

        if(par["relevel"].as<bool>())
        {
          // TODO: use closest int?
          Eigen::VectorXi labels=D.col(idx_field).cast<int>();
          std::set<int> levels;
          for(int i=0;i<labels.size();++i)
          {
            levels.insert(labels(i));
          }
          if(verbose)
            std::cout<<"Found "<<levels.size()<<" levels"<<std::endl;
          std::map<int,int> relabel;
          int lv=0;
          for(auto l:levels)
            relabel[l]=lv++;

          for(int i=0;i<labels.size();++i)
          {
            field(i)=relabel[labels(i)];
          }
        }

        double v_min,v_max;

        if(par.count("min"))
          v_min=par["min"].as<double>();
        else
          v_min=field.minCoeff();

        if(par.count("max"))
          v_max=par["max"].as<double>();
        else
          v_max=field.maxCoeff();

        Eigen::MatrixXd C;
        // select field
        if(verbose) 
          std::cout<<"Field:"<<idx_field<<" "<<field.minCoeff()<<" "<<field.maxCoeff()<<std::endl;
        
        if(par["cmap"].as<std::string>()=="jet")
          igl::colormap(igl::COLOR_MAP_TYPE_JET,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="inferno")
          igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="magma")
          igl::colormap(igl::COLOR_MAP_TYPE_MAGMA,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="parula")
          igl::colormap(igl::COLOR_MAP_TYPE_PARULA,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="plasma")
          igl::colormap(igl::COLOR_MAP_TYPE_PLASMA,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="viridis")
          igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS,field,v_min,v_max,C);
        else if(par["cmap"].as<std::string>()=="turbo")
          igl::colormap(igl::COLOR_MAP_TYPE_TURBO,field,v_min,v_max,C);

        // calculate per-face color for flat coloring
        Eigen::MatrixXd CF;
        if(flat_color)
        {
          std::cout<<"Using flat coloring"<<std::endl;
          CF=Eigen::MatrixXd::NullaryExpr(F.rows(),3, [&C,&F](auto r,auto c) {
            return (C(F(r,0),c)+C(F(r,1),c)+C(F(r,2),c))/3;
          });
          er.set_colors(CF);
        } else {
          er.set_colors(C);
        }
      } else {
        // using uniform colour
        er.set_uniform_color(Eigen::RowVector3d(1.0,0.0,0.0));
      }

      // the  viewport
      int width=par["width"].as<int>();
      int height=par["height"].as<int>();

      Eigen::Matrix3d rot_matrix;

      auto  make_rot_matrix = [](auto rx,auto ry,auto rz) {
        return  Eigen::AngleAxisd( rx*M_PI/180.0, Eigen::Vector3d::UnitX())
              * Eigen::AngleAxisd( ry*M_PI/180.0, Eigen::Vector3d::UnitY())
              * Eigen::AngleAxisd( rz*M_PI/180.0, Eigen::Vector3d::UnitZ());
      };
      
      if(par.count("view"))
      { // front back top bottom left right
        if(par["view"].as<std::string>()=="top")
          rot_matrix =  make_rot_matrix(0,0,0);  
        else if(par["view"].as<std::string>()=="bottom")
          rot_matrix =  make_rot_matrix(0,180,0);  
        else if(par["view"].as<std::string>()=="front")
          rot_matrix =  make_rot_matrix(90,0,180);  
        else if(par["view"].as<std::string>()=="back")
          rot_matrix =  make_rot_matrix(-90,0,0);  
        else if(par["view"].as<std::string>()=="left")
          rot_matrix =  make_rot_matrix(0,90,90);  
        else if(par["view"].as<std::string>()=="rigth")
          rot_matrix =  make_rot_matrix(0,-90,-90);  
        else 
        {
          std::cerr<<"Unknown view:"<<par["view"].as<std::string>()<<std::endl;
          return 1;
        }
      } else {      
        rot_matrix =  make_rot_matrix(par["rx"].as<double>(),par["ry"].as<double>(),par["rz"].as<double>());
      }

      er.set_zoom(par["zoom"].as<double>());
      er.set_rot(rot_matrix);
      er.set_orthographic(par["ortho"].as<bool>());

      // render view using embree
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(width, height);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(width, height);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(width, height);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(width, height);

      er.render_buffer(R,G,B,A);

      igl::png::writePNG(R,G,B,A,par["output"].as<std::string>());

    } else {
      std::cerr<<"Error reading ply:"<<argv[1]<<std::endl;
      return 1;
    }
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }

  return 0;
}