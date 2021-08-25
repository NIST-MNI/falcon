#include <iostream>
#include <set>


#include <igl/read_triangle_mesh.h>
#include <igl/jet.h>
#include <igl/parula.h>

#ifdef STB_PNG
#include <writePNG.h>
#else
#include <igl/png/writePNG.h>
#endif

#include <Eigen/Geometry> 

// embree
#include <igl/embree/EmbreeRenderer.h>
#include <igl/embree/reorient_facets_raycast.h>

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
    ("RGB",           "Use RGB values from mesh, of present",  cxxopts::value<bool>()->default_value("false"))

    ("w,width",       "Frame width ",  cxxopts::value<int>()->default_value("1280"))
    ("h,height",      "Frame height ", cxxopts::value<int>()->default_value("960"))
    ("z,zoom",        "Zoom ", cxxopts::value<double>()->default_value("1.0"))

    ("rx",        "Rotation around X ", cxxopts::value<double>()->default_value("0.0"))
    ("ry",        "Rotation around Y ", cxxopts::value<double>()->default_value("0.0"))
    ("rz",        "Rotation around Z ", cxxopts::value<double>()->default_value("0.0"))

    ("view",      "preset view: front back top bottom left right",     cxxopts::value<std::string>())
    ("six",       "Display all six views on one rendering, multiply output size X*col x Y*6/col",     cxxopts::value<bool>()->default_value("false"))
    ("cols",      "Number of columns for six view",     cxxopts::value<int>()->default_value("2"))

    ("flat",    "Flat color per face (average of corners)",  cxxopts::value<bool>()->default_value("false"))
    ("ortho",   "Use orthographic projection",  cxxopts::value<bool>()->default_value("false"))

    // colormaps
    ("cmap", "Color map: inferno,jet,magma,parula,plasma,viridis,turbo ", cxxopts::value<std::string>()->default_value("jet"))
    ("min",  "Min value ", cxxopts::value<double>())
    ("max",  "Max value ", cxxopts::value<double>())
    ("relevel", "Relevel discrete labels befor coloring ", cxxopts::value<bool>()->default_value("false"))

    ("r",  "flat color r ", cxxopts::value<double>())
    ("g",  "flat color g ", cxxopts::value<double>())
    ("b",  "flat color b ", cxxopts::value<double>())

    ("fix", "Fix normals (long!) ", cxxopts::value<bool>()->default_value("false"))
    ("double", "Double sided triangles ", cxxopts::value<bool>()->default_value("false"))


    ("help", "Print help")
  ;
  
  options.parse_positional({"mesh","output"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();
  if(par.count("mesh") && par.count("output"))
  {
    Eigen::MatrixXd V,N,UV,VD,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> headerV;
    std::vector<std::string> headerF, headerE, comments;
 
    if(igl::readPLY(par["mesh"].as<std::string>(), V, F, E, N, UV, VD, headerV,
              FD,headerF, ED,headerE, comments))
    {
      const size_t k = 5;
      int idx_field = -1; ;

      if(verbose) {
        std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        
        std::cout << "Vertex Data: " << VD.rows() << "x"<< VD.cols() << std::endl;

        if(!headerV.empty())
        {
          std::cout<<"\t";
          for(auto h:headerV)
            std::cout<<h<<"\t";

          std::cout<<std::endl;
        }

        std::cout << "Face Data: " << FD.rows() << "x"<< FD.cols() << std::endl;
        if(!headerF.empty())
        {
          std::cout<<"\t";
          for(auto h:headerF)
            std::cout<<h<<"\t";

          std::cout<<std::endl;
        }

        std::cout << "Edge Data: " << ED.rows() << "x"<< ED.cols() << std::endl;
        if(!headerE.empty())
        {
          std::cout<<"\t";
          for(auto h:headerE)
            std::cout<<h<<"\t";

          std::cout<<std::endl;
        }
      }

      if(par.count("fix"))
      {
        //need to fix normals
        Eigen::VectorXi I;
        igl::embree::reorient_facets_raycast(V,F,F,I);
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
          Eigen::MatrixXd D_(V.rows(), VD.cols()+Dcsv.cols());
          D_ << VD,Dcsv;
          for(auto const &h:csv_header)
            headerV.push_back(h);
          VD = D_;
        }
      }

      if(par.count("field"))
      {
        auto idx_field_=std::find(headerV.begin(),headerV.end(),par["field"].as<std::string>() );
        if(idx_field_==headerV.end())
        {
          std::cerr<<"Can't find field:\""<<par["field"].as<std::string>()<<"\""<<std::endl;
          return 1;
        }
        idx_field=idx_field_-headerV.begin();
      } else if(par.count("n") && par["n"].as<int>()>=0 && par["n"].as<int>()<headerV.size()) {
        idx_field=par["n"].as<int>();
      } 

      bool flat_color  =par["flat"].as<bool>();
      bool rgb_color   =par["RGB"].as<bool>();
      bool orthographic=par["ortho"].as<bool>();

      // embree object
      igl::embree::EmbreeRenderer er;
      er.set_mesh(V,F,true);
      er.set_double_sided(par["double"].as<bool>());

      if(idx_field < VD.cols() && idx_field>=0)
      {
        Eigen::VectorXd field=VD.col(idx_field);

        if(par["relevel"].as<bool>())
        {
          // TODO: use closest int?
          Eigen::VectorXi labels=VD.col(idx_field).cast<int>();
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
          CF = Eigen::MatrixXd::NullaryExpr(F.rows(),3, [&C,&F](auto r,auto c) {
            return (C(F(r,0),c)+C(F(r,1),c)+C(F(r,2),c))/3;
          });
          er.set_colors(CF);
        } 
        else {
          er.set_colors(C);
        }
      } else if(rgb_color) {
          std::cout<<"Using RGB per-face coloring from mesh"<<std::endl;
          auto red_idx=std::find(headerF.begin(),headerF.end(),"red");
          auto green_idx=std::find(headerF.begin(),headerF.end(),"green");
          auto blue_idx=std::find(headerF.begin(),headerF.end(),"blue");
          if(red_idx==headerF.end())
            std::cerr<<"Can't find red channel"<<std::endl;
          else if(green_idx==headerF.end())
            std::cerr<<"Can't find green channel"<<std::endl;
          else if(blue_idx==headerF.end())
            std::cerr<<"Can't find blue channel"<<std::endl;
          else
          {
            Eigen::MatrixXd CF(F.rows(),3);
            CF <<  FD.col(red_idx-headerF.begin()),
                   FD.col(green_idx-headerF.begin()),
                   FD.col(blue_idx-headerF.begin());
            CF/=255.0;
            er.set_colors(CF);

          }
      } else if(par.count("r") && par.count("g") && par.count("b")){
        // using uniform colour
        er.set_colors(Eigen::RowVector3d(par["r"].as<double>(),par["g"].as<double>(),par["b"].as<double>()));
      } else {
        // using uniform colour
        er.set_colors(Eigen::RowVector3d(1,0.6,0.6));
      }

      // the  viewport
      int width=par["width"].as<int>();
      int height=par["height"].as<int>();

      Eigen::Matrix3d rot_matrix;

      auto  make_rot_matrix = [](auto rx,auto ry,auto rz)->Eigen::Matrix3d {
        Eigen::Matrix3d r;
        r =  Eigen::AngleAxisd( rz*M_PI/180.0, Eigen::Vector3d::UnitZ())
           * Eigen::AngleAxisd( ry*M_PI/180.0, Eigen::Vector3d::UnitY())
           * Eigen::AngleAxisd( rx*M_PI/180.0, Eigen::Vector3d::UnitX());
        return r;
      };
      er.set_zoom(par["zoom"].as<double>());
      er.set_orthographic(par["ortho"].as<bool>());
      
      if(par["six"].as<bool>())
      {
        int cols=par["cols"].as<int>();
        int rows=6/cols;

        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(width*cols, height*rows);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(width*cols, height*rows);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(width*cols, height*rows);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(width*cols, height*rows);

        std::vector< Eigen::Matrix3d> views {
          make_rot_matrix(0,0,0),   make_rot_matrix(0,180,0),
          make_rot_matrix(90,0,180),make_rot_matrix(-90,0,0),
          make_rot_matrix(0,90,90), make_rot_matrix(0,-90,-90)
        };
        for(int i=0;i<views.size();++i)
        {
          er.set_rot(views[i]);
          // render view using embree
          Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> r(width, height);
          Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> g(width, height);
          Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> b(width, height);
          Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> a(width, height);
          er.render_buffer(r,g,b,a);
          R.block( (i%cols)*width,(i/cols)*height,width,height )=r;
          G.block( (i%cols)*width,(i/cols)*height,width,height )=g;
          B.block( (i%cols)*width,(i/cols)*height,width,height )=b;
          A.block( (i%cols)*width,(i/cols)*height,width,height )=a;
        }

        igl::png::writePNG(R,G,B,A,par["output"].as<std::string>());
        
      } else {
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

        er.set_rot(rot_matrix);

        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(width, height);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(width, height);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(width, height);
        Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(width, height);

        // render view using embree
        er.render_buffer(R,G,B,A);

        igl::png::writePNG(R,G,B,A,par["output"].as<std::string>());
      }
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
