#include <igl/read_triangle_mesh.h>
#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <iostream>
#include <igl/jet.h>
#include <igl/parula.h>

#include <Eigen/Sparse>

// Eigen solver
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

//UI
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include "igl/readPLY.h"
#include "csv/readCSV.h"


#include "cxxopts.hpp"


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "View mesh fields");

  options
      .positional_help("<mesh> [csv]")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("m,mesh",    "Input mesh ", cxxopts::value<std::string>())
    ("c,csv",     "Input csv ", cxxopts::value<std::string>())
    ("help", "Print help")
  ;
  
  options.parse_positional({"mesh", "csv"});
  auto par = options.parse(argc, argv);

  if(par.count("mesh"))
  {
    Eigen::MatrixXd V,N,UV,VD,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> headerV;
    std::vector<std::string> headerF,headerE,comments;

    //face colors
    Eigen::MatrixXd FC;

    if(igl::readPLY(par["mesh"].as<std::string>(), V, F, E, N, UV, VD, headerV,
                  FD,headerF, ED,headerE, comments))
    {
      const size_t k = 5;
      static int idx_field = 0;
      static bool show_mesh_rgb=false;

      if(argc>2)
        idx_field=atoi(argv[2]);

      if(idx_field>=VD.cols()||idx_field<0) idx_field=0;

      std::cout << "Vertices:   " << V.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Faces:      " << F.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Vertex Data:" << VD.rows() << "x"<< VD.cols() << std::endl;

      // read RGB if present
      {
        auto idx_red=std::find(headerF.begin(),headerF.end(),"red");
        auto idx_green=std::find(headerF.begin(),headerF.end(),"green");
        auto idx_blue=std::find(headerF.begin(),headerF.end(),"blue");

        if(idx_red!=headerF.end()&&idx_green!=headerF.end()&& idx_blue!=headerF.end())
        {
          FC.resize(F.rows(),3);
          FC << FD.col(idx_red-headerF.begin()),
                FD.col(idx_green-headerF.begin()),
                FD.col(idx_blue-headerF.begin());
          FC/=255.0;
          show_mesh_rgb=true;
        }
      }

      if(par.count("csv"))
      {
        Eigen::MatrixXd Dcsv;
        std::vector<std::string> csv_header;
        std::cout<<"Reading from :"<<par["csv"].as<std::string>()<<std::endl;
        igl::readCSV(par["csv"].as<std::string>(), Dcsv, csv_header);
        if(Dcsv.cols()>0)
        {
          std::cout << "CSV Data:     " << Dcsv.rows() << "x"<< Dcsv.cols() << std::endl;

          Eigen::MatrixXd D_(V.rows(), VD.cols()+Dcsv.cols());
          D_ << VD,Dcsv;
          for(auto const &h:csv_header)
            headerV.push_back(h);
          VD = D_;
        }
      }

      if(!headerV.empty())
      {
        std::cout<<"\t";
        for(auto h:headerV)
          std::cout<<h<<"\t";

        std::cout<<std::endl;
      }

      Eigen::MatrixXd U;
      Eigen::VectorXd S;

      // Plot the mesh
      igl::opengl::glfw::Viewer viewer;
      viewer.data().set_mesh(V, F);
      viewer.data().set_face_based(true);

      Eigen::MatrixXd field_range(VD.cols(),2);
      double field_range_current[2];

      field_range.col(0) = VD.colwise().minCoeff();
      field_range.col(1) = VD.colwise().maxCoeff();


      auto update_field=[&]() {
        if(show_mesh_rgb)
        {
          viewer.data().set_colors(FC);
        } else if(idx_field<VD.cols() && idx_field>=0) {
            Eigen::MatrixXd C;
            Eigen::MatrixXd _D=VD.col(idx_field);
            igl::colormap(igl::COLOR_MAP_TYPE_JET, _D, field_range(idx_field,0), field_range(idx_field,1), C);
            viewer.data().set_colors(C);
            field_range_current[0] = field_range(idx_field,0);
            field_range_current[1] = field_range(idx_field,1);
        }
      };

      update_field();

      // Attach a custom menu
      igl::opengl::glfw::imgui::ImGuiMenu menu;

      // Customize default menu
      menu.callback_draw_viewer_menu = [&]()
      {
        menu.draw_viewer_menu(); // draw default menu, comment to overwrite

        if(FC.cols()!=0) // add option to show RGB field
        {
          if(ImGui::Checkbox("Mesh RGB",&show_mesh_rgb))
          {
            update_field();
          }
        }

        if(ImGui::Combo("Measure", &idx_field, headerV))
        {
          update_field();
        }

        if(idx_field<VD.cols() && idx_field>=0)
        {
          if(ImGui::InputDouble("min",&field_range_current[0]))
          {
            field_range(idx_field,0)=field_range_current[0];
            update_field();
          }

          if(ImGui::InputDouble("max",&field_range_current[1]))
          {
            field_range(idx_field,1)=field_range_current[1];
            update_field();
          }
        }
      };

      viewer.plugins.push_back(&menu);

      viewer.launch();
    } else {
      std::cerr<<"Error reding ply:"<<argv[1]<<std::endl;
      return 1;
    }
  } else {
    std::cerr<<"Usage:"<<argv[0]<<" input.off/ply [measurements field no]"<<std::endl;
    return 1;
  }

  return 0;
}
