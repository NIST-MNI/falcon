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
    ("m,mesh", "Input mesh ", cxxopts::value<std::string>())
    ("c,csv", "Input csv ", cxxopts::value<std::string>())
    ("help", "Print help")
  ;
  
  options.parse_positional({"mesh", "csv"});
  auto par = options.parse(argc, argv);

  if(par.count("mesh"))
  {
    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header;
    std::vector<std::string> headerF,headerE,comments;

    if(igl::readPLY(par["mesh"].as<std::string>(), V, F, E, N, UV, D, header,
                  FD,headerF, ED,headerE, comments))
    {
      const size_t k = 5;
      static int idx_field = 0;

      if(argc>2)
        idx_field=atoi(argv[2]);

      if(idx_field>=D.cols()||idx_field<0) idx_field=0;

      std::cout << "Vertices: " << V.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;

      if(par.count("csv"))
      {
        Eigen::MatrixXd Dcsv;
        std::vector<std::string> csv_header;
        std::cout<<"Reading from :"<<par["csv"].as<std::string>()<<std::endl;
        igl::readCSV(par["csv"].as<std::string>(), Dcsv, csv_header);
        if(Dcsv.cols()>0)
        {
          std::cout << "CSV Data:     " << Dcsv.rows() << "x"<< Dcsv.cols() << std::endl;

          Eigen::MatrixXd D_(V.rows(), D.cols()+Dcsv.cols());
          D_ << D,Dcsv;
          for(auto const &h:csv_header)
            header.push_back(h);
          D = D_;
        }
      }

      if(!header.empty())
      {
        std::cout<<"Data:";
        for(auto h:header)
          std::cout<<h<<"\t";

        std::cout<<std::endl;
      }


      Eigen::MatrixXd U;
      Eigen::VectorXd S;

      // Plot the mesh
      igl::opengl::glfw::Viewer viewer;
      viewer.data().set_mesh(V, F);
      viewer.data().set_face_based(true);

      if(idx_field<D.cols())
      {
          std::cout<<"Field:"<<idx_field<<" "<<D.col(idx_field).minCoeff()<<" "<<D.col(idx_field).maxCoeff()<<std::endl;

          Eigen::MatrixXd C;
          Eigen::MatrixXd _D=D.col(idx_field);
          igl::colormap(igl::COLOR_MAP_TYPE_JET, _D , _D.minCoeff(), _D.maxCoeff(), C);

          viewer.data().set_colors(C);
      }


      // Attach a custom menu
      igl::opengl::glfw::imgui::ImGuiMenu menu;

      // Customize default menu
      menu.callback_draw_viewer_menu = [&]()
      {
        menu.draw_viewer_menu(); // draw default menu, comment to overwrite

        if(ImGui::Combo("Measure", &idx_field, header))
        {
          Eigen::MatrixXd C;
          Eigen::MatrixXd _D=D.col(idx_field);
          igl::colormap(igl::COLOR_MAP_TYPE_JET, _D , _D.minCoeff(), _D.maxCoeff(), C);
          viewer.data().set_colors(C);
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
