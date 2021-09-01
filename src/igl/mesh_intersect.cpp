#include <iostream>
#include <unistd.h>
#include <chrono>
#include <numeric>
#include "Guigue2003_tri_tri_intersect.h"

#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/AABB.h>

#include "cxxopts.hpp"

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Resample field");

  options
      .positional_help("<source> <output>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("s,source",  "Source mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output path ",   cxxopts::value<std::string>())

    ("help", "Print help") ;
  
  options.parse_positional({"source", "output"});
  auto par = options.parse(argc, argv);

  if( par.count("source") && 
      par.count("output") )
  {

    Eigen::MatrixXd V,N,UV,D,FD,ED;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header;
    std::vector<std::string> headerF,headerE;
    std::vector<std::string> comments;

    double alpha=0.05;

    if(igl::readPLY(par["source"].as<std::string>(), V, F, E, N, UV, D, header, FD,headerF,ED,headerE,comments ))
    {
      if(par["verbose"].as<bool>()) {
        std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
        std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
        std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
        std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
        std::cout<<"header:";

        for(auto h:header)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

      // igl::AABB<MatrixXd,3> tree;
      // tree.init(V,T);
      
      // intersect with a triangle along X-Y plane at Z=mean
      Eigen::RowVectorXd lo=V.colwise().minCoeff();
      Eigen::RowVectorXd hi=V.colwise().maxCoeff();
      Eigen::RowVectorXd cnt=(lo+hi)/2.0;

      std::cout<<"Lo:"<<lo<<" Hi:"<<hi<<std::endl;

      Eigen::Matrix<double,3,3,Eigen::ColMajor> sect_tri;

      sect_tri<< lo(0),lo(1),cnt(2),
                 lo(0),2*hi(1)-lo(1),cnt(2),
                 2*hi(0)-lo(0),lo(1),cnt(2);

      std::vector<Eigen::Matrix<double,2,3,Eigen::ColMajor> > edges;
      int coplanar_ctr=0;
      // go over all triangles and check intersections
      auto tick = std::chrono::steady_clock::now();
      for(int i=0;i<F.rows();++i)
      {
        int coplanar=0;
        Eigen::Matrix<double,2,3,Eigen::ColMajor> edge;
        Eigen::Matrix<double,3,3,Eigen::ColMajor> tri;
        tri << V.row(F(i,0)),
               V.row(F(i,1)),
               V.row(F(i,2));
        
        if(tri_tri_intersection_test_3d(
           tri.data(),tri.data()+3,tri.data()+6,
           sect_tri.data(),sect_tri.data()+3,sect_tri.data()+6,&coplanar,edge.data(),edge.data()+3))
        {
          if(!coplanar)
            edges.push_back(edge);
          else
            coplanar_ctr++; // HACK: add all three edges?
        }
      }
      double duration=std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count();
      std::cout<<"Found "<<edges.size()<<" intersections"<<" and "<< coplanar_ctr << " coplanar"<<std::endl;
      std::cout<<"Took:"<<duration<<" ms"<<std::endl;

      // crude self-intersection check
      std::cout<<"Checking for self-intersections"<<std::endl;
      Eigen::VectorXi intersect=Eigen::VectorXi::Zero(F.rows());

      tick = std::chrono::steady_clock::now();
      for(int i=0;i<F.rows();++i)
      {
        if(intersect(i)) continue;
        
        Eigen::Matrix<double,3,3,Eigen::ColMajor> tri1;
        tri1 << V.row(F(i,0)),
                V.row(F(i,1)),
                V.row(F(i,2));
        for(int j=i+1;j<F.rows();++j)
        {
          if(intersect(i)) continue;
          int coplanar=0;
          Eigen::Matrix<double,2,3,Eigen::ColMajor> edge;
          Eigen::Matrix<double,3,3,Eigen::ColMajor> tri2;
          tri2 << V.row(F(i,0)),
                  V.row(F(i,1)),
                  V.row(F(i,2));

          if(tri_tri_intersection_test_3d(
            tri1.data(),tri1.data()+3,tri1.data()+6,
            tri2.data(),tri2.data()+3,tri2.data()+6,&coplanar,edge.data(),edge.data()+3))
          {
            // if(!coplanar)
            //   edges.push_back(edge);
            // else
            //   coplanar_ctr++; // HACK: add all three edges?
            intersect(i)=1;
            intersect(j)=1;
          }
        }
      }
      duration=std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count();
      std::cout<<"Found "<<intersect.sum()<<" intersections"<<std::endl;
      std::cout<<"Took:"<<duration<<" ms"<<std::endl;

      
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
