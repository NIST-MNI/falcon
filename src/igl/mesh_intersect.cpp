#include <iostream>
#include <unistd.h>
#include <chrono>
#include <numeric>
#include "Guigue2003_tri_tri_intersect.h"

#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/AABB.h>
#include <igl/tinyply.h>

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

        std::cout << "FData:     " << FD.rows() << "x"<< FD.cols() << std::endl;
        std::cout<<"headerF:";

        for(auto h:headerF)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

      #if 0
      // intersect with a triangle along X-Y plane at Z=mean
      Eigen::RowVectorXd lo=V.colwise().minCoeff();
      Eigen::RowVectorXd hi=V.colwise().maxCoeff();
      Eigen::RowVectorXd cnt=(lo+hi)/2.0;

      std::cout<<"Lo:"<<lo<<" Hi:"<<hi<<std::endl;

      Eigen::Matrix<double,3,3,Eigen::RowMajor> sect_tri;

      sect_tri<< lo(0),lo(1),cnt(2),
                 lo(0),2*hi(1)-lo(1),cnt(2),
                 2*hi(0)-lo(0),lo(1),cnt(2);


      std::vector<Eigen::Matrix<double,1,3,Eigen::RowMajor> > edges;

      int coplanar_ctr=0;
      // go over all triangles and check intersections
      auto tick = std::chrono::steady_clock::now();
      for(int i=0; i<F.rows(); ++i)
      {
        bool coplanar=false;

        Eigen::Matrix<double,1,3,Eigen::RowMajor> edge1,edge2;
        
        if(igl::tri_tri_intersection_test_3d(
           V.row(F(i,0)),  V.row(F(i,1)),  V.row(F(i,2)),
           sect_tri.row(0),sect_tri.row(1),sect_tri.row(2),
           coplanar,
           edge1, edge2))
        {
          if(!coplanar)
          {

            edges.push_back(edge1);
            edges.push_back(edge2);
          }
          else
            coplanar_ctr++; // HACK: add all three edges?
        }
      }
      double duration=std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count();
      std::cout<<"Found "<<edges.size()/2<<" intersections"<<" and "<< coplanar_ctr << " coplanar"<<std::endl;
      std::cout<<"Took:"<<duration<<" ms"<<std::endl;
      #endif //

      #if 0
      // save to a file
      // using .ply file
      {
        tinyply::PlyFile file;
        std::vector<double> _v;
        std::vector<int>    _e;
        int ec=0;
        //TODO: remove duplicated vertex

        for(auto e:edges)
        { 
          for(int j=0;j<3;++j)
            _v.push_back(e(j));
          _e.push_back(ec++);
        }

        file.add_properties_to_element("vertex",{"x","y","z"},
          tinyply::Type::FLOAT64, edges.size(),
          reinterpret_cast<uint8_t*>(&_v[0]), tinyply::Type::INVALID,0);

        file.add_properties_to_element("edge",{"vertex1","vertex2"},
          tinyply::Type::INT32, edges.size()/2,
          reinterpret_cast<uint8_t*>(&_e[0]), tinyply::Type::INVALID,0);

        std::filebuf fb;
        fb.open("section.ply",std::ios::out|std::ios::binary);

        std::ostream stream(&fb);
         
        file.write(stream,true);
      }
      #endif


      auto tick = std::chrono::steady_clock::now();

      using AABBTree = igl::AABB<Eigen::MatrixXd,3>;
      // using AABB for self intersection check

      AABBTree tree;
      tree.init(V,F);
      
      std::cout<<"Checking for self-intersections"<<std::endl;
      Eigen::VectorXi intersect=Eigen::VectorXi::Zero(F.rows());


      for(int i=0; i<F.rows(); ++i)
      {
        int _inter=-1;
        int depth=0;
        if( intersect(i) ) continue;
        using BBOX=Eigen::AlignedBox<AABBTree::Scalar,3>;

        BBOX tri_box;

        for(int j=0;j<3;++j)
          tri_box.extend( V.row( F(i,j) ).transpose() );

        auto adjacent_faces = [](auto A,auto B)->bool {
          for(int i=0;i<3;++i)
            for(int j=0;j<3;j++)
            {
              if(A(i)==B(j))
                return true;
            }
          return false;
        };
        
        // find leaf nodes containing intersecting tri_box
        std::function<bool(const AABBTree &,int)> check_intersect;
        check_intersect = [&](const AABBTree &t,int d) -> bool
        {
          if(t.is_leaf()) //check for the actual intersection
          {
            if(t.m_primitive==i) //itself
               return false;
            if(adjacent_faces(F.row(i), F.row(t.m_primitive)) )
               return false;

            bool coplanar=false;
            Eigen::Matrix<double,1,3,Eigen::RowMajor> edge1,edge2;

            if(igl::tri_tri_intersection_test_3d(
              V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),
              V.row(F(t.m_primitive,0)),V.row(F(t.m_primitive,1)),V.row(F(t.m_primitive,2)),
              coplanar,
              edge1,edge2))
            {
              if(!coplanar)
              {
                _inter=t.m_primitive;
                intersect(i)=1;
                intersect(t.m_primitive)=1;
                depth=d;
                return true;
              }
            }
          } else {
            if(t.m_box.intersects(tri_box)) 
              return check_intersect(*t.m_left,d+1) ||
                     check_intersect(*t.m_right,d+1);
          }
          return false;
        };

        // actual check
        bool rr=check_intersect(tree,0);
      }

      double duration=std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count();
      std::cout<<"Found "<<intersect.sum()<<" intersections"<<std::endl;
      std::cout<<"Took:"<<duration<<" ms"<<std::endl;
      // saving as FACE data
      headerF.push_back("red");
      headerF.push_back("green");
      headerF.push_back("blue");
      Eigen::MatrixXd FD1(F.rows(),FD.cols()+3);
      intersect.array() += 1.0;
      intersect.array() *= 127;

      FD1<<FD,intersect.cast<double>(),intersect.cast<double>(),intersect.cast<double>();

      igl::writePLY(par["output"].as<std::string>(), V, F, E, N, UV, D, header, FD1,headerF,ED,headerE,comments,igl::FileEncoding::Binary );      
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
