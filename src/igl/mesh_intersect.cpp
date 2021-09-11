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

        std::cout << "FData:     " << FD.rows() << "x"<< FD.cols() << std::endl;
        std::cout<<"headerF:";

        for(auto h:headerF)
          std::cout<<h<<"\t";
        std::cout << std::endl;
      }

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
      for(int i=0; i<F.rows(); ++i)
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

      using AABBTree = igl::AABB<Eigen::MatrixXd,3>;
      // using AABB for self intersection check

      AABBTree tree;
      tree.init(V,F);
      
      std::cout<<"Checking for self-intersections"<<std::endl;
      Eigen::VectorXi intersect=Eigen::VectorXi::Zero(F.rows());

      tick = std::chrono::steady_clock::now();
      for(int i=0; i<F.rows(); ++i)
      {
        int _inter=-1;
        int depth=0;
        if( intersect(i) ) continue;
        using BBOX=Eigen::AlignedBox<AABBTree::Scalar,3>;

        BBOX tri_box;

        for(int j=0;j<3;++j)
          tri_box.extend( V.row( F(i,j) ).transpose() );

        if(i<10)
          std::cout << i<<"("
                    << tri_box.corner(BBOX::CornerType::TopLeftCeil).transpose()
                    << "):(" 
                    << tri_box.corner(BBOX::CornerType::BottomRightFloor).transpose() 
                    << ") ";

        Eigen::Matrix<double, 3,3, Eigen::ColMajor> tri1;
        tri1 << V.row(F(i,0)),
                V.row(F(i,1)),
                V.row(F(i,2));


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

            int coplanar=0;
            Eigen::Matrix<double,2,3,Eigen::ColMajor> edge;
            Eigen::Matrix<double,3,3,Eigen::ColMajor> tri2;

            tri2 << V.row(F(t.m_primitive,0)),
                    V.row(F(t.m_primitive,1)),
                    V.row(F(t.m_primitive,2));

            if(tri_tri_intersection_test_3d(
              tri1.data(),tri1.data()+3,tri1.data()+6,
              tri2.data(),tri2.data()+3,tri2.data()+6,
              &coplanar,
              edge.data(),edge.data()+3))
            {
              //if(!coplanar)
              //{
                _inter=t.m_primitive;
                intersect(i)=1;
                intersect(t.m_primitive)=1;
                depth=d;
                return true;
              //}
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
        if(i<10)
          std::cout<<_inter<<" "<<depth<<std::endl;
      }
      duration=std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count();
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
