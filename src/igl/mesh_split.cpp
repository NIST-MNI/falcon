#include <iostream>
#include <unistd.h>

#include "igl/adjacency_matrix.h"
#include "igl/vertex_components.h"
#include "igl/readPLY.h"
#include "igl/writePLY.h"

#include "cxxopts.hpp"

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Resample field");

  options
      .positional_help("<source> <output printf format>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("s,source",  "Source mesh ",   cxxopts::value<std::string>())
    ("o,output",  "Output printf format ", cxxopts::value<std::string>())

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

      Eigen::SparseMatrix<typename Eigen::MatrixXi::Scalar> A;
      Eigen::VectorXi C, counts;

      igl::adjacency_matrix(F, A );
      igl::vertex_components(A, C, counts );

      // use component of the 0th vertex of the face
      Eigen::VectorXi FC = Eigen::VectorXi::NullaryExpr( F.rows(), [&] (Eigen::Index row) { return C( F( row, 0) ); } );

      // use component of the 0th vertex of the edge
      Eigen::VectorXi EC;

      if( E.rows()>0 )
       EC = Eigen::VectorXi::NullaryExpr( E.rows(), [&] (Eigen::Index row) {return C( E( row, 0)); } );

      if(par["verbose"].as<bool>()) std::cout << "Connected components:" << counts.rows() << std::endl;

      for(typename Eigen::VectorXi::Scalar i=0; i<counts.rows(); ++i)
      {
        Eigen::MatrixXd V1,N1,UV1,D1,FD1,ED1;
        Eigen::MatrixXi F1,E1;

        std::vector<std::string> header1F,header1E;

        struct collect_idx {
          Eigen::VectorXi idx;
          Eigen::VectorXi rev_idx;
          Eigen::Index val;
          Eigen::Index next;

          // called for the first coefficient
          void init(const int& value, Eigen::Index r, Eigen::Index c)
          {
            if(value==val)
            {
              idx[next] = r;
              rev_idx[r] = next;
              ++next;
            }
          }

          // called for all other coefficients
          void operator() (const int& value, Eigen::Index r, Eigen::Index c)
          {
            init(value,r,c);
          }

          collect_idx(Eigen::Index max_count,Eigen::Index v): 
            idx(max_count), rev_idx(max_count), val(v), next(0)
          { 
            rev_idx.setConstant(-1); 
          }
        };

        // build indexes
        collect_idx idx_v(V.rows(),i);
        C.visit(idx_v);

        collect_idx idx_f(F.rows(),i);
        FC.visit(idx_f);

        // remap vertices
        V1 = Eigen::MatrixXd::NullaryExpr(counts(i), V.cols(), [&] (Eigen::Index row, Eigen::Index col)     { return V(idx_v.idx[row], col);} );

        // remap faces
        F1 = Eigen::MatrixXi::NullaryExpr(idx_f.next, F.cols(), [&] (Eigen::Index row, Eigen::Index col)    { return  idx_v.rev_idx(F(idx_f.idx[row], col));} );

        if(E.rows()>0)
        {
          collect_idx idx_e(E.rows(),i);
          EC.visit(idx_e);
          E1 = Eigen::MatrixXi::NullaryExpr(idx_e.next, E.cols(), [&] (Eigen::Index row, Eigen::Index col) {return idx_v.rev_idx(E(idx_e.idx[row], col));} );
        } else {
          std::cout<< std::endl;
        }

        if(N.rows()==V.rows())
          N1 =  Eigen::MatrixXd::NullaryExpr(counts(i), N.cols(), [&] (Eigen::Index row, Eigen::Index col) {return N(idx_v.idx[row], col);} );
        if(UV.rows()==V.rows())
          UV1 = Eigen::MatrixXd::NullaryExpr(counts(i), UV.cols(), [&] (Eigen::Index row, Eigen::Index col) {return UV(idx_v.idx[row], col);} );
        if(D.rows()==V.rows())
          D1 = Eigen::MatrixXd::NullaryExpr(counts(i), D.cols(), [&] (Eigen::Index row, Eigen::Index col) {return D(idx_v.idx[row], col);} );

        char tmp[1024];
        snprintf(tmp, 1023, par["output"].as<std::string>().c_str(), i);
        if(par["verbose"].as<bool>()) {
          std::cout << tmp << std::endl;
          std::cout << "  Vertices: " << V1.rows() << "x"<< V1.cols() << std::endl;
          std::cout << "  Faces:    " << F1.rows() << "x"<< F1.cols() << std::endl;
          std::cout << "  Data:     " << D1.rows() << "x"<< D1.cols() << std::endl;
          std::cout << "  Edges:    " << E1.rows() << "x"<< E1.cols() << std::endl;
        }

        if(!igl::writePLY(tmp, V1, F1, E1, N1, UV1, D1, header, FD1, header1F, ED1, header1E, comments, igl::FileEncoding::Binary ))
        {
          std::cerr<<"Error writing:"<<tmp<<std::endl;
          return 1;
        }
      }
             
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
