#include <iostream>

#include "igl/adjacency_matrix.h"
#include "igl/vertex_components.h"
#include "igl/readPLY.h"
#include "igl/writePLY.h"


int main(int argc, char *argv[])
{
  if(argc>2)
  {
    Eigen::MatrixXd V,N,UV,D;
    Eigen::MatrixXi F,E;
    std::vector<std::string> header;

    double alpha=0.05;

    if(igl::readPLY(argv[1], V, F, E, N, UV, D, header))
    {
      size_t k = 5;

      if(argc>3) k=atoi(argv[3]);

      std::cout << "Vertices: " << V.rows() << "x"<< V.cols() << std::endl;
      std::cout << "Faces:    " << F.rows() << "x"<< F.cols() << std::endl;
      std::cout << "Data:     " << D.rows() << "x"<< D.cols() << std::endl;
      std::cout << "Edges:    " << E.rows() << "x"<< E.cols() << std::endl;
      std::cout<<"header:";

      for(auto h:header)
        std::cout<<h<<"\t";
      std::cout << std::endl;

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

      std::cout << "Connected components:" << counts.rows() << std::endl;

      for(typename Eigen::VectorXi::Scalar i=0; i<counts.rows(); ++i)
      {
        Eigen::MatrixXd V1,N1,UV1,D1;
        Eigen::MatrixXi F1,E1;

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
        V1 = Eigen::MatrixXd::NullaryExpr(counts(i), V.cols(), [&] (Eigen::Index row, Eigen::Index col) { return V(idx_v.idx[row], col);} );

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
        sprintf(tmp, argv[2], i);
        std::cout << tmp << std::endl;
        std::cout << "  Vertices: " << V1.rows() << "x"<< V1.cols() << std::endl;
        std::cout << "  Faces:    " << F1.rows() << "x"<< F1.cols() << std::endl;
        std::cout << "  Data:     " << D1.rows() << "x"<< D1.cols() << std::endl;
        std::cout << "  Edges:    " << E1.rows() << "x"<< E1.cols() << std::endl;


        igl::writePLY(tmp, V1, F1, E1, N1, UV1, D1, header );
      }
      
       
    } else {
      std::cerr<<"Error reding ply:"<<argv[1]<<std::endl;
      return 1;
    }
  } else {
    std::cerr<<"Usage:"<<argv[0]<<" input.ply output_base_%d.ply"<<std::endl;
    return 1;
  }
  return 0;
}
