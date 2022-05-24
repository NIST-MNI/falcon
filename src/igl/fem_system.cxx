#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fem_system.h"

// Sparse solvers
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>

//#include <igl/readMSH.h>
//#include <igl/writeMSH.h>
//#include "tri_to_tet.h"

#define USE_SVD 0
//#include "mesh_util.h"



// print out gradient only when MESH_DEBUG is defined
const int __debug_tet=2684719; //2341252; // 2152426; //2341252; // 2152426; //;




void compute_gradient(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G)
{
// from Simnibs:
// ''' G calculates the gradient of a function in each tetrahedra
// The way it works: The operator has 2 parts
// G = T^{-1}A
// A is a projection matrix
// A = [-1, 1, 0, 0]
//     [-1, 0, 1, 0]
//     [-1, 0, 0, 1]
// And T is the transfomation to baricentric coordinates
//
// here the output G is 4*Tet,3 with colum correspoding to x,y,z

    // WARNING: output is th * 4 x XYZ ( simnibs output is th x XYZ x 4 ? )
    G.resize(Tet.rows()*4, 3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,4> A;
        A << -1, 1, 0, 0,
             -1, 0, 1, 0,
             -1, 0, 0, 1;

        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,4> _G;

        T << ( X.row(Tet(i,1)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,2)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,3)) - X.row(Tet(i,0)) );

        Eigen::PartialPivLU<Eigen::Matrix<double,3,3> > lu(T);
        _G = lu.solve(A);
        //_G = T.colPivHouseholderQr().solve(A);

        // TODO: check if this is what we want
        G.template block<4,3>(i*4, 0) = _G.transpose();
        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == __debug_tet ) {
            std::cout<<"Debug:Tetrahedral gradient:" << debug_tet << std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A="      << A << std::endl;
            std::cout<<"T="      << T << std::endl;
            std::cout<<"G="      <<_G << std::endl;
        }
        #endif //MESH_DEBUG
    }       
}

void compute_gradient_svd(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G,
                            double epsilon)
{
// from Simnibs:
// ''' G calculates the gradient of a function in each tetrahedra
// The way it works: The operator has 2 parts
// G = T^{-1}A
// A is a projection matrix
// A = [-1, 1, 0, 0]
//     [-1, 0, 1, 0]
//     [-1, 0, 0, 1]
// And T is the transfomation to baricentric coordinates
//
// here the output G is 4*Tet,3 with colum correspoding to x,y,z

    // WARNING: output is th * 4 x XYZ ( simnibs output is th x XYZ x 4 ? )
    G.resize(Tet.rows()*4, 3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,4> A;
        A << -1, 1, 0, 0,
             -1, 0, 1, 0,
             -1, 0, 0, 1;

        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,4> _G;

        T << ( X.row(Tet(i,1)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,2)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,3)) - X.row(Tet(i,0)) );

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(epsilon);
        _G = svd.solve(A);

        // TODO: check if this is what we want
        G.template block<4,3>(i*4, 0) = _G.transpose();
        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == __debug_tet ) {
            std::cout<<"Debug:Tetrahedral gradient:" << debug_tet << std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A="      << A << std::endl;
            std::cout<<"T="      << T << std::endl;
            std::cout<<"G="      <<_G << std::endl;
        }
        #endif //MESH_DEBUG
    }       
}



void compute_field_gradient(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                      const Eigen::VectorXd &Field,
                      Eigen::MatrixXd &G)
{
    // output is th * 4 * XYZ 
    G.resize(Tet.rows(),3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,1>  A;
        Eigen::Matrix<double,3,1> _G;

        T << (X.row(Tet(i,1)) - X.row(Tet(i,0))),
             (X.row(Tet(i,2)) - X.row(Tet(i,0))),
             (X.row(Tet(i,3)) - X.row(Tet(i,0)));

        A << Field(Tet(i,1)) - Field(Tet(i,0)),
             Field(Tet(i,2)) - Field(Tet(i,0)),
             Field(Tet(i,3)) - Field(Tet(i,0));

        //_G = T.colPivHouseholderQr().solve(A);
        Eigen::PartialPivLU<Eigen::Matrix<double,3,3> > lu(T);
        _G = lu.solve(A);

        // TODO: check if this is what we want
        G.row(i) = _G;

        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == debug_tet ) {
            std::cout<<"Debug:field gradient:"<<debug_tet<<std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A=" <<A<<std::endl;
            std::cout<<"T=" <<T<<std::endl;
            std::cout<<"G=" <<_G<<std::endl;
        }
        #endif        
    }       
}


void compute_field_gradient_svd(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                      const Eigen::VectorXd &Field,
                      Eigen::MatrixXd &G,
                      double epsilon)
{
    // output is th * 4 * XYZ 
    G.resize(Tet.rows(),3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,1>  A;
        Eigen::Matrix<double,3,1> _G;

        T << (X.row(Tet(i,1)) - X.row(Tet(i,0))),
             (X.row(Tet(i,2)) - X.row(Tet(i,0))),
             (X.row(Tet(i,3)) - X.row(Tet(i,0)));

        A << Field(Tet(i,1)) - Field(Tet(i,0)),
             Field(Tet(i,2)) - Field(Tet(i,0)),
             Field(Tet(i,3)) - Field(Tet(i,0));

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(epsilon);
        _G = svd.solve(A);

        // TODO: check if this is what we want
        G.row(i) = _G;

        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == debug_tet ) {
            std::cout<<"Debug:field gradient:"<<debug_tet<<std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A=" <<A<<std::endl;
            std::cout<<"T=" <<T<<std::endl;
            std::cout<<"G=" <<_G<<std::endl;
        }
        #endif        
    }       
}


void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxb, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyb,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzb)
{
    int total_voxels=si*sj*sk;
    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> Dx_,Dy_,Dz_;
    std::vector<triplet> Dx_2,Dy_2,Dz_2;
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };
    Dx_.reserve(total_voxels*2);
    Dy_.reserve(total_voxels*2);
    Dz_.reserve(total_voxels*2);

    Dx_2.reserve(total_voxels*2);
    Dy_2.reserve(total_voxels*2);
    Dz_2.reserve(total_voxels*2);

    for(Eigen::Index k=0;k<sk;++k) 
        for(Eigen::Index j=0;j<sj;++j)
            for(Eigen::Index i=0;i<si;++i)
            {
                if(k>0)      Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), -1));
                Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));

                if(j>0)      Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), -1));
                Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));
                
                if(i>0)      Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), -1));
                Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));


                Dx_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(k<(sk-1)) Dx_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 1));

                Dy_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(j<(sj-1)) Dy_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 1));
                
                Dz_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(i<(si-1)) Dz_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 1));
            }
    Dxf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxb.setFromTriplets(Dx_.begin(), Dx_.end());
    Dyb.setFromTriplets(Dy_.begin(), Dy_.end());
    Dzb.setFromTriplets(Dz_.begin(), Dz_.end());

    Dxf.setFromTriplets(Dx_2.begin(), Dx_2.end());
    Dyf.setFromTriplets(Dy_2.begin(), Dy_2.end());
    Dzf.setFromTriplets(Dz_2.begin(), Dz_2.end());

    Dxf.makeCompressed();
    Dyf.makeCompressed();
    Dzf.makeCompressed();
    Dxb.makeCompressed();
    Dyb.makeCompressed();
    Dzb.makeCompressed();
}

void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzc)
{
    // centered version
    int total_voxels=si*sj*sk;
    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> Dx_,Dy_,Dz_;
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };
    Dx_.reserve(total_voxels*2);
    Dy_.reserve(total_voxels*2);
    Dz_.reserve(total_voxels*2);

    for(Eigen::Index k=0;k<sk;++k) 
        for(Eigen::Index j=0;j<sj;++j)
            for(Eigen::Index i=0;i<si;++i)
            {
                if(k>0)      Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), -0.5));
                if(k<(sk-1)) Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 0.5));

                if(j>0)      Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), -0.5));
                if(j<(sj-1)) Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 0.5));
                
                if(i>0)      Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), -0.5));
                if(i<(si-1)) Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 0.5));

            }
    Dxc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxc.setFromTriplets(Dx_.begin(), Dx_.end());
    Dyc.setFromTriplets(Dy_.begin(), Dy_.end());
    Dzc.setFromTriplets(Dz_.begin(), Dz_.end());

    Dxc.makeCompressed();
    Dyc.makeCompressed();
    Dzc.makeCompressed();
}


void create_laplacian_matrix(int si,int sj,int sk, Eigen::SparseMatrix<double, Eigen::RowMajor> &A, bool s27)
{
    // now we can determine 2nd order differential equation
    int total_voxels=si*sj*sk;

    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> A_;
    
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };

    if(s27)
    {
        auto ijk_hit = [&](int i,int j,int k) -> bool
        {
            return i>=0  && j>=0 && k>=0 && 
                    i<si && j<sj && k<sk;
        };

        auto ijk_27stencil = [&](int i,int j,int k) -> double
        {
            switch(abs(i)+abs(j)+abs(k))
            {
                case 0:return -88.0/26;
                case 1:return   6.0/26;
                case 2:return   3.0/26;
                case 3:return   2.0/26;
                default: return 0.0;
            }

        };

        A_.reserve(total_voxels*27);
        // using 27 point stencil 
        // see https://en.wikipedia.org/wiki/Discrete_Laplace_operator
        for(Eigen::Index k=0;k<sk;++k) 
            for(Eigen::Index j=0;j<sj;++j)
                for(Eigen::Index i=0;i<si;++i)
                {
                    // fill out stencil:
                    for(int ii=-1;ii<2;++ii)
                        for(int jj=-1;jj<2;++jj)
                            for(int kk=-1;kk<2;++kk)
                                if(ijk_hit(i+ii,j+jj,k+kk))
                                    A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+ii,j+jj,k+kk),  
                                        ijk_27stencil(ii,jj,kk) ));
                }
    } else {
        A_.reserve(total_voxels*7);
        // using 7 point stencil for now 
        // see https://en.wikipedia.org/wiki/Discrete_Laplace_operator
        for(Eigen::Index k=0;k<sk;++k) 
            for(Eigen::Index j=0;j<sj;++j)
                for(Eigen::Index i=0;i<si;++i)
                {
                    // fill out stencil:
                    A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k),  -6.0));

                    if(i>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), 1.0));
                    if(j>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), 1.0));
                    if(k>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), 1.0));

                    if(i<(si-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 1.0));
                    if(j<(sj-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 1.0));
                    if(k<(sk-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 1.0));
                }
    }
    A=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    A.setFromTriplets(A_.begin(), A_.end());
    A.makeCompressed();
}

