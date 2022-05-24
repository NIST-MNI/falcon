#ifndef __FEM_SYSTEM_H__
#define __FEM_SYSTEM_H__

#include <string>
#include <chrono>
#include <map>
#include <set>

#include "minc_volume.h"

// Eigen
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>


// MSH IO
#include <igl/readMSH.h>
#include <igl/writeMSH.h>

//#include "tri_to_tet.h"

#define USE_SVD 0
//#include "mesh_util.h"

// openmp stuff
// #ifdef USE_OPENMP
// #include <omp.h>
// #else
// #define omp_get_num_threads() 1
// #define omp_get_thread_num() 0
// #define omp_get_max_threads() 1
// #endif


template<typename DerivedX,
         typename DerivedT
        > 
void create_baricentric_matrix(
    const Eigen::MatrixBase<DerivedX> &X,
    const Eigen::MatrixBase<DerivedT> &Tet,
    std::vector<Eigen::Matrix<typename DerivedX::Scalar,3,3> > &XYZ_to_BAR
)
{
    XYZ_to_BAR.reserve(Tet.rows());
    assert(Tet.cols()==4);
    assert(X.cols()==3);

    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<typename DerivedX::Scalar,3,3> bar_to_xyz;

        bar_to_xyz << 
            (X.row(Tet(i,0))- X.row(Tet(i,3))),
            (X.row(Tet(i,1))- X.row(Tet(i,3))),
            (X.row(Tet(i,2))- X.row(Tet(i,3)));
        
        XYZ_to_BAR.push_back( bar_to_xyz.transpose().inverse());
    }
}


template<class Image> void mesh_to_image(
    const Eigen::MatrixXd &C,
    const Eigen::MatrixXi &Tet,
    const std::vector<Eigen::Matrix3d> &IDX_to_BAR,
    const Eigen::MatrixXd &field,
    minc_volume & out_vol,
    minc_volume & out_norm
)
{
    // linerly interpolate field from volumetric mesh
    // to a regularly sampled volume

    assert( IDX_to_BAR.size()==Tet.rows() );
    assert( C.cols()==3 );
    double eps = 1e-9;
    Eigen::RowVector3i  image_size = out_vol.dims;
 
    // here C is already continious index, rather then real XYZ coordinates
    using VoxelType = double;
    assert(field.cols() == out_vol.n_comp);

    // create a bounding box for each tetrahedra
    // using integer voxel coordinates
    Eigen::Matrix<int,Eigen::Dynamic,3> low=
             Eigen::Matrix<int, Eigen::Dynamic, 3>::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (int)std::floor(std::min( { C(Tet(i,0), j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j)  })) ;
                    } );

    Eigen::Matrix<int,Eigen::Dynamic,3> high=
             Eigen::Matrix<int, Eigen::Dynamic, 3>::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (int)std::ceil(std::max( { C(Tet(i,0),j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j) })) ;
                    } );
    
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        const Eigen::Matrix3d& idx_to_bar=IDX_to_BAR[i];
        Eigen::Vector3d idx;
        Eigen::Vector4d bar;
        Eigen::Matrix<double,Eigen::Dynamic,4> tet_fields(field.cols(), 4);

        // initialize by column
        tet_fields << field.row(Tet(i, 0)).transpose(),
                      field.row(Tet(i, 1)).transpose(),
                      field.row(Tet(i, 2)).transpose(),
                      field.row(Tet(i, 3)).transpose();  

        for(int x=low(i,0)-1; x<=high(i,0) && x<image_size[0]; ++x)
        {
            idx(0)=(double)x - C(Tet(i,3),0);
            for(int y=low(i,1)-1; y<=high(i,1) && y<image_size[1]; ++y)
            {
                idx(1)=(double)y - C(Tet(i,3),1);
                for(int z=low(i,2)-1; z<=high(i,2) && z<image_size[2]; ++z)
                {
                    idx(2) = (double)z - C(Tet(i,3),2);
                    bar.head<3>() = idx_to_bar * idx;
                    bar(3) = 1.0 - bar.head<3>().sum();

                    if( ( bar.array() >= -eps     ).all() && 
                        ( bar.array() <=(1.0+eps) ).all() )
                    {
                        Eigen::RowVector3i voxelIndex(x, y, z);

                        Eigen::RowVectorXd interpolated_value = tet_fields*bar;
                        // summarize fields
                        interpolated_value += out_vol.sample_nn_vec( voxelIndex ); 
                        out_vol.set_voxel_vec(voxelIndex,interpolated_value);

                        // add the normalization factor
                        Eigen::RowVectorXd norm_value(field.cols());norm_value.setOnes();
                        norm_value += out_norm.sample_nn_vec( voxelIndex ); 
                        out_norm.set_voxel_vec(voxelIndex, norm_value );
                    }
                }
            }
        }
    }
}


template<class Image, typename Derived> 
void mesh_to_image_nn(
    const Eigen::MatrixXd &C,
    const Eigen::MatrixXi &Tet,
    const std::vector<Eigen::Matrix3d> &IDX_to_BAR,
    const Eigen::PlainObjectBase<Derived> &tetField,
    minc_volume & out_vol
)
{
    // set voxels inside each tetrahedra to the same value
    // here C is already continious index, rather then real XYZ coordinates
    assert( IDX_to_BAR.size() == Tet.rows() );
    assert( C.cols() == 3 );
    assert( Tet.cols() == 4 );
    assert( tetField.rows() == Tet.rows());
    assert( tetField.cols() == out_vol.n_comp);

    using Scalar=double;

    double eps = 1e-9;
    Eigen::RowVector3i  image_size = out_vol.dims;
 
    using VoxelType = double;

    // create a bounding box for each tetrahedra
    // using integer voxel coordinates
    Eigen::Matrix<int,Eigen::Dynamic,3> low=
             Eigen::Matrix<int, Eigen::Dynamic, 3>::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (int)std::floor(std::min( { C(Tet(i,0), j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j)  })) ;
                    } );

    Eigen::Matrix<int,Eigen::Dynamic,3> high=
             Eigen::Matrix<int, Eigen::Dynamic, 3>::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (int)std::ceil(std::max( { C(Tet(i,0),j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j) })) ;
                    } );
    
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        const Eigen::Matrix3d& idx_to_bar = IDX_to_BAR[i];
        Eigen::Vector3d idx;
        Eigen::Vector4d bar;

        // TODO: figure out if it's possible to have seams
        for(int x=low(i,0)-1; x<=high(i,0) && x<image_size[0]; ++x)
        {
            idx(0) = (double)x - C(Tet(i,3),0);
            for(int y=low(i,1)-1; y<=high(i,1) && y<image_size[1]; ++y)
            {
                idx(1) = (double)y - C(Tet(i,3),1);
                for(int z=low(i,2)-1; z<=high(i,2) && z<image_size[2]; ++z)
                {
                    idx(2) = (double)z - C(Tet(i,3),2);
                    bar.head<3>() = idx_to_bar * idx;
                    bar(3) = 1.0 - bar.head<3>().sum();

                    if( ( bar.array() >= -eps     ).all() && 
                        ( bar.array() <=(1.0+eps) ).all() )
                    {
                        Eigen::RowVector3i voxelIndex(x, y, z);
                        out_vol.set_voxel_vec(voxelIndex,tetField.row(i));
                    }
                }
            }
        }
    }
}



template< 
        typename Image,
        typename DerivedC,
        typename DerivedT,
        typename DerivedF
        > 
void image_to_mesh_majority(
    const Eigen::MatrixBase<DerivedC> &C,
    const Eigen::MatrixBase<DerivedT> &Tet,
    const std::vector<Eigen::Matrix<typename DerivedC::Scalar,3,3>> &IDX_to_BAR,
    Eigen::PlainObjectBase<DerivedF>   &tetField,
    const minc_volume & volume,
    int n_cls,
    typename DerivedF::Scalar          default_value=0
)
{
    // sample voxels inside each tetrahedra and get majority vote
    // assign to the the output
    // here C is already continious index, rather then real XYZ coordinates
    assert( IDX_to_BAR.size() == Tet.rows() );
    assert( C.cols() == 3 );
    assert( Tet.cols() == 4 );
    assert( tetField.rows() == Tet.rows());
    //assert( tetField.cols() == image->GetNumberOfComponentsPerPixel());
    assert( tetField.cols() == 1);

    using Index  = typename DerivedT::Scalar;
    using Scalar = typename DerivedC::Scalar;
    using Matrix3x3 = Eigen::Matrix<typename DerivedC::Scalar,3,3>;
    using IndexVec = Eigen::Matrix<Index,Eigen::Dynamic,3>;
    Scalar eps = 1e-9;

    // image index bounds 
    Eigen::RowVector3i  i2 = volume.dims;
    Eigen::RowVector3i  i1(0,0,0);
    using VoxelType = double;

    // create a bounding box for each tetrahedra
    // using integer voxel coordinates
    IndexVec low=
             IndexVec::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (Index)std::floor(std::min( { C(Tet(i,0), j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j)  })) ;
                    } );

    IndexVec high=
             IndexVec::NullaryExpr( Tet.rows(), 3,
                    [&C,&Tet](auto i, auto j)->double {
                        return (Index)std::ceil(std::max( { C(Tet(i,0),j), C(Tet(i,1),j), C(Tet(i,2),j), C(Tet(i,3),j) })) ;
                    } );
    using  _PixelType = double;
    
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::VectorXi cls_counts=Eigen::VectorXi::Zero(n_cls+1);
        const Matrix3x3& idx_to_bar = IDX_to_BAR[i];
        Eigen::Matrix<Scalar,3,1> idx;
        Eigen::Matrix<Scalar,4,1> bar;

        
        if( (high(i,0)-low(i,0)<=1.0) && (high(i,1)-low(i,1)<=1.0) && (high(i,2)-low(i,2)<=1.0) )
        {
            // whole tetrahedra fits inside one voxel
            // so, just sample at barycenter
            idx = (C.row(Tet(i,0))+C.row(Tet(i,1))+C.row(Tet(i,2))+C.row(Tet(i,3)))/4;
            if(idx(0)>=i1[0] && idx(1)>=i1[1] && idx(2)>=i1[2] &&
               idx(0)< i2[0] && idx(1) <i2[2] && idx(2) <i2[2] )
            {
                Eigen::RowVector3i voxelIndex( (int)std::floor(idx(0)+0.5), (int)std::floor(idx(1)+0.5), (int)std::floor(idx(2)+0.5));

                int pix = (int)volume.sample_nn(voxelIndex);
                cls_counts(pix)++;
            }
        } else {
            for(Index x=low(i,0)-1; x<=high(i,0); ++x)
            {   
                if(x<i1[0]||x>=i2[0]) continue;

                idx(0) = (Scalar)x - C(Tet(i,3),0);
                for(Index y=low(i,1)-1; y<=high(i,1); ++y)
                {
                    if(y<i1[1]||y>=i2[1]) continue;

                    idx(1) = (Scalar)y - C(Tet(i,3),1);
                    for(Index z=low(i,2)-1; z<=high(i,2); ++z)
                    {
                        if(z<i1[2]||z>=i2[2]) continue;

                        idx(2) = (Scalar)z - C(Tet(i,3),2);
                        bar.template head<3>() = idx_to_bar * idx;
                        bar(3) = 1.0 - bar.template head<3>().sum();

                        if( ( bar.array() >= -eps     ).all() && 
                            ( bar.array() <=(1.0+eps) ).all() )
                        {
                            Eigen::RowVector3i voxelIndex(x,y,z);
                            int pix = (int)volume.sample_nn(voxelIndex);
                            cls_counts(pix)++;
                        }
                    }
                }
            }
        }
        // determine majority
        if(cls_counts.sum()>0)
        {
            int maj_count = 0;
            int maj_idx = 0;
            maj_count = cls_counts.maxCoeff(&maj_idx);
            // TODO: store the count?
            tetField(i) = maj_idx;
        } else {
            tetField(i) = default_value;
        }
    }
}

void compute_gradient(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G);

void compute_gradient_svd(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G,
                            double epsilon);

void compute_field_gradient(const Eigen::MatrixXd &X,
                            const Eigen::MatrixXi &Tet,
                            const Eigen::VectorXd &Field,
                                  Eigen::MatrixXd &G);

void compute_field_gradient_svd(const Eigen::MatrixXd &X,
                            const Eigen::MatrixXi &Tet,
                            const Eigen::VectorXd &Field,
                                  Eigen::MatrixXd &G,
                                  double epsilon);


void create_laplacian_matrix(int si,int sj,int sk, Eigen::SparseMatrix<double, Eigen::RowMajor> &A, bool s27=false);

// void create_grad_matrixes(int si,int sj,int sk, 
//     Eigen::SparseMatrix<double, Eigen::RowMajor> &Dx, 
//     Eigen::SparseMatrix<double, Eigen::RowMajor> &Dy,
//     Eigen::SparseMatrix<double, Eigen::RowMajor> &Dz);


void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxb, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyb,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzb);


void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzc);




#endif //__FEM_SYSTEM_H__
