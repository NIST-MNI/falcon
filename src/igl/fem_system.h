#ifndef __FEM_SYSTEM_H__
#define __FEM_SYSTEM_H__

#include <string>
#include <chrono>
#include <map>
#include <set>


// Eigen
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/LU>

// Sparse solvers
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>

// MSH IO
#include <igl/readMSH.h>
#include <igl/writeMSH.h>

#include "tri_to_tet.h"

#define USE_SVD 0
#include "mesh_util.h"

// openmp stuff
#ifdef USE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


template<class Interpolator, class Image> void sample_vector_field_v2(
    typename Interpolator::Pointer interpolator,
    const Eigen::Matrix<double,4,4> &coil_mat,
    typename Image::Pointer image,
    const Eigen::MatrixXd &coords,
    Eigen::MatrixXd &field,
    const Eigen::VectorXd & default_value,
    bool flip
    )
{
    using VoxelType = typename  Image::PixelType;
    assert( field.cols() == image->GetNumberOfComponentsPerPixel() );
    Eigen::Matrix<double,4,4> inv_coil_mat = coil_mat.inverse(); // .inverse();
    Eigen::Matrix<double,3,3> vect_rotate=coil_mat.topLeftCorner<3,3>();

    if(flip) {
        Eigen::Matrix<double,3,3> flip_mx;
        flip_mx<< -1,0,0,
                   0,-1,0,
                   0,0,1;
        vect_rotate=flip_mx*flip_mx*flip_mx;
    }

    size_t out=0;
    for(size_t i=0; i<coords.rows(); ++i)
    {
        // vertex coordinate, not needed?
        Eigen::Vector3d vertex = inv_coil_mat.topLeftCorner<3,3>()*coords.row(i).transpose()+inv_coil_mat.topRightCorner<3,1>();

        // convert to ITK point
        itk::Point<double,3> itkvertex(vertex.data());
        // apply inverse transform , to move it into coil space
        itk::ContinuousIndex<double, 3 > img_idx;

        //TODO: figure out if we need to adapt for LPS->RAS system here
        
        if( image->TransformPhysicalPointToContinuousIndex(itkvertex, img_idx) )
        {
            using  _MatrixType = Eigen::Matrix<typename Image::InternalPixelType, 3 , 1 >;
            using  _MapType = Eigen::Map<const _MatrixType>;
            VoxelType vox = interpolator->EvaluateAtContinuousIndex(img_idx);
            Eigen::Vector3d fld=_MapType( vox.GetDataPointer() , 3, 1 ).template cast<double>();

            field.row(i) = vect_rotate * fld;
        } else {
            // out of bounds
            ++out;
            field.row(i) = default_value;
        }
    }
    std::cout<<"Out:"<< out <<" of "<< coords.rows()<<std::endl;
}

template<class Image> void xyz_to_index(
    typename Image::Pointer image,
    const Eigen::MatrixXd &X,
    Eigen::MatrixXd &C
    )
{
    using PointValueType=typename Image::PointValueType;
    //assert(C.rows()==X.rows());
    assert(C.cols()==3);
    C.resize(X.rows(),3);
    for(int i=0; i<X.rows(); ++i)
    {
        // convert to ITK point
        typename Image::PointType itkvertex;
        Eigen::Map< Eigen::Matrix<PointValueType,3,1> >(itkvertex.GetDataPointer(),3,1) = X.row(i);
        typename itk::ContinuousIndex<PointValueType, 3> img_idx;

        // TODO: check if we are inside image here?
        image->TransformPhysicalPointToContinuousIndex(itkvertex, img_idx);

        C.row(i) = Eigen::Map<const Eigen::Vector3d>( img_idx.GetDataPointer() , 3, 1 );
    }
}

template<typename Image,
        typename Derived1,
        typename XT> 
void index_to_xyz(
    typename Image::Pointer image,
    const Eigen::MatrixBase<Derived1> &C,
    Eigen::Matrix<XT,-1,-1> &X
    )
{
    using PointValueType=typename Image::PointValueType;

    assert(C.cols()==3);
    X.resize(C.rows(),3);

    for(int i=0; i<C.rows(); ++i)
    {
        // convert to ITK point
        typename Image::PointType itkvertex;
        typename itk::ContinuousIndex<PointValueType, 3> img_idx;

        Eigen::Map< Eigen::Matrix<PointValueType,3,1> >(img_idx.GetDataPointer(),3,1) = C.row(i).template cast<PointValueType>();

        // TODO: check if we are inside image here?
        image->TransformContinuousIndexToPhysicalPoint(img_idx, itkvertex);

        X.row(i) = Eigen::Map<const Eigen::Matrix<PointValueType,3,1> >(itkvertex.GetDataPointer(),3,1).template cast<XT>();
    }
}

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




template<class Image> 
void define_vec_image_bbox(
    const Eigen::MatrixXd &X,
    typename Image::Pointer &image,
    int ncomp=3,
    double step=2.0,
    double border=20.0)
{
    assert(X.cols()==3);
    // create empty image (filled with zeros)
    Eigen::VectorXd bbox_bottom=X.colwise().minCoeff();
    Eigen::VectorXd bbox_top=X.colwise().maxCoeff();

    image=Image::New();

    typename Image::IndexType start;
    typename Image::SizeType size;
    typename Image::PointType org;
    typename Image::SpacingType spc;
    typename Image::RegionType region;
    typename Image::DirectionType identity;

    typename Image::PixelType default_pixel;
    default_pixel.SetSize(ncomp);
    default_pixel.Fill(0.0);

    identity.SetIdentity();
    spc.Fill(step);
    
    for(int j=0;j<3;j++)
    {
        // add border to make zero-padding
        size[j]  = ceil( (bbox_top[j]-bbox_bottom[j]+border*2)/step );
        org[j]   = bbox_bottom[j]-border;
        start[j] = 0;
    }
    region.SetSize  (size);
    region.SetIndex (start);

    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spc );
    image->SetOrigin( org );
    image->SetDirection(identity);
    image->SetNumberOfComponentsPerPixel(ncomp);
}

template<class Image> void define_image_bbox(
    const Eigen::MatrixXd &X,
    typename Image::Pointer &image,
    double step=2.0,
    double border=20.0)
{
    assert(X.cols()==3);
    // create empty image (filled with zeros)
    Eigen::VectorXd bbox_bottom=X.colwise().minCoeff();
    Eigen::VectorXd bbox_top=X.colwise().maxCoeff();

    image = Image::New();

    typename Image::IndexType start;
    typename Image::SizeType size;
    typename Image::PointType org;
    typename Image::SpacingType spc;
    typename Image::RegionType region;
    typename Image::DirectionType identity;

    identity.SetIdentity();
    spc.Fill(step);
    
    for(int j=0;j<3;j++)
    {
        // add border to make zero-padding
        size[j]  = ceil( (bbox_top[j]-bbox_bottom[j]+border*2)/step );
        org[j]   = bbox_bottom[j]-border;
        start[j] = 0;
    }
    region.SetSize  (size);
    region.SetIndex (start);

    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spc );
    image->SetOrigin( org );
    image->SetDirection(identity);
}




template<class ImageRef,class Image> void define_vec_image_byref(
    ImageRef* ref,
    typename Image::Pointer &image,
    int ncomp=3)
{
    // create empty image (filled with zeros)

    image = Image::New();

    typename ImageRef::IndexType start;
    typename ImageRef::SizeType size;
    typename ImageRef::PointType org;
    typename ImageRef::SpacingType spc;
    typename ImageRef::RegionType region;
    typename ImageRef::DirectionType dir;

    // typename Image::PixelType default_pixel;
    // default_pixel.SetSize(ncomp);
    // default_pixel.Fill(0.0);

    region=ref->GetLargestPossibleRegion ();
    spc=ref->GetSpacing( );
    org=ref->GetOrigin( );
    dir=ref->GetDirection();

    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spc );
    image->SetOrigin( org );
    image->SetDirection( dir );
    image->SetNumberOfComponentsPerPixel(ncomp);
}

template<class ImageRef,class Image> void define_image_byref(
    ImageRef* ref,
    typename Image::Pointer &image)
{
    image=Image::New();

    typename ImageRef::IndexType start;
    typename ImageRef::SizeType size;
    typename ImageRef::PointType org;
    typename ImageRef::SpacingType spc;
    typename ImageRef::RegionType region;
    typename ImageRef::DirectionType dir;

    region=ref->GetLargestPossibleRegion ();
    spc=ref->GetSpacing( );
    org=ref->GetOrigin( );
    dir=ref->GetDirection();

    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spc );
    image->SetOrigin( org );
    image->SetDirection( dir );
}


// normalize voxels in image by norm , where it's !=0
template<class Image> void normalize_vec_image(
    typename Image::Pointer image,
    typename Image::Pointer norm )
{
    using ImageIterator=typename itk::ImageRegionIterator<Image>;
    using PixelType=typename Image::PixelType; //will be VariableVector
    int ncomp=image->GetNumberOfComponentsPerPixel();
    ImageIterator it1(image,image->GetLargestPossibleRegion());
    ImageIterator it2(norm,norm->GetLargestPossibleRegion());
    for(it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd()&&!it2.IsAtEnd(); ++it1,++it2)
    {

        PixelType n=it2.Get();
        PixelType o=it1.Get();
        for(int i=0;i<ncomp;++i)
        {
            if(n[i]>0.0)
            {
                o[i]/=n[i];
            }
        }
        it1.Set(o);
    }
}

template<class Image> void normalize_image(
    typename Image::Pointer image,
    typename Image::Pointer norm )
{
    using ImageIterator=typename itk::ImageRegionIterator<Image>;
    using PixelType=typename Image::PixelType; //will be VariableVector
    int ncomp=image->GetNumberOfComponentsPerPixel();
    ImageIterator it1(image,image->GetLargestPossibleRegion());
    ImageIterator it2(norm,norm->GetLargestPossibleRegion());

    for(it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd()&&!it2.IsAtEnd(); ++it1,++it2)
    {

        PixelType n=it2.Get();
        PixelType o=it1.Get();

        for(int i=0;i<ncomp;++i)
        {
            if(n[i]>0.0)
            {
                o[i]/=n[i];
            }
        }
        if(n>0.0)
            it1.Set(o/n);
        else
            it1.Set(o);
    }
}



template<class Image1,class Image2> 
typename Image2::Pointer convert_4d_image_type(typename Image1::Pointer img)
{
    using Image3D=itk::Image<double,3>;
    using ExtractFilterType = typename itk::ExtractImageFilter<Image1, Image3D>;
    using ComposeFilterType = typename itk::ComposeImageFilter<Image3D, Image2>;

	using OrientFilterType=itk::OrientImageFilter<Image3D, Image3D> ;

    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    typename ComposeFilterType::Pointer composeFilter = ComposeFilterType::New();

    // TODO: needed for NIFTI input ?
	typename OrientFilterType::Pointer orienter = OrientFilterType::New();
    // 

    extractFilter->SetDirectionCollapseToSubmatrix();
    typename Image1::RegionType inputRegion = img->GetBufferedRegion();
    typename Image1::SizeType   size = inputRegion.GetSize();

    assert(size[3]==3); // 3D vector
    size[3] = 0; // we extract along t direction
    typename Image1::IndexType start = inputRegion.GetIndex();
    extractFilter->SetInput(img);

    for(int i=0;i<3;i++)
    {
        typename Image1::RegionType desiredRegion;
        start[3]=i;
        desiredRegion.SetSize(size);  
        desiredRegion.SetIndex(start);

        extractFilter->SetExtractionRegion(desiredRegion);

    	// orienter->SetInput(extractFilter->GetOutput());
        // orienter->UseImageDirectionOff();
        // orienter->SetGivenCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI);
        // orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI);
    	// orienter->Update();
        // Image3D::Pointer slice=orienter->GetOutput();

        extractFilter->Update();
        Image3D::Pointer slice=extractFilter->GetOutput();
        slice->DisconnectPipeline();

        composeFilter->SetInput(i, slice);
    }
    composeFilter->Update();
    return composeFilter->GetOutput();
}


template<class Image> void mesh_to_image(
    const Eigen::MatrixXd &C,
    const Eigen::MatrixXi &Tet,
    const std::vector<Eigen::Matrix3d> &IDX_to_BAR,
    const Eigen::MatrixXd &field,
    typename Image::Pointer image,
    typename Image::Pointer norm
)
{
    // linerly interpolate field from volumetric mesh
    // to a regularly sampled volume

    assert( IDX_to_BAR.size()==Tet.rows() );
    assert( C.cols()==3 );
    double eps = 1e-9;
    typename Image::SizeType image_size = image->GetBufferedRegion().GetSize();
 
    // here C is already continious index, rather then real XYZ coordinates
    using VoxelType = typename  Image::PixelType;
    assert(field.cols() == image->GetNumberOfComponentsPerPixel());

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
                        typename Image::IndexType voxelIndex = {{x, y, z}};
                        using  _InternalPixelType = typename Image::InternalPixelType;
                        using  _PixelType = typename Image::PixelType;
                        using  _MatrixType = Eigen::Matrix<_InternalPixelType, -1 , 1 >;
                        using  _MapTypeC  = Eigen::Map<const _MatrixType>;
                        using  _MapType  = Eigen::Map<_MatrixType>;

                        Eigen::VectorXd interpolated_value = tet_fields*bar;
                        // summarize fields
                        interpolated_value += _MapTypeC( image->GetPixel( voxelIndex ).GetDataPointer(), field.cols(), 1); 
                        _PixelType out_pixel(interpolated_value.data(), field.cols());
                        image->SetPixel(voxelIndex, out_pixel );

                        // add the normalization factor
                        Eigen::VectorXd norm_value(field.cols());norm_value.setOnes();
                        norm_value += _MapTypeC( norm->GetPixel( voxelIndex ).GetDataPointer(), field.cols(), 1); 
                        _PixelType out_pixel_norm(norm_value.data(), field.cols());
                        norm->SetPixel(voxelIndex, out_pixel_norm );
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
    typename Image::Pointer image
)
{
    // set voxels inside each tetrahedra to the same value
    // here C is already continious index, rather then real XYZ coordinates
    assert( IDX_to_BAR.size() == Tet.rows() );
    assert( C.cols() == 3 );
    assert( Tet.cols() == 4 );
    assert( tetField.rows() == Tet.rows());
    assert( tetField.cols() == image->GetNumberOfComponentsPerPixel());

    using Scalar=typename Derived::Scalar;

    double eps = 1e-9;
    typename Image::SizeType image_size = image->GetBufferedRegion().GetSize();
 
    using VoxelType = typename  Image::PixelType;

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
                        typename Image::IndexType voxelIndex = {{x, y, z}};
                        using  _InternalPixelType = typename Image::InternalPixelType;
                        using  _PixelType = typename Image::PixelType;
                        using  _MatrixType = Eigen::Matrix<_InternalPixelType, -1 , 1 >;
                        using  _MapTypeC  = Eigen::Map<const _MatrixType>;
                        using  _MapType  = Eigen::Map<_MatrixType>;

                        _MatrixType out = tetField.row(i).template cast<_InternalPixelType>();
                        _PixelType out_pixel(out.data(), tetField.cols());

                        image->SetPixel(voxelIndex, out_pixel );
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
    typename Image::Pointer            image,
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
    typename  Image::RegionType region;
    region = image->GetBufferedRegion(); 
    typename  Image::IndexType i1 = region.GetIndex();
    typename  Image::IndexType i2 = region.GetUpperIndex();

    using VoxelType = typename  Image::PixelType;

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
    using  _PixelType = typename Image::PixelType;
    
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
                typename Image::IndexType voxelIndex = {{ (Index)std::floor(idx(0)+0.5), (Index)std::floor(idx(1)+0.5), (Index)std::floor(idx(2)+0.5)}};
                _PixelType pix = image->GetPixel(voxelIndex);
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
                            typename Image::IndexType voxelIndex = {{x, y, z}};

                            _PixelType pix = image->GetPixel(voxelIndex);
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

template<class Image, 
         typename DerivedC,
         typename DerivedF> 
void sample_image_nn(
    const Eigen::MatrixBase<DerivedC> &C,
    typename Image::Pointer image,
    Eigen::MatrixBase<DerivedF> &vField,
    typename DerivedF::Scalar default_value=0
)
{
    // simple nearest neigbor sampler
    // here C is already continious index, rather then real XYZ coordinates
    assert( C.cols() == 3 );
    assert( vField.rows() == C.rows());
    assert( vField.cols() == 1);

    using Scalar = typename DerivedF::Scalar;

    typename Image::SizeType image_size = image->GetBufferedRegion().GetSize();
 
    using VoxelType = typename  Image::PixelType;

    for(size_t i=0; i<C.rows(); ++i)
    {
        Eigen::Vector3d idx;
        typename Image::IndexType voxelIndex = {{std::floor(C(i,0)+0.5), std::floor(C(i,1)+0.5), std::floor(C(i,2)+0.5)}};
        if(voxelIndex[0]<image_size[0] && voxelIndex[1]<image_size[1] && voxelIndex[2]<image_size[2] &&
           voxelIndex[0]>=0 && voxelIndex[1]>=0 && voxelIndex[2]>=0)
            vField(i) = image->GetPixel(voxelIndex);
        else
            vField(i) = default_value;
    }
}


template< class Interpolator, 
          class Image,
          typename DerivedC,
          typename DerivedF > 
void sample_image(
    Interpolator*                     interpolator,
    Image*                            image,
    const Eigen::MatrixBase<DerivedC> &C,
    Eigen::PlainObjectBase<DerivedF>  &vField,
    typename DerivedF::Scalar         default_value=0
    )
{
    assert( vField.rows()==C.rows() );
    assert( vField.cols()==1 );

    using VoxelType = typename Image::PixelType;
    using Scalar = typename DerivedF::Scalar;
    using Coord  = typename Interpolator::CoordRepType;
    typename Image::RegionType region = image->GetLargestPossibleRegion(); 

    typename Image::IndexType i1 = region.GetIndex();
    typename Image::IndexType i2 = region.GetUpperIndex();

    for(size_t i=0; i<C.rows(); ++i)
    {
        Eigen::Matrix<Coord,1,3> _ci = C.row(i).template cast<Coord>();
        // convert to ITK point
        typename Interpolator::ContinuousIndexType  ci( _ci.data() );

        if( _ci[0] <i2[0] && _ci[1] <i2[1] && _ci[2] <i2[2] &&
            _ci[0]>=i1[0] && _ci[1]>=i1[1] && _ci[2]>=i1[2] )
            vField(i) = static_cast<Scalar>(interpolator->EvaluateAtContinuousIndex(ci));
        else
            vField(i) = default_value;
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

void print_image_info(itk::ImageBase<3>* ref);


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
