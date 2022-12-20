#ifndef MINC_VOLUME_H__
#define MINC_VOLUME_H__

#include <iostream>

#include "minc2-simple.h"


struct minc_volume {
  Eigen::Matrix3d      dir;//direction cosines
  Eigen::Matrix3d      idir;//direction cosines

  Eigen::RowVector3d  start;//(_dims[0].start, _dims[1].start, _dims[2].start);
  Eigen::RowVector3d  step;//(_dims[0].step,  _dims[1].step,  _dims[2].step);
  Eigen::RowVector3i  dims;
  Eigen::RowVector3i  stride;
  Eigen::Array<double,-1,-1,Eigen::RowMajor> volume;

  int                 n_comp;

  Eigen::RowVector3d voxel_to_world(const Eigen::RowVector3d& ijk) const 
  {
    return ((ijk.array()*step.array()).matrix()+start) * dir;
  }

  Eigen::RowVector3d world_to_voxel(const Eigen::RowVector3d& xyz) const
  {
    return ((xyz*idir)-start).array()/step.array();
  }

  double sample_nn(const Eigen::RowVector3d& ijk) const
  {
    return this->volume( (ijk.array().round().cast<int>()*this->stride.array()).sum(),0 );
  }

  double set_voxel(const Eigen::RowVector3d& ijk,double val)
  {
    return this->volume( (ijk.array().round().cast<int>()*this->stride.array()).sum() )=val;
  }

  const Eigen::RowVectorXd sample_nn_vec(const Eigen::RowVector3d& ijk) const
  {
    return this->volume.row( (ijk.array().round().cast<int>()*this->stride.array()).sum() );
  }

  const Eigen::RowVectorXd set_voxel_vec(const Eigen::RowVector3d& ijk,const Eigen::RowVectorXd& val)
  {
    return this->volume.row((ijk.array().round().cast<int>()*this->stride.array()).sum() )=val;
  }


  double sample_interpolate(const Eigen::RowVector3d& ijk) const
  {
    using Array3i=Eigen::Array<int,1,3>;
    using Array3d=Eigen::Array<double,1,3>;

    Array3i ijk_  = ijk.array().floor().cast<int>();
    Array3d ijk_d = ijk.array()-ijk_.cast<double>();

    Array3i ijk_h(ijk_(0)<(this->dims[0]-1)?ijk_(0)+1:ijk_(0), 
                  ijk_(1)<(this->dims[1]-1)?ijk_(1)+1:ijk_(1), 
                  ijk_(2)<(this->dims[2]-1)?ijk_(2)+1:ijk_(2));
    

    Array3i ijk_0(ijk_h(0),ijk_(1),ijk_(2));
    Array3i ijk_1(ijk_(0),ijk_h(1),ijk_(2));
    Array3i ijk_2(ijk_(0),ijk_(1),ijk_h(2));

    Array3i ijk_3(ijk_h(0),ijk_(1),ijk_h(2));
    Array3i ijk_4(ijk_(0),ijk_h(1),ijk_h(2));
    Array3i ijk_5(ijk_h(0),ijk_h(1),ijk_(2));
    Array3i ijk_6(ijk_h(0),ijk_h(1),ijk_h(2));

    //trilinear intrpolation
    return     (1.0-ijk_d[0])*(1.0-ijk_d[1])*(1.0-ijk_d[2])*this->volume( (ijk_*this->stride.array()).sum() )+
    
                ijk_d[0]*(1.0-ijk_d[1])*(1.0-ijk_d[2])* this->volume( (ijk_0*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*ijk_d[1]*(1.0-ijk_d[2])* this->volume( (ijk_1*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*(1.0-ijk_d[1])*ijk_d[2]* this->volume( (ijk_2*this->stride.array()).sum() )+
                
                ijk_d[0]*(1.0-ijk_d[1])*ijk_d[2]*this->volume( (ijk_3*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*ijk_d[1]*ijk_d[2]*this->volume( (ijk_4*this->stride.array()).sum() )+
                ijk_d[0]*ijk_d[1]*(1.0-ijk_d[2])*this->volume( (ijk_5*this->stride.array()).sum() )+
                
                ijk_d[0]*ijk_d[1]*ijk_d[2]*      this->volume( (ijk_6*this->stride.array()).sum() );
  }

  Eigen::ArrayXd sample_interpolate_vec(const Eigen::RowVector3d& ijk) const
  {
    using Array3i=Eigen::Array<int,1,3>;
    using Array3d=Eigen::Array<double,1,3>;

    Array3i ijk_  = ijk.array().floor().cast<int>();
    Array3d ijk_d = ijk.array()-ijk_.cast<double>();

    Array3i ijk_h(ijk_(0)<(this->dims[0]-1)?ijk_(0)+1:ijk_(0), 
                  ijk_(1)<(this->dims[1]-1)?ijk_(1)+1:ijk_(1), 
                  ijk_(2)<(this->dims[2]-1)?ijk_(2)+1:ijk_(2));
    

    Array3i ijk_0(ijk_h(0),ijk_(1),ijk_(2));
    Array3i ijk_1(ijk_(0),ijk_h(1),ijk_(2));
    Array3i ijk_2(ijk_(0),ijk_(1),ijk_h(2));

    Array3i ijk_3(ijk_h(0),ijk_(1),ijk_h(2));
    Array3i ijk_4(ijk_(0),ijk_h(1),ijk_h(2));
    Array3i ijk_5(ijk_h(0),ijk_h(1),ijk_(2));
    Array3i ijk_6(ijk_h(0),ijk_h(1),ijk_h(2));

    //trilinear intrpolation
    return     (1.0-ijk_d[0])*(1.0-ijk_d[1])*(1.0-ijk_d[2])*this->volume.row( (ijk_*this->stride.array()).sum() )+
    
                ijk_d[0]*(1.0-ijk_d[1])*(1.0-ijk_d[2])* this->volume.row( (ijk_0*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*ijk_d[1]*(1.0-ijk_d[2])* this->volume.row( (ijk_1*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*(1.0-ijk_d[1])*ijk_d[2]* this->volume.row( (ijk_2*this->stride.array()).sum() )+
                
                ijk_d[0]*(1.0-ijk_d[1])*ijk_d[2]*this->volume.row( (ijk_3*this->stride.array()).sum() )+
                (1.0-ijk_d[0])*ijk_d[1]*ijk_d[2]*this->volume.row( (ijk_4*this->stride.array()).sum() )+
                ijk_d[0]*ijk_d[1]*(1.0-ijk_d[2])*this->volume.row( (ijk_5*this->stride.array()).sum() )+
                
                ijk_d[0]*ijk_d[1]*ijk_d[2]*      this->volume.row( (ijk_6*this->stride.array()).sum() );
  }
};


inline bool define_similar(minc_volume& vol, const minc_volume& like )
{
    vol.n_comp=like.n_comp;

    vol.dims = like.dims;
    vol.stride = like.stride;
    vol.dir = like.dir;
    vol.idir = like.idir;

    vol.volume.resize(vol.dims.array().prod(),vol.n_comp);
    
    return true;
}

inline bool load_volume(const char*minc_file, minc_volume& vol, bool initialize_only=false )
{
    struct minc2_dimension * _dims;
    int ndim;
    // load minc file
    minc2_file_handle h = minc2_allocate0();

    if(minc2_open(h, minc_file)!=MINC2_SUCCESS)
    {
        std::cerr << "Can't open " << minc_file  << " for reading" << std::endl;
        minc2_free(h);
        return false;
    }
    minc2_setup_standard_order(h);
    minc2_ndim(h, &ndim);

  
    if(ndim<3||ndim>4)
    {
      std::cerr<<"Can read only 3D or 4D volumes, "<<minc_file<<" is "<<ndim<<std::endl;
      minc2_close(h);
      minc2_free(h);
      return false;
    }

    minc2_get_representation_dimensions(h, &_dims);

    vol.n_comp=1;
    if(ndim==4)
      vol.n_comp=_dims[0].length;

    vol.dims << _dims[ndim-3].length, _dims[ndim-2].length, _dims[ndim-1].length;
    vol.volume.resize(vol.dims.array().prod(),vol.n_comp);
    
    if(initialize_only)
    {
        vol.volume=0.0;
    } else {
        if( minc2_load_complete_volume(h, vol.volume.data(), MINC2_DOUBLE) != MINC2_SUCCESS )
        {
            std::cerr << "Error reading data from minc file" << std::endl;
            minc2_close(h);
            minc2_free(h);
            return false;
        }
    }

    vol.start<<_dims[ndim-3].start, _dims[ndim-2].start, _dims[ndim-1].start;
    vol.step <<_dims[ndim-3].step,  _dims[ndim-2].step,  _dims[ndim-1].step;

    vol.dir << _dims[ndim-3].dir_cos[0] , _dims[ndim-3].dir_cos[1] , _dims[ndim-3].dir_cos[2] ,
               _dims[ndim-2].dir_cos[0] , _dims[ndim-2].dir_cos[1] , _dims[ndim-2].dir_cos[2] ,
               _dims[ndim-1].dir_cos[0] , _dims[ndim-1].dir_cos[1] , _dims[ndim-1].dir_cos[2];

    vol.idir = vol.dir.inverse();
    vol.stride = Eigen::RowVector3i(1, vol.dims(0), vol.dims(0)*vol.dims(1));

    minc2_close(h);
    minc2_free(h);
    
    return true;
}

inline bool save_volume(const char*minc_file, const char *minc_file_ref, minc_volume& vol, bool binarize=true)
{
    struct minc2_dimension * _dims;
    struct minc2_dimension * store_dims;
    int ndim;
    // load minc file
    minc2_file_handle h = minc2_allocate0();
    minc2_file_handle o=minc2_allocate0();

    if(minc2_open(h, minc_file_ref)!=MINC2_SUCCESS)
    {
        std::cerr << "Can't open " << minc_file  << " for reading" << std::endl;
        minc2_free(h);
        return false;
    }

    minc2_ndim(h, &ndim);

    if(ndim!=3)
    {
      std::cerr<<"Can write only 3D , "<<minc_file<<" is "<<ndim<<std::endl;
      minc2_close(h);
      minc2_free(h);
      return false;
    }

    minc2_get_store_dimensions(h,&store_dims);
    minc2_define(o,store_dims,binarize?MINC2_UBYTE:MINC2_USHORT, MINC2_DOUBLE);
    minc2_create(o,minc_file);

    minc2_close(h);
    minc2_free(h);
    minc2_setup_standard_order(o);

    if(minc2_save_complete_volume(o,vol.volume.data(), MINC2_DOUBLE)!=MINC2_SUCCESS)
    {
      std::cerr << "Error writing data to minc file" << std::endl;
      return false;
    }

    minc2_close(o);
    minc2_free(o);
    
    return true;
}


template <typename Derived>
Eigen::Matrix<typename Derived::Scalar,-1,1> sample_values_nn(const minc_volume& vol, const Eigen::MatrixXd &C, typename Derived::Scalar _default=0.0 ) 
// C - ijk coordinates
// O - values
{
  Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
  Eigen::RowVector3d hi(vol.dims[0]-1,vol.dims[1]-1,vol.dims[2]-1);

  return Eigen::VectorXd::NullaryExpr(C.rows(),[&](auto i) {
    return (C.row(i).array()>=lo.array() && C.row(i).array()<hi.array()).all() ? 
      static_cast<typename Derived::Scalar>(vol.sample_nn(C.row(i))) : _default;
    });
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar,-1,1> sample_values(const minc_volume& vol, const Eigen::MatrixXd &C, typename Derived::Scalar _default=0.0 )
// C - ijk coordinates
// O - values
{
  Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
  Eigen::RowVector3d hi(vol.dims[0]-1,vol.dims[1]-1,vol.dims[2]-1);

  return Eigen::Matrix<typename Derived::Scalar,-1,1>::NullaryExpr(C.rows(),[&](auto i) 
    { return (C.row(i).array()>=lo.array() && C.row(i).array()<hi.array()).all()? 
      static_cast<typename Derived::Scalar>(vol.sample_interpolate(C.row(i))):
      _default;
    });
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar,-1,-1> sample_values_vec(const minc_volume& vol, const Eigen::MatrixXd &C, 
  Eigen::Matrix<typename Derived::Scalar,1,-1> _default )
// C - ijk coordinates
// O - values
{
  Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
  Eigen::RowVector3d hi(vol.dims[0]-1,vol.dims[1]-1,vol.dims[2]-1);

  // return Eigen::Matrix<typename Derived::Scalar,-1,-1>::NullaryExpr(C.rows(),vol.n_comp,[&](auto i) 
  //   { return (C.row(i).array()>lo.array() && C.row(i).array()<hi.array()).all()? 
  //     static_cast<typename Derived::Scalar>(vol.sample_interpolate(C.row(i))):
  //     _default;
  //   });
  Eigen::Matrix<typename Derived::Scalar,-1,-1> ret(C.rows(),vol.n_comp);
  for(int i=0;i<C.rows();++i)
  {
    if((C.row(i).array()>lo.array() && C.row(i).array()<hi.array()).all())
      ret.row(i)=vol.sample_interpolate_vec(C.row(i));
    else
      ret.row(i)=_default;
  }
  return ret;
}

#endif //MINC_VOLUME_H__