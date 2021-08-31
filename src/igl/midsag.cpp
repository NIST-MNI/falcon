#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/sum.h>
#include <igl/AABB.h>
#include <igl/signed_distance.h>

#include <igl/writePLY.h>
#include <unistd.h>

#include "cxxopts.hpp"

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
    return this->volume( (ijk.array().round().cast<int>()*this->stride.array()).sum() );
  }

  const Eigen::ArrayXd sample_nn_vec(const Eigen::RowVector3d& ijk) const
  {
    return this->volume.row( (ijk.array().round().cast<int>()*this->stride.array()).sum() );
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

bool load_volume(const char*minc_file, minc_volume& vol )
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

    if( minc2_load_complete_volume(h, vol.volume.data(), MINC2_DOUBLE) != MINC2_SUCCESS )
    {
        std::cerr << "Error reading data from minc file" << std::endl;
        minc2_close(h);
        minc2_free(h);
        return false;
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

bool save_volume(const char*minc_file, const char *minc_file_ref, minc_volume& vol, bool binarize=true)
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


void gen_plane(const Eigen::MatrixXd & edges, double step, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
  // generate uniformly sampled triangular mesh between 4 edges, using only 3 points from edges
  Eigen::RowVector3d origin = edges.row(0);
  Eigen::RowVector3d si = (edges.row(1)-edges.row(0));
  Eigen::RowVector3d sj = (edges.row(2)-edges.row(0));

  int li=std::ceil(si.norm()/step);
  int lj=std::ceil(sj.norm()/step);

  si/=li;
  sj/=lj;

  V.resize((li+1)*(lj+1),3); 
  
  for(int i=0,c=0;i<=li; ++i)
    for(int j=0;  j<=lj; ++j,++c)
    {
      V.row(c) = origin+si*i+sj*j;
    }

  // generate faces
  F.resize((li)*(lj)*2,3);

  for(int i=0,c=0; i<li; ++i)
    for(int j=0; j<lj; ++j)
    {
      F.row(c++) << j+i*(lj+1),  j+1+i*(lj+1),    j+(i+1)*(lj+1);
      F.row(c++) << j+1+i*(lj+1),j+1+(i+1)*(lj+1),j+(i+1)*(lj+1);
    }
}


void sample_values_nn(const minc_volume& vol, const Eigen::MatrixXd &C, Eigen::VectorXd &O ) 
// C - ijk coordinates
// O - values
{
  O=Eigen::VectorXd::NullaryExpr(C.rows(),[&](auto i) {return vol.sample_nn(C.row(i));});
}

template <typename Derived>
Eigen::Matrix<typename Derived::Scalar,-1,1> sample_values(const minc_volume& vol, const Eigen::MatrixXd &C, typename Derived::Scalar _default=0.0 )
// C - ijk coordinates
// O - values
{
  Eigen::RowVector3d lo = Eigen::RowVector3d::Zero();
  Eigen::RowVector3d hi(vol.dims[0]-1,vol.dims[1]-1,vol.dims[2]-1);

  return Eigen::Matrix<typename Derived::Scalar,-1,1>::NullaryExpr(C.rows(),[&](auto i) 
    { return (C.row(i).array()>lo.array() && C.row(i).array()<hi.array()).all()? 
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

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Make mid-sagittal separation");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",         cxxopts::value<bool>()->default_value("false"))
    ("i,input", "Source volume ",           cxxopts::value<std::string>())

    ("g,grad", "distance gradient, 4D volume ",         cxxopts::value<std::string>())

    ("o,output", "Output mesh ",             cxxopts::value<std::string>())
    ("l,left", "Left side of the mask ",     cxxopts::value<std::string>())
    ("r,right", "Right side of the mask ",   cxxopts::value<std::string>())
    ("c,central", "Central strip,(experimental) ",   cxxopts::value<std::string>())

    ("clobber", "Clobber output file ",      cxxopts::value<bool>()->default_value("false"))

    ("step", "Step size",                    cxxopts::value<double>()->default_value("0.5"))
    ("ftol", "Stopping criteria",            cxxopts::value<double>()->default_value("1e-5"))
    ("iter", "Maximum number of iterations", cxxopts::value<int>()->default_value("1000"))

    ("help", "Print help") ;
  
  options.parse_positional({"input","grad"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();

  if( par.count("input") && 
      par.count("grad")  )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }
    minc_volume vol,grad_vol;
    if( !load_volume(par["input"].as<std::string>().c_str(),vol ) )
    {
        std::cerr<<"Error loading:"<<par["input"].as<std::string>().c_str()<<std::endl;
        return 1;
    }
    if( !load_volume(par["grad"].as<std::string>().c_str(),grad_vol ) )
    {
        std::cerr<<"Error loading:"<<par["grad"].as<std::string>().c_str()<<std::endl;
        return 1;
    }

    Eigen::Matrix3d edges(3,3);

    edges << Eigen::RowVector3d(vol.dims[0]/2.0,0.0,0.0),
             Eigen::RowVector3d(vol.dims[0]/2.0,vol.dims[1],0.0),
             Eigen::RowVector3d(vol.dims[0]/2.0,0.0,vol.dims[2]-0.0);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // TODO: find the plane with minimal intersect first ?
    gen_plane(edges, 2.0, V, F);

    // optimize cut surface
    // fake normals
    Eigen::MatrixXd N(V.rows(),3);
    
    int    iter = par["iter"].as<int>();
    double step = par["step"].as<double>();
    double ftol = par["ftol"].as<double>();

    Eigen::MatrixXd grad(V.rows(),3);
    Eigen::VectorXd mask(V.rows());

    // smooth using neighborhood averaging 
    // create smoothing matrix
    Eigen::SparseMatrix<double> nei;
    {
      Eigen::VectorXd Asum;
      igl::adjacency_matrix(F,nei);

      // remove edge vertices (that has less then 5 neighbors)
      // to avoid shrinking towards the center
      igl::sum(nei,1,Asum);

      Eigen::VectorXd amask( ((Asum.array()>4).cast<double>()).matrix());
      nei = nei*amask.asDiagonal();
      
      // create smoothing matrix
      Eigen::SparseMatrix<double> ones(Eigen::VectorXd::Ones(V.rows()).asDiagonal());
      nei += ones;

      igl::sum(nei,1,Asum);

      Asum = Asum.cwiseInverse();
      nei = nei*Asum.asDiagonal();

      nei.makeCompressed();
    }

    Eigen::VectorXd cost=Eigen::VectorXd::Zero(iter);
    int    sum_iter=20;
    double iter_change=1e7;    
    double prev_cost_moving_sum=0.0;
    double center_dist=1.0;

    for(int i=0;i<iter;++i) {
      double cost_moving_sum=0.0;

      igl::per_vertex_normals(V,F,N);

      grad = sample_values_vec<Eigen::VectorXd>(grad_vol, V, Eigen::RowVector3d::Zero());

      mask=sample_values<Eigen::VectorXd>(vol, V);
      cost(i)=mask.mean();

      if(i>sum_iter)
      {
        cost_moving_sum=cost.block(i-sum_iter,0,i,1).mean();
        if(i>(sum_iter+1))
          iter_change = prev_cost_moving_sum - cost_moving_sum;

        prev_cost_moving_sum = cost_moving_sum;
      }
      if(verbose)
        std::cout<<"Iter:"<<i<<" cost fun change:"<<iter_change;

      if(iter_change < ftol)
      {
        if(verbose)
          std::cout<< " converged due to ftol"<<std::endl;
        break;
      }

      // project on Normal
      // HACK
      for(int j=0;j<grad.rows();++j)
        grad.row(j).stableNormalize();

      Eigen::VectorXd update_step=(N.array()*grad.array()).matrix().rowwise().sum(); //calculate projection

      // apply mask (?)
      // mask.array() *= -1.0;
      // mask.array() += 1.0;
      //update_step.array()*=mask.array();

      // TODO: smooth
      update_step *= step;
      Eigen::MatrixXd U=N;

      // Eigen doens't do broadcast?
      U.col(0).array() *= update_step.array();
      U.col(1).array() *= update_step.array();
      U.col(2).array() *= update_step.array();

      if(verbose)
        std::cout<<" Update:"<<U.rowwise().norm().mean()<<std::endl;

      V += U;

      // smooth solution
      V.col(0) = V.col(0).transpose()*nei;
      V.col(1) = V.col(1).transpose()*nei;
      V.col(2) = V.col(2).transpose()*nei;
    }

    if(par.count("left")||par.count("right")||par.count("central"))
    {
      // need to split the voxel mask
      minc_volume left_vol  = vol;
      minc_volume right_vol = vol;
      minc_volume center_vol= vol;

      right_vol.volume.setZero();
      left_vol.volume.setZero();
      center_vol.volume.setZero();

      igl::AABB<Eigen::MatrixXd,3> tree;
      igl::FastWindingNumberBVH fwn_bvh;

      tree.init(V,F);
      Eigen::MatrixXd FN(F.rows(),3);
      igl::per_face_normals(V,F,FN);
      
      if(verbose)
        std::cout<<"Using surface distance for splitting volume..."<<std::endl;
      
      for(int i=0;i<vol.dims(0);++i)
      {
        for(int j=0;j<vol.dims(1);++j)
          for(int k=0;k<vol.dims(2);++k)
          {
            double proj=0.0;
            int idx=i*vol.stride(0)+j*vol.stride(1)+k*vol.stride(2);
            Eigen::RowVector3d ijk(i,j,k);
            Eigen::VectorXi I(1);
            Eigen::VectorXd sqrD(1);
            Eigen::MatrixXd C(1,3);
            
            if(i<(tree.m_box.corner(Eigen::AlignedBox<double,3>::BottomLeftCeil)(0)-center_dist) )
            {
              proj=100.0;
              sqrD(0)=(i-(tree.m_box.corner(Eigen::AlignedBox<double,3>::BottomLeftCeil)(0)));
              sqrD(0)*=sqrD(0);
            } else if (i>(tree.m_box.corner(Eigen::AlignedBox<double,3>::BottomRightCeil)(0)+center_dist) ) {
              proj=-100.0;
              sqrD(0)=(i-(tree.m_box.corner(Eigen::AlignedBox<double,3>::BottomRightCeil)(0)));
              sqrD(0)*=sqrD(0);
            } else {
              tree.squared_distance(V,F,ijk,sqrD,I,C);
              proj=(FN.row(I(0)).array()*(ijk-C).array()).sum();
            }

            // HACK : due to triangle orintation the Normals are pointing left
            if(proj<0.0) 
            {
              right_vol.volume(idx)=vol.volume(idx);
              left_vol.volume(idx)=0;
            } else {
              left_vol.volume(idx)=vol.volume(idx);
              right_vol.volume(idx)=0;
            }
            if(sqrD(0)<(center_dist*center_dist))
            {
              center_vol.volume(idx)=vol.volume(idx);
            } else {
              center_vol.volume(idx)=0.0;
            }
          }
          if(verbose)
            std::cout<<"."<<std::flush;
      }
      if(verbose)
        std::cout<<std::endl;
      
      if(par.count("left"))
        save_volume(par["left"].as<std::string>().c_str(),par["input"].as<std::string>().c_str(),left_vol,true);
      
      if(par.count("right"))
        save_volume(par["right"].as<std::string>().c_str(),par["input"].as<std::string>().c_str(),right_vol,true);

      if(par.count("central"))
        save_volume(par["central"].as<std::string>().c_str(),par["input"].as<std::string>().c_str(),center_vol,true);
    }

    if(par.count("output"))
    {
      // convert from voxel to world coordinates
      Eigen::MatrixXd W(V.rows(),3);
      for(int i=0;i<V.rows();++i)
        W.row(i)=vol.voxel_to_world(V.row(i));

      igl::writePLY(par["output"].as<std::string>(), W, F );
    // igl::writePLY(par["output"].as<std::string>(),V,F, 
    //                       Eigen::MatrixXi(0,0), Eigen::MatrixXd(0,0), Eigen::MatrixXd(0,0), 
    //                       mask , {"mask"},
    //                      {"midsag"});
    }

  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
