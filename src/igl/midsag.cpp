#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/sum.h>
#include <igl/AABB.h>
#include <igl/signed_distance.h>

#include <igl/read_triangle_mesh.h>
#include <igl/writePLY.h>
#include <unistd.h>

#include "cxxopts.hpp"
#include "minc_volume.h"

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


int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Make mid-sagittal separation");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",         cxxopts::value<bool>()->default_value("false"))
    ("i,input", "Source volume ",           cxxopts::value<std::string>())
    ("mesh", "Initial mesh, if absent generate plane ",cxxopts::value<std::string>())
    ("g,grad", "distance gradient, 4D volume ",         cxxopts::value<std::string>())
    ("mask", "Non moving mask ",   cxxopts::value<std::string>())

    ("o,output", "Output mesh ",             cxxopts::value<std::string>())
    ("l,left", "Left side of the mask ",     cxxopts::value<std::string>())
    ("r,right", "Right side of the mask ",   cxxopts::value<std::string>())
    ("c,central", "Central strip,(experimental) ",   cxxopts::value<std::string>())


    ("clobber", "Clobber output file ",      cxxopts::value<bool>()->default_value("false"))

    ("step", "Step size",                    cxxopts::value<double>()->default_value("0.5"))
    ("ftol", "Stopping criteria",            cxxopts::value<double>()->default_value("1e-7"))
    ("iter", "Maximum number of iterations", cxxopts::value<int>()->default_value("1000"))
    ("dist", "Cut distance",                 cxxopts::value<double>()->default_value("0.5"))

    ("help", "Print help") ;
  
  options.parse_positional({"input"});
  auto par = options.parse(argc, argv);
  bool verbose=par["verbose"].as<bool>();

  if( par.count("input") )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }
    minc_volume vol, grad_vol,mask_vol;
    if( !load_volume(par["input"].as<std::string>().c_str(),vol ) )
    {
        std::cerr<<"Error loading:"<<par["input"].as<std::string>().c_str()<<std::endl;
        return 1;
    }
    
    if(par.count("grad") && !load_volume(par["grad"].as<std::string>().c_str(),grad_vol ) )
    {
        std::cerr<<"Error loading:"<<par["grad"].as<std::string>().c_str()<<std::endl;
        return 1;
    }

    if(par.count("mask") && !load_volume(par["mask"].as<std::string>().c_str(),mask_vol ) )
    {
        std::cerr<<"Error loading:"<<par["mask"].as<std::string>().c_str()<<std::endl;
        return 1;
    }


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if( par.count("mesh") )
    {
      Eigen::MatrixXd W;
      if(!igl::read_triangle_mesh(par["mesh"].as<std::string>(),W,F))
      { 
        std::cerr<<"Error loading:"<<par["mesh"].as<std::string>().c_str()<<std::endl;
        return 1;
      }
      V.resize(W.rows(),3);
      // convert to voxel coordinates
      for(int i=0;i<V.rows();++i)
        V.row(i)=vol.world_to_voxel(W.row(i));

    } else {
      // TODO: find the plane with minimal intersect first ?
      Eigen::Matrix3d edges(3,3);

      // determine voxel coodinates of 0,0,0
      Eigen::RowVector3d center_ijk=vol.world_to_voxel(Eigen::RowVector3d::Zero());

      edges << Eigen::RowVector3d(center_ijk[0],0.0,0.0),
              Eigen::RowVector3d(center_ijk[0],vol.dims[1],0.0),
              Eigen::RowVector3d(center_ijk[0],0.0,vol.dims[2]); 
      gen_plane(edges, 2.0, V, F);
    }

    int    iter = par["iter"].as<int>();
    double step = par["step"].as<double>();
    double ftol = par["ftol"].as<double>();
    double center_dist=par["dist"].as<double>();

    // A hack, without the gradients just dump surface into output
    if(par.count("grad"))
    {

      // optimize cut surface
      // fake normals
      Eigen::MatrixXd N(V.rows(),3);
      

      Eigen::MatrixXd grad(V.rows(),3);
      Eigen::VectorXd mask_intersect(V.rows());
      Eigen::VectorXd vmask;

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
        
        // do not smooth over masked area
        if(par.count("mask")) {
          vmask = ((sample_values<Eigen::VectorXd>(mask_vol, V).array()<1.0).cast<double>()).matrix();
          nei = nei*vmask.asDiagonal();
        }

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

      for(int i=0;i<iter;++i) {
        double cost_moving_sum=0.0;

        igl::per_vertex_normals(V,F,N);

        grad = sample_values_vec<Eigen::VectorXd>(grad_vol, V, Eigen::RowVector3d::Zero());

        mask_intersect = sample_values<Eigen::VectorXd>(vol, V);
        cost(i)=mask_intersect.mean();

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
        if(par.count("mask"))
          update_step.array()*=vmask.array();

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
    } else {
      std::cout<<"Warning: Can't deform surface without gradients"<<std::endl;
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
            if(sqrD(0)>=(center_dist*center_dist)) 
            {
              if(proj<0.0) 
              {
                right_vol.volume(idx)=vol.volume(idx);
                left_vol.volume(idx)=0;
              } else {
                left_vol.volume(idx)=vol.volume(idx);
                right_vol.volume(idx)=0;
              }
            } else {
              left_vol.volume(idx)=0;
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
    }

  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
