#include <unistd.h>

#include <igl/read_triangle_mesh.h>
#include <iostream>

#include "igl/readPLY.h"
#include "igl/writePLY.h"
#include "nanoflann.hpp"

#include "csv/readCSV.h"
#include "csv/writeCSV.h"

#include "cxxopts.hpp"
#include <algorithm>


bool extract_psi_the(const std::vector<std::string> &header, 
                     const Eigen::MatrixXd& D, 
                           Eigen::MatrixXd& out_psi_the )
{
  auto idx_psi=std::find(header.begin(),header.end(),"psi");
  auto idx_the=std::find(header.begin(),header.end(),"the");
  if(idx_psi!=header.end() && idx_the!=header.end())
  {
    size_t _idx_psi=idx_psi-header.begin();
    size_t _idx_the=idx_the-header.begin();

    out_psi_the.resize(D.rows(),2);

    out_psi_the<<D.col(_idx_psi),D.col(_idx_the);
    return true;
  } else {
    return false;
  }
}

void sph_to_xyz(const Eigen::MatrixXd& sph, 
                      Eigen::MatrixXd& xyz )
{

  xyz.resize(sph.rows(), 3);

  xyz << sph.col(1).array().sin() * sph.col(0).array().cos(), 
         sph.col(1).array().sin() * sph.col(0).array().sin(),
         sph.col(1).array().cos() ;

}


template <typename Distance,
          typename DerivedC, 
          typename DerivedD, 
          typename Aggregate
          > 
  void resample_field(
    const Eigen::PlainObjectBase<DerivedC> & coord_src,
    const Eigen::PlainObjectBase<DerivedC> & coord_trg,
    const Eigen::PlainObjectBase<DerivedD> & in_data,
    Eigen::PlainObjectBase<DerivedD>       & out_data,
    Aggregate agg,
    const int knn = 1,
    const int leaf_max_size = 10
    )
{
  typedef nanoflann::KDTreeEigenMatrixAdaptor< typename Eigen::PlainObjectBase<DerivedC>,3, Distance > kd_tree_t;
    kd_tree_t idx(coord_src.cols() , coord_src, leaf_max_size );
    idx.index->buildIndex();

    //Eigen::MatrixXd oD(sph2.rows(),D.cols());
    out_data.resize(coord_trg.rows(), in_data.cols());

    // for each vertex in the referece surface search for the k closest in the source surface
    for(size_t i=0; i<coord_trg.rows(); ++i)
    {
        Eigen::Array< typename Eigen::PlainObjectBase<DerivedC>::Scalar,-1,1 > row_sph = coord_trg.row(i);

        // do a knn search
        const size_t num_results = knn;
        Eigen::Array< size_t,-1,1 > ret_indexes( num_results ) ;
        Eigen::Array< typename Eigen::PlainObjectBase<DerivedD>::Scalar,-1, 1 > out_dists_sqr( num_results );

        nanoflann::KNNResultSet< typename Eigen::PlainObjectBase<DerivedD>::Scalar > resultSet(num_results);
        resultSet.init( ret_indexes.data(), out_dists_sqr.data() );

        // check that .data() is correct
        assert(row_sph.size() == coord_src.cols());
        assert(row_sph.innerStride() == 1);

        idx.index->findNeighbors(resultSet, row_sph.data(), nanoflann::SearchParams());

        out_data.row(i) = agg(
          Eigen::Array<typename DerivedD::Scalar,-1,-1>::NullaryExpr(num_results, in_data.cols(), 
            [&] (Eigen::Index r, Eigen::Index c) { return in_data(ret_indexes(r), c);}),
            out_dists_sqr );
    }

}

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], "Resample field");

  options
      .positional_help("<source> <target>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",     cxxopts::value<bool>()->default_value("false"))
    ("s,source", "Source mesh ",        cxxopts::value<std::string>())
    ("t,target", "Target mesh ",        cxxopts::value<std::string>())
    ("i,input",  "Input  field (csv) ", cxxopts::value<std::string>())
    ("chess",    "Generate input data as chessboard ", cxxopts::value<bool>()->default_value("false"))
    ("o,output", "Output field (csv) ", cxxopts::value<std::string>())
    ("k,knn",    "number of nearest neighbors ", cxxopts::value<int>()->default_value("1"))

    ("nearest",  "Use nearest-neighbor ", cxxopts::value<bool>()->default_value("false"))
    ("weighted", "Use weighted average ", cxxopts::value<bool>()->default_value("false"))
    ("invexp",   "Use inverse-exp weighted average ", cxxopts::value<bool>()->default_value("false"))

    ("majority", "Use majority voting, mostly for labelled data", cxxopts::value<bool>()->default_value("false"))
    ("majority_weighted", "Use weighted majority voting ", cxxopts::value<bool>()->default_value("false"))
    ("majority_invexp",   "Use inverse-exp weighted majority voting ", cxxopts::value<bool>()->default_value("false"))

    ("SO3",      "Use SO3 metric (angular distance)", cxxopts::value<bool>()->default_value("false"))

    ("clobber", "Clobber output file ", cxxopts::value<bool>()->default_value("false"))

    ("help", "Print help") ;
  
  options.parse_positional({"source", "target" });
  auto par = options.parse(argc, argv);

  if( par.count("source") && 
      par.count("target") && 
      par.count("output") )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }

    Eigen::MatrixXd V1,N1,UV1,D1;
    Eigen::MatrixXi F1,E1;
    std::vector<std::string> header1;

    Eigen::MatrixXd V2,N2,UV2,D2;
    Eigen::MatrixXi F2,E2;
    std::vector<std::string> header2;

    Eigen::MatrixXd D;
    std::vector<std::string> header;

    if(! igl::readPLY(par["source"].as<std::string>(), V1, F1, E1, N1, UV1, D1, header1))
    {
      std::cerr<<"Error reding ply:"<<par["source"].as<std::string>()<<std::endl;
      return 1;
    }

    if(! igl::readPLY(par["target"].as<std::string>(), V2, F2, E2, N2, UV2, D2, header2) )
    {
      std::cerr<<"Error reding ply:"<<par["target"].as<std::string>()<<std::endl;
      return 1;
    }

    if(par["verbose"].as<bool>())
    {
      std::cout << "Source Mesh 1:" << argv[1] << std::endl;
      std::cout << " Vertices: " << V1.rows() << "x"<< V1.cols() << std::endl;
      std::cout << " Faces:    " << F1.rows() << "x"<< F1.cols() << std::endl;
      std::cout << " Data:     " << D1.rows() << "x"<< D1.cols() << std::endl;
      std::cout << " Header:   ";
      for(auto i: header1)
          std::cout<<i<<"\t";
      std::cout<<std::endl;

      std::cout << "Reference Mesh 2:" << argv[2] << std::endl;
      std::cout << " Vertices: " << V2.rows() << "x"<< V2.cols() << std::endl;
      std::cout << " Faces:    " << F2.rows() << "x"<< F2.cols() << std::endl;
      std::cout << " Data:     " << D2.rows() << "x"<< D2.cols() << std::endl;
      std::cout << " Header:   ";
      for(auto i: header2)
          std::cout<<i<<"\t";
      std::cout<<std::endl;
    }
  
    Eigen::MatrixXd  PT1, PT2;
    Eigen::MatrixXd  sph1, sph2;

    // convert spherical coordinates to x,y,z on sphere 
    if( extract_psi_the(header1, D1, PT1) &&
        extract_psi_the(header2, D2, PT2) )
    {
      int knn=10;
      sph_to_xyz(PT1, sph1);
      sph_to_xyz(PT2, sph2);

      if( par["chess"].as<bool>() ) {
        // generate chessboard patern at 10deg
        double freq = 10.0;
        D = ((PT1.array()*freq).sin()*0.5+0.5).round().rowwise().prod().matrix() ;
        header = {"chess"};
      } else if(par.count("input") ) {
        if(!igl::readCSV(par["input"].as<std::string>(), D, header))
        {
          std::cerr<<"Error reding csv:"<<par["input"].as<std::string>()<<std::endl;
          return 1;
        }
        if(D.rows()!=sph1.rows())
        {
          std::cerr<<"Unexpected number of rows:"<< D.rows() <<" instead of:"<< sph1.rows() <<std::endl;
          return 1;
        }
      } else {
        // use fields from the input mesh
        std::vector<int> out_idx;
        for(int i=0;i<header1.size();++i)
        {
          if( header1[i]!="psi" && 
              header1[i]!="the" )
            out_idx.push_back(i);
        }
        header.resize(out_idx.size());
        D.resize(PT1.rows(),out_idx.size());
        for(int i=0;i<out_idx.size();++i)
        {
          header[i]=header1[out_idx[i]];
          D.col(i)=PT1.col(out_idx[i]);
        }
      }
      if(D.cols()==0)
      {
        std::cerr << "Empty input data, probably something is wrong "<<std::endl;
        return 1;
      }

      Eigen::MatrixXd oD; 

      // single nearest neighbor, the simplest case, ignore dist completely
      auto nearest_neighbor_agg = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
           return data.row(0);
      };

      // simple inv-distance-weighted average
      auto weighted_avg = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
          Eigen::ArrayXd w = (dist+1e-6).cwiseInverse();
          w/=w.sum();
          return (data*w).colwise().sum();
      };

      auto invexp_weighted_avg = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
          Eigen::ArrayXd w = Eigen::exp(-1.0*dist);
          w/=w.sum();
          return (data*w).colwise().sum();
      };
      
      // simple majority voting, works for single dimension data only
      auto majority_voting = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
        Eigen::VectorXd r=Eigen::VectorXd::Zero(data.cols());
        for(int c=0;c<data.cols();++c) 
        {
            std::unordered_map<double, int> cnts;
            int max_vote_count=0;
            double max_vote=data(0,c);
            for(int i=0;i<data.rows();++i)
            {
              auto t = cnts.insert( std::make_pair(data(i,c),0) );
              t.first->second++;
              if(t.first->second>max_vote_count)
              {
                max_vote_count=t.first->second;
                max_vote=data(i,c);
              }
            }
            r[c]=max_vote;
        }
        return r;
      };

      // weighted majority voting, works for single dimension data only
      auto majority_weighted_voting = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
          Eigen::ArrayXd w = (dist+1e-6).cwiseInverse();
          w/=w.sum();

          Eigen::VectorXd r=Eigen::VectorXd::Zero(data.cols());
          for(int c=0;c<data.cols();++c) 
          {
            std::unordered_map<double,double> cnts;
            double max_vote_weight=0;
            double max_vote=data(0,c);
            for(int i=0;i<data.rows();++i)
            {
              auto t = cnts.insert( std::make_pair(data(i,c), 0.0) );
              t.first->second+=w(i,0);
              if(t.first->second>max_vote_weight)
              {
                max_vote_weight=t.first->second;
                max_vote=data(i,c);
              }
            }
            r[c]=max_vote;
          }
          return r;
      };

      // inverse-exponential weighted majority voting, works for single dimension data only
      auto majority_invexp_weighted_voting = [] (const Eigen::ArrayXd& data, const Eigen::ArrayXd& dist ) -> Eigen::VectorXd {
          Eigen::ArrayXd w = Eigen::exp(-1.0*dist);
          w/=w.sum();
          Eigen::VectorXd r=Eigen::VectorXd::Zero(data.cols());

          for(int c=0;c<data.cols();++c) 
          {
            std::unordered_map<double,double> cnts;
            double max_vote_weight=0;
            double max_vote=data(0,c);
            for(int i=0;i<data.rows();++i)
            {
              auto t = cnts.insert( std::make_pair(data(i,c), 0.0) );
              t.first->second+=w(i,0);
              if(t.first->second>max_vote_weight)
              {
                max_vote_weight=t.first->second;
                max_vote=data(i,c);
              }
            }
            r[c]=max_vote;
          }
          return r;
      };

      // TODO: make templated version with double and int
      // TODO: allow specifying different algorithm for different columns
      if( par["SO3"].as<bool>() )
      {
        if( par["nearest"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, nearest_neighbor_agg, knn);
        } else if( par["weighted"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, weighted_avg, knn);
        } else if( par["invexp"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, invexp_weighted_avg, knn);
        } else if( par["majority"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, majority_voting, knn);
        } else if( par["majority_weighted"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, majority_weighted_voting, knn);
        } else if( par["majority_invexp"].as<bool>() ) {
          resample_field<nanoflann::metric_SO3>(sph1, sph2, D, oD, majority_invexp_weighted_voting, knn);
        }        
        else {
          std::cerr<< " Need to specify --weighted --invexp or --nearest" << std::endl;
          return 1;
        }
      } else {
        if( par["nearest"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, nearest_neighbor_agg, knn);
        } else if( par["weighted"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, weighted_avg, knn);
        } else if( par["invexp"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, invexp_weighted_avg, knn);
        } else if( par["majority"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, majority_voting, knn);
        } else if( par["majority_weighted"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, majority_weighted_voting, knn);
        } else if( par["majority_invexp"].as<bool>() ) {
          resample_field<nanoflann::metric_L2>(sph1, sph2, D, oD, majority_invexp_weighted_voting, knn);
        } else {
          std::cerr<< " Need to specify --weighted --invexp or --nearest" << std::endl;
          return 1;
        }
      }

      // aggregate information from k neighbors: average, make weighted average, take median, mode (for labels ) or highest probability (for labels)
      if(igl::check_ext(par["output"].as<std::string>(),"ply"))
      {
        Eigen::MatrixXd ND2(V2.rows(), D2.cols()+oD.cols() );
        ND2 << D2,oD;

        for(auto const &h:header)
          header2.push_back(h);

        igl::writePLY(par["output"].as<std::string>(), V2, F2, E2, N2, UV2, ND2, header2 );
      }
      else
      {
        igl::writeCSV(par["output"].as<std::string>(), oD, header); 
      }
    }
  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
