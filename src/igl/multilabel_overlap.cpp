#include <unistd.h>
#include <iostream>
#include <map>

#include "igl/readPLY.h"
#include "igl/writePLY.h"

#include "csv/readCSV.h"
#include "csv/writeCSV.h"

#include "cxxopts.hpp"
#include <algorithm>

template<typename Derived> 
void unique_values_map(
          const Eigen::DenseBase<Derived> & iL, 
          std::map<typename Derived::Scalar,int> & lmap)
{
  int lc=0;
  for(int i=0;i<iL.rows();++i)
    for(int j=0;j<iL.cols();++j)
    {
      auto in = lmap.insert(std::pair<int,int>( iL(i,j), lc )); //TODO: remove BG?
      if(in.second) //new label
        ++lc;
    }
}


template<typename Derived> 
void remap(
          Eigen::DenseBase<Derived> & iL,
          const std::map<typename Derived::Scalar, typename Derived::Scalar> & lmap )
{
  for(int i=0;i<iL.rows();++i)
    for(int j=0;j<iL.cols();++j)
    {
      auto q=lmap.find(iL(i,j));
      iL(i,j) = q->second;
    }
}



template <typename DerivedA, 
          typename DerivedB
          > 
bool extract_column(const std::string &column,
                    const std::vector<std::string> &header, 
                    const Eigen::PlainObjectBase<DerivedA>& in, 
                          Eigen::PlainObjectBase<DerivedB>& out )
{
  size_t _idx=0;
  if(column != "" )
  {
    auto idx = std::find(header.begin(), header.end(), column);
    if(idx == header.end())
      return false;
    _idx = idx-header.begin();
  }
  out = in.col(_idx).template cast<typename DerivedB::Scalar>();

  return true;
}

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0], 
  "Calculate per-vertex label overlap metric\n"
  "\n"
  "Implements vertex-wise Generalized Tanimoto coefficient (GTC)\n"
  "\n"
  " based on :  William R. Crum, Oscar Camara, and Derek L. G. Hill\n"
  "\"Generalized Overlap Measures for Evaluation and Validation in Medical Image Analysis\" \n"
  "IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 25, NO. 11, NOVEMBER 2006\n"
  "http://dx.doi.org/10.1109/TMI.2006.880587"
  );

  options
      .positional_help("<input>")
      .show_positional_help();
  
  options.add_options()
    ("v,verbose", "Verbose output",      cxxopts::value<bool>()->default_value("false"))
    ("c,column", "Column name, default use first available",
                                         cxxopts::value<std::string>()->default_value(""))
    ("l,list",   "list of input files",  cxxopts::value<std::string>())
    ("o,output", "Output file (csv) ",   cxxopts::value<std::string>())
    ("clobber",  "Clobber output file ", cxxopts::value<bool>()->default_value("false"))
    ("input",    "Input files (.csv,.ply) ",  cxxopts::value<std::vector<std::string>>())
    ("help", "Print help") ;
  
  options.parse_positional({"input"});
  auto par = options.parse(argc, argv);
  std::vector<std::string> inputs;

  if(  par.count("output") && 
      (par.count("input")>0 || par.count("list")>0) )
  {
    if ( !par["clobber"].as<bool>() && 
         !access( par["output"].as<std::string>().c_str(), F_OK)) {
      std::cerr << par["output"].as<std::string>()<<" Exists!"<<std::endl;
      return 1;
    }

    if(par.count("input")>0)
    {
      inputs = par["input"].as< std::vector<std::string> >();
    } else {
      std::ifstream in(par["list"].as<std::string>());
      if(!in)
      {
        std::cerr<<"Can't read list file:"<<par["list"].as<std::string>()<<std::endl;
        return 1;
      }
      std::string ln;
      while(std::getline(in, ln))
      {
        if(ln.size() > 0)
            inputs.push_back(ln);
      }
    }
    std::vector< Eigen::VectorXi > input_labels;
    int nfiles=inputs.size();

    for(const auto &f: inputs)
    {
      Eigen::VectorXi _l;
      Eigen::MatrixXd _D;
      std::vector<std::string> _header;

      if(!igl::readCSV(f, _D, _header))
      {
        std::cerr<<"Error reding csv:"<< f <<std::endl;
        return 1;
      }
      if(!extract_column(par["column"].as<std::string>(), _header, _D,_l))
      {
        std::cerr<<"Error finding column:"<<par["column"].as<std::string>()<<" in:"<<f<<std::endl;
        return 1;
      }
      input_labels.push_back(_l);
      if(!input_labels.empty())
      {
        if(_l.size()!=input_labels[0].size())
        {
          std::cerr<<"Inconsistent number of rows in:" << f <<std::endl;
          return 1;
        }
      }
    }

    Eigen::MatrixXi iL(input_labels[0].size(),nfiles);
    for(int i=0;i<nfiles;++i)
    {
      iL.col(i) = input_labels[i];
    }

    std::map<int,int> lmap, ilmap;
    unique_values_map(iL, lmap);

    if(par["verbose"].as<bool>())
      std::cout<<"Found:"<<lmap.size()<<std::endl;

    for(auto i:lmap)
      ilmap.insert(std::pair<int,int>(i.second, i.first));
    
    // remap 
    remap(iL, lmap);
    Eigen::VectorXi ml  = Eigen::VectorXi::Zero(iL.rows());
    Eigen::VectorXd ovl = Eigen::VectorXd::Zero(iL.rows());

    double tU=0.0;
    double tI=0.0;
    //TODO: recode this to use vectorized operations
    for(int i=0; i<iL.rows(); ++i)
    {
      Eigen::VectorXd hist = Eigen::VectorXd::Zero(lmap.size());

      //calculate histogram
      int max_count=0;
      for(int j=0; j<nfiles; ++j)
      {
        int _l    = iL(i, j);
        hist(_l) += 1.0;

        if( hist(_l) > max_count )
        {
          max_count = hist(_l);
          ml(i) = _l;
        }
      }

      //discrete version of the formula 4 from http://dx.doi.org/10.1109/TMI.2006.880587
      double I=0.0;
      double U=0.0;

      for(int j=0; j<hist.size(); ++j) {
        // group-wise union
        U += nfiles*(nfiles-1)/2.0 - (nfiles-hist(j))*(nfiles-hist(j)-1)/2.0;

        // group-wise intersection: val[j]+val[j]-1+.....1 
        I += hist(j)*(hist(j)-1)/2;
      }
      tU+=U;//TODO add label weights
      tI+=I;

      if(U>0.0)
        ovl(i) = I/U;
    }

    std::cout.precision(10);
    if(par["verbose"].as<bool>())
    {
      std::cout<<"Global GTC :";
    }
    std::cout<<tI/tU<<std::endl;

    std::vector<std::string> header{"majority", "overlap"};

    Eigen::MatrixXd oD(iL.rows(),2);
    // map back 
    remap(ml,ilmap);

    oD << ml.cast<double>() , ovl ; // , mc.cast<double>()

    igl::writeCSV(par["output"].as<std::string>(), oD, header); 

  } else {
    std::cerr << options.help({"", "Group"}) << std::endl;
    return 1;
  }
  return 0;
}
