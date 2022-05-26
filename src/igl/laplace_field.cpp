#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#include <string>
#include <chrono>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cxxopts.hpp"

#include "minc_volume.h"
#include "fem_system.h"

// Sparse solvers
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>




int main(int argc,char **argv)
{
    cxxopts::Options options(argv[0], "Run laplace equation solver example");

    options
        .positional_help("<in> <out> ")
        .show_positional_help();
    
    options.add_options()
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
        ("s27", "27 point stencil", cxxopts::value<bool>()->default_value("false"))
        ("grad", "grad operator", cxxopts::value<bool>()->default_value("false"))
        ("in",  "Input volume",  cxxopts::value<std::string>())
        ("out", "Output volume", cxxopts::value<std::string>())
        ("droptol", "Preconditioning tolerance ",           cxxopts::value<double>()->default_value("1e-3"))
        ("h,help", "Print help")
    ;
    
    options.parse_positional({"in","out"});
    auto par = options.parse(argc, argv);

    if(par.count("in") && par.count("out") && !par.count("help"))
    {

#ifdef HAVE_OPENMP
        std::cout<<"Using OpenMP, max number of threads="<<omp_get_max_threads()<<std::endl;
#endif
        
        // ITK setup
        // using PixelType = double;
        // using ImageType = itk::Image<PixelType, 3>;
        // using LabelImageType = itk::Image<unsigned char, 3>;
        auto tick = std::chrono::steady_clock::now();

        minc_volume image;
        if( !load_volume(par["in"].as<std::string>().c_str(),image, false ) )
        {
            std::cerr<<"Error loading:"<<par["in"].as<std::string>().c_str()<<std::endl;
            return 1;
        }
        std::cout<<"image.min="<<image.volume.minCoeff()
                 <<" image.max="<<image.volume.maxCoeff()
                 <<std::endl;

        std::cout<<"Image.dims="<<image.dims(0)<<","<<image.dims(1)<<","<<image.dims(2)<<std::endl;

        minc_volume out_image;
        define_similar(out_image,image);
        out_image.volume = 0.0;        
        
        using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        int total_voxels=image.dims.prod();

        SparseMatrix A;
        if(par["grad"].as<bool>())
        {
            // using gradient operators
            // FOR DEBUGGING ONLY!
            SparseMatrix Dxf,Dyf,Dzf;
            SparseMatrix Dxb,Dyb,Dzb;
            create_grad_matrixes(image.dims(0), image.dims(1), image.dims(2), Dxf,Dxb,Dyf,Dyb,Dzf,Dzb);
            A=Dxf*Dxb+Dyf*Dyb+Dzf*Dzb;
        } else {
            // directly generate laplacian 
            create_laplacian_matrix(image.dims(0), image.dims(1), image.dims(2), A, par["s27"].as<bool>());
        }
        tick = std::chrono::steady_clock::now();
        ////

        std::cout << "Created sparse matrix:" << std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count() << " ms" << std::endl;
        std::cout << "Total number of voxels: "<< total_voxels<<std::endl;
        std::cout << "Number of nonzero elements in A:" << A.nonZeros()<<std::endl;

        auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
        {
            return i+j*image.dims(0)+k*image.dims(0)*image.dims(1);
        };

        tick = std::chrono::steady_clock::now();
        // setting the dirichlet boundary condition
        using triplet = Eigen::Triplet<double>;
        std::vector<triplet> B_;
        
        for(Eigen::Index k=0;k<image.dims(2);++k) 
            for(Eigen::Index j=0;j<image.dims(1);++j)
                for(Eigen::Index i=0;i<image.dims(0);++i)
                {
                    if(image.volume(ijk_to_idx(i,j,k))>0)
                    {
                        double val=0.0;
                        Eigen::Index idx = ijk_to_idx(i,j,k);
                        
                        if(image.volume(ijk_to_idx(i,j,k))>1) val = 1.0;
 
                        B_.push_back(triplet( idx, 0, val ));

                        // set BC:
                        A.row(idx) *= 0.0;
                        A.coeffRef(idx, idx) = 1.0;
                    }
                }
        Eigen::SparseMatrix<double, Eigen::ColMajor> B(total_voxels, 1);
        B.setFromTriplets(B_.begin(), B_.end());
        B.makeCompressed();

        Eigen::VectorXd dB=Eigen::MatrixXd(B);

        std::cout << "Set boundary condition:" << std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count() << " ms" << std::endl;
        std::cout << "Number of nonzero elements in B:"<<B.nonZeros()<<std::endl;

        Eigen::VectorXd solution;

        // solving....
        Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double > > solver; // seem to be the fastest single threaded
        solver.preconditioner().setDroptol( par["droptol"].as<double>() ); // preconditioner works slower, but solution is quite fast after

        tick = std::chrono::steady_clock::now();
        solver.compute(A);
        std::cout << "Solver preconditioner:" << std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count() << " ms" << std::endl;
        tick = std::chrono::steady_clock::now();

        if(solver.info()!= Eigen::Success) {
            std::cerr<<"Solving failed"<<std::endl;
            return -1;
        } else {
            solution = solver.solve(dB); //TODO: use B ?
        }
        std::cout << "Solving:         " << std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count() << " ms" << std::endl;
        std::cout << "#iterations:     " << solver.iterations() << std::endl;
        std::cout << "estimated error: " << solver.error()      << std::endl;

        for(Eigen::Index k=0;k<(image.dims(2));++k) 
            for(Eigen::Index j=0;j<(image.dims(1));++j)
                for(Eigen::Index i=0;i<(image.dims(0));++i)
                {
                    out_image.volume(ijk_to_idx(i,j,k)) = solution(ijk_to_idx(i,j,k));
                }
        

        if(par.count("out"))
            save_volume(par["out"].as<std::string>().c_str(),
                        par["in"].as<std::string>().c_str(),
                        out_image,false);

        return 0;
    } else {
        std::cout<<options.help()<<std::endl;
        return 1;
    }
}