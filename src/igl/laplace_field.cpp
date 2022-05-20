#include <string>
#include <chrono>

#include "cxxopts.hpp"

//#include "fem_system.h"

// unsupported
#include <unsupported/Eigen/CXX11/Tensor>


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
        ("h,help", "Print help")
    ;
    
    options.parse_positional({"in","out"});
    auto par = options.parse(argc, argv);

    if(par.count("in") && par.count("out") && !par.count("help"))
    {

#ifdef USE_OPENMP
        std::cout<<"Using OpenMP, max number of threads="<<omp_get_max_threads()<<std::endl;
#endif
        
        // ITK setup
        // using PixelType = double;
        // using ImageType = itk::Image<PixelType, 3>;
        // using LabelImageType = itk::Image<unsigned char, 3>;
        // auto tick = std::chrono::steady_clock::now();

        // using ReaderType = itk::ImageFileReader<LabelImageType>;
        // ReaderType::Pointer reader = ReaderType::New();
        // reader->SetFileName(par["in"].as<std::string>());
        // reader->Update();
        // LabelImageType::Pointer image = reader->GetOutput();
        // LabelImageType::RegionType region = image->GetLargestPossibleRegion();

        // using LabelTensor = Eigen::Tensor<unsigned char,3,Eigen::ColMajor>;
        // using Tensor = Eigen::Tensor<double,3,Eigen::ColMajor>;

        // LabelTensor in_vol=Eigen::TensorMap< LabelTensor>(image->GetBufferPointer(),
        //     region.GetSize(2), region.GetSize(1), region.GetSize(0));
    
        // ImageType::Pointer out_image;
        // define_image_byref<LabelImageType,ImageType>(image, out_image);

        // out_image->Allocate();
        // out_image->FillBuffer(0.0);

        // Eigen::TensorMap<Tensor> out_vol(out_image->GetBufferPointer(),
        //     region.GetSize(2), region.GetSize(1), region.GetSize(0));
        
        
        using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        int total_voxels=region.GetSize(2)*region.GetSize(1)*region.GetSize(0);
        SparseMatrix A;
        if(par["grad"].as<bool>())
        {
            // using gradient operators
            // FOR DEBUGGING ONLY!
            SparseMatrix Dxf,Dyf,Dzf;
            SparseMatrix Dxb,Dyb,Dzb;
            create_grad_matrixes(region.GetSize(2), region.GetSize(1), region.GetSize(0), Dxf,Dxb,Dyf,Dyb,Dzf,Dzb);
            A=Dxf*Dxb+Dyf*Dyb+Dzf*Dzb;
        } else {
            create_laplacian_matrix(region.GetSize(2), region.GetSize(1), region.GetSize(0), A, par["s27"].as<bool>());
        }
        tick = std::chrono::steady_clock::now();
        ////

        std::cout << "Created sparse matrix:" << std::chrono::duration <double, std::milli> (std::chrono::steady_clock::now() - tick).count() << " ms" << std::endl;
        std::cout << "Total number of voxels: "<< total_voxels<<std::endl;
        std::cout << "Number of nonzero elements in A:" << A.nonZeros()<<std::endl;

        auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
        {
            return k+j*region.GetSize(0)+i*region.GetSize(0)*region.GetSize(1);
        };

        tick = std::chrono::steady_clock::now();
        // setting the dirichlet boundary condition
        using triplet = Eigen::Triplet<double>;
        std::vector<triplet> B_;

        for(Eigen::Index k=0;k<region.GetSize(0);++k) 
            for(Eigen::Index j=0;j<region.GetSize(1);++j)
                for(Eigen::Index i=0;i<region.GetSize(2);++i)
                {
                    if(in_vol(i,j,k)>0)
                    {
                        double val=0.0;
                        Eigen::Index idx = ijk_to_idx(i,j,k);
                        if(in_vol(i,j,k)>1) val = 1.0;
 
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
        solver.preconditioner().setDroptol(1e-4); // preconditioner works slower, but solution is quite fast after

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

        for(Eigen::Index k=0;k<(region.GetSize(2));++k) 
            for(Eigen::Index j=0;j<(region.GetSize(1));++j)
                for(Eigen::Index i=0;i<(region.GetSize(0));++i)
                {
                    out_vol(k,j,i) = solution(ijk_to_idx(i,j,k));
                }
        
        // using WriterType = itk::ImageFileWriter<ImageType>;
        // WriterType::Pointer writer = WriterType::New();
        // writer->SetFileName(par["out"].as<std::string>());
        // writer->SetInput(out_image);
        // writer->Update();

        return 0;
    } else {
        std::cout<<options.help()<<std::endl;
        return 1;
    }
}