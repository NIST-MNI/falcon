#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>


#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <igl/principal_curvature.h>


#include "depth_potential.h"

#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>


#include  <chrono>
#include  <iostream>


bool depth_potential(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,double alpha,Eigen::VectorXd &dp_out)
{
    // implemeting depth potential calculation from 
    // http://dx.doi.org/10.1016/j.media.2008.09.001
    // Maxime Boucher  Sue Whitesides Alan Evans 
    // "Depth potential function for folding pattern representation, registration and analysis"
    //

    #if 1
    // Compute curvature directions via quadratic fitting
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;

    igl::principal_curvature(V,F, PD1, PD2, PV1, PV2);

    // mean curvature
    Eigen::RowVectorXd curvature = (0.5*(PV1+PV2));
    #endif

    Eigen::SparseMatrix<double> L, femB, invFemB, LB;
    Eigen::VectorXd femBS;

    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI,femB);

    #if 0
    // compute descrete curvature
    igl::invert_diag(femB,invFemB);
    Eigen::MatrixXd HN = -invFemB*(L*V);
    Eigen::RowVectorXd curvature = (HN.rowwise().norm()); //up to sign
    #endif

    igl::sum(femB,1,femBS);
    
    LB =  Eigen::SparseMatrix<double>(femBS.asDiagonal());
    LB.makeCompressed();

    L *= -0.5; // convert to notation from 
    Eigen::VectorXd  B = (curvature - curvature*LB/LB.sum())*LB*-2.0;
    Eigen::SparseMatrix<double>  M = LB*alpha + L; 

    //Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double > > solver; 
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double > > solver;
    //Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double >  > solver;
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper,Eigen::IncompleteLUT<double > > solver;
    solver.preconditioner().setDroptol(0.001);

    solver.compute(M);

    if(solver.info()!= Eigen::Success) {
        std::cerr<<"Solving failed"<<std::endl;
        return false;
    } else {

        dp_out=solver.solve(B);
        return true;
    }
}



bool depth_potential_benchmark(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,double alpha,Eigen::VectorXd &dp_out)
{
    #if 1
    // Compute curvature directions via quadric fitting
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;

    igl::principal_curvature(V,F, PD1, PD2, PV1, PV2);

    // mean curvature
    Eigen::RowVectorXd curvature = (0.5*(PV1+PV2));
    #endif

    Eigen::SparseMatrix<double> L,femB,invFemB,LB;
    Eigen::VectorXd femBS;

    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI,femB);

    #if 0
    // compute descrete curvature
    igl::invert_diag(femB,invFemB);
    Eigen::MatrixXd HN = -invFemB*(L*V);
    Eigen::RowVectorXd curvature = (HN.rowwise().norm()); //up to sign
    #endif
    igl::sum(femB,1,femBS);

    LB =  Eigen::SparseMatrix<double>(femBS.asDiagonal());
    LB.makeCompressed();

    L *= -0.5;
    Eigen::VectorXd  B = (curvature - curvature*LB/LB.sum())*LB*-2.0;
    Eigen::SparseMatrix<double>  M = LB*alpha + L;
 
    auto elapsed = [](const std::chrono::time_point<std::chrono::steady_clock> &a,
                      const std::chrono::time_point<std::chrono::steady_clock> &b)
              {return std::chrono::duration<double, std::milli>(a-b).count();};

    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);
        std::cout << "CG IdentityPreconditioner: "<<std::endl;
        auto start = std::chrono::steady_clock::now();

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> solver;

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() <<
                     ", time: "<< elapsed(std::chrono::steady_clock::now(),start) <<" "<<x.sum()<<std::endl;

    }
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);
        std::cout << "BiCGSTAB IdentityPreconditioner: "<<std::endl;
        auto start = std::chrono::steady_clock::now();

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() <<
                     ", time: "<< elapsed(std::chrono::steady_clock::now(),start) <<" "<<x.sum()<<std::endl;

    }
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);

        std::cout << "GMRES IdentityPreconditioner: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() <<
                     ", time: "<< elapsed(std::chrono::steady_clock::now(),start) <<" "<<x.sum()<<std::endl<<std::flush;

    }
#if 0
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);

        std::cout << "DGMRES IdentityPreconditioner: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;

        solver.compute(M);
        auto x = solver.solve(B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() <<
                     ", time: "<<elapsed(std::chrono::steady_clock::now(),start)<<" "<<x.sum()<<std::endl<<std::flush;
    }
#endif
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);
        
        std::cout << "CG ILUT: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteLUT<double > > solver;
        solver.preconditioner().setDroptol(0.001);

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << 
                     ", time: "<<elapsed(std::chrono::steady_clock::now(),start)<<" "<<x.sum()<<std::endl<<std::flush;
    }
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);

        std::cout << "BiCGSTAB ILUT: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double >> solver;
        solver.preconditioner().setDroptol(0.001);

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << 
                     ", time: "<<elapsed(std::chrono::steady_clock::now(),start)<<" "<<x.sum()<<std::endl;
    }
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);

        std::cout << "GMRES ILUT: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double > > solver;
        solver.preconditioner().setDroptol(0.001);

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << 
                     ", time: "<<elapsed(std::chrono::steady_clock::now(),start)<<" "<<x.sum()<<std::endl;
    }
    {
        Eigen::SparseMatrix<double> _M(M);
        Eigen::VectorXd  _B(B);

        std::cout << "DGMRES ILUT: "<<std::endl;
        auto start = std::chrono::steady_clock::now();
        Eigen::DGMRES<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double > > solver;
        solver.preconditioner().setDroptol(0.001);

        solver.compute(_M);
        auto x = solver.solve(_B);
        std::cout << "  #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << 
                     ", time: "<<elapsed(std::chrono::steady_clock::now(),start)<<" "<<x.sum()<<std::endl;

        dp_out=x;
    }


    return true;
}
