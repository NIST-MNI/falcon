#include "fem_system.h"


#include <igl/readMSH.h>
#include <igl/writeMSH.h>
#include "tri_to_tet.h"

#define USE_SVD 0
#include "mesh_util.h"


// HACK
#ifdef USE_HYPRE
#include <stdio.h>
#include "simnibs/_solver.h"
#endif


// print out gradient only when MESH_DEBUG is defined
const int __debug_tet=2684719; //2341252; // 2152426; //2341252; // 2152426; //;




void compute_gradient(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G)
{
// from Simnibs:
// ''' G calculates the gradient of a function in each tetrahedra
// The way it works: The operator has 2 parts
// G = T^{-1}A
// A is a projection matrix
// A = [-1, 1, 0, 0]
//     [-1, 0, 1, 0]
//     [-1, 0, 0, 1]
// And T is the transfomation to baricentric coordinates
//
// here the output G is 4*Tet,3 with colum correspoding to x,y,z

    // WARNING: output is th * 4 x XYZ ( simnibs output is th x XYZ x 4 ? )
    G.resize(Tet.rows()*4, 3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,4> A;
        A << -1, 1, 0, 0,
             -1, 0, 1, 0,
             -1, 0, 0, 1;

        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,4> _G;

        T << ( X.row(Tet(i,1)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,2)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,3)) - X.row(Tet(i,0)) );

        Eigen::PartialPivLU<Eigen::Matrix<double,3,3> > lu(T);
        _G = lu.solve(A);
        //_G = T.colPivHouseholderQr().solve(A);

        // TODO: check if this is what we want
        G.template block<4,3>(i*4, 0) = _G.transpose();
        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == __debug_tet ) {
            std::cout<<"Debug:Tetrahedral gradient:" << debug_tet << std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A="      << A << std::endl;
            std::cout<<"T="      << T << std::endl;
            std::cout<<"G="      <<_G << std::endl;
        }
        #endif //MESH_DEBUG
    }       
}

void compute_gradient_svd(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                            Eigen::MatrixXd &G,
                            double epsilon)
{
// from Simnibs:
// ''' G calculates the gradient of a function in each tetrahedra
// The way it works: The operator has 2 parts
// G = T^{-1}A
// A is a projection matrix
// A = [-1, 1, 0, 0]
//     [-1, 0, 1, 0]
//     [-1, 0, 0, 1]
// And T is the transfomation to baricentric coordinates
//
// here the output G is 4*Tet,3 with colum correspoding to x,y,z

    // WARNING: output is th * 4 x XYZ ( simnibs output is th x XYZ x 4 ? )
    G.resize(Tet.rows()*4, 3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,4> A;
        A << -1, 1, 0, 0,
             -1, 0, 1, 0,
             -1, 0, 0, 1;

        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,4> _G;

        T << ( X.row(Tet(i,1)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,2)) - X.row(Tet(i,0)) ),
             ( X.row(Tet(i,3)) - X.row(Tet(i,0)) );

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(epsilon);
        _G = svd.solve(A);

        // TODO: check if this is what we want
        G.template block<4,3>(i*4, 0) = _G.transpose();
        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == __debug_tet ) {
            std::cout<<"Debug:Tetrahedral gradient:" << debug_tet << std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A="      << A << std::endl;
            std::cout<<"T="      << T << std::endl;
            std::cout<<"G="      <<_G << std::endl;
        }
        #endif //MESH_DEBUG
    }       
}



void compute_field_gradient(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                      const Eigen::VectorXd &Field,
                      Eigen::MatrixXd &G)
{
    // output is th * 4 * XYZ 
    G.resize(Tet.rows(),3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,1>  A;
        Eigen::Matrix<double,3,1> _G;

        T << (X.row(Tet(i,1)) - X.row(Tet(i,0))),
             (X.row(Tet(i,2)) - X.row(Tet(i,0))),
             (X.row(Tet(i,3)) - X.row(Tet(i,0)));

        A << Field(Tet(i,1)) - Field(Tet(i,0)),
             Field(Tet(i,2)) - Field(Tet(i,0)),
             Field(Tet(i,3)) - Field(Tet(i,0));

        //_G = T.colPivHouseholderQr().solve(A);
        Eigen::PartialPivLU<Eigen::Matrix<double,3,3> > lu(T);
        _G = lu.solve(A);

        // TODO: check if this is what we want
        G.row(i) = _G;

        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == debug_tet ) {
            std::cout<<"Debug:field gradient:"<<debug_tet<<std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A=" <<A<<std::endl;
            std::cout<<"T=" <<T<<std::endl;
            std::cout<<"G=" <<_G<<std::endl;
        }
        #endif        
    }       
}


void compute_field_gradient_svd(const Eigen::MatrixXd &X,
                      const Eigen::MatrixXi &Tet,
                      const Eigen::VectorXd &Field,
                      Eigen::MatrixXd &G,
                      double epsilon)
{
    // output is th * 4 * XYZ 
    G.resize(Tet.rows(),3);

    // TODO: paralellize?
    for(size_t i=0; i<Tet.rows(); ++i)
    {
        Eigen::Matrix<double,3,3>  T;
        Eigen::Matrix<double,3,1>  A;
        Eigen::Matrix<double,3,1> _G;

        T << (X.row(Tet(i,1)) - X.row(Tet(i,0))),
             (X.row(Tet(i,2)) - X.row(Tet(i,0))),
             (X.row(Tet(i,3)) - X.row(Tet(i,0)));

        A << Field(Tet(i,1)) - Field(Tet(i,0)),
             Field(Tet(i,2)) - Field(Tet(i,0)),
             Field(Tet(i,3)) - Field(Tet(i,0));

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
        svd.setThreshold(epsilon);
        _G = svd.solve(A);

        // TODO: check if this is what we want
        G.row(i) = _G;

        //DEBUG:
        #ifdef MESH_DEBUG
        if( i == debug_tet ) {
            std::cout<<"Debug:field gradient:"<<debug_tet<<std::endl;
            std::cout<<"det(T)=" << T.determinant() << std::endl;
            std::cout<<"A=" <<A<<std::endl;
            std::cout<<"T=" <<T<<std::endl;
            std::cout<<"G=" <<_G<<std::endl;
        }
        #endif        
    }       
}


FemSystem::FemSystem()
{
    cond_table={
        {0, 0.126}, //WM
        {1, 0.275}, //GM
        {2, 1.654}, //CSF
        {3, 0.010}, //Bone
        {4, 0.465}, //Scalp
        {5, 0.5  }, //Eye_balls
        {6, 0.008}, //Compact_bone
        {7, 0.025}, //Spongy_bone
        {8, 0.6  }, //Blood
        {9, 0.16 }, //Muscle
        {99, 29.4}, //Electrode_rubber
        {499, 1.0}  //Saline
    };
};

void FemSystem::load_mesh(const std::string &mesh_file)
{
    igl::readMSH( mesh_file, X, Tri, Tet, TriTag, TetTag, XFields, XF, EFields, TriF, TetF);
}

void FemSystem::prepare_mesh(void)
{
    build_tri_to_tet_index(Tet, Tri, tri_to_tet);
    separate_by_tag(Tet, TetTag, tag_to_tet, Tet_by_Tag, Edge_by_Tag);

    TetBC.resize(Tet.rows(),3);
    mesh_barycenter(X, Tet, TetBC);
    vols.resize(Tet.rows());
    tet_volumes(X, Tet, vols);

    // prepare gradients
    if(grad_epsilon>0.0)
        compute_gradient_svd(X, Tet, G,grad_epsilon);
    else
        compute_gradient(X, Tet, G);

    // initalize conductivity
    // this is scalar case
    // TODO: implement directional conductivity
    conductivity = Eigen::MatrixXd::NullaryExpr( Tet.rows(), 3,
            [this](auto i, auto j)->double {
                auto it=this->cond_table.find(this->TetTag(i)-1);
                if(it!=this->cond_table.end())
                    return it->second;
                else
                    return 0.0; // TODO: abort?
                } );

    A.resize(X.rows(), X.rows());
    vGc.resize(G.rows(), G.cols());

    // this works only for scalar conductivities
    // TODO: implement directional conductivity
    // vGc = vols[:, None, None]*G*cond[:, None, None]
    for(size_t i=0; i<G.rows(); ++i)
        vGc.row(i) = vols(i/4) * G.row(i) * conductivity(i/4,0); // TODO: implement non-uniform conductivity

    // initializing main stiffness matrix
    std::vector<triplet> A_;

    A_.reserve(6*2*Tet.rows());

    for(size_t i=0; i<4; ++i)
        for(size_t j=i; j<4; ++j)
            for(size_t k=0; k<Tet.rows(); ++k)
            {
                //Kg = (vGc[:, i, :]*G[:, j, :]).sum(axis=1)
                int r = Tet(k, i);
                int c = Tet(k, j);

                double Kg = vGc.row( k*4+i ).dot( G.row( k*4+j ) );
                // A
                A_.push_back(triplet(r, c, Kg));

                // A.T 
                if(i!=j)
                    A_.push_back(triplet(c, r, Kg));
            }
    A.setFromTriplets(A_.begin(), A_.end());
    A /= scaling_factor;  // # * 1e6 from the gradiend operator, 1e-9 from the volume
}


void FemSystem::setup_volumes(BaseImageType* reference)
{
    if(reference!=nullptr)
    {
        define_vec_image_byref<BaseImageType, ImageType>(reference, ref_out_image);
    } else {
        define_vec_image_bbox<ImageType>(X, ref_out_image);
    }

    ref_out_image->Allocate();
    C.resize(X.rows(), X.cols());
    xyz_to_index<ImageType>(ref_out_image, X, C);
}


TmsFemSolution * FemSystem::solve_tms(
    Eigen::MatrixXd dAdt)
{
    // allocate solution
    TmsFemSolution *solution = new TmsFemSolution(this);
    
    solution->dAdt=dAdt;
    //
    solution->sigma_dadt = solution->dAdt.cwiseProduct(conductivity).eval();

    std::vector<triplet> b_;
    b_.reserve( Tet.rows()*4 );
    for(size_t j=0; j<Tet.rows(); ++j)
    {
        for(size_t i=0; i<Tet.cols(); ++i)
        {
            double elm_node_integral = solution->sigma_dadt.row(j).dot(G.row(j*4 + i)) * vols(j)*-1.0;
            elm_node_integral /= (scaling_factor*scaling_factor);

            b_.push_back( triplet(Tet(j, i), 0, elm_node_integral) );
            //B(Tet(j, i)) += elm_node_integral;
        }
    }

    // this unrolls it into weighet sum of elm_node_integral
    // using a trick that sparse matrix with repeated entries sums them up
    Eigen::SparseMatrix<double> b(X.rows(), 1);
    b.setFromTriplets(b_.begin(), b_.end());
    b.makeCompressed();

    Eigen::VectorXd B=Eigen::MatrixXd(b);

    // impose Dirichlet boundary condition
    Eigen::Index lowest_element=0;
    // find the lowest Z element
    X.col(2).minCoeff(&lowest_element);

    // applying dirichlet boundary condition for TMS
    // by "grounding" the lowest element
    // and forcing it's potential to 0
    B(lowest_element) = 0.0;
    A.row(lowest_element) *= 0.0;
    A.coeffRef(lowest_element, lowest_element) = 1.0;
    A.makeCompressed();

    // solution will be put here
    solution->VFieldX.resize(B.rows());

    #ifndef USE_HYPRE //HACK
    
    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double > > solver; // seem to be the fastest single threaded
    //Eigen::BiCGSTAB<SparseMatrix > solver; // slow without preconditioner

    //Eigen::DGMRES<SparseMatrix, Eigen::IncompleteLUT<double > > solver;  // need "unsupported" part of Eigen
    //Eigen::ConjugateGradient<SparseMatrix,Eigen::Lower|Eigen::Upper, Eigen::IncompleteLUT<double> > solver; // really slow with preconditioner

    // ConjugateGradient for some reason doens't quote work with IncompleteLUT
    // OpenMP accelerates CG stronger then BiCGSTAB
    // Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower|Eigen::Upper> solver;

    // these parameters will increase solution time for no reason
    // solver.setMaxIterations( 1000 );
    // solver.setTolerance(1e-20); 

    solver.preconditioner().setDroptol(1e-4); // preconditioner works slower, but solution is quite fast after
    solver.compute(A);

    if(solver.info()!= Eigen::Success) {
        //std::cerr<<"Solving failed"<<std::endl;
        //return -1;
        //TODO: throw an exception
        abort();
    } else {
        solution->VFieldX = solver.solve(B);
    }
    #else // HACK: using Hypre library, with the same parameters as SimNIBS
    // solving using interface from SimNIBS
    const char *petsc_argv[] = { "-ksp_type", "cg",
        "-ksp_rtol","1e-10","-pc_type","hypre",
        "-pc_hypre_type","boomeramg",
        "-pc_hypre_boomeramg_coarsen_type","HMIS", NULL };
    int petsc_argc=sizeof(petsc_argv)/sizeof(petsc_argv[0])-1;
    _petsc_initialize();
    KSP ksp;
    PetscErrorCode err;

    // CONVERT Eigen RowMajor data to something that petsc understands
    int N=A.nonZeros();
    double *data=A.valuePtr();
    int *indptr=A.outerIndexPtr();
    int *indices=A.innerIndexPtr();

    err=_petsc_prepare_ksp(petsc_argc,(char **)petsc_argv, A.rows(), indptr, indices, data , stdout, &ksp);
    _print_ksp_info(ksp, stdout);
    // solve

    Eigen::VectorXd VFieldX(B.rows());
    err = _petsc_solve_with_ksp(ksp, A.rows(), B.data(), stdout, VFieldX.data());

    solution->VFieldX=VFieldX;
    _dealloc(ksp);
    _petsc_finalize();
    #endif

    //std::cout << "estimated error: " << solver.error()      << std::endl;

    // DEBUG: calculate error
    solution->rel_error = solver.error(); // (A*solution->VFieldX - B).norm()/B.norm();

    return solution;
}


TmsFemSolution * FemSystem::solve_tms(
        ImageType::Pointer coil_field,
        double dIDt,
        const Eigen::Matrix4d& coil_pos)
{

    //typedef itk::NearestNeighborInterpolateImageFunction< image3d, double >  InterpolatorType;
    using InterpolatorType = itk::LinearInterpolateImageFunction< ImageType, double >;
    //using InterpolatorType = itk::VariableVectorBSplineInterpolate< ImageType, double >;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    //interpolator->SetSplineOrder(1);
    interpolator->SetInputImage(coil_field);

    Eigen::MatrixXd dAdt_X(X.rows(),3);

    sample_vector_field_v2<InterpolatorType, ImageType>(interpolator,
            coil_pos, coil_field, X, dAdt_X, Eigen::Vector3d::Zero(), false);

    // project to tetrahedra by averaging corners
    Eigen::MatrixXd dAdt = Eigen::MatrixXd::NullaryExpr(Tet.rows(), 3,
            [&dAdt_X, this](auto i,auto j)->double {
                return ( dAdt_X(this->Tet(i,0),j) + dAdt_X(this->Tet(i,1),j) + 
                         dAdt_X(this->Tet(i,2),j) + dAdt_X(this->Tet(i,3),j) )/4;} );
    dAdt *= dIDt;

    return this->solve_tms(dAdt);
}


TmsFemSolution::TmsFemSolution(FemSystem* sys):fem(sys),rel_error(-1.0)
{
}


void TmsFemSolution::save_result_mesh(const std::string &fname, bool save_e_field, bool save_v_field )
{
    std::vector<std::string> XFields;
    std::vector<std::string> EFields;
    std::vector<Eigen::MatrixXd> XF,TriF,TetF;

    if(save_e_field)
    {
        assert( this->EFieldTet.rows()==fem->Tet.rows() );
        assert( this->EFieldTet.cols()==3 );

        EFields.push_back("E");
        TetF.push_back(this->EFieldTet);
        Eigen::MatrixXd EFieldTri(fem->Tri.rows(), 3);
        propagate_on_tri(fem->tri_to_tet, this->EFieldTet, EFieldTri);
        TriF.push_back(EFieldTri);
    }

    if(save_v_field)
    {
        assert( this->VFieldX.rows()==fem->X.rows() );
        assert( this->VFieldX.cols()==1 );

        XFields.push_back("v");
        XF.push_back( this->VFieldX );
    }

    igl::writeMSH(fname,
        fem->X, fem->Tri, fem->Tet, fem->TriTag, fem->TetTag, 
        XFields, XF, EFields, TriF, TetF);        
}


void TmsFemSolution::save_result_image(ImageType::Pointer img,
    const std::string &fname)
{
    using WriterType = itk::ImageFileWriter<ImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(fname);
    writer->SetInput(img);
    writer->Update();
}


void TmsFemSolution::calculate_E_field(void)
{
    assert( this->fem!=nullptr );
    assert( this->fem->X.rows()==this->VFieldX.rows() );
    assert( this->VFieldX.cols()==1 );

    Eigen::MatrixXd Vgrad(this->fem->Tet.rows(), 3 );

    if(this->fem->grad_epsilon>0.0)
        compute_field_gradient_svd(this->fem->X, this->fem->Tet, VFieldX, Vgrad, this->fem->grad_epsilon);
    else
        compute_field_gradient(this->fem->X, this->fem->Tet, VFieldX, Vgrad);

    Vgrad*= this->fem->scaling_factor;

    this->EFieldTet = ( this->dAdt + Vgrad )*-1.0;
}


void TmsFemSolution::map_E_field_to_volume(void)
{
    // initialize image
    if(fem->ref_out_image.IsNotNull()) {
        define_vec_image_byref<ImageType, ImageType>(fem->ref_out_image, this->E_volume);
        define_vec_image_byref<ImageType, ImageType>(fem->ref_out_image, this->norm_volume);

        this->E_volume->Allocate();
        this->norm_volume->Allocate();
 
        ImageType::PixelType default_pixel;
        default_pixel.SetSize(3);
        default_pixel.Fill(0.0);

        this->E_volume->FillBuffer(default_pixel);
        this->norm_volume->FillBuffer(default_pixel);
    }
    // need to interpolate on volume
    // iterate over all tags
    for(auto e_it:this->fem->Edge_by_Tag)
    {
        int _tag = e_it.first;
        int _tag_n = this->fem->tag_to_tet.count(_tag);
        auto t_it = this->fem->Tet_by_Tag.find(_tag);

        std::unordered_set<int> &_edge = e_it.second;
        Eigen::MatrixXi &_tet = t_it->second;

        Eigen::MatrixXd _TetBC(    _tag_n, this->fem->TetBC.cols());
        Eigen::MatrixXd _vols(     _tag_n, this->fem->vols.cols());
        Eigen::MatrixXd _EFieldTet(_tag_n, this->EFieldTet.cols());

        extract_tag(_tag, this->fem->tag_to_tet, this->fem->TetBC,     _TetBC);
        extract_tag(_tag, this->fem->tag_to_tet, this->fem->vols,      _vols);
        extract_tag(_tag, this->fem->tag_to_tet, this->EFieldTet, _EFieldTet);

        Eigen::MatrixXd fldX(this->fem->X.rows(), this->EFieldTet.cols());
        fldX.setZero();

        tet_to_node_superconvergent(this->fem->X, _tet, _TetBC, _vols, _EFieldTet, _edge, fldX);
        //tet_to_node_average(X, _tet, _conductivity, condX);

        std::vector<Eigen::Matrix3d> idx_to_bar;
        create_baricentric_matrix(this->fem->C, _tet, idx_to_bar); 
        mesh_to_image<ImageType>(this->fem->C, _tet, idx_to_bar, fldX, this->E_volume, this->norm_volume);
    }

    // normalize image, TODO: maybe get rid of this
    normalize_vec_image<ImageType>(this->E_volume, this->norm_volume);
}


void print_image_info(itk::ImageBase<3>* ref)
{
    using Image=itk::ImageBase<3>;
    Image::PointType org=ref->GetOrigin( );
    Image::SpacingType spc=ref->GetSpacing( );
    Image::RegionType region=ref->GetLargestPossibleRegion ();
    Image::DirectionType dir=ref->GetDirection();
    Image::RegionType::SizeType size=region.GetSize();
    // check if the image was flipped
    std::vector<double> orig_spacing;
    itk::ImageBase<3>::DirectionType orig_dir;

    std::cout<<"sz:"<<size<<" sp:"<<spc<<" o:"<< org <<std::endl;
}


void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxb, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyb,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzf,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzb)
{
    int total_voxels=si*sj*sk;
    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> Dx_,Dy_,Dz_;
    std::vector<triplet> Dx_2,Dy_2,Dz_2;
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };
    Dx_.reserve(total_voxels*2);
    Dy_.reserve(total_voxels*2);
    Dz_.reserve(total_voxels*2);

    Dx_2.reserve(total_voxels*2);
    Dy_2.reserve(total_voxels*2);
    Dz_2.reserve(total_voxels*2);

    for(Eigen::Index k=0;k<sk;++k) 
        for(Eigen::Index j=0;j<sj;++j)
            for(Eigen::Index i=0;i<si;++i)
            {
                if(k>0)      Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), -1));
                Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));

                if(j>0)      Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), -1));
                Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));
                
                if(i>0)      Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), -1));
                Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), 1));


                Dx_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(k<(sk-1)) Dx_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 1));

                Dy_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(j<(sj-1)) Dy_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 1));
                
                Dz_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k), -1));
                if(i<(si-1)) Dz_2.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 1));
            }
    Dxf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzf=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzb=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxb.setFromTriplets(Dx_.begin(), Dx_.end());
    Dyb.setFromTriplets(Dy_.begin(), Dy_.end());
    Dzb.setFromTriplets(Dz_.begin(), Dz_.end());

    Dxf.setFromTriplets(Dx_2.begin(), Dx_2.end());
    Dyf.setFromTriplets(Dy_2.begin(), Dy_2.end());
    Dzf.setFromTriplets(Dz_2.begin(), Dz_2.end());

    Dxf.makeCompressed();
    Dyf.makeCompressed();
    Dzf.makeCompressed();
    Dxb.makeCompressed();
    Dyb.makeCompressed();
    Dzb.makeCompressed();
}

void create_grad_matrixes(int si,int sj,int sk, 
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dxc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dyc,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &Dzc)
{
    // centered version
    int total_voxels=si*sj*sk;
    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> Dx_,Dy_,Dz_;
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };
    Dx_.reserve(total_voxels*2);
    Dy_.reserve(total_voxels*2);
    Dz_.reserve(total_voxels*2);

    for(Eigen::Index k=0;k<sk;++k) 
        for(Eigen::Index j=0;j<sj;++j)
            for(Eigen::Index i=0;i<si;++i)
            {
                if(k>0)      Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), -0.5));
                if(k<(sk-1)) Dx_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 0.5));

                if(j>0)      Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), -0.5));
                if(j<(sj-1)) Dy_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 0.5));
                
                if(i>0)      Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), -0.5));
                if(i<(si-1)) Dz_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 0.5));

            }
    Dxc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dyc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    Dzc=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);

    Dxc.setFromTriplets(Dx_.begin(), Dx_.end());
    Dyc.setFromTriplets(Dy_.begin(), Dy_.end());
    Dzc.setFromTriplets(Dz_.begin(), Dz_.end());

    Dxc.makeCompressed();
    Dyc.makeCompressed();
    Dzc.makeCompressed();
}


void create_laplacian_matrix(int si,int sj,int sk, Eigen::SparseMatrix<double, Eigen::RowMajor> &A, bool s27)
{
    // now we can determine 2nd order differential equation
    int total_voxels=si*sj*sk;

    using triplet = Eigen::Triplet<double>;
    std::vector<triplet> A_;
    
    auto ijk_to_idx = [&](auto i,auto j,auto k) -> Eigen::Index 
    {
        return k+j*sk+i*sk*sj;
    };

    if(s27)
    {
        auto ijk_hit = [&](int i,int j,int k) -> bool
        {
            return i>=0  && j>=0 && k>=0 && 
                    i<si && j<sj && k<sk;
        };

        auto ijk_27stencil = [&](int i,int j,int k) -> double
        {
            switch(abs(i)+abs(j)+abs(k))
            {
                case 0:return -88.0/26;
                case 1:return   6.0/26;
                case 2:return   3.0/26;
                case 3:return   2.0/26;
                default: return 0.0;
            }

        };

        A_.reserve(total_voxels*27);
        // using 27 point stencil 
        // see https://en.wikipedia.org/wiki/Discrete_Laplace_operator
        for(Eigen::Index k=0;k<sk;++k) 
            for(Eigen::Index j=0;j<sj;++j)
                for(Eigen::Index i=0;i<si;++i)
                {
                    // fill out stencil:
                    for(int ii=-1;ii<2;++ii)
                        for(int jj=-1;jj<2;++jj)
                            for(int kk=-1;kk<2;++kk)
                                if(ijk_hit(i+ii,j+jj,k+kk))
                                    A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+ii,j+jj,k+kk),  
                                        ijk_27stencil(ii,jj,kk) ));
                }
    } else {
        A_.reserve(total_voxels*7);
        // using 7 point stencil for now 
        // see https://en.wikipedia.org/wiki/Discrete_Laplace_operator
        for(Eigen::Index k=0;k<sk;++k) 
            for(Eigen::Index j=0;j<sj;++j)
                for(Eigen::Index i=0;i<si;++i)
                {
                    // fill out stencil:
                    A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k),  -6.0));

                    if(i>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i-1, j, k), 1.0));
                    if(j>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j-1, k), 1.0));
                    if(k>0) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k-1), 1.0));

                    if(i<(si-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i+1, j, k), 1.0));
                    if(j<(sj-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j+1, k), 1.0));
                    if(k<(sk-1)) A_.push_back(triplet(ijk_to_idx(i,j,k), ijk_to_idx(i, j, k+1), 1.0));
                }
    }
    A=Eigen::SparseMatrix<double, Eigen::RowMajor>(total_voxels, total_voxels);
    A.setFromTriplets(A_.begin(), A_.end());
    A.makeCompressed();
}

