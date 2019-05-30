/* Filename:     nifti1_kunio_aregister.c
 * Description:  affine registration function
 * Author:       Kunio Nakamura, knakamura@mrs.mni.mcgill.ca
 * Date:         February 24, 2012
 *
 *
 * MAJOR VERSION: 1
 * MINOR VERSION: 1
 * MICRO VERSION: 0
 *
 * VERSION INFO 1.0.0, Kunio Nakamura, knakamura@mrs.mni.mcgill.ca
 * -initial version for keeping track
 *
 * VERSION INFO 1.1.0, Kunio Nakamura, knakamura@mrs.mni.mcgill.ca
 * -initial version for keeping track
 *
 */

/* definition of affine parameter vector
 *
 *    0 = registration error
 *    1-3   rotations in about x,y,z
 *    4-6   translations in x,y,z
 *    7     global scaling
 *    8-10  individual scaling
 *    11-13 shear
 *    14-16 pivot (center)
 *
 *
 * -qform and sform are not used during registration
 * -but affine matrix can be saved in qform or sform
 *
 */

#ifndef _FALCON_AREGISTER_C_
#define _FALCON_AREGISTER_C_

#include "falcon.h"


/* global variables for register */
nifti_image *g_niik_aregister_movimg=NULL;
nifti_image *g_imgs[12]; /* global variable for registration:
			  *  [0] = refimg, [1] = movimg,
			  *  [2] = refseg, [3] = movseg  */

double *g_ref_data=NULL;     /* reference image's sampled data */
int g_reg_method=-1;         /* registration methods (NIIK_REGISTER_*) */
double g_sample = -1;        /* sampling distance */
int *g_fixpar=NULL;          /* fixed parameters -corresponds to g_affpar */
int g_invmat=0;              /* inverse matrix -useful for inverse registration of symmetric registration */
double *g_affpar=NULL;       /* affine parameters
			      * [0] = flag
			      * [1-3] = 3 rotations
			      * [4-6] = 3 translations
			      * [7-9] = 3 scaling
			      * [10-12] = 3 shearing
			      * [13-15] = 3 centr of rotation (pivot) */
nmi_obj *g_nmiobj = NULL;    /* variable for normalized mutual info */
int g_iter = 0;              /* iteration */
int g_verbose_errfunc=1;     /* error functions' verbose */
niikmat *g_niik_aregister_imat=NULL; /* initial (pre) matrix */
double apodization_thresh=2; /* apodization threshold */


/* FUNCTIONS FOR GLOBAL VARIABLES */
int niik_image_aregister_set_g_img(nifti_image *img,int idx) {
  g_imgs[idx]=img;
  return 1;
}
int niik_image_aregister_set_g_reg_method(int method) {
  g_reg_method=method;
  return 1;
}
int niik_image_aregister_set_g_img_to_null() {
  int n;
  for(n=0; n<12; n++) g_imgs[n]=NULL;
  return 1;
}
int niik_image_aregister_set_g_invmat(int imat) {
  g_invmat=imat;
  return 1;
}
int niik_image_aregister_get_g_invmat() {
  return g_invmat;
}
int niik_image_aregister_set_g_verbose_errfunc(int v) {
  g_verbose_errfunc=v;
  return 1;
}

/* INTERNAL FUNCTIONS */
double niik_image_aregister_obj_func(double *v); /* objective cost function */
double *niik_image_aregister_obj_func_set_ref_data(nifti_image *refimg, nifti_image *maskimg, double sample); /* set reference data */
double niik_image_aregister_obj_func_apodization(nifti_image *refimg,niikpt p,double thresh);


/********************************************************
 *
 * display functions
 *
 *
 *********************************************************/

char *niik_aregister_method_string( int reg_method ) {
  switch(reg_method) {
  case NIIK_REGISTER_SSQ:
    return "Sum of squared difference";
  case NIIK_REGISTER_CC:
    return "Correlation coefficient";
  case NIIK_REGISTER_NMI:
    return "Normalized mutual information";
  default:
    return "ERROR: unknown registration method";
  }
  return "ERROR: unknown registration method";
}

int niik_aregister_display_affine( double *affpar )
/* affine parameters are defined at the top of this document */
{
  if(affpar==NULL) {
    fprintf(stderr,"ERROR: affpar is null\n");
    return 0;
  }
  fprintf(stdout,"      r  = %7.3f %7.3f %7.3f \n",affpar[1],affpar[2],affpar[3]); /* in degrees */
  fprintf(stdout,"      t  = %7.3f %7.3f %7.3f \n",affpar[ 4],affpar[ 5],affpar[ 6]);
  fprintf(stdout,"      s  = %7.3f %7.3f %7.3f   includes global scaling of %-7.3f\n",affpar[ 7]*affpar[8],affpar[7]*affpar[ 9],affpar[7]*affpar[10],affpar[7]);
  fprintf(stdout,"      sh = %7.3f %7.3f %7.3f \n",affpar[11],affpar[12],affpar[13]);
  fprintf(stdout,"      c  = %7.2f %7.2f %7.2f \n",affpar[14],affpar[15],affpar[16]);
  return 1;
}


/***************************************************
 *
 * GENERAL REGISTRATION MATRIX
 *
 ***************************************************/

niikmat *niik_aregister_matrix_from_affpar(double *affpar) {
  niikmat *afmat=NULL;
  afmat = niikmat_identity(4,4);
  if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
    fprintf(stderr,"ERROR: niik_aregister_matrix_from_affpar_update\n");
    afmat = niikmat_free(afmat);
    return NULL;
  }
  return afmat;
}

int niik_aregister_matrix_from_affpar_update( niikmat *afmat, double *affpar ) {
  if(afmat ==NULL) {
    fprintf(stderr,"ERROR: afmat is null\n");
    return 0;
  }
  if(affpar==NULL) {
    fprintf(stderr,"ERROR: affpar is null\n");
    return 0;
  }
  if(afmat->row!=4) {
    fprintf(stderr,"ERROR: matrix is not 4-by-4 (%i,%i)\n",afmat->row,afmat->col);
    return 0;
  }
  if(afmat->col!=4) {
    fprintf(stderr,"ERROR: matrix is not 4-by-4 (%i,%i)\n",afmat->row,afmat->col);
    return 0;
  }
  if(!niikmat_affine_matrix(afmat,
                            affpar[ 1],affpar[ 2],affpar[ 3],
                            affpar[ 4],affpar[ 5],affpar[ 6],
                            affpar[ 7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                            affpar[11],affpar[12],affpar[13],
                            affpar[14],affpar[15],affpar[16])) {
    fprintf(stderr,"ERROR: niikmat_affine_matrix\n");
    return 0;
  }
  return 1;
}


/***** end of general functions *****/



/********************************************************
 *
 * generic registration method
 *
 * -refimg   reference image
 * -regseg   reference image's ROI
 * -movimg   moving image
 * -movseg   moving image'as mask, currently not used
 * -affpar   initial affine parameters
 * -perturb  perturbation to the affpar
 *           if zero, then that parameter is fixed
 * -register_method   CC, NMI, ...
 * -sample   sampling distance in mm
 * -filFWHM  gaussian blurring's FWHM
 *
 *
 * see also
 *
 *   double nifti_k_register_obj_func(double *v);
 *
 *
 *********************************************************/

int niik_image_aregister(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,nifti_image *movseg,
                         double *affpar,double *perturb,int register_method,double sample,double filFWHM)
/*
 * nifti_k_affine_register
 *   refimg     -reference target image
 *   refseg     -reference target mask image
 *   movimg     -moving image
 *   movseg     -moving mask image (currently not used)
 *   affpar     -initial affine registration parameters
 *   perturb    -amount of perturbation for each affine parameter
 *              -if zero, then fixed
 *   register_method
 *              -registration method
 *   sample     -sampling distance in mm, smaller -> longer (probably more accurate)
 *   filFWHM    -filtering FWHM
 *
 * 2012-02-26 Kunio
 * -this function became a wrapper
 * -imat for initial matrix is new
 */
{
  if(!niik_image_aregister_imat(refimg,refseg,movimg,movseg,NULL,affpar,perturb,register_method,sample,filFWHM)) {
    fprintf(stderr,"[niik_image_aregister] ERROR: niik_image_aregister_imat\n");
    return 0;
  }
  return 1;
}

int niik_image_aregister_imat(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,nifti_image *movseg,niikmat *imat,
                              double *affpar,double *perturb,int register_method,double sample,double filFWHM)
/*
 * nifti_k_affine_register
 *   refimg     -reference target image
 *   refseg     -reference target mask image
 *   movimg     -moving image
 *   movseg     -moving mask image (currently not used)
 *   imat       -initial matrix
 *   affpar     -initial affine registration parameters
 *   perturb    -amount of perturbation for each affine parameter
 *              -if zero, then fixed
 *   register_method
 *              -registration method
 *   sample     -sampling distance in mm, smaller -> longer (probably more accurate)
 *   filFWHM    -filtering FWHM
 *
 */
{

  char fcname[64]="niik_image_aregister_imat";
  int
  ndim=18,
  flag_keep_gimgs=0,
  maxiter=1e5,
  i,j,m,n,nvar,
  verbose=1;
  niikmat *p;
  double
  plim=1e-7,
  tol=1e-4,(* pfn)();

  if(verbose>1) fprintf(stdout,"  niik_affine_register top\n");
  /* make sure images are present */
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is a null pointer\n");
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"ERROR: movimg is a null pointer\n");
    return 0;
  }
  if(affpar==NULL) {
    fprintf(stderr,"ERROR: affpar is a null pointer\n");
    return 0;
  }
  if(refseg!=NULL) {
    if(refseg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"ERROR: refseg is not NIFTI_TYPE_UINT8\n");
      return 0;
    }
  }
  if(movseg!=NULL) {
    if(movseg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"ERROR: movseg is not NIFTI_TYPE_UINT8\n");
      return 0;
    }
  }
  if(verbose>1) fprintf(stdout,"  niik_affine_register nvar\n");
  for(n=1,nvar=0; n<=16; n++)
    if(fabs(perturb[n])>plim) nvar++;
  /* check number of variables */
  if(nvar==0) {
    fprintf(stderr,"ERROR: no variable \n");
    return 0;
  }
  /* check nmi object */
  if(register_method==NIIK_REGISTER_NMI)  {
    if(verbose>1) fprintf(stdout,"  niik_affine_register check g_nmiobj\n");
    if(g_nmiobj==NULL) {
      fprintf(stderr,"[%s] ERROR: g_nmiobj is undefined\n",fcname);
      return 0;
    }
  }

  if(verbose>=1) {
    fprintf(stdout,"[%s] info\n",fcname);
    fprintf(stdout,"    Gaussian blur   %7.2f mm (FWHM)\n",filFWHM);
    fprintf(stdout,"    sampling dist   %7.4f\n",sample);
    fprintf(stdout,"    DOF             %2i: ",nvar);
    if(fabs(perturb[ 1])>plim) fprintf(stdout,"Rx ");
    if(fabs(perturb[ 2])>plim) fprintf(stdout,"Ry ");
    if(fabs(perturb[ 3])>plim) fprintf(stdout,"Rz ");
    if(fabs(perturb[ 4])>plim) fprintf(stdout,"Tx ");
    if(fabs(perturb[ 5])>plim) fprintf(stdout,"Ty ");
    if(fabs(perturb[ 6])>plim) fprintf(stdout,"Tz ");
    if(fabs(perturb[ 7])>plim) fprintf(stdout,"Sg ");
    if(fabs(perturb[ 8])>plim) fprintf(stdout,"Sx ");
    if(fabs(perturb[ 9])>plim) fprintf(stdout,"Sy ");
    if(fabs(perturb[10])>plim) fprintf(stdout,"Sz ");
    if(fabs(perturb[11])>plim) fprintf(stdout,"Kx ");
    if(fabs(perturb[12])>plim) fprintf(stdout,"Ky ");
    if(fabs(perturb[13])>plim) fprintf(stdout,"Kz ");
    fprintf(stdout,"\n");
    if(imat==NULL) {
      fprintf(stdout,"    no initial matrix\n");
    } else {
      fprintf(stdout,"    using initial matrix\n");
    }
    fprintf(stdout,"    cost function   %s\n",niik_aregister_method_string(register_method));
    if(refseg==NULL) {
      fprintf(stdout,"    no mask\n");
    } else {
      fprintf(stdout,"    mask            %s\n",refseg->fname);
    }
    if(register_method==NIIK_REGISTER_NMI)  {
      fprintf(stdout,"    NMI hist        %5.1f %5.1f %6.1f %-5i (ref)\n",
              g_nmiobj->hmin[0],g_nmiobj->hran[0],g_nmiobj->hmax[0],g_nmiobj->hnum[0]);
      fprintf(stdout,"                    %5.1f %5.1f %6.1f %-5i (mov)\n",
              g_nmiobj->hmin[1],g_nmiobj->hran[1],g_nmiobj->hmax[1],g_nmiobj->hnum[1]);
    }
  } /* show info */


  /*
   * copy images to global variables
   * gaussian filter
   * convert to double
   */

  if(g_imgs[0]==NULL) {
    if(filFWHM<=0) {
      if(verbose>1) fprintf(stdout,"    niik_image_copy\n");
      if((g_imgs[0]=niik_image_copy(refimg))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy\n");
        return 0;
      }
    } else {
      if(verbose>1) fprintf(stdout,"    niik_image_filter_gaussian %i %5.1f [ref]\n",(int)floor(filFWHM),filFWHM);
      if((g_imgs[0]=niik_image_filter_gaussian(refimg,floor(filFWHM),filFWHM))==NULL) {
        fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
        return 0;
      }
    }
  } else {
    flag_keep_gimgs += 1;
    if(verbose)
      fprintf(stdout,"      using pre-computed refimg\n");
  }

  if(g_imgs[1]==NULL) {
    if(filFWHM<=0) {
      if((g_imgs[1]=niik_image_copy(movimg))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy\n");
        return 0;
      }
    } else {
      if(verbose>1) fprintf(stdout,"    niik_image_filter_gaussian %i %5.1f [mov]\n",(int)floor(filFWHM),filFWHM);
      if((g_imgs[1]=niik_image_filter_gaussian(movimg,floor(filFWHM),filFWHM))==NULL) {
        fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
        return 0;
      }
    }
  } else {
    flag_keep_gimgs += 2;
    if(verbose)
      fprintf(stdout,"      using pre-computed movimg\n");
  }

  for(n=0; n<2; n++) {
    if(verbose>1) fprintf(stdout,"    niik_image_type_convert %i\n",n+1);
    if(!niik_image_type_convert(g_imgs[n],NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"ERROR: niik_image_convert %i\n",n);
      return 0;
    }
  }

  /* set other global variables */
  g_imgs[2] = refseg;
  g_imgs[3] = NULL;

  if(0) {
    fprintf(stdout,"writing tmp imgs\n");
    niik_image_write("tmp_frefimg.nii.gz",g_imgs[0]);
    niik_image_write("tmp_fmovimg.nii.gz",g_imgs[1]);
    if(g_imgs[2]!=NULL) niik_image_write("tmp_frefseg.nii.gz",g_imgs[2]);
  }

  /* update the global variables */
  g_reg_method = register_method;
  g_sample = sample;
  g_affpar = affpar;

  /* pre-set the reference data because this is not moving
   * -g_imgs[0] (reference) data is no longer needed */
  if(verbose>1) fprintf(stdout,"  niik_image_aregister_obj_func_set_ref_data\n");
  if((g_ref_data = niik_image_aregister_obj_func_set_ref_data(g_imgs[0],g_imgs[2],sample))==NULL) {
    fprintf(stderr,"ERROR: niik_image_aregister_obj_func_set_ref_data\n");
    return 0;
  }
  if(!(flag_keep_gimgs%2)) {
    /*if(verbose){ fprintf(stdout,"-v1 (niik_image_aregister) free g_imgs[0]->data\n"); }
    free(g_imgs[0]->data);
    g_imgs[0]->data=NULL;*/
  }

  /* the sampling space info */
  n =
    (int)floor(refimg->dx*(refimg->nx-1) / sample + 1) *
    (int)floor(refimg->dy*(refimg->ny-1) / sample + 1) *
    (int)floor(refimg->dz*(refimg->nz-1) / sample + 1);
  if(verbose>1) fprintf(stdout,"  n = %i | %i %i %i\n",n,
                          (int)floor(refimg->dx*refimg->nx / sample + 1),
                          (int)floor(refimg->dy*refimg->ny / sample + 1),
                          (int)floor(refimg->dz*refimg->nz / sample + 1) );
  for(i=j=0; i<n; i++) {
    if(!niik_check_double_problem(g_ref_data[i])) j++;
  }
  if(verbose)
    fprintf(stdout,"    mask size         %i\n",j);

  if(g_fixpar==NULL)  g_fixpar = (int *)calloc(16,sizeof(int));
  for(n=1; n<=16; n++) {
    if(fabs(perturb[n])<=1e-7)  /* no perturbation -> don't change */
      g_fixpar[n] = 1;
    else
      g_fixpar[n] = 0;
  }

  /* prepare for simplex optimization variables */
  if(verbose) fprintf(stdout,"    initialize\n");
  niik_image_aregister_obj_func(NULL);
  pfn=niik_image_aregister_obj_func;
  g_niik_aregister_imat = imat;

  /* update with initial values */
  if(verbose) fprintf(stdout,"    initialize NM matrix\n");
  p=niikmat_init(ndim+1,ndim);
  for(m=0; m<=ndim; m++) {
    for(n=1; n<17; n++) {
      p->m[m][n] = affpar[n];
    }
  }
  for(n=1,m=0; n<17; n++) {
    p->m[m][n] += (fabs(perturb[n])>plim) * p->m[m][n] * niik_get_rand() * 0.01;
  }
  for(m=1; m<17; m++) {
    if(fabs(perturb[m])>plim) {
      p->m[m][m] += perturb[m];
    }
    /*else {
      for(n=0;n<=16;n++) {
      p->m[m][n] += (niik_get_rand()-0.5)*2.0 * perturb[n]; } }*/
  }
  for(m=17; m<=ndim; m++) {
    for(n=0; n<=16; n++) {
      p->m[m][n] += (niik_get_rand()-0.5)*2.0 * perturb[n];
    }
  }

  /* for debug
  fprintf(stdout,"\naffpar = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7.2f ",affpar[m]); }
  fprintf(stdout,"\npeturb = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7.2f ",perturb[m]); }
  fprintf(stdout,"\nfixpar = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7i ",g_fixpar[m]); }
  fprintf(stdout,"\n");
  niikmat_display(p); */

  tol = 1e-4;
  maxiter=100;
  if(verbose>=1) fprintf(stdout,"[%s] start nelder mead\n",fcname);
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
    return 0;
  }

  /* try again with the best parameters */
  if(verbose) {
    fprintf(stdout,"    improved  %9.5f  %i\n",tol,g_iter);
    niik_aregister_display_affine(p->m[0]);
  }
  for(n=0; n<=16; n++) {
    p->m[ndim][n]=affpar[n];
  }
  for(m=1; m<ndim; m++) {
    for(n=1; n<=16; n++) {
      p->m[m][n] += (niik_get_rand()-0.5)*0.5 * perturb[n];
    }
  }

  tol = 1e-6;
  maxiter=2000;
  if(verbose) fprintf(stdout,"    start nelder mead\n");
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
    return 0;
  }

  /* put the registration paramters into the output vector */
  if(verbose>1) fprintf(stdout,"    output vector\n");
  for(n=0; n<=16; n++) {
    affpar[n] = p->m[0][n];
  }
  niikmat_free(p);

  /* update error (saved at index=0) */
  affpar[0] = tol;
  if(verbose>=1) {
    fprintf(stdout,"[%s] optimized %9.5f  %i\n",fcname,tol,g_iter);
    niik_aregister_display_affine(affpar);
  }

  /* free all variables */
  free(g_ref_data);
  g_ref_data = NULL;

  if(verbose>2) fprintf(stdout,"-v2 nifti_image_free\n");
  if(!(flag_keep_gimgs%2)) {
    if(verbose>2) fprintf(stdout,"-v2 (niik_image_aregister) nifti_image_free refimg\n");
    g_imgs[0] = niik_image_free(g_imgs[0]);
  }

  if(!(flag_keep_gimgs/2)) {
    if(verbose>2) fprintf(stdout,"-v2 (niik_image_aregister) nifti_image_free movimg\n");
    g_imgs[1] = niik_image_free(g_imgs[1]);
  }

  free(g_fixpar);
  g_fixpar = NULL;

  if(verbose>=3) niik_fc_display(fcname,0);
  return 1;
} /* niik_affine_register */



/************************************************************
 *
 * registration's objective function
 *
 * -should be used in all affine registration here
 * -need these global variables:
 * --g_imgs[0]     reference image
 * --g_imgs[1]     moving image
 * --g_imgs[2]     reference's mask  (dont need this?)
 * --g_ref_data    sampled reference data array
 * --g_sample      sample distance
 * --g_reg_method  registration method (e.g., CC, NMI ...)
 * --g_affpar      affine parameters (rx, ry, rz, tx, ty, tz, sx, sy, sz, kx, ky, kz, px, py, pz) <- 15 of them
 * --g_fixpar      degrees of freedom (list of fixed parameters)
 * --g_nmiobj      normalized mutual information object (i.e. histogram info)
 * (--g_dof         degrees of freedom -> to be changed to an array of fixed parameters)
 *
 ************************************************************/

double niik_image_aregister_obj_func(double *v) {
  static double val_optim=1e5,init_wsum=0;
  static niikpt fovref,fovmov;
  static double *mov_data=NULL;
  static unsigned char *mov_mask=NULL;
  niikmat
  *afmat = NULL;
  niikpt
  p,q;
  int
  ni[5],nj[2],
  m,n,
  xmax,ymax,zmax,narea,nsize,
  num;
  double
  fi[5],
  dd[5],
  errval,errval2,
  cc,nmi,
  dint1,dint2,
  affpar[25],
  wval,
  xsum,ysum,xssq,yssq,xysum,wsum;
  int verbose=g_verbose_errfunc;
  verbose=1;

  afmat=niikmat_init(4,4);
  p.w=q.w=0;

  /* initialize when inputs are null */
  if(v==NULL) {
    if(verbose>1) fprintf(stdout,"-d (niik_image_aregister) initialization\n");
    fovref.x = g_imgs[0]->pixdim[1] * (g_imgs[0]->dim[1]-1);
    fovref.y = g_imgs[0]->pixdim[2] * (g_imgs[0]->dim[2]-1);
    fovref.z = g_imgs[0]->pixdim[3] * (g_imgs[0]->dim[3]-1);
    if(verbose>1) fprintf(stdout,"  ref fov:   %7.2f %7.2f %7.2f\n",fovref.x,fovref.y,fovref.z);
    fovmov.x = g_imgs[1]->pixdim[1] * (g_imgs[1]->dim[1]-1);
    fovmov.y = g_imgs[1]->pixdim[2] * (g_imgs[1]->dim[2]-1);
    fovmov.z = g_imgs[1]->pixdim[3] * (g_imgs[1]->dim[3]-1);
    if(verbose>1) fprintf(stdout,"  mov fov:   %7.2f %7.2f %7.2f\n",fovmov.x,fovmov.y,fovmov.z);
    val_optim = 1e6;
    g_iter = 0;
    if(mov_data != NULL) {
      free(mov_data);
      free(mov_mask);
      mov_data = NULL;
      mov_mask = NULL;
    }
    if(verbose>1) fprintf(stdout,"niik_aregister_obj_func initialization done\n");
    return 0;
  } /* initialization */


  /* create a matrix */
  if(verbose>1) fprintf(stdout,"-d (niik_aregister_obj_func) create matrix\n");
  for(n=0; n<=17; n++) {
    affpar[n]=0;
    if(g_fixpar[n]) affpar[n] = g_affpar[n];
    else            affpar[n] = v[n];
  }
  if(verbose>1) {
    fprintf(stdout,"-d (niik_aregister_obj_func) affpar: ");
    niik_display_double_vector(affpar,17);
  }

  /* updating matrix in the objective function */
  if(verbose>1) fprintf(stdout,"-d (niik_aregister_obj_func) create matrix from affpar\n");
  if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
    fprintf(stderr,"ERROR: niik_aregister_matrix_from_affpar_update\n");
    return NIIKMAX;
  }

  if(g_niik_aregister_imat!=NULL) {
    if((afmat=niikmat_multiply_free12(g_niik_aregister_imat,afmat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_multiply_free12\n");
      return NIIKMAX;
    }
  }

  if(verbose>2) {
    /*fprintf(stdout,"  method = %i\n",g_reg_method);*/
    fprintf(stdout,"%4i | %5.1f %5.1f %5.1f | %3.0f %3.0f %3.0f | %4.2f %4.2f %4.2f | %4.2f %4.2f %4.2f | %3.0f %3.0f %3.0f \n",g_iter,
            affpar[1],affpar[2],affpar[3],
            affpar[4],affpar[5],affpar[6],
            affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
            affpar[11],affpar[12],affpar[13],
            affpar[14],affpar[15],affpar[16]);
    if(verbose>2) niikmat_display(afmat);
  }

  g_iter++;

  /* normally inverse b/c g_invmat is zero by default */
  if(!g_invmat) {
    if(!niikmat_inverse_update(afmat)) {
      fprintf(stderr,"ERROR: niikmat_inverse_update\n");
      return NIIKMAX;
    }
  }
  if(verbose>2) niikmat_display(afmat);

  /* initialize variables */
  num = 0;
  xsum = ysum = xysum = xssq = yssq = wsum = 0;

  xmax = (int)floor(fovref.x / g_sample + 1);
  ymax = (int)floor(fovref.y / g_sample + 1);
  zmax = (int)floor(fovref.z / g_sample + 1);
  narea = xmax * ymax;
  nsize = narea * zmax;

  switch(g_reg_method) {

  case NIIK_REGISTER_SAD:      /* Sum of absolute difference  */
    return (double)NIIKMAX;
    return xsum;

  case NIIK_REGISTER_SSQ:      /* Sum of squared difference */
    return (double)NIIKMAX;
    return xssq;

  case NIIK_REGISTER_NMI:    /* NMI registration */

    if(mov_data==NULL) mov_data = (double *)calloc(nsize*zmax,sizeof(double));
    if(mov_mask==NULL) mov_mask = (unsigned char *)calloc(nsize*zmax,sizeof(char));
    wsum = wval = 0;

    #pragma omp parallel for private(p,q) reduction(+:wsum) reduction(+:wval)
    for(m=0; m<nsize; m++) {
      p.x = (m%xmax)        * g_sample;
      p.y = ((m/xmax)%ymax) * g_sample;
      p.z = (m/narea)       * g_sample;
      p.w = q.w = 0;

      mov_data[m] = 0;
      mov_mask[m] = 0;
      if(niik_check_double_problem(g_ref_data[m])) continue;

      /* transform and check with fov, each direction separately */
      q.x=afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0)         continue;
      if(q.x > fovmov.x) continue;
      q.y=afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0)         continue;
      if(q.y > fovmov.y) continue;
      q.z=afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0)         continue;
      if(q.z > fovmov.z) continue;

      /* dint1 = niik_image_interpolate_3d_linear(g_imgs[0],p);*/
      mov_data[m] = niik_image_interpolate_float_image_3d_linear(g_imgs[1],q);
      mov_mask[m] = 1;
    } /* omp loop */

    /* initialize histograms */
    for(n=0; n<2; n++)
      for(m=0; m<g_nmiobj->hnum[1]; m++)
        g_nmiobj->h1[n][m]=0;
    for(n=0; n<g_nmiobj->hnum[0]; n++)
      for(m=0; m<g_nmiobj->hnum[1]; m++)
        g_nmiobj->h2[n][m]=0;
    wval = wsum = 0;

    for(n=0,num=0; n<nsize; n++) {
      if(mov_mask[n]==0) continue;
      /*if(niik_check_double_problem(  mov_data[n])) continue;*/
      if(niik_check_double_problem(g_ref_data[n])) continue;
      wsum += 1;
      num++;

      /* histogram positions */
      fi[0] = g_ref_data[n] / g_nmiobj->hran[0];
      fi[1] =   mov_data[n] / g_nmiobj->hran[1];

      if(0) {
        for(m=0; m<2; m++) {
          ni[m] = floor(fi[m]+0.5);
          if(ni[m]<0) ni[m]=0;
          else if(ni[m]>=g_nmiobj->hnum[m]) ni[m]=g_nmiobj->hnum[m]-1;
          g_nmiobj->h1[m][ni[m]] += 1.0;
        }
        g_nmiobj->h2[ni[0]][ni[1]] += 1.0;
      }

      else {
        /* 1d histogram -weighted histogram similar to Jenkinson 2002 NeuroImage */
        for(m=0; m<2; m++) {
          ni[m] = floor(fi[m]);
          nj[m] = ni[m]+1;
          dd[m] = fi[m] - ni[m];
          if(ni[m]>=g_nmiobj->hnum[m]) ni[m]=-1;
          if(nj[m]>=g_nmiobj->hnum[m]) nj[m]=-1;
          if(ni[m]>=0) g_nmiobj->h1[m][ni[m]] += 1.0 - dd[m];
          if(nj[m]>=0) g_nmiobj->h1[m][nj[m]] += dd[m];
        }
        /* 2d histogram -weighted histogram */
        if(ni[0]>=0) {
          if(ni[1]>=0) {
            g_nmiobj->h2[ni[0]][ni[1]] += (1.0-dd[0]) * (1.0-dd[1]);
          }
          if(nj[1]>=0) {
            g_nmiobj->h2[ni[0]][nj[1]] += (1.0-dd[0]) * dd[1];
          }
        }
        if(nj[0]>=0) {
          if(ni[1]>=0) {
            g_nmiobj->h2[nj[0]][ni[1]] += dd[0] * (1.0-dd[1]);
          }
          if(nj[1]>=0) {
            g_nmiobj->h2[nj[0]][nj[1]] += dd[0] * dd[1];
          }
        }
      }

      /*fprintf(stdout,"%6.2f %6.2f %6.2f %6.2f  ->  %8.3f %8.3f %8.3f %8.3f -> %8.4f\n",
        dd[0],dd[1],dd[2],dd[3],
        dd[0]*dd[1],dd[0]*dd[3],dd[2]*dd[1],dd[2]*dd[3],
        dd[0]*dd[1]+dd[0]*dd[3]+dd[2]*dd[1]+dd[2]*dd[3]);*/
    } /* each sample point */

    /* entropy calculation */
    fi[0]=fi[1]=0;
    for(n=0; n<2; n++)
      for(m=0; m<g_nmiobj->hnum[1]; m++) {
        g_nmiobj->h1[n][m]/=num;
        if(g_nmiobj->h1[n][m]<1e-5) continue;
        fi[0]-=g_nmiobj->h1[n][m]*log10(g_nmiobj->h1[n][m]);
      }

    for(n=0; n<g_nmiobj->hnum[0]; n++)
      for(m=0; m<g_nmiobj->hnum[1]; m++) {
        g_nmiobj->h2[n][m]/=num;
        if(g_nmiobj->h2[n][m]<1e-5) continue;
        fi[1]-=g_nmiobj->h2[n][m]*log10(g_nmiobj->h2[n][m]);
      }

    /* nmi calculation */
    nmi=fi[0]/(fi[1]+1e-6);

    if(g_iter<=1) {
      init_wsum = wsum;
    }

    /*errval =  2.0 * NIIK_Heaviside(0.3 - wsum / num,0.1);*/
    errval = errval2 = 2.0 * NIIK_Heaviside(0.75-wsum/init_wsum,0.05);
    errval += (2.0 - nmi);

    if(verbose>=2) {
      fprintf(stdout,"e %6.3f | nmi %6.3f %6.3f %6.3f | %4i | %8i | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.3f %5.3f %5.3f | %5.2f %5.2f %5.2f | %5.1f %5.1f %5.1f\n",
              errval2,nmi,
              fi[0],fi[1],
              g_iter,num,
              affpar[1],affpar[2],affpar[3],
              affpar[4],affpar[5],affpar[6],
              affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
              affpar[11],affpar[12],affpar[13],
              affpar[14],affpar[15],affpar[16]);
    }

    if(val_optim > errval) {
      val_optim = errval;
      if(verbose>=1) {
        fprintf(stdout,"e %6.3f | nmi %6.3f %6.3f %6.3f | %4i | %8i | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.3f %5.3f %5.3f | %5.2f %5.2f %5.2f | %5.1f %5.1f %5.1f *\n",
                errval2,nmi,fi[0],fi[1],g_iter,num,
                affpar[1],affpar[2],affpar[3],
                affpar[4],affpar[5],affpar[6],
                affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                affpar[11],affpar[12],affpar[13],
                affpar[14],affpar[15],affpar[16]);
      }
    }

    niikmat_free(afmat);
    return errval;


  case NIIK_REGISTER_CC:    /* CC registration, aka normalized correlation (Jenkinson 2002 NeuroImage) */

    num = 0;

    #pragma omp parallel for private(p,q,dint1,dint2,wval) reduction(+:num) reduction(+:xsum) reduction(+:ysum) reduction(+:xssq) reduction(+:yssq) reduction(+:xysum) reduction(+:wsum)
    for(m=0; m<nsize; m++) {
      p.x = (m%xmax)        * g_sample;
      p.y = ((m/xmax)%ymax) * g_sample;
      p.z = floor(m/narea)  * g_sample;
      p.w = q.w = 0;

      if(niik_check_double_problem(g_ref_data[m])) continue;
      num++;

      /* transform and check with fov, each direction separately */
      q.x = afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0) continue;
      if(q.x > fovmov.x) continue;

      q.y = afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0) continue;
      if(q.y > fovmov.y) continue;

      q.z = afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0) continue;
      if(q.z > fovmov.z) continue;

      dint1 = g_ref_data[m];
      dint2 = niik_image_interpolate_float_image_3d_linear(g_imgs[1],q);

      /* apodization
       * as in Jenkinson 2002 NeuroImage */
      wval = 1; /*niik_image_aregister_obj_func_apodization(g_imgs[1],q,apodization_thresh);  */
      if(fabs(wval)<1e-5) continue;

      xsum  += dint1 * wval;
      ysum  += dint2 * wval;
      xssq  += dint1*dint1 * wval;
      yssq  += dint2*dint2 * wval;
      xysum += dint1*dint2 * wval;
      wsum  += wval;
    } /* each voxel */

    if(verbose>2) {
      fprintf(stdout,"  xsum  = %12.6f\n",xsum);
      fprintf(stdout,"  ysum  = %12.6f\n",ysum);
      fprintf(stdout,"  xssq  = %12.6f\n",xssq);
      fprintf(stdout,"  yssq  = %12.6f\n",yssq);
      fprintf(stdout,"  xysum = %12.6f\n",xysum);
      fprintf(stdout,"  wsum  = %12.6f\n",wsum);
      fprintf(stdout,"  nsum  = %12i\n",num);
    }

    cc = (wsum * xysum - xsum * ysum) / sqrt(wsum*xssq-xsum*xsum) / sqrt(wsum*yssq-ysum*ysum);
    errval =  3.0 * NIIK_Heaviside(0.3 - wsum / num,0.2);
    errval += (1.0 - cc);

    if(verbose>=2) {
      fprintf(stdout,"   err %-6.4f cc %-6.4f | %4i | %6i | %6.3f %6.3f %6.3f | %6.2f %6.2f %6.2f | %5.3f %5.3f %5.3f | %6.3f %6.3f %6.3f | %5.1f %5.1f %5.1f\n",
              3.0*NIIK_Heaviside(0.3 - wsum / num,0.2),cc,g_iter,(int)wsum,
              affpar[1],affpar[2],affpar[3],
              affpar[4],affpar[5],affpar[6],
              affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
              affpar[11],affpar[12],affpar[13],
              affpar[14],affpar[15],affpar[16]);
    }

    if(val_optim > errval) {
      val_optim = errval;
      if(verbose>=1) {
        fprintf(stdout,"   err %-6.4f cc %-6.4f | %4i | %6i | %6.3f %6.3f %6.3f | %6.2f %6.2f %6.2f | %5.3f %5.3f %5.3f | %6.3f %6.3f %6.3f | %5.1f %5.1f %5.1f *\n",
                3.0*NIIK_Heaviside(0.3 - wsum / num,0.2),cc,g_iter,(int)wsum,
                affpar[1],affpar[2],affpar[3],
                affpar[4],affpar[5],affpar[6],
                affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                affpar[11],affpar[12],affpar[13],
                affpar[14],affpar[15],affpar[16]);
      }
    }
    niikmat_free(afmat);
    return errval;

  default:
    fprintf(stderr,"ERROR: unkown registration type, %s\n",niik_aregister_method_string(g_reg_method));
    return NIIKMAX;
  }

  return NIIKMAX;
}


/*
 * apodization
 *   Reference = Jenkinson 2002 NeuroImage
 *
 * -refimg is the reference image for finding the borders
 * -p is the point, has a unit of mm, and to be compared with ijk*pixdim[?]
 *  and does not involve qform
 * -thresh is the threshold distance
 *
*/


double niik_image_aregister_obj_func_apodization(nifti_image *refimg,niikpt p,double thresh) {
  double dval;
  if(thresh<=1e-4) return 1;
  if(p.x<=0) return 0;
  if(p.y<=0) return 0;
  if(p.z<=0) return 0;
  if(p.x >= (refimg->nx-1)*(refimg->dx)) return 0;
  if(p.y >= (refimg->ny-1)*(refimg->dy)) return 0;
  if(p.z >= (refimg->nz-1)*(refimg->dz)) return 0;
  dval  = (p.x>=thresh)?1:p.x/thresh;
  dval *= (p.y>=thresh)?1:p.y/thresh;
  dval *= (p.z>=thresh)?1:p.z/thresh;
  p.x = (refimg->nx-1)*(refimg->dx) - p.x;
  p.y = (refimg->ny-1)*(refimg->dy) - p.y;
  p.z = (refimg->nz-1)*(refimg->dz) - p.z;
  dval *= (p.x>=thresh)?1:p.x/thresh;
  dval *= (p.y>=thresh)?1:p.y/thresh;
  dval *= (p.z>=thresh)?1:p.z/thresh;
  return dval;
}



/********************************************************
 *
 * niik_image_aregister_obj_func_set_ref_data
 *
 * -interpolate reference image for faster registration
 *
 *
 ********************************************************/

double *niik_image_aregister_obj_func_set_ref_data(nifti_image *refimg, nifti_image *maskimg, double sample) {
  int
  verbose=0,
  i,j,k,
  masksize=0,
  xmax,ymax,zmax,
  m,nsize;
  niikpt fovref,p;
  double *dimg;
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is a null pointer\n");
    return NULL;
  }
  if(refimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"ERROR: refimg is not a float type\n");
    return NULL;
  }
  if(maskimg!=NULL) {
    if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"ERROR: maskimg is not a uint8 type\n");
      return NULL;
    }
  }
  fovref.x = refimg->pixdim[1] * (refimg->dim[1]-1);
  fovref.y = refimg->pixdim[2] * (refimg->dim[2]-1);
  fovref.z = refimg->pixdim[3] * (refimg->dim[3]-1);
  xmax = (int)floor(fovref.x / sample + 1);
  ymax = (int)floor(fovref.y / sample + 1);
  zmax = (int)floor(fovref.z / sample + 1);
  /* fprintf(stdout,"xyz = %i %i %i \n",xmax,ymax,zmax); */
  /* calculate the size of output vector
   * -this loop must be the same as used in the actual
   *  niik_image_aregister_obj_func */
  nsize = xmax*ymax*zmax;
  if(verbose) fprintf(stdout,"  n = %i | %i %i %i\n",nsize,xmax,ymax,zmax);
  dimg = (double *)calloc(nsize,sizeof(double));
  p.w=0;
  if(maskimg!=NULL) { /* with mask image */
    for(k=m=0; k<zmax; k++) {
      p.z = k*sample;
      for(j=0; j<ymax; j++) {
        p.y = j*sample;
        for(i=0; i<xmax; m++,i++) {
          p.x = i*sample;
          if(niik_image_interpolate_uint8_image_3d_nn(maskimg,p.x,p.y,p.z)==0) {
            dimg[m] = NIIKMAX;
          } /* set it to a crazy value */
          else {
            masksize++;
            dimg[m] = niik_image_interpolate_float_image_3d_linear(refimg,p);
          }
        }
      }
    }
    /*fprintf(stdout,"  masksize = %i / %i %.4f | ds = %.4f\n",masksize,xmax*ymax*zmax,(double)masksize/nsize,sample);*/
  } /* with mask image */
  else { /* no mask image */
    for(k=m=0; k<zmax; k++) {
      p.z = k*sample;
      for(j=0; j<ymax; j++) {
        p.y = j*sample;
        for(i=0; i<xmax; m++,i++) {
          p.x = i*sample;
          dimg[m] = niik_image_interpolate_float_image_3d_linear(refimg,p);
        }
      }
    }
  } /* no mask image */
  return dimg;
} /* niik_image_aregister_obj_func_set_ref_data */



/***** end of generic registration *****/














/************************************************************
 *
 * functions for nmi_obj
 *
 *   nmi_obj *nifti_k_register_nmi_obj_alloc();
 *   void nifti_k_register_nmi_obj_free(nmi_obj *obj);
 *   nmi_obj *nifti_k_register_nmi_obj_init(double min1,double max1,double min2,double max2,int num1,int num2);
 *   int nifti_k_register_nmi_obj_update(nmi_obj *obj,double min1,double max1,double min2,double max2,int num1,int num2);
 *   int nifti_k_register_nmi_obj_update_var(double min1,double max1,double min2,double max2,int num1,int num2);
 *   int nifti_k_register_nmi_obj_display_var();
 *
 ************************************************************/

nmi_obj *niik_aregister_nmi_obj_alloc() {
  nmi_obj *obj;
  obj = (nmi_obj *)calloc(1,sizeof(nmi_obj));
  obj->nmi = 0;
  obj->hmin[0]=obj->hmin[1]=0;
  obj->hmax[0]=obj->hmax[1]=0;
  obj->hran[0]=obj->hran[1]=0;
  obj->h1=obj->h2=NULL;
  obj->hnum[0]=obj->hnum[1]=0;
  return obj;
}

void niik_aregister_nmi_obj_free(nmi_obj *obj) {
  int n;
  if(obj==NULL) return;
  if(obj->h1!=NULL) {
    free(obj->h1[0]);
    free(obj->h1[1]);
    free(obj->h1);
  }
  if(obj->h2!=NULL) {
    for(n=0; n<obj->hnum[0]; n++) {
      free(obj->h2[n]);
    }
    free(obj->h2);
  }
  free(obj);
}

nmi_obj *niik_aregister_nmi_obj_init(double min1,double max1,double min2,double max2,int num1,int num2) {
  nmi_obj *obj;
  int n;
  obj = niik_aregister_nmi_obj_alloc();
  obj->hnum[0]=num1;
  obj->hnum[1]=num2;
  obj->hmin[0]=min1;
  obj->hmin[1]=min2;
  obj->hmax[0]=max1;
  obj->hmax[1]=max2;
  obj->hran[0]=(max1-min1)/(num1-1);
  obj->hran[1]=(max2-min2)/(num2-1);
  obj->h1=(double **)calloc(   2,sizeof(double *));
  obj->h2=(double **)calloc(num1,sizeof(double *));
  for(n=0; n<   2; n++) {
    obj->h1[n]=(double *)calloc(num2,sizeof(double));
  }
  for(n=0; n<num1; n++) {
    obj->h2[n]=(double *)calloc(num2,sizeof(double));
  }
  return obj;
}

int niik_aregister_nmi_obj_update(nmi_obj *obj,double min1,double max1,double min2,double max2,int num1,int num2) {
  int n;
  if(obj==NULL) return 0;
  if(obj->h1!=NULL) {
    free(obj->h1[0]);
    free(obj->h1[1]);
    free(obj->h1);
  }
  if(obj->h2!=NULL) {
    for(n=0; n<obj->hnum[0]; n++) {
      free(obj->h2[n]);
    }
    free(obj->h2);
  }
  obj->hnum[0]=num1;
  obj->hnum[1]=num2;
  obj->hmin[0]=min1;
  obj->hmin[1]=min2;
  obj->hmax[0]=max1;
  obj->hmax[1]=max2;
  obj->hran[0]=(max1-min1)/(num1-1);
  obj->hran[1]=(max2-min2)/(num2-1);
  obj->h1=(double **)calloc(   2,sizeof(double *));
  obj->h2=(double **)calloc(num1,sizeof(double *));
  for(n=0; n<   2; n++) {
    obj->h1[n]=(double *)calloc(num2,sizeof(double));
  }
  for(n=0; n<num1; n++) {
    obj->h2[n]=(double *)calloc(num2,sizeof(double));
  }
  return 1;
}

int niik_aregister_nmi_obj_update_var(double min1,double max1,double min2,double max2,int num1,int num2) {
  if(g_nmiobj==NULL) {
    g_nmiobj=niik_aregister_nmi_obj_init(min1,max1,min2,max2,num1,num2);
    return 1;
  }
  if(!niik_aregister_nmi_obj_update(g_nmiobj,min1,max1,min2,max2,num1,num2)) {
    fprintf(stderr,"ERROR: niik_aregister_nmi_obj_update\n");
    return 0;
  }
  return 1;
}

int niik_aregister_nmi_obj_display_var() {
  char fcname[64]="niik_aregister_nmi_obj_display_var";
  if(g_nmiobj==NULL) {
    fprintf(stdout,"[%s] ERROR: g_nmiobj is null\n",fcname);
    return 0;
  }
  if(!niik_aregister_nmi_obj_display(g_nmiobj)) {
    fprintf(stderr,"[%s] ERROR: niik_aregister_nmi_obj_display\n",fcname);
    return 0;
  }
  return 1;
}

int niik_aregister_nmi_obj_display(nmi_obj *nmiobj) {
  char fcname[64]="niik_aregister_nmi_obj_display";
  if(nmiobj==NULL) {
    fprintf(stdout,"[%s] ERROR: nmiobj is null\n",fcname);
    return 0;
  }
  fprintf(stdout,"    movimg:  %5.2f %5.2f %5.2f %5i\n",nmiobj->hmin[0],nmiobj->hran[0],nmiobj->hmax[0],nmiobj->hnum[0]);
  fprintf(stdout,"    refimg:  %5.2f %5.2f %5.2f %5i\n",nmiobj->hmin[1],nmiobj->hran[1],nmiobj->hmax[1],nmiobj->hnum[1]);
  return 1;
}

nmi_obj *niik_aregister_g_nmi_obj_get() {
  return g_nmiobj;
}
int niik_aregister_g_nmi_obj_set(nmi_obj *nmiobj) {
  g_nmiobj=nmiobj;
  return 1;
}

/***** end of nmi object ******/








/***************************************************
 *
 * STANDARD SPACE REGISTRATION FUNCTION
 *
 * int niik_aregister_align_mni(nifti_image *mni_img, nifti_image *mni_seg, nifti_image *img, nifti_image *seg,
 *     double *affpar, int cost_method, double filFWHM, double sample);
 *
 *
 ***************************************************/

niikpt g_niik_aregister_align_mni_maximize_symmetry_center;

int niik_aregister_align_mni(nifti_image *mni_img, nifti_image *mni_seg, nifti_image *img, nifti_image *seg, double *affpar,
                             int cost_method, double filFWHM, double sample)
/* main function for starndard space affine registration
 *
 * 1. symmetry maximization: determines ry and rz with a image's center
 * 2. translation (tx,ty,tz) and adjust the pivot point (for better registration, dof)
 * 3. global scaling
 * 4. rigid-body registartion with fixed ry, rz, and sg
 * 5. scaling and shearing with fixed rx,ry,rz,sg,sx,sy,sz
 */

{
  const char *fcname="niik_aregister_align_mni";
  nifti_image
  *tmpseg=NULL,
   *tmp_mni_seg=NULL,
    *tmpimg=NULL;
  int
  m,n,
  gs_num=6,
  nmi_num[2],
  ndim,nm_maxiter;      /* symmetry-maximization function */
  niikpt
  imgctr,mnictr;
  niikmat
  *gs_par=NULL,
   *afmat=NULL,
    *p;
  double
  thresh=-1,
  *daffpar=NULL,
   step,
   nmi_max[2],nmi_p99[2],
   nmi_min[2],nmi_p01[2],
   *gs_list=NULL,     /* global scaling optimization */
    gs_limit=0.2,      /* global scaling optimization */
    tol=1e-5,
    (* pfn)();         /* symmetry-maximization function */
  char
  fname[512];


  niik_fc_display(fcname,1);
  NIIK_RET0((mni_img==NULL),fcname,"mni_img is a null pointer");
  NIIK_RET0((    img==NULL),fcname,"img is a null pointer");
  NIIK_RET0(( affpar==NULL),fcname,"affpar is a null pointer");
  fprintf(stdout,"    cost function              %s\n",niik_aregister_method_string(cost_method));
  fprintf(stdout,"    base sampling distance     %8.3f\n",sample);
  fprintf(stdout,"    base filter FWHM           %8.3f\n",filFWHM);
  fprintf(stdout,"    initial parameters\n");
  NIIK_RET0((!niik_aregister_display_affine(affpar)),fcname,"niik_aregister_display_affine");

  if(mni_seg!=NULL) {
    if(!niik_image_type_convert(mni_seg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"ERROR: niik_image_type_convert(mni_seg,NIFTI_TYPE_UINT8)\n");
      return 0;
    }
  }
  if(!niik_image_type_convert(mni_img,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(mni_img,NIFTI_TYPE_FLOAT32)\n");
    return 0;
  }
  if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)\n");
    return 0;
  }
  if(seg!=NULL) {
    if(!niik_image_type_convert(seg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"ERROR: niik_image_type_convert(seg,NIFTI_TYPE_UINT8)\n");
      return 0;
    }
    tmp_mni_seg=mni_seg;
  } else {
    tmp_mni_seg=NULL;
    fprintf(stdout,"[%s]   not using MNI mask\n",fcname);
  }

  /* check brain mask */
  if(seg!=NULL) {
    if(niik_image_count_mask(seg)==0) {
      fprintf(stderr,"ERROR: no volume in the mask\n");
      return 0;
    }
  }
  if(mni_seg!=NULL) {
    if(niik_image_count_mask(mni_seg)==0) {
      fprintf(stderr,"ERROR: no volume in the MNI mask\n");
      return 0;
    }
  }

  if(seg==NULL || mni_seg==NULL) { /* if at least one of the masks is missing */
    NIIK_RET0((!niik_image_thresh_otsu(img,NULL,&thresh)),fcname,"niik_image_thresh_otsu");
    fprintf(stdout,"[%s] using threshold (%.2f) for center calculation\n",fcname,thresh);
    NIIK_RET0(((tmpseg=niik_image_threshold_new(img,thresh))==NULL),fcname,"niik_image_threshold_new");
    imgctr=niikpt_image_get_centroid(img,tmpseg);
    tmpseg=niik_image_free(tmpseg);
    if(niik_check_double_problem(imgctr.x)) {
      fprintf(stderr,"ERROR: niikpt_image_get_centroid\n");
      return 0;
    }
    fprintf(stdout,"\timg center:   %7.3f %7.3f %7.3f    %7.3f %7.3f %7.3f\n",
            imgctr.x,imgctr.y,imgctr.z,
            imgctr.x/img->dx,imgctr.y/img->dy,imgctr.z/img->dz);
    mnictr=niikpt_image_get_centroid(mni_img,NULL);
    if(niik_check_double_problem(mnictr.x)) {
      fprintf(stderr,"ERROR: niikpt_image_get_centroid\n");
      return 0;
    }
    fprintf(stdout,"\tmni center:   %7.3f %7.3f %7.3f    %7.3f %7.3f %7.3f\n",
            mnictr.x,mnictr.y,mnictr.z,
            mnictr.x/mni_img->dx,mnictr.y/mni_img->dy,mnictr.z/mni_img->dz);
  } else { /* both masks exist */
    imgctr=niikpt_image_get_centroid(img,seg);
    if(niik_check_double_problem(imgctr.x)) {
      fprintf(stderr,"ERROR: niikpt_image_get_centroid\n");
      return 0;
    }
    fprintf(stdout,"\timg center:   %7.3f %7.3f %7.3f\n",imgctr.x,imgctr.y,imgctr.z);
    mnictr=niikpt_image_get_centroid(mni_img,mni_seg);
    if(niik_check_double_problem(mnictr.x)) {
      fprintf(stderr,"ERROR: niikpt_image_get_centroid\n");
      return 0;
    }
    fprintf(stdout,"\tmni center:   %7.3f %7.3f %7.3f\n",mnictr.x,mnictr.y,mnictr.z);
  }

  NIIK_RET0(((afmat=niikmat_init(4,4))==NULL),fcname,"niikmat_init");
  NIIK_RET0(((daffpar = (double *)calloc(20,sizeof(double)))==NULL),fcname,"calloc");


  /*
   *  1. symmetry maximization
   */
  fprintf(stdout,"    1. symmetry maximization\n");
  g_niik_aregister_align_mni_maximize_symmetry_center = imgctr;
  niik_aregister_align_mni_maximize_symmetry(NULL);
  ndim=12;
  p=niikmat_init(ndim+1,ndim);
  for(m=1; m<ndim+1; m++) {
    for(n=0; n<ndim; n++) {
      p->m[m][n] = niik_get_rand() * 60.0 - 30.0;
    }
  }
  pfn=niik_aregister_align_mni_maximize_symmetry;
  nm_maxiter = 700;
  if(filFWHM<=0) {
    fprintf(stdout,"\tno gaussian filter\n");
    NIIK_RET0(((g_niik_aregister_movimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  } else {
    fprintf(stdout,"\tapply gaussian filter %5.2f\n",filFWHM);
    NIIK_RET0(((g_niik_aregister_movimg=niik_image_filter_gaussian(img,floor(filFWHM+1.5),filFWHM))==NULL),fcname,
              "niik_image_filter_gaussian");
  }

  fprintf(stdout,"    optimization (symmetry-maximization)\n");
  NIIK_RET0((!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&nm_maxiter)),fcname,"niik_nelder_mead");

  fprintf(stdout,"     symmetry-maximization: r = %7.3f %7.3f %7.3f | iter %i | error %9.5f\n",0.0f,p->m[0][0],p->m[0][1],nm_maxiter,tol);
  affpar[2]=p->m[0][0];
  affpar[3]=p->m[0][1];
  p = niikmat_free(p);

  if(0) {
    niik_aregister_matrix_from_affpar_update(afmat,affpar);
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    sprintf(fname,"tmp_reg1.nii.gz");
    fprintf(stdout,"    transform image for checking: %s\n",fname);
    if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      exit(0);
    }
    niik_image_write(fname,tmpimg);
    tmpimg = niik_image_free(tmpimg);
    exit(0);
  }


  /* NMI object
     nmi_num[0]=nmi_num[1]=32;
     nmi_num[0]*=niik_image_get_percentile(    img,NULL,0.99)/niik_image_get_percentile(    img,NULL,0.9)/2.0;
     nmi_num[1]*=niik_image_get_percentile(mni_img,NULL,0.99)/niik_image_get_percentile(mni_img,NULL,0.9)/2.0;
     nmi_max[0] = niik_image_get_percentile(    img,NULL,0.99);
     nmi_max[1] = niik_image_get_percentile(mni_img,NULL,0.99); */

  nmi_num[0] = nmi_num[1] = 32;
  nmi_max[0] = niik_image_get_max(    img,NULL);
  nmi_max[1] = niik_image_get_max(mni_img,NULL);
  nmi_min[0] = niik_image_get_min(    img,NULL);
  nmi_min[1] = niik_image_get_min(mni_img,NULL);
  nmi_p99[0] = niik_image_get_percentile(    img,NULL,0.99);
  nmi_p99[1] = niik_image_get_percentile(mni_img,NULL,0.99);
  nmi_p01[0] = niik_image_get_percentile(    img,NULL,0.01);
  nmi_p01[1] = niik_image_get_percentile(mni_img,NULL,0.01);
  nmi_num[0]*=(nmi_max[0]-nmi_min[0])/(nmi_p99[0]-nmi_p01[0])/2.0;
  nmi_num[1]*=(nmi_max[1]-nmi_min[1])/(nmi_p99[1]-nmi_p01[1])/2.0;
  NIIK_RET0(((g_nmiobj=niik_aregister_nmi_obj_init(nmi_min[1],nmi_max[1],
                       nmi_min[0],nmi_max[0],
                       nmi_num[1],nmi_num[0]))==NULL),
            fcname,
            "niik_aregister_nmi_obj_init");

  /*
   * 2. translation and adjust pivot point
   */

  /* parameters are not initialized
   * -update the pivot point  */
  fprintf(stdout,"    2. translation based on centers\n");
  if(affpar[0]==0.0) {
    /*affpar[14]=(imgctr.x+mnictr.x)/2.0;
      affpar[15]=(imgctr.y+mnictr.y)/2.0;
      affpar[16]=(imgctr.z+mnictr.z)/2.0;*/
    affpar[14]=imgctr.x;
    affpar[15]=imgctr.y;
    affpar[16]=imgctr.z;
  }
  /* translation */
  affpar[4]=(mnictr.x-imgctr.x);
  affpar[5]=(mnictr.y-imgctr.y);
  affpar[6]=(mnictr.z-imgctr.z);
  fprintf(stdout,"  image center-based translation  %7.3f %7.3f %7.3f\n",affpar[4],affpar[5],affpar[6]);

  if(0) {
    niik_aregister_matrix_from_affpar_update(afmat,affpar);
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    sprintf(fname,"tmp_reg2.nii.gz");
    fprintf(stdout,"    transform image for checking: %s\n",fname);
    if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      exit(0);
    }
    niik_image_write(fname,tmpimg);
    tmpimg = niik_image_free(tmpimg);
    exit(0);
  }

  if(!niik_aregister_display_affine(affpar)) {
    fprintf(stderr,"ERROR: niik_aregister_display_affine\n");
    return 0;
  }


  /*
   * 2a. rotation and translation
   *
   * 2013-01-21, Kunio
   * -added to do z-translation and x-rotation so that scaling can start with better initial registration
   */
  fprintf(stdout,"    2a. rotation and translation\n");
  for(m=0; m<20; m++) {
    daffpar[m]=0;
  }
  daffpar[1]=10;
  daffpar[6]=10;
  if(!niik_image_aregister(mni_img,tmp_mni_seg,g_niik_aregister_movimg,seg,affpar,daffpar,cost_method,sample,filFWHM)) {
    fprintf(stderr,"ERROR: niik_image_aregister\n");
    return 0;
  }

  if(0) {
    niik_aregister_matrix_from_affpar_update(afmat,affpar);
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    sprintf(fname,"tmp_reg2a.nii.gz");
    fprintf(stdout,"    transform image for checking: %s\n",fname);
    if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      exit(0);
    }
    niik_image_write(fname,tmpimg);
    tmpimg = niik_image_free(tmpimg);
    exit(0);
  }



  /*
   * 3. global scaling
   */

  fprintf(stdout,"    3. global scaling\n");
  gs_num=2;
  gs_limit=0.2;
  step = gs_limit / (gs_num-1.0) / (1.0+1.0/(gs_num-1));
  gs_num=gs_num*2+1;
  gs_list = (double *)calloc(gs_num,sizeof(double));
  fprintf(stdout,"\t  step = %8.4f\n",step);
  for(n=0; n<gs_num; n++) {
    gs_list[n] = (n-gs_num/2)*step + 1.0;
    /* fprintf(stdout,"\t  %i  gs = %8.4f\n",n,gs_list[n]); */
  }
  gs_par=niikmat_init(gs_num,20);
  g_imgs[0]=g_imgs[1]=g_imgs[2]=g_imgs[3]=NULL;

  /* blur images for faster processing */
  if((g_imgs[0]=niik_image_filter_gaussian(mni_img,filFWHM*2,filFWHM*3))==NULL) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
    return 0;
  }
  if((g_imgs[1]=niik_image_filter_gaussian(    img,filFWHM*2,filFWHM*2))==NULL) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
    return 0;
  }
  /* global scaling testing */
  for(n=0; n<gs_num; n++) {
    for(m=0; m<17; m++) {
      daffpar[m]=0;
      gs_par->m[n][m]=affpar[m];
    }
    daffpar[1]=30;
    daffpar[4]=10;
    daffpar[5]=10;
    daffpar[6]=10;
    gs_par->m[n][7]=gs_list[n];
    fprintf(stdout,"\n\ttesting global scaling = %5.2f\n",gs_list[n]);
    if(mni_seg==NULL || seg==NULL) {
      fprintf(stdout,"[%s] not using mask\n",fcname);
      if(!niik_image_aregister(mni_img,NULL,img,NULL,gs_par->m[n],daffpar,cost_method,sample*2.2,filFWHM*3)) {
        fprintf(stderr,"ERROR: niik_image_aregister\n");
        return 0;
      }
    } else {
      fprintf(stdout,"[%s] using both masks\n",fcname);
      if(!niik_image_aregister(mni_img,mni_seg,img,seg,gs_par->m[n],daffpar,cost_method,sample*2.2,filFWHM*3)) {
        fprintf(stderr,"ERROR: niik_image_aregister\n");
        return 0;
      }
    }

    if(0) {
      niik_aregister_matrix_from_affpar_update(afmat,gs_par->m[n]);
      niikmat_display(afmat);
      niik_aregister_display_affine(gs_par->m[n]);
      sprintf(fname,"tmp%i.nii.gz",n+1);
      fprintf(stdout,"    transform image for checking: %s\n",fname);
      if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
        fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
        exit(0);
      }
      niik_image_write(fname,tmpimg);
      tmpimg = niik_image_free(tmpimg);
    }

  } /* global scaling testing */
  g_imgs[1] = niik_image_free(g_imgs[1]);

  for(n=m=0; n<gs_num; n++) {
    if(gs_par->m[n][0]<gs_par->m[m][0]) {
      m=n;
    }
  }
  fprintf(stdout,"\t   gs optim: %i / %i at %6.3f with %8.4f\n",m,gs_num,gs_list[m],gs_par->m[m][0]);
  free(gs_list);
  for(n=0; n<17; n++) {
    affpar[n] = gs_par->m[m][n];
  }
  niik_aregister_display_affine(affpar);
  gs_par = niikmat_free(gs_par);

  if(0) {
    niik_aregister_matrix_from_affpar_update(afmat,affpar);
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    sprintf(fname,"tmp_reg3.nii.gz");
    fprintf(stdout,"    transform image for checking: %s\n",fname);
    if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      exit(0);
    }
    niik_image_write(fname,tmpimg);
    tmpimg = niik_image_free(tmpimg);
    exit(0);
  }

  /*
   * 4. rotation and translation
   */
  fprintf(stdout,"    4. rotation and translations\n");
  for(m=0; m<20; m++) {
    daffpar[m]=0;
  }
  daffpar[1]=10;
  daffpar[4]=3;
  daffpar[5]=3;
  daffpar[6]=3;
  if(!niik_image_aregister(mni_img,tmp_mni_seg,g_niik_aregister_movimg,seg,affpar,daffpar,cost_method,sample,filFWHM)) {
    fprintf(stderr,"ERROR: niik_image_aregister\n");
    return 0;
  }

  if(0) {
    niik_aregister_matrix_from_affpar_update(afmat,affpar);
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    sprintf(fname,"tmp_reg4.nii.gz");
    fprintf(stdout,"    transform image for checking: %s\n",fname);
    if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      exit(0);
    }
    niik_image_write(fname,tmpimg);
    tmpimg = niik_image_free(tmpimg);
  }


  /*
   * 5. scaling and shearing
   */
  fprintf(stdout,"    5. scaling and shearing\n");
  for(n=8; n<=13; n++) {
    for(m=0; m<20; m++) {
      daffpar[m]=0;
    }
    daffpar[n]=0.1;
    NIIK_RET0((!niik_image_aregister(mni_img,tmp_mni_seg,g_niik_aregister_movimg,seg,affpar,daffpar,cost_method,sample,filFWHM)),
              fcname,"niik_image_aregister");
    if(0) {
      niik_aregister_matrix_from_affpar_update(afmat,affpar);
      niikmat_display(afmat);
      niik_aregister_display_affine(affpar);
      sprintf(fname,"tmp_reg5.nii.gz");
      fprintf(stdout,"    transform image for checking: %s\n",fname);
      if((tmpimg = niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
        fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
        exit(0);
      }
      niik_image_write(fname,tmpimg);
      tmpimg = niik_image_free(tmpimg);
    }
  } /* scaling / skewing */




  /*
   * 6. each parameter
   */
  fprintf(stdout,"    6. each parameter\n");
  for(n=1; n<=13; n++) {
    for(m=0; m<20; m++) {
      daffpar[m]=0;
    }
    daffpar[n]=0.1;
    NIIK_RET0((!niik_image_aregister(mni_img,tmp_mni_seg,g_niik_aregister_movimg,seg,affpar,daffpar,cost_method,sample,filFWHM)),
              fcname,"niik_image_aregister");
  } /* each parameter */


  g_niik_aregister_movimg = niik_image_free(g_niik_aregister_movimg);
  g_imgs[0] = niik_image_free(g_imgs[0]);
  free(daffpar);
  return 1;
}

int niik_aregister_align_mni_predefined_imgs(nifti_image *img, nifti_image *seg, double *affpar,
    int cost_method, double filFWHM, double sample) {
  nifti_image *mni_img=NULL,*mni_seg=NULL;
  char *FSLDIR=NULL,fname[4096],fcname[512]="niik_aregister_align_mni_predefined_imgs";
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv FSLDIR\n",fcname);
    fprintf(stderr,"  export FSLDIR=/usr/local/fsl\n");
    fprintf(stderr,"  export FSLDIR=/lab1/fsl/4.1.7\n");
    return 0;
  }
  /* read mni images */
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
  fprintf(stdout,"[%s] reading mni img   %s\n",fcname,fname);
  if((mni_img=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
    return 0;
  }
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
  fprintf(stdout,"[%s] reading mni seg   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
    return 0;
  }
  if(!niik_aregister_align_mni(mni_img,mni_seg,img,seg,affpar,cost_method,filFWHM,sample)) {
    fprintf(stderr,"ERROR: niik_aregister_align_mni\n");
    return 0;
  }
  mni_img=niik_image_free(mni_img);
  mni_seg=niik_image_free(mni_seg);
  return 1;
}


int niik_aregister_align_mni_default1(nifti_image *img, nifti_image *seg, double *affpar) {
  nifti_image *mni_img=NULL,*mni_seg=NULL;
  char *FSLDIR,fname[4096];
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[niikmath] ERROR: please setenv FSLDIR\n");
    return 0;
  }
  /* read mni images */
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
  fprintf(stdout,"[niikmath] reading mni img   %s\n",fname);
  if((mni_img=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
  fprintf(stdout,"[niikmath] reading mni seg   %s\n",fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  if(!niik_aregister_align_mni(mni_img,mni_seg,img,seg,affpar,NIIK_REGISTER_CC,2,2)) {
    fprintf(stderr,"ERROR: niik_aregister_align_mni\n");
    return 0;
  }
  mni_img=niik_image_free(mni_img);
  mni_seg=niik_image_free(mni_seg);
  return 1;
}


/*
 * function for calculating 2 rotation based on symmetry
 *
 * -matrix is xyz_to_ijk ( rotate ( translate(center) ) )
 * -then we can compare positive and negative sides along x-axis
 */

double niik_aregister_align_mni_maximize_symmetry(double *v) {
  niikmat *afmat=NULL;
  double
  err=0,
  phi,
  dphi = 0.05,
  rad,
  drad = 0.4,
  dist,
  dd = 1.0;
  niikpt
  pp[3];
  int
  num=0,
  verbose=niik_verbose();
  /* initialization */
  if(v==NULL) {
    g_iter = 0;
    return 0;
  }
  g_iter++;
  if(verbose) fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) start\n");
  if((afmat = niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    fprintf(stderr,"ERROR: niik_aregister_align_mni_maximize_symmetry\n");
    exit(0);
  }
  if(verbose>2) {
    fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) identity\n");
    niikmat_display(afmat);
  }
  niikmat_rotate_matrix_update(afmat,0,-v[0],-v[1]);
  if(verbose>2) {
    fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) rotate\n");
    niikmat_display(afmat);
  }
  niikmat_multiply_mat2_free1(niikmat_translate_matrix(g_niik_aregister_align_mni_maximize_symmetry_center.x,
                              g_niik_aregister_align_mni_maximize_symmetry_center.y,
                              g_niik_aregister_align_mni_maximize_symmetry_center.z),
                              afmat);
  if(verbose>2) {
    fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) translate\n");
    niikmat_display(afmat);
  }
  niikmat_multiply_mat2_free1(niikmat_scale_matrix   (1.0/g_niik_aregister_movimg->dx,
                              1.0/g_niik_aregister_movimg->dy,
                              1.0/g_niik_aregister_movimg->dz),
                              afmat);
  if(verbose>2) {
    fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) scale (pixel size)\n");
    niikmat_display(afmat);
  }
  if(verbose>1) {
    fprintf(stdout,"-d (niik_aregister_align_mni_maximize_symmetry) R = %7.3f %7.3f\n",v[0],v[1]);
    niikmat_display(afmat);
  }
  phi = 0;
  rad = 0;
  for(;;) {
    phi += dphi;
    rad += drad;
    pp[0].x =0;
    pp[0].y = cos(phi) * rad;
    pp[0].z = sin(phi) * rad;
    pp[1] = niikpt_affine_transform(afmat,pp[0]);
    if(verbose>2) {
      fprintf(stdout,"\npp0  %6.1f %6.1f %6.1f\n",pp[0].x,pp[0].y,pp[0].z);
      fprintf(stdout,"pp1  %6.1f %6.1f %6.1f\n",pp[1].x,pp[1].y,pp[1].z);
    }
    if(pp[1].x<0) break;
    if(pp[1].y<0) break;
    if(pp[1].z<0) break;
    if(pp[1].x>g_niik_aregister_movimg->dim[1]) break;
    if(pp[1].y>g_niik_aregister_movimg->dim[2]) break;
    if(pp[1].z>g_niik_aregister_movimg->dim[3]) break;
    for(dist=0; dist<100; dist+=dd) {
      /* Right direction */
      pp[0].x = dist;
      pp[1] = niikpt_affine_transform(afmat,pp[0]);
      if(verbose>3) fprintf(stdout,"  %6.1f %6.1f %6.1f\n",pp[1].x,pp[1].y,pp[1].z);
      if(pp[1].x<0) continue;
      if(pp[1].y<0) continue;
      if(pp[1].z<0) continue;
      if(pp[1].x >= g_niik_aregister_movimg->nx) continue;
      if(pp[1].y >= g_niik_aregister_movimg->ny) continue;
      if(pp[1].z >= g_niik_aregister_movimg->nz) continue;
      /* Left direction */
      pp[0].x = -dist;
      pp[2] = niikpt_affine_transform(afmat,pp[0]);
      if(verbose>2) fprintf(stdout,"  %6.1f %6.1f %6.1f\n",pp[2].x,pp[2].y,pp[2].z);
      if(pp[2].x<0) continue;
      if(pp[2].y<0) continue;
      if(pp[2].z<0) continue;
      if(pp[2].x >= g_niik_aregister_movimg->nx) continue;
      if(pp[2].y >= g_niik_aregister_movimg->ny) continue;
      if(pp[2].z >= g_niik_aregister_movimg->nz) continue;
      /* interpolate */
      pp[1].w = niik_image_interpolate_float_image_3d_linear_ijk(g_niik_aregister_movimg,pp[1]);
      pp[2].w = niik_image_interpolate_float_image_3d_linear_ijk(g_niik_aregister_movimg,pp[2]);
      if(verbose>3) {
        fprintf(stdout,"  %6.1f %6.1f %6.1f | %5.1f\n",pp[1].x,pp[1].y,pp[1].z,pp[1].w);
        fprintf(stdout,"  %6.1f %6.1f %6.1f | %5.1f\n",pp[2].x,pp[2].y,pp[2].z,pp[2].w);
      }
      /* absolute difference is the cost function */
      err+=fabs(pp[1].w - pp[2].w);
      num++;
    } /* symmetric interpolation about mid-plane */
    if(rad>250) break;
  } /* spiral loop on the mid-plane */
  if(verbose>1) fprintf(stdout,"  error = %9.3f   num = %8i   r[1] = %6.1f r[2] = %6.1f \n",
                          err/num,num,
                          v[0],v[1]);
  if(num==0) exit(0);
  niikmat_free(afmat);
  return err / num;
}


#endif /* _FALCON_AREGISTER_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/