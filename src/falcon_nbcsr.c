/* Filename:     nifti1_kunio_aregister.c
 * Description:  affine registration function
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 * Reference:    Chen JT et al 2008 ISMRM; Nakamura et al 2011 NeuroImage;
 *
 * R1. rough initial standard space registration
 * R2. 6-dof (R+T) head registration
 * R3. 12-dof full affine brain registration
 * R4. 12-dof full affine non-brain registration
 * R5. 6-dof brain registration with fixed scaling/skewing
 *
 *
 *
 * int niik_aregister_nbcr_display_var(double *daffpar);
 * int niik_aregister_nbcr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg,
 *     double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag);
 * int niik_aregister_nbcsr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg,
 *     double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag);
 *
 */

#ifndef _FALCON_NBCSR_C_
#define _FALCON_NBCSR_C_

#include "falcon.h"

int g_FALCON_nbcsr_debug=0;

void nifti1_kunio_nbcsr_turn_on_debug()  {
  g_FALCON_nbcsr_debug=1;
}
void nifti1_kunio_nbcsr_turn_off_debug() {
  g_FALCON_nbcsr_debug=0;
}

int niik_aregister_nbcr_display_var(double *daffpar);

int niik_aregister_nbcr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg, double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag)
/* niik_aregister_nbcr
 * -function to do image registration of two images
 * -movseg can be null, but more robust with movseg
 * -returns 1 for success and zero for failure
 * -arguments include:
 *   refimg    reference (target) gray scale image
 *   refseg    reference (target) mask image --usually brain or closed brain
 *   movimg    moving grayscale image
 *   movseg    moving mask image
 *   affpar    affine parameters (replaced on output), see nifti1_kunio_aregister.c for description
 *   cost_method    cost function type, see nifti1_kunio.h (NIIK_REGISTER_NMI NIIK_REGISTER_CC)
 *   filFWHM   basic filter size (FWHM)
 *   sample    basic sampling distance
 */
{

  char fcname[32]="niik_aregister_nbcr";
  nifti_image
  *tmpimg=NULL,*tmpimglist[4],
   *blurimgs[2],
   *mni_roi=NULL,
    *mni_img=NULL,
     *mni_seg=NULL;
  niikmat *rpar=NULL;
  double
  *daffpar=NULL;
  int
  verbose=1,
  R4_nlevel=1,
  *R4_maxiter=NULL,
   *R4_nseed=NULL,
    m,n;
  niikvec
  *R4_FWHM=NULL,
   *R4_delta=NULL,
    *R4_tol=NULL;
  niikmat
  *R4_daffpar=NULL,
   *R4_seed=NULL,
    *afmat=NULL;
  const char *NIIKDIR=NULL;
  char
  *FSLDIR=NULL,
   fname[4096];
  unsigned char
  *bimg,*bseg;
  nmi_obj *nmiobj=NULL;
  niikvec *idaffpar;

  if(verbose>=1) niik_fc_display(fcname,1);
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(refseg==NULL) {
    fprintf(stderr,"[%s] ERROR: refseg is null\n",fcname);
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"[%s] ERROR: movimg is null\n",fcname);
    return 0;
  }
  if(affpar==NULL) {
    fprintf(stderr,"[%s] ERROR: affpar is null\n",fcname);
    return 0;
  }
  if(refseg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"[%s] ERROR: refseg->datatype != NIFTI_TYPE_UINT8\n",fcname);
    return 0;
  }
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv FSLDIR\n",fcname);
    return 0;
  }
  if((NIIKDIR = get_NIIKDIR())==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv NIIKDIR\n",fcname);
    return 0;
  }
  sprintf(fname,"%s/data/NBCSR/MNI152_T1_2mm_NBCSR_ROI.nii.gz",NIIKDIR);
  if((mni_roi=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
    return 0;
  }

  if((daffpar=(double *)calloc(21,sizeof(double)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc\n",fcname);
    return 0;
  }
  idaffpar=niikvec_init(21);
  rpar=niikmat_init(6,21);
  rpar->m[0][14]=rpar->m[1][14]=affpar[14];
  rpar->m[0][15]=rpar->m[1][15]=affpar[15];
  rpar->m[0][16]=rpar->m[1][16]=affpar[16];

  niik_image_aregister_set_g_verbose_errfunc(0);

  if(cost_method==NIIK_REGISTER_NMI) {
    /* nmi object */
    fprintf(stdout,"[%s] crate NMI object\n",fcname);
    if((nmiobj = niik_aregister2_nmi_obj_init(niik_image_get_min(refimg,NULL),niik_image_get_max(refimg,NULL),
                 niik_image_get_min(movimg,NULL),niik_image_get_max(movimg,NULL),32,32))==NULL) {
      fprintf(stderr,"ERROR: niik_aregister2_nmi_obj_init\n");
      return 0;
    }
    niik_image_aregister_set_nmi_obj(nmiobj);
    niik_aregister_g_nmi_obj_set(nmiobj);
  }


  /*
   *
   * R1 standard space registration
   *
   */
  fprintf(stdout,"\n\n");
  fprintf(stdout,"[%s] R1 standard space registration\n",fcname);
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
  fprintf(stdout,"[%s] reading mni img   %s\n",fcname,fname);
  if((mni_img=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
    return 0;
  }
  if(refimg->sform_code == NIFTI_XFORM_MNI_152) {
    fprintf(stdout,"[niik_aregister_nbcr] using sform_code %s:\n",nifti_xform_string(refimg->sform_code));
    /*niikmat_display(niikmat_mat44_matrix(refimg->sto_xyz));*/
    afmat=niikmat_mat44_matrix(refimg->sto_xyz);
    niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(1.0/refimg->dx,1.0/refimg->dy,1.0/refimg->dz));
    niikmat_multiply_mat2_free1(niikmat_mat44_matrix(mni_img->sto_ijk),afmat);
    niikmat_display(afmat);
    /*if(!niikmat_decompose_affine(afmat,rpar->m[0],12)){
      fprintf(stderr,"ERROR: niikmat_decompose_affine\n");
      return 0; }*/
  } /* using previously calculated sform */
  else {
    sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
    fprintf(stdout,"[niik_aregister_nbcr] reading mni seg   %s\n",fname);
    if((mni_seg=nifti_image_read(fname,1))==NULL) {
      fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
      return 0;
    }
    rpar->m[0][7]=rpar->m[0][8]=rpar->m[0][9]=rpar->m[0][10]=1;
    /* nmi object */
    if(!niik_aregister_nmi_obj_update_var(niik_image_get_min(refimg,NULL),niik_image_get_max(refimg,NULL),
                                          niik_image_get_min(movimg,NULL),niik_image_get_max(movimg,NULL),32,32)) {
      fprintf(stderr,"[%s] ERROR: niik_aregister_nmi_obj_update_var\n",fcname);
      return 0;
    }
    if(!niik_aregister_align_mni(mni_img,mni_seg,refimg,refseg,rpar->m[0],NIIK_REGISTER_NMI,6,3)) {
      fprintf(stderr,"[%s] ERROR: niik_aregister_align_mni\n",fcname);
      return 0;
    }
    mni_seg=niik_image_free(mni_seg);
    afmat=niik_aregister_matrix_from_affpar(rpar->m[0]);
    fprintf(stdout,"[niik_aregister_nbcr] R1 affine parameters\n");
    niik_aregister_display_affine( rpar->m[0] );
  } /* calculate sform matrix */
  mni_img=niik_image_free(mni_img);

  fprintf(stdout,"[niik_aregister_nbcr] transform crop ROI\n");
  niikmat_inverse_update(afmat);
  if(!niik_image_affine_transform_3d_update(mni_roi,refimg,afmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
    return 0;
  }
  if(!niik_image_type_convert(mni_roi,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }

  if(verbose/10 || g_FALCON_nbcsr_debug) {
    fprintf(stdout,"[niik_aregister_nbcr] writing output  tmpnbcsr_mni_roi.nii.gz\n");
    niik_image_write("tmpnbcsr_mni_roi.nii.gz",mni_roi);
  }


  /*
   *
   * R2 whole head 6-dof registration
   *
   */
  fprintf(stdout,"\n[niik_aregister_nbcr] R2 whole head 6-dof registration\n");
  for(n=0; n<21; n++) {
    rpar->m[1][n]=affpar[n];
    daffpar[n]=0;
  }
  daffpar[1] =
    daffpar[2] =
      daffpar[3] = 30;
  daffpar[4] =
    daffpar[5] =
      daffpar[6] = 20;
  niik_image_aregister_set_g_img_to_null();
  fprintf(stdout,"                      delta    %8.4f\n",sample*2.2);
  fprintf(stdout,"                      FWHM     %8.4f\n",filFWHM*3);
  fprintf(stdout,"                      cost     %s\n",niik_aregister_method_string(cost_method));
  for(n=0; n<17; n++) {
    idaffpar->v[n]=daffpar[n];
  }
  niik_aregister_nbcr_display_var(daffpar);

  if(verbose) fprintf(stdout,"[niik_aregister_nbcr] gaussian filter\n");
  if((blurimgs[0] = niik_image_filter_gaussian(refimg,filFWHM*2,filFWHM))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
    return 0;
  }
  if((blurimgs[1] = niik_image_filter_gaussian(movimg,filFWHM*2,filFWHM))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
    return 0;
  }

  if(1) {
    if(verbose) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister2_test2\n");
    if(!niik_image_aregister2_test2(blurimgs[0],NULL,blurimgs[1],NULL,rpar->m[1],cost_method,3,filFWHM,100,2,sample*2.2,idaffpar,nmiobj)) {
      fprintf(stderr,"[%s] ERROR: niik_image_aregister2_test2\n",fcname);
      return 0;
    }
  } else {
    if(verbose) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister\n");
    if(!niik_image_aregister(blurimgs[0],NULL,blurimgs[1],NULL,rpar->m[1],daffpar,cost_method,sample*3.2,filFWHM)) {
      fprintf(stderr,"[%s] ERROR: niik_image_aregister\n",fcname);
      return 0;
    }
  }

  fprintf(stdout,"[niik_aregister_nbcr] R2 affine parameters\n");
  niik_aregister_display_affine( rpar->m[1] );

  if(verbose/10 || g_FALCON_nbcsr_debug) {
    niik_aregister_matrix_from_affpar_update(afmat,rpar->m[1]);
    if(niik_image_aregister_get_g_invmat()) niikmat_inverse_update(afmat);
    tmpimg=niik_image_affine_transform_3d(movimg,refimg,afmat,NIIK_INTERP_LINEAR);
    tmpimglist[0]=tmpimg;
    tmpimglist[1]=refimg;
    fprintf(stdout,"\t  writing output  tmpnbcsr_mni_R2.nii.gz\n");
    niik_image_combine_and_write("tmpnbcsr_mni_R2.nii.gz",tmpimglist,2,'t',140);
    /*niik_image_write("tmpnbcsr_mni_R2.nii.gz",tmpimg); */
    tmpimg=niik_image_free(tmpimg);
  }


  /*
   *
   * R3 brain full affine 12-dof registration
   *
   */
  fprintf(stdout,"\n[niik_aregister_nbcr] R3 brain full affine 12-dof registration\n");
  if(cost_method==NIIK_REGISTER_NMI) {
    niik_aregister_g_nmi_obj_set(nmiobj);
  }

  for(n=0; n<21; n++) {
    rpar->m[2][n]=rpar->m[1][n];
    daffpar[n]=0;
  }
  daffpar[1] =
    daffpar[2] =
      daffpar[3] = 10;
  daffpar[4] =
    daffpar[5] =
      daffpar[6] = 10;
  if(!rigid_flag) {
    daffpar[7] =
      daffpar[8] =
        daffpar[9] =
          daffpar[10] = 0.05;
    daffpar[11] =
      daffpar[12] =
        daffpar[13] = 0.05;
  }
  niik_image_aregister_set_g_img_to_null();
  fprintf(stdout,"                      delta    %8.4f\n",sample*1.2);
  fprintf(stdout,"                      FWHM     %8.4f\n",filFWHM*2.0);
  fprintf(stdout,"                      cost     %s\n",niik_aregister_method_string(cost_method));
  niik_aregister_nbcr_display_var(daffpar);

  if(verbose) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister\n");
  if(!niik_image_aregister(blurimgs[0],refseg,blurimgs[1],NULL,rpar->m[2],daffpar,cost_method,sample*2.2,filFWHM)) {
    fprintf(stderr,"[%s] ERROR: niik_image_aregister\n",fcname);
    return 0;
  }

  fprintf(stdout,"[niik_aregister_nbcr] R3 affine parameters\n");
  niik_aregister_display_affine( rpar->m[2] );

  if(verbose/10 || g_FALCON_nbcsr_debug) {
    niik_aregister_matrix_from_affpar_update(afmat,rpar->m[2]);
    if(niik_image_aregister_get_g_invmat()) niikmat_inverse_update(afmat);
    tmpimg=tmpimglist[0]=niik_image_affine_transform_3d(movimg,refimg,afmat,NIIK_INTERP_LINEAR);
    fprintf(stdout,"\t  writing output  tmpnbcsr_mni_R3.nii.gz\n");
    niik_image_combine_and_write("tmpnbcsr_mni_R3.nii.gz",tmpimglist,2,'t',140);
    /*niik_image_write("tmpnbcsr_mni_R3.nii.gz",tmpimg); */
    tmpimg=niik_image_free(tmpimg);
  }


  /*
   *
   * R4 non-brain full affine 12-dof registration
   *
   */
  fprintf(stdout,"\n[niik_aregister_nbcr] R4 non-brain full affine 12-dof registration\n");

  if(rigid_flag) {
    for(n=0; n<21; n++) {
      rpar->m[3][n]=rpar->m[2][n];
    }
  }

  else {
    /* prepare crop ROI */
    if(verbose) fprintf(stdout,"[niik_aregister_nbcr] prepare non-brain mask \n");
    fprintf(stdout,"    prepare cropped ROI\n");
    bimg = mni_roi -> data;
    bseg = refseg -> data;
    for(n=0; n<refseg->nvox; n++) {
      if(bseg[n]==0 && bimg[n]>0)
        bimg[n] = 1;
      else
        bimg[n] = 0;
    }
    if(verbose/10 || g_FALCON_nbcsr_debug) {
      fprintf(stdout,"[niik_aregister_nbcr] writing output  tmpnbcsr_mni_roi.nii.gz\n");
      niik_image_write("tmpnbcsr_mni_roi2.nii.gz",mni_roi);
    }

    for(n=0; n<21; n++) {
      rpar->m[3][n]=rpar->m[2][n];
      daffpar[n]=0;
    }
    daffpar[1] = daffpar[2] = daffpar[3] = 5;
    daffpar[4] = daffpar[5] = daffpar[6] = 5;
    if(!rigid_flag) {
      daffpar[ 7] = 0;
      daffpar[ 8] = daffpar[ 9] = daffpar[10] = 0.05;
      daffpar[11] = daffpar[12] = daffpar[13] = 0.05;
    }
    niik_image_aregister_set_g_img_to_null();
    fprintf(stdout,"                      delta    %8.4f\n",sample/1.2);
    fprintf(stdout,"                      FWHM     %8.4f\n",filFWHM/2.2);
    fprintf(stdout,"                      cost     %s\n",niik_aregister_method_string(cost_method));
    niik_aregister_nbcr_display_var(daffpar);

    if(0) {
      R4_nlevel = 9;
      R4_FWHM = niikvec_init(R4_nlevel);
      R4_delta = niikvec_init(R4_nlevel);
      R4_tol = niikvec_init(R4_nlevel);
      R4_daffpar = niikmat_init(R4_nlevel,25);
      R4_maxiter = (int *)calloc(R4_nlevel,sizeof(int));
      R4_nseed   = (int *)calloc(R4_nlevel,sizeof(int));
      R4_nseed[0] = 10;
      R4_seed    = niikmat_init(R4_nseed[0],25);
      for(n=0; n<R4_nseed[0]; n++) {
        for(m=0; m<17; m++) {
          R4_seed->m[n][m] = rpar->m[3][m];
        }
        if(n==0) continue;
        R4_seed->m[n][1] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][2] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][3] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][4] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][5] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][6] += 5 * (niik_get_rand()-0.5);
        R4_seed->m[n][7] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][8] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][9] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][10] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][11] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][12] += 0.02 * (niik_get_rand()-0.5);
        R4_seed->m[n][13] += 0.02 * (niik_get_rand()-0.5);
      }
      for(n=0; n<R4_nlevel; n++) {
        R4_daffpar->m[n][1] =
          R4_daffpar->m[n][2] =
            R4_daffpar->m[n][3] = 5;
        R4_daffpar->m[n][4] =
          R4_daffpar->m[n][5] =
            R4_daffpar->m[n][6] = 5;
      }
      R4_daffpar->m[0][ 7]=0.05;
      R4_daffpar->m[1][ 8]=0.05;
      R4_daffpar->m[2][ 9]=0.05;
      R4_daffpar->m[3][10]=0.05;
      R4_daffpar->m[4][11]=0.05;
      R4_daffpar->m[5][12]=0.05;
      R4_daffpar->m[6][13]=0.05;
      R4_daffpar->m[7][8]=R4_daffpar->m[7][9]=R4_daffpar->m[7][10]=
          R4_daffpar->m[7][11]=R4_daffpar->m[7][12]=R4_daffpar->m[7][13]=0.02;
      R4_daffpar->m[8][8]=R4_daffpar->m[8][9]=R4_daffpar->m[8][10]=
          R4_daffpar->m[8][11]=R4_daffpar->m[8][12]=R4_daffpar->m[8][13]=0.02;
      for(n=1; n<R4_nlevel-1; n++)  {
        R4_nseed[n] = 5;
      }
      for(n=0; n<R4_nlevel-1; n++)  {
        R4_maxiter[n] = 36;
        R4_delta->v[n]=sample*2.0;
        R4_FWHM->v[n]=filFWHM*2.0;
      }
      R4_nseed[n] = 2;
      R4_maxiter[n] = 50;
      R4_delta->v[n]=sample;
      R4_FWHM->v[n]=filFWHM/1.2;
      if(!niik_image_aregister2_multilevel(blurimgs[0],mni_roi,blurimgs[1],NULL,rpar->m[3],cost_method,nmiobj,
                                           R4_nlevel,R4_FWHM,R4_delta,R4_tol,R4_nseed,R4_maxiter,R4_daffpar,R4_seed)) {
        fprintf(stderr,"[%s] ERROR: niik_image_aregister2_multilevel\n",fcname);
        return 0;
      }
      free(R4_nseed);
      free(R4_maxiter);
      R4_FWHM=niikvec_free(R4_FWHM);
      R4_delta=niikvec_free(R4_delta);
      R4_daffpar=niikmat_free(R4_daffpar);
      R4_seed=niikmat_free(R4_seed);
    } else if(0) {
      for(n=0; n<17; n++) {
        idaffpar->v[n]=daffpar[n];
      }
      if(verbose>=1) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister2_test2\n");
      if(!niik_image_aregister2_test2(blurimgs[0],mni_roi,blurimgs[1],NULL,rpar->m[3],cost_method,3,filFWHM/1.2,150,2,sample,idaffpar,nmiobj)) {
        fprintf(stderr,"[%s] ERROR: niik_image_aregister2_test2\n",fcname);
        return 0;
      }
    } else {
      if(verbose) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister\n");
      if(!niik_image_aregister(blurimgs[0],mni_roi,blurimgs[1],NULL,rpar->m[3],daffpar,cost_method,sample,filFWHM)) {
        fprintf(stderr,"[%s] ERROR: niik_image_aregister\n",fcname);
        return 0;
      }
    }

    mni_roi=niik_image_free(mni_roi);
    fprintf(stdout,"[niik_aregister_nbcr] R4 affine parameters\n");
    niik_aregister_display_affine( rpar->m[3] );

    if(verbose/10 || g_FALCON_nbcsr_debug) {
      niik_aregister_matrix_from_affpar_update(afmat,rpar->m[3]);
      if(niik_image_aregister_get_g_invmat()) niikmat_inverse_update(afmat);
      tmpimg=tmpimglist[0]=niik_image_affine_transform_3d(movimg,refimg,afmat,NIIK_INTERP_LINEAR);
      fprintf(stdout,"\t  writing output  tmpnbcsr_mni_R4.nii.gz\n");
      niik_image_combine_and_write("tmpnbcsr_mni_R4.nii.gz",tmpimglist,2,'t',140);
      /*niik_image_write("tmpnbcsr_mni_R4.nii.gz",tmpimg); */
      tmpimg=niik_image_free(tmpimg);
    }
  }

  /*
   *
   * R5 brain 6-dof registration with fixed scaling and skewing
   *
   */
  fprintf(stdout,"\n[niik_aregister_nbcr] R5 brain full affine 6-dof registration\n");
  for(n=0; n<21; n++) {
    rpar->m[4][n]=rpar->m[3][n];
    daffpar[n]=0;
  }
  daffpar[1] =
    daffpar[2] =
      daffpar[3] = 2;
  daffpar[4] =
    daffpar[5] =
      daffpar[6] = 2;
  niik_image_aregister_set_g_img_to_null();
  fprintf(stdout,"                      delta    %8.4f\n",sample);
  fprintf(stdout,"                      FWHM     %8.4f\n",filFWHM);
  fprintf(stdout,"                      cost     %s\n",niik_aregister_method_string(cost_method));
  niik_aregister_nbcr_display_var(daffpar);
  if(verbose>=1) fprintf(stdout,"[niik_aregister_nbcr] niik_image_aregister\n");
  if(!niik_image_aregister(blurimgs[0],refseg,blurimgs[1],NULL,rpar->m[4],daffpar,cost_method,sample,0)) {
    fprintf(stderr,"[%s] ERROR: niik_image_aregister\n",fcname);
    return 0;
  }
  fprintf(stdout,"[niik_aregister_nbcr] R5 affine parameters\n");
  niik_aregister_display_affine( rpar->m[4] );

  if(verbose/10 || g_FALCON_nbcsr_debug) {
    niik_aregister_matrix_from_affpar_update(afmat,rpar->m[4]);
    if(niik_image_aregister_get_g_invmat()) niikmat_inverse_update(afmat);
    tmpimg=tmpimglist[0]=niik_image_affine_transform_3d(movimg,refimg,afmat,NIIK_INTERP_LINEAR);
    fprintf(stdout,"\t  writing output  tmpnbcsr_mni_R5.nii.gz\n");
    /*niik_image_write("tmpnbcsr_mni_R5.nii.gz",tmpimg); */
    niik_image_combine_and_write("tmpnbcsr_mni_R5.nii.gz",tmpimglist,2,'t',140);
    tmpimg=niik_image_free(tmpimg);
  }

  for(n=0; n<17; n++) {
    affpar[n] = rpar->m[4][n];
  }

  niik_aregister_g_nmi_obj_set(NULL);
  /* free memory */
  rpar=niikmat_free(rpar);
  blurimgs[0]=niik_image_free(blurimgs[0]);
  blurimgs[1]=niik_image_free(blurimgs[1]);
  return 1;
} /* niik_aregister_nbcsr */



int niik_aregister_nbcr_display_var(double *daffpar) {
  int n,nvar;
  double lim=1e-5;
  for(n=1,nvar=0; n<=16; n++)
    if(fabs(daffpar[n])>lim) nvar++;
  fprintf(stdout,"                      var      %i ",nvar);
  if(fabs(daffpar[1])>lim) fprintf(stdout,"Rx ");
  if(fabs(daffpar[2])>lim) fprintf(stdout,"Ry ");
  if(fabs(daffpar[3])>lim) fprintf(stdout,"Rz ");
  if(fabs(daffpar[4])>lim) fprintf(stdout,"Tx ");
  if(fabs(daffpar[5])>lim) fprintf(stdout,"Ty ");
  if(fabs(daffpar[6])>lim) fprintf(stdout,"Tz ");
  if(fabs(daffpar[7])>lim) fprintf(stdout,"Sg ");
  if(fabs(daffpar[8 ])>lim) fprintf(stdout,"Sx ");
  if(fabs(daffpar[9 ])>lim) fprintf(stdout,"Sy ");
  if(fabs(daffpar[10])>lim) fprintf(stdout,"Sz ");
  if(fabs(daffpar[11])>lim) fprintf(stdout,"Kx ");
  if(fabs(daffpar[12])>lim) fprintf(stdout,"Ky ");
  if(fabs(daffpar[13])>lim) fprintf(stdout,"Kz ");
  fprintf(stdout,"\n");
  return 1;
}


int niik_aregister_nbcsr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg, double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag)
/* niik_aregister_nbcsr
 * -function to do symmetric image registration of two images
 * -rest is the same as niik_aregister_nbcr
 */

{
  char fcname[64]="niik_aregister_nbcsr";
  nifti_image
  *tmpimg,*tmplist[2],
  *tmpseg=NULL;
  double
  *fwdpar,
  *invpar;
  niikmat
  *fwdmat=NULL,
   *invmat=NULL,
    *afmat=NULL;
  niikpt
  pctr[3];
  int
  dof=12,
  debug=0,
  n;

  fwdpar=(double *)calloc(21,sizeof(double));
  invpar=(double *)calloc(21,sizeof(double));

  for(n=0; n<3; n++) pctr[n]=niikpt_zero();

  /*
   * initial parameters
   */
  if(rigid_flag) dof=6;
  if(affpar[0]>0) {
    for(n=0; n<17; n++) {
      fwdpar[n]=invpar[n]=affpar[n];
    }
  } /* use current affpar as initial values */
  else {
    pctr[1] = niikpt_image_get_centroid(refimg,NULL);
    pctr[2] = niikpt_image_get_centroid(movimg,NULL);
    pctr[0] = niikpt_avg(pctr[1],pctr[2]);
    fwdpar[14] = pctr[0].x;
    fwdpar[15] = pctr[0].y;
    fwdpar[16] = pctr[0].z;
    fprintf(stdout,"[niik_aregister_nbcsr] centroid [ref] %9.4f %9.4f %9.4f\n",pctr[1].x,pctr[1].y,pctr[1].z);
    fprintf(stdout,"                       centroid [mov] %9.4f %9.4f %9.4f\n",pctr[2].x,pctr[2].y,pctr[2].z);
    fprintf(stdout,"                       centroid       %9.4f %9.4f %9.4f\n",pctr[0].x,pctr[0].y,pctr[0].z);
    for(n=7; n<=10; n++)
      fwdpar[n]=invpar[n]=1;
    fwdpar[4]=invpar[4]=pctr[1].x-pctr[2].x;
    fwdpar[5]=invpar[5]=pctr[1].y-pctr[2].y;
    fwdpar[6]=invpar[6]=pctr[1].z-pctr[2].z;
    fprintf(stdout,"[niik_aregister_nbcsr] preset translation   %9.4f %9.4f %9.4f\n",fwdpar[4],fwdpar[5],fwdpar[6]);
  } /* calculate initial values */


  /*
   * forward registration
   */
  if(rigid_flag) fprintf(stdout,"[niik_aregister_nbcsr] forward direction (rigid)\n");
  else           fprintf(stdout,"[niik_aregister_nbcsr] forward direction\n");
  if(!niik_aregister_nbcr(refimg,refseg,movimg,movseg,fwdpar,cost_method,filFWHM,sample,rigid_flag)) {
    fprintf(stderr,"[%s] ERROR: niik_aregister_nbcr\n",fcname);
    return 0;
  }
  fwdmat=niik_aregister_matrix_from_affpar(fwdpar);

  if(debug) {
    afmat=niik_aregister_matrix_from_affpar(fwdpar);
    fprintf(stdout,"[niik_aregister_nbcsr] transform image (debug)\n");
    tmplist[1]=niik_image_copy(movimg);
    niik_image_affine_transform_3d_update(tmplist[1],refimg,afmat,NIIK_INTERP_LINEAR);
    tmplist[0]=refimg;
    tmpimg=niik_image_combine(tmplist,2,4,140);
    niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8);
    niik_image_write("tmp_nbcsr_fwd.nii.gz",tmpimg);
    tmpimg=niik_image_free(tmpimg);
    tmplist[1]=niik_image_free(tmpimg);
    afmat = niikmat_free(afmat);
  }


  /*
   * -if moving image's mask is missing,
   *  create one
   */
  niik_image_aregister_set_g_invmat(1);
  if(movseg==NULL) {
    fprintf(stdout,"[niik_aregister_nbcsr] creating moving mask image\n");
    afmat=niik_aregister_matrix_from_affpar(fwdpar);
    niikmat_inverse_update(afmat);
    if((tmpseg=niik_image_affine_transform_3d(refseg,movimg,afmat,NIIK_INTERP_NN))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d\n",fcname);
      return 0;
    }
    afmat=niikmat_free(afmat);
  } else {
    fprintf(stdout,"[niik_aregister_nbcsr] using moving mask image\n");
    tmpseg=movseg;
  }

  /*
   * apply the inverse registration
   */
  invpar[14] = pctr[0].x;
  invpar[15] = pctr[0].y;
  invpar[16] = pctr[0].z;
  fprintf(stdout,"[niik_aregister_nbcsr] inverse direction\n");
  if(!niik_aregister_nbcr(movimg,tmpseg,refimg,refseg,invpar,cost_method,filFWHM,sample,rigid_flag)) {
    fprintf(stderr,"[%s] ERROR: niik_aregister_nbcr\n",fcname);
    return 0;
  }
  if(movseg==NULL)
    tmpseg=niik_image_free(tmpseg);
  invmat=niik_aregister_matrix_from_affpar(invpar);

  if(debug) {
    fprintf(stdout,"[niik_aregister_nbcsr] transform image (debug)\n");
    afmat=niik_aregister_matrix_from_affpar(invpar);
    tmplist[1]=niik_image_copy(movimg);
    niik_image_affine_transform_3d_update(tmplist[1],refimg,afmat,NIIK_INTERP_LINEAR);
    tmplist[0]=refimg;
    tmpimg=niik_image_combine(tmplist,2,4,140);
    niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8);
    niik_image_write("tmp_nbcsr_inv.nii.gz",tmpimg);
    tmpimg=niik_image_free(tmpimg);
    tmplist[1]=niik_image_free(tmpimg);
    afmat = niikmat_free(afmat);
  }

  if(0) {
    fwdpar[14]=invpar[14]=pctr[0].x;
    fwdpar[15]=invpar[15]=pctr[0].y;
    fwdpar[16]=invpar[16]=pctr[0].z;
    fprintf(stdout,"[niik_aregister_nbcsr] decompose fwd parameters\n");
    niikmat_decompose_affine(fwdmat,fwdpar,dof);
    fprintf(stdout,"[niik_aregister_nbcsr] decompose inv parameters\n");
    niikmat_decompose_affine(invmat,invpar,dof);
  }


  /*
   * calculate symmetric registration parameters
   */
  fprintf(stdout,"[niik_aregister_nbcsr] forward parameters\n");
  niik_aregister_display_affine(fwdpar);
  fprintf(stdout,"[niik_aregister_nbcsr] inverse parameters\n");
  niik_aregister_display_affine(invpar);
  for(n=1; n<17; n++) {
    affpar[n] = (fwdpar[n] + invpar[n]) * 0.5;
  }
  fprintf(stdout,"[niik_aregister_nbcsr] averaged parameters\n");
  niik_aregister_display_affine(affpar);

  /* reset/free things */
  niik_image_aregister_set_g_invmat(0);
  free(fwdpar);
  free(invpar);
  return 1;
}

#endif /* _FALCON_NBCSR_C_ */
