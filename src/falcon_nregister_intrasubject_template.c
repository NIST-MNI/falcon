/* Filename:     nifti1_kunio_nregister_intrasubject_template.c
 * Description:  nonlinear registration functions for intra-subject template
 * Author:       Kunio Nakamura
 * Date:         November 9, 2012
 */

#ifndef _FALCON_NREGISTER_INTRASUBJECT_TEMPLATE_C_
#define _FALCON_NREGISTER_INTRASUBJECT_TEMPLATE_C_

#include "falcon.h"

int niik_image_nregister_average_with_fov(nifti_image *refimg,nifti_image **warpimg,nifti_image **imglist,niikmat **matlist,int num,int interp);
int niik_image_nregister_intrasubject_template(nifti_image *refimg,nifti_image *maskimg,
    nifti_image **warpimglist,nifti_image **imglist,niikmat **matlist,int num,
    niikvec *gFWHM,niikvec *wFWHM,
    niikvec *maxiter,
    int nlevel);


int niik_image_nregister_intrasubject_template_test1(nifti_image *refimg,nifti_image *maskimg,
    nifti_image **warpimglist,nifti_image **imglist,niikmat **matlist,int num)
/* using pre-defined values */
{

  char fcname[128]="niik_image_nregister_intrasubject_template_test1";
  niikvec *gFWHM,*wFWHM,*maxiter;
  int nlevel=2;
  gFWHM=niikvec_init(nlevel);
  wFWHM=niikvec_init(nlevel);
  maxiter=niikvec_init(nlevel);
  gFWHM  ->v[0]=1.5;
  gFWHM->v[1]=1.0;
  wFWHM  ->v[0]=2.5;
  wFWHM->v[1]=1.5;
  maxiter->v[0]=5;
  maxiter->v[1]=5;
  if(!niik_image_nregister_intrasubject_template(refimg,maskimg,warpimglist,imglist,matlist,num,gFWHM,wFWHM,maxiter,nlevel)) {
    fprintf(stderr,"[%s] ERROR: niik_image_nregister_intrasubject_template\n",fcname);
    return 0;
  }
  gFWHM=niikvec_free(gFWHM);
  wFWHM=niikvec_free(wFWHM);
  maxiter=niikvec_free(maxiter);
  return 1;
} /* niik_image_nregister_intrasubject_template_test1 */


int niik_image_nregister_intrasubject_template(nifti_image *refimg,nifti_image *maskimg,
    nifti_image **warpimglist,nifti_image **imglist,niikmat **matlist,int num,
    niikvec *gFWHM,niikvec *wFWHM,
    niikvec *maxiter,
    int nlevel) {

  char
  fname[4096],
        fcname[64]="niik_image_nregister_intrasubject_template";
  int
  n,level,iter,
  verbose=1;
  niikmat *imat=NULL;
  nifti_image
  *gradimg=NULL;

  niik_fc_display(fcname,1);
  /* CHECK INPUTS */
  if(verbose>=2) fprintf(stdout,"[%s] checking inputs\n",fcname);
  if(gFWHM==NULL) {
    fprintf(stderr,"[%s] ERROR: gFWHM is null\n",fcname);
    return 0;
  }
  if(wFWHM==NULL) {
    fprintf(stderr,"[%s] ERROR: wFWHM is null\n",fcname);
    return 0;
  }
  if(maxiter==NULL) {
    fprintf(stderr,"[%s] ERROR: maxiter is null\n",fcname);
    return 0;
  }
  if(warpimglist==NULL) {
    fprintf(stderr,"[%s] ERROR: warpimg is null\n",fcname);
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(verbose>=2) fprintf(stdout,"[%s]   checking warpimg\n",fcname);
  for(n=0; n<num; n++) {
    if(warpimglist[n]==NULL) {
      fprintf(stderr,"[%s] ERROR: warpimglist[%i] is null\n",fcname,n);
      return 0;
    }
  }
  if(verbose>=2) fprintf(stdout,"[%s]   checking imglist\n",fcname);
  for(n=0; n<num; n++) {
    if(imglist[n]==NULL) {
      fprintf(stderr,"[%s] ERROR: imglist[%i] is null\n",fcname,n);
      return 0;
    }
  }
  if(verbose>=2) fprintf(stdout,"[%s]   checking matlist\n",fcname);
  for(n=0; n<num; n++) {
    if(matlist[n]==NULL) {
      fprintf(stderr,"[%s] ERROR: imglist[%i] is null\n",fcname,n);
      return 0;
    }
  }

  if(verbose>=2) fprintf(stdout,"[%s]   checking warp dimensions\n",fcname);
  for(n=0; n<num; n++) {
    if(warpimglist[n]->nx!=refimg->nx) {
      fprintf(stderr,"[%s] ERROR: x-dimension for warpimg[%i] and refimg, %i %i\n",fcname,n,warpimglist[n]->nx,refimg->nx);
      return 0;
    }
    if(warpimglist[n]->ny!=refimg->ny) {
      fprintf(stderr,"[%s] ERROR: x-dimension for warpimg[%i] and refimg, %i %i\n",fcname,n,warpimglist[n]->ny,refimg->ny);
      return 0;
    }
    if(warpimglist[n]->nz!=refimg->nz) {
      fprintf(stderr,"[%s] ERROR: x-dimension for warpimg[%i] and refimg, %i %i\n",fcname,n,warpimglist[n]->nz,refimg->nz);
      return 0;
    }
    if(warpimglist[n]->nu!=3) {
      fprintf(stderr,"[%s] ERROR: warp image[%i] nu is not 3, %i\n",fcname,n,warpimglist[n]->nu);
      return 0;
    }
    if(warpimglist[n]->datatype!=NIFTI_TYPE_FLOAT32) {
      fprintf(stderr,"[%s] ERROR: warp image[%i] datatype is not float32, %s\n",fcname,n,nifti_datatype_string(warpimglist[n]->datatype));
      return 0;
    }
  }
  if(niik_image_cmp_dim(refimg,maskimg)!=0) {
    fprintf(stderr,"[%s] ERROR: image dimension for refimg and maskimg\n",fcname);
    return 0;
  }

  imat=niikmat_identity(4,4);

  for(level=0; level<nlevel; level++) {
    if(verbose>=1) fprintf(stdout,"[%s] level %i gFWHM %6.2f \n",fcname,level,gFWHM->v[level]);
    if((gradimg=niik_image_gauss_sobel_filters_with_mag_single_output(refimg,gFWHM->v[level]))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag_single_output\n",fcname);
      return 0;
    }

    for(iter=0; iter<(int)floor(0.5+maxiter->v[level]); iter++) {
      for(n=0; n<num; n++) {
        if(!niik_image_nregister_demons_update(warpimglist[n],gradimg,refimg,maskimg,refimg,imat,imglist[n],matlist[n])) {
          fprintf(stderr,"[%s] ERROR: niik_image_nregister_demons_update %i\n",fcname,n);
          return 0;
        }
        if(!niik_image_filter_gaussian_update(warpimglist[n],wFWHM->v[level]*2,wFWHM->v[level])) {
          fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian_update, %i\n",fcname,n);
          return 0;
        }
        if(verbose>=10) {
          fprintf(stdout,"[%s] writing output image\n",fcname);
          sprintf(fname,"tmp_nregimg_warp%i_img%i.nii.gz",iter+1,n+1);
          niik_image_write(fname,warpimglist[n]);
          /*sprintf(fname,"tmp_nregimg_diff%i_img%i.nii.gz",iter+1,n+1);  niik_image_write(fname,refimg);*/
        }
      } /* each time point */
    } /* iteration */

    if(verbose>=10) {
      for(n=0; n<num; n++) {
        if(!niik_image_nregister_intrasubject_apply_warp(warpimglist[n],refimg,imglist[n],matlist[n],NIIK_INTERP_BSPLINE)) {
          fprintf(stderr,"[%s] ERROR: niik_image_nregister_intrasubject_apply_warp\n",fcname);
          return 0;
        }
        sprintf(fname,"tmp_nregimg_img%i.nii.gz",n+1);
        niik_image_write(fname,refimg);
      } /* each image */
    } /* writing warped image */

    fprintf(stdout,"[%s] making a new reference image\n",fcname);
    if(!niik_image_nregister_average_with_fov(refimg,warpimglist,imglist,matlist,num,NIIK_INTERP_BSPLINE)) {
      fprintf(stderr,"[%s] ERROR: niik_image_nregister_average_with_fov\n",fcname);
      return 0;
    }
    if(verbose>=10) {
      sprintf(fname,"tmp_nregimg_refimg%i.nii.gz",level+1);
      niik_image_write(fname,refimg);
    }

    gradimg=niik_image_free(gradimg);
  } /* level */

  imat=niikmat_free(imat);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_nregister_intrasubject */


int niik_image_nregister_average_with_fov(nifti_image *refimg,nifti_image **warpimg,nifti_image **imglist,niikmat **matlist,int num,int interp)
/*
 * average images with field-of-view weighting
 *
 */
{
  char fcname[128]="niik_image_nregister_average_with_fov";
  nifti_image
  *fovimg=NULL,
   *tmpimg=NULL;
  niikmat
  *ww=NULL,
   *afmat=NULL;
  niikpt
  pfov,p,q;
  int
  m,n,i,j,k;
  double *dimg=NULL,*dfov=NULL;
  if( refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return 0;
  }
  if(matlist==NULL) {
    fprintf(stderr,"[%s] ERROR: matlist is null\n",fcname);
    return 0;
  }
  /*interp=NIIK_INTERP_LINEAR;*/
  fprintf(stdout,"[%s] %s \n",fcname,niik_interpolate_string(interp));
  dimg = (double *)calloc(refimg->nvox,sizeof(double));
  dfov = (double *)calloc(refimg->nvox,sizeof(double));
  for(n=0; n<num; n++) {
    if((afmat = niikmat_inverse(matlist[n]))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse, %i \n",n);
      return 0;
    }
    if((fovimg = niik_image_copy(imglist[n]))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(!niik_image_one(fovimg)) {
      fprintf(stderr,"ERROR: niik_image_one\n");
      return 0;
    }
    if(interp==NIIK_INTERP_BSPLINE) {
      if((tmpimg = niik_image_copy(imglist[n]))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy\n");
        return 0;
      }
      if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT32)) {
        fprintf(stderr,"ERROR: niik_image_type_convert, %i\n",n);
        return 0;
      }
      if(!niik_image_interpolate_convert_3d_bspline_coeff(tmpimg)) {
        fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff, %i\n",n);
        return 0;
      }
    } else {
      tmpimg = imglist[n];
    }
    pfov.x=(imglist[n]->nx-1)*imglist[n]->dx;
    pfov.y=(imglist[n]->ny-1)*imglist[n]->dy;
    pfov.z=(imglist[n]->nz-1)*imglist[n]->dz;
    ww=niikmat_init(refimg->nz,3);
    #pragma omp parallel for private(i,j,p,m,q)
    for(k=0; k<refimg->nz; k++) {
      p.w=0;
      p.z=k*refimg->dz;
      m=k*refimg->nx*refimg->ny;
      for(j=0; j<refimg->ny; j++) {
        p.y=j*refimg->dy;
        for(i=0; i<refimg->nx; m++,i++) {
          p.x=i*refimg->dx;
          niik_image_interpolate_3d_xyz_update(warpimg[n],p,NIIK_INTERP_LINEAR,ww->m[k]);
          q.x=p.x+ww->m[k][0];
          q.y=p.y+ww->m[k][1];
          q.z=p.z+ww->m[k][2];
          q = niikpt_affine_transform(afmat,q);
          if(interp==NIIK_INTERP_BSPLINE) {
            if     (q.x<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.y<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.z<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.x>pfov.x-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.y>pfov.y-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.z>pfov.z-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else {
              dimg[m] += niik_image_interpolate_3d(tmpimg,q,NIIK_INTERP_BSPLINE);
            }
          } else {
            dimg[m] += niik_image_interpolate_3d(tmpimg,q,interp);
          }
          dfov[m] += niik_image_interpolate_3d(fovimg,q,NIIK_INTERP_LINEAR);
          /*if(i==80 && j==138 && k==18) fprintf(stdout,"[%3i %3i %3i, %3i] %8.4f %5.2f  | %6.1f %6.1f %6.1f\n",
                  i,j,k,n,dimg[m],dfov[m],q.x,q.y,q.z);*/
        }
      }
    }
    if(interp==NIIK_INTERP_BSPLINE) {
      tmpimg=niik_image_free(tmpimg);
    }
    fovimg = niik_image_free(fovimg);
    afmat = niikmat_free(afmat);
    ww=niikmat_free(ww);
  } /* each image */
  for(i=0; i<refimg->nvox; i++) {
    if(dfov[i]<0.25) dimg[i]=0;
    else dimg[i] = dimg[i] / dfov[i];
  }
  if(!niik_image_set_voxels_from_double_vector(refimg,dfov)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  /*niik_image_write("tmp_fov.nii.gz",refimg);*/
  if(!niik_image_set_voxels_from_double_vector(refimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(dimg);
  free(dfov);
  return 1;
} /* niik_image_nregister_average_with_fov */


/*
 *
 * Convet between niikmat and nifti_image
 *
 * -using intent_code
 *
 */

nifti_image *niikmat_convert_to_image(niikmat *afmat) {
  nifti_image *out=NULL;
  char fcname[32]="niikmat_convert_to_image";
  int i,j,k;
  if(afmat==NULL) {
    fprintf(stderr,"[%s] ERROR: afmat is null\n",fcname);
    return NULL;
  }
  if((out=niik_image_init(1,1,1,1,afmat->row*afmat->col,1,1,
                          1,1,1,1,1,1,1,NIFTI_TYPE_FLOAT64))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_init\n",fcname);
    return NULL;
  }
  out->intent_code=NIFTI_INTENT_GENMATRIX;
  out->intent_p1=afmat->col;
  out->intent_p2=afmat->row;
  out->dim[0]=out->ndim=5;
  for(i=k=0; i<afmat->col; i++) {
    for(j=0; j<afmat->row; k++,j++) {
      niik_image_set_voxel(out,k,afmat->m[i][j]);
    }
  }
  return out;
} /* niikmat_convert_to_image */

niikmat *niikmat_convert_from_image(nifti_image *matimg) {
  niikmat *out=NULL;
  char fcname[32]="niikmat_convert_from_image";
  int i,j,k;
  if(matimg==NULL) {
    fprintf(stderr,"[%s] ERROR: matimg is null\n",fcname);
    return NULL;
  }
  if(matimg->intent_code!=NIFTI_INTENT_GENMATRIX) {
    fprintf(stderr,"[%s] ERROR: intent_code is not NIFTI_INTENT_GENMATRIX (%s)\n",
            fcname,nifti_intent_string(matimg->intent_code));
    return NULL;
  }
  if((out=niikmat_init((int)floor(0.5+matimg->intent_p1),
                       (int)floor(0.5+matimg->intent_p2)))==NULL) {
    fprintf(stderr,"[%s] ERROR: niikmat_init\n",fcname);
    return NULL;
  }
  for(i=k=0; i<out->col; i++) {
    for(j=0; j<out->row; k++,j++) {
      out->m[i][j]=niik_image_get_voxel(matimg,k);
    }
  }
  return out;
} /* niikmat_convert_from_image */




#endif /* _FALCON_NREGISTER_INTRASUBJECT_TEMPLATE_C_ */
