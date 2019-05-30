/* Filename:     nifti1_kunio_nregister_bspline.c
 * Description:  bspline-based nonlinear registration functions
 * Author:       Kunio Nakamura
 * Date:         April 9, 2012
 *
 * Parent:       nifti1_kunio_nregister.c
 *
 * Functions:
   int    niik_image_nregister_bspline
   int    niik_image_nregister_bspline_update_warp_coeff
   double niik_image_nregister_bspline_obj_func
 *
 *
 */

#ifndef _FALCON_NREGISTER_BSPLINE_C_
#define _FALCON_NREGISTER_BSPLINE_C_

#include "falcon.h"

int niik_image_nregister_bspline_update_warp_coeff(nifti_image *img,nifti_image *refimg,double delta,double regularization_size);


/**************************************************************************************
 *
 * FUNCTIONS FOR NONLINEAR REGISTRATION
 *
 **************************************************************************************/

double niik_image_nregister_bspline_obj_func(nifti_image *refimg,nifti_image *refseg,
    nifti_image *movimg,nifti_image *movseg,
    nifti_image *wbsimg,
    nifti_image *refgradimg,nifti_image *movgradimg,
    double reg_dist, double grad_weight,
    int xmin,int ymin,int zmin,
    int xmax,int ymax,int zmax)
/* objective function for niik_image_nregister_bspline
 * -refimg and refseg are the reference (target) image and its mask
 * -movimg and movseg are the moving image and its mask
 *   -both masks are not currently ued
 * -xmin,ymin,zmin ...  are the inclusive cooordinates in refime's space
 *  where the similarity is computed.
 */
{
  double
  vr,vm,
  *vrlist,*vmlist,
  gr,gm,gmax=0,gv,
        rsum,msum,csum,rssq,mssq,
        E=0,C=0,G=0;
  niikpt
  fovM,
  plim,
  pt,pt2,qt;
  int
  i,j,k,n,
  ii,jj,kk,nn,
  xdim,ydim,zdim,xydim,area,size,
  verbose=0,
  nvox=0;
  char fcname[64]="niik_image_nregister_bspline_obj_func]";
  unsigned char
  *bimg;
  float
  *fimg,
  *fx,*fy,*fz,*fv;
  niikmat *bt;

  /************************
   * check inputs
   ************************/
  if(verbose)
    fprintf(stdout,"[%s] start %4i %4i %4i -> %4i %4i %4i\n",fcname,
            xmin,ymin,zmin,xmax,ymax,zmax);
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"ERROR: movimg is null\n");
    return 0;
  }
  if(refseg==NULL) {
    fprintf(stderr,"ERROR: refseg is null\n");
    return 0;
  }
  if(verbose)  {
    if(refgradimg==NULL && movgradimg==NULL) {
      fprintf(stdout,"[niik_image_nregister_bspline_obj_func] no gradient images\n");
    } else if(refgradimg!=NULL && movgradimg!=NULL) {
      fprintf(stdout,"[niik_image_nregister_bspline_obj_func] using gradient images\n");
    } else {
      fprintf(stdout,"[niik_image_nregister_bspline_obj_func] only one gradient images\n");
      return 0;
    }
  }

  if(wbsimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[niik_image_nregister_bspline_obj_func] ERROR: wbsimg is not NIFTI_TYPE_FLOAT32\n");
    exit(0);
  }

  if(refimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[niik_image_nregister_bspline_obj_func] ERROR: refimg is not NIFTI_TYPE_FLOAT32\n");
    exit(0);
  }
  fimg = (float *)refimg->data;
  if(refseg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"[niik_image_nregister_bspline_obj_func] ERROR: refseg is not NIFTI_TYPE_UINT8\n");
    exit(0);
  }
  bimg = (unsigned char *)refseg->data;

  /*fovR.x = refimg->nx * refimg->dx;
    fovR.y = refimg->ny * refimg->dy;
    fovR.z = refimg->nz * refimg->dz;*/
  fovM.x = movimg->nx * movimg->dx;
  fovM.y = movimg->ny * movimg->dy;
  fovM.z = movimg->nz * movimg->dz;

  plim.x = -(1.0-reg_dist)*refimg->dx;
  plim.y = -(1.0-reg_dist)*refimg->dy;
  plim.z = -(1.0-reg_dist)*refimg->dz;

  xdim = (xmax-xmin) + 1;
  ydim = (ymax-ymin) + 1;
  zdim = (zmax-zmin) + 1;
  xydim = xdim*ydim;
  size = xdim * ydim * zdim;
  area = refimg->nx * refimg->ny;

  if(verbose) fprintf(stdout,"[%s] dim = %4i %4i %4i\n",fcname,xdim,ydim,zdim);

  fv = fx = (float *)calloc(size*3,sizeof(float));
  fy = fx + size;
  fz = fy + size;

  rsum=msum=csum=rssq=mssq=0;
  if(verbose)
    fprintf(stdout,"[niik_image_nregister_bspline_obj_func] main loop\n");

  bt=niikmat_init(zmax+1,5);

  if(1) {
    pt.w=0;
    for(k=zmin,kk=0,n=0; k<=zmax; kk++,k++) {
      pt.z  = k * refimg->dz;
      pt2.z = pt.z / wbsimg->dz;
      for(j=ymin,jj=0; j<=ymax; jj++,j++) {
        pt.y  = j * refimg->dy;
        pt2.y = pt.y / wbsimg->dy;
        for(i=xmin,ii=0,nn=xmin+j*refimg->nx+k*area; i<=xmax; nn++,ii++,i++,n++) {
          pt.x  = i * refimg->dx;
          pt2.x = pt.x / wbsimg->dx;
          /* GET DISPLACEMENT VECTOR */
          if(!niik_image_interpolate_3d_bspline_ijk_update(wbsimg,pt2,bt->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_ijk_update\n");
            exit(0);
          }
          /* UPDATE POSITION COORDINATES */
          qt.x = pt.x + bt->m[k][0];
          qt.y = pt.y + bt->m[k][1];
          qt.z = pt.z + bt->m[k][2];
          qt.w = 0;
          fx[n]=bt->m[k][0];
          fy[n]=bt->m[k][1];
          fz[n]=bt->m[k][2];
          /* check for positive jacobian */
          if(ii>0) {
            if((fx[n]-fx[n-1])<=plim.x) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f x\n",ii,jj,kk,fx[n-1],fx[n]);
              E+=100;
            }
          }
          if(jj>0) {
            if((fy[n]-fy[n-xdim])<=plim.y) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f y\n",ii,jj,kk,fy[n-xdim],fy[n]);
              E+=100;
            }
          }
          if(kk>0) {
            if((fz[n]-fz[n-xydim])<=plim.z) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f z\n",ii,jj,kk,fz[n-xydim],fz[n]);
              E+=100;
            }
          }
          if(bimg[nn]==0) continue;
          if(qt.x>=fovM.x) continue;
          if(qt.y>=fovM.y) continue;
          if(qt.z>=fovM.z) continue;
          vr = fimg[nn]; /*niik_image_interpolate_3d(refimg,pt,NIIK_INTERP_LINEAR);*/
          vm = niik_image_interpolate_3d(movimg,qt,NIIK_INTERP_LINEAR);
          rsum += vr;
          msum += vm;
          rssq += vr*vr;
          mssq += vm*vm;
          csum += vr*vm;
          nvox++;
          /* sanity check */
          if(!(vr>-1e35 && vr<1e35)) {
            fprintf(stderr,"\n\n  vr is un-real\n");
            fprintf(stdout,"  nn     %i\n",nn);
            fprintf(stdout,"  vr     %12.3f\n",vr);
            fprintf(stdout,"  ijk    %3i %3i %3i\n",i,j,k);
            exit(0);
          }
          /*if(vr<-2 || vm<-2) {
            fprintf(stderr,"\n\n  vm/vr is un-real\n");
            fprintf(stdout,"  nn     %i\n",nn);
            fprintf(stdout,"  vr     %12.3f\n",vr);
            fprintf(stdout,"  vm     %12.3f\n",vm);
            fprintf(stdout,"  ijk    %3i %3i %3i\n",i,j,k);
            fprintf(stdout,"  p      %7.2f %7.2f %7.2f\n",pt.x,pt.y,pt.z);
            fprintf(stdout,"  q      %7.2f %7.2f %7.2f\n",qt.x,qt.y,qt.z);
            return 1e10; }*/
          if(refgradimg!=NULL) {
            gr = niik_image_interpolate_3d(refgradimg,pt,NIIK_INTERP_LINEAR);
            gm = niik_image_interpolate_3d(movgradimg,pt,NIIK_INTERP_LINEAR);
            gv = sqrt(gr*gm);
            gmax = NIIK_IMAX(gv,gmax);
            G += gv;
          }
          if(verbose>1)
            fprintf(stdout,"%7.3f %7.3f %7.3f [%5.0f]  ->  %7.3f %7.3f %7.3f [%5.0f]  |  %8.4f %8.4f %8.4f\n",
                    pt.x,pt.y,pt.z,vr,  qt.x,qt.y,qt.z,vm,  bt->m[k][0],bt->m[k][1],bt->m[k][2]);
        }
      }
    }
  }

  else { /* testing */
    vrlist=(double *)calloc(size,sizeof(double));
    vmlist=(double *)calloc(size,sizeof(double));

    #pragma omp parallel for private(i,j,pt,pt2,qt,n,nn)
    for(k=0; k<zdim; k++) {
      pt.z=(k+zmin)*refimg->dz;
      pt2.z=pt.z/wbsimg->dz;
      for(j=0; j<ydim; j++) {
        pt.y=(j+ymin)*refimg->dy;
        pt2.y=pt.y/wbsimg->dy;
        nn=xmin+(j+ymin)*refimg->nx+(k+zmin)*area;
        for(i=0; i<xdim; i++,nn++) {
          pt.x=(i+xmin)*refimg->dx;
          pt2.x=pt.x/wbsimg->dx;
          /* GET DISPLACEMENT VECTOR */
          if(!niik_image_interpolate_3d_bspline_ijk_update(wbsimg,pt2,bt->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_ijk_update\n");
            exit(0);
          }
          /* UPDATE POSITION COORDINATES */
          qt.x = pt.x + bt->m[k][0];
          qt.y = pt.y + bt->m[k][1];
          qt.z = pt.z + bt->m[k][2];
          n=i+j*xdim+k*xydim;
          fx[n]=bt->m[k][0];
          fy[n]=bt->m[k][1];
          fz[n]=bt->m[k][2];
          vrlist[n]=-1e11;
          vmlist[n]=-1e11;
          /* check for positive jacobian */
          if(i>0) {
            if((fx[n]-fx[n-1])<=plim.x) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f x\n",i,j,k,fx[n-1],fx[n]);
              E+=100;
            }
          }
          if(j>0) {
            if((fy[n]-fy[n-xdim])<=plim.y) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f y\n",i,j,k,fy[n-xdim],fy[n]);
              E+=100;
            }
          }
          if(k>0) {
            if((fz[n]-fz[n-xydim])<=plim.z) {
              if(verbose) fprintf(stdout,"\t%3i %3i %3i   %9.4f %9.4f z\n",i,j,k,fz[n-xydim],fz[n]);
              E+=100;
            }
          }
          if(bimg[nn]==0) continue;
          if(qt.x>=fovM.x) continue;
          if(qt.y>=fovM.y) continue;
          if(qt.z>=fovM.z) continue;
          vrlist[n] = fimg[nn]; /*niik_image_interpolate_3d(refimg,pt,NIIK_INTERP_LINEAR);*/
          vmlist[n] = niik_image_interpolate_3d(movimg,qt,NIIK_INTERP_NN);
        }
      }
    }

    for(n=0,nvox=0; n<size; n++) {
      if(vrlist[n]<-1e10) continue;
      if(vmlist[n]<-1e10) continue;
      rsum += vrlist[n];
      msum += vmlist[n];
      rssq += vrlist[n]*vrlist[n];
      mssq += vmlist[n]*vmlist[n];
      csum += vrlist[n]*vmlist[n];
      nvox++;
    }
    free(vrlist);
    free(vmlist);
  }

  free(fv);
  bt=niikmat_free(bt);

  /* correlation coefficient calculation */
  if     (fabs(nvox*rssq-rsum*rsum)<0.1) {
    C = -1;
  } else if(fabs(nvox*mssq-msum*msum)<0.1) {
    C = -1;
  } else {
    C = (nvox*csum - rsum * msum) / sqrt(nvox*rssq-rsum*rsum) / sqrt(nvox*mssq-msum*msum);
  }

  /* sanity check */
  if(isnan(C)) {
    fprintf(stderr,"\n\n  C is nan\n");
    fprintf(stdout,"  nvox   %i\n",nvox);
    fprintf(stdout,"  csum   %12.3f\n",csum);
    fprintf(stdout,"  rsum   %12.3f\n",rsum);
    fprintf(stdout,"  msum   %12.3f\n",msum);
    fprintf(stdout,"  mssq   %12.3f\n",mssq);
    fprintf(stdout,"  rssq   %12.3f\n",rssq);
    fprintf(stderr,"  C = %8.6lf / SQRT(%8.6lf) / SQRT(%8.6lf)\n",nvox*csum-rsum*msum,
            nvox*rssq-rsum*rsum, nvox*mssq-msum*msum);
    return -1e9;
  } else if(!(C>-1e35 && C<1e35)) {
    fprintf(stderr,"\n\n  C is un-real\n");
    fprintf(stdout,"  nvox   %i\n",nvox);
    fprintf(stdout,"  csum   %12.3f\n",csum);
    fprintf(stdout,"  rsum   %12.3f\n",rsum);
    fprintf(stdout,"  msum   %12.3f\n",msum);
    fprintf(stdout,"  mssq   %12.3f\n",mssq);
    fprintf(stdout,"  rssq   %12.3f\n",rssq);
    return -1e9;
  }
  if(isnan(E)) {
    fprintf(stderr,"\n\n  E is nan\n");
    return -1e9;
  } else if(!(E>-1e35 && E<1e35)) {
    fprintf(stderr,"\n\n  E is un-real\n");
    return -1e9;
  }

  if(nvox*8 < size) return 0;
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_obj_func] %i / %i  %8.5f\n",nvox,size,C);

  if(refgradimg!=NULL) G = G / gmax / nvox * grad_weight;

  if(verbose) {
    fprintf(stdout,"[niik_image_nregister_bspline_obj_func] outputs\n");
    fprintf(stdout,"  C %-12.5f G %-12.5f E %-12.5f %s\n",C,G,E,wbsimg->fname);
  }

  return C - E + G;
} /* niik_image_nregister_bspline_obj_func */



int niik_image_nregister_bspline(nifti_image *refimg,nifti_image *refseg,
                                 nifti_image *movimg,nifti_image *movseg,
                                 nifti_image *wbsimg,
                                 int flag_sobel,
                                 double *idelta,double *wdelta,double *filFWHM,int *maxiter,double *maxdfm,
                                 int nlevel,int localsize,int ndir,double reg_dist,double grad_weight,double dfm_fc)
/* nonlinear registration function
 * -based on b-spline
 * -refimg is the target image
 * -refseg is target's ROI (optional)
 * -movimg is the moving image
 * -movseg is the moving image's ROI (NOT USED)
 * -wbsimg is replaced with warp bspline coefficient image
 * -flag_sobel is the flag for using sobel filter images
 *    if 1, use 3-direction for deformation along the gradient direction
 *    if 2, use sobel filtered gradient image for a part of registration cost function
 * -idelta is image spacing for each level
 * -wdelta is the warp spacing for each level
 * -maxiter is the maximum iteration for each level
 * -maxdfm is the max deformation for each level
 * -nlevel is the number of levels
 * -localsize is the size of local ROI
 * -ndir is the number of warp directions to test
 * -reg_dist is the regularization distance factor (probably 0.1-0.5, never below 0 or above 1)
 * -dfm_fc is the regularization weighting factor usually quite small
 */

{

  nifti_image
  *tmpimglist[4],
  **sobelimg=NULL,
    **warplist=NULL,
      *crefimg=NULL,
       *crefseg=NULL,
        *cmovimg=NULL;
  niikpt
  cw,
  pt,st,
  *dirlist=NULL;
  niikvec
  *num_mask_vec=NULL;
  char
  fname[512];
  unsigned char
  *bimg;
  int
  m,n,
  i,j,k,p,
  ii,jj,kk,nn,mm,
  xlo,xhi,ylo,yhi,zlo,zhi,
  nvox3d1,
  verbose=1,
  ndfm,nskip,
  iter,
  level;
  float
  *fs[3];
  float
  *wxv[20],*wyv[20],*wzv[20];
  double
  *C,Cprev=0,Ccurr=0,
     dfm;
  struct tm *stm; /* showing times */
  time_t ctm;
  char tmstr[256];
  char fcname[64]="niik_image_nregister_bspline";

  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline] start\n");
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  /* if(refseg==NULL) { fprintf(stderr,"ERROR: refseg is null\n"); return 0; } */
  if(movimg==NULL) {
    fprintf(stderr,"ERROR: movimg is null\n");
    return 0;
  }
  /* if(movseg==NULL) { fprintf(stderr,"ERROR: movseg is null\n"); return 0; } */
  if(wbsimg==NULL) {
    fprintf(stderr,"ERROR: wbsimg is null\n");
    return 0;
  }
  if(idelta==NULL) {
    fprintf(stderr,"ERROR: idelta is null\n");
    return 0;
  }
  if(wdelta==NULL) {
    fprintf(stderr,"ERROR: wdelta is null\n");
    return 0;
  }
  if(maxdfm==NULL) {
    fprintf(stderr,"ERROR: maxdfm is null\n");
    return 0;
  }
  if(filFWHM==NULL) {
    fprintf(stderr,"ERROR: filFWHM is null\n");
    return 0;
  }
  if(maxiter==NULL) {
    fprintf(stderr,"ERROR: maxiter is null\n");
    return 0;
  }
  if(reg_dist<0) {
    fprintf(stderr,"ERROR: reg_dist is below zero %9.5f\n",reg_dist);
    return 0;
  }
  if(reg_dist>1) {
    fprintf(stderr,"ERROR: reg_dist is above one %9.5f\n",reg_dist);
    return 0;
  }

  if(!niik_image_type_convert(wbsimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(wbsimg,NIFTI_TYPE_FLOAT32)\n");
    return 0;
  }
  if(!niik_image_type_convert(movimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(movimg,NIFTI_TYPE_FLOAT32)\n");
    return 0;
  }
  if(!niik_image_type_convert(refimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(refimg,NIFTI_TYPE_FLOAT32)\n");
    return 0;
  }
  if(refseg!=NULL) {
    if(!niik_image_type_convert(refseg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"ERROR: niik_image_type_convert(refseg,NIFTI_TYPE_UINT8)\n");
      return 0;
    }
  }
  if(movseg!=NULL) {
    if(!niik_image_type_convert(movseg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"ERROR: niik_image_type_convert(movseg,NIFTI_TYPE_UINT8)\n");
      return 0;
    }
  }

  if(verbose>1 && 0) {
    fprintf(stdout,"\twriting initial warp\n");
    if(!niik_image_write("tmp_nregimg_iwbsp.nii.gz",wbsimg)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      return 0;
    }
  }

  if(flag_sobel==1) ndir=3;

  fprintf(stdout,"[%s] bspline-based nonlinear registration\n",fcname);
  fprintf(stdout,"  ref img      %s\n",refimg->fname);
  if(refseg!=NULL) fprintf(stdout,"  ref mask     %s\n",refseg->fname);
  fprintf(stdout,"  mov img      %s\n",movimg->fname);
  if(movseg!=NULL) fprintf(stdout,"  mov mask     %s\n",movseg->fname);
  fprintf(stdout,"  img delta    ");
  for(n=0; n<nlevel; n++) fprintf(stdout,"%7.3f ",idelta[n]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  warp delta   ");
  for(n=0; n<nlevel; n++) fprintf(stdout,"%7.3f ",wdelta[n]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  FWHM         ");
  for(n=0; n<nlevel; n++) fprintf(stdout,"%7.3f ",filFWHM[n]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  max deform   ");
  for(n=0; n<nlevel; n++) fprintf(stdout,"%7.3f ",maxdfm[n]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  max iter     ");
  for(n=0; n<nlevel; n++) fprintf(stdout,"%i ",maxiter[n]);
  fprintf(stdout,"\n");
  if(flag_sobel==1)
    fprintf(stdout,"  #dir         %i along gradient direction\n",ndir);
  else if(flag_sobel==2) {
    fprintf(stdout,"  #dir         %i\n",ndir);
    fprintf(stdout,"               using gradient in cost function\n");
  } else
    fprintf(stdout,"  #dir         %i\n",ndir);
  fprintf(stdout,"  local size   %i\n",localsize);
  fprintf(stdout,"  regularization size  %8.4f\n",reg_dist);

  if(verbose>1) fprintf(stdout,"[%s] scale\n",fcname);
  if(!niik_image_iscale(refimg,3,niik_image_get_percentile(refimg,NULL,0.95),0,1000)) {
    fprintf(stderr,"[%s] ERROR: niik_image_iscale\n",fcname);
    return 0;
  }
  if(!niik_image_iscale(movimg,3,niik_image_get_percentile(movimg,NULL,0.95),0,1000)) {
    fprintf(stderr,"[%s] ERROR: niik_image_iscale\n",fcname);
    return 0;
  }

  if(verbose>1) fprintf(stdout,"[%s] define dir vec\n",fcname);
  C=(double *)calloc(ndir,sizeof(double));
  dirlist=(niikpt *)calloc(ndir,sizeof(niikpt));
  for(n=0; n<ndir; n++) dirlist[n]=niikpt_zero();
  switch(ndir) {
  case 4:  /* kind of random --i made them up */
    dirlist[1].x = 0.9;
    dirlist[1].y = 0.3;
    dirlist[1].z = -0.2;
    dirlist[2].x =-0.2;
    dirlist[2].y = 0.8;
    dirlist[2].z =  0.3;
    dirlist[3].x = 0.3;
    dirlist[3].y =-0.2;
    dirlist[3].z =  0.7;
    break;
  case 5: /* tetrahedron */
    n=1;
    dirlist[n].x = dirlist[n].y = dirlist[n].z =  1.0;
    n++;
    dirlist[n].x = dirlist[n].y = -1;
    dirlist[n].z =  1.0;
    n++;
    dirlist[n].x = dirlist[n].z = -1;
    dirlist[n].y =  1.0;
    n++;
    dirlist[n].x = 1.0;
    dirlist[n].y = dirlist[n].z = -1.0;
    n++;
    break;
  case 7: /* standard set */
    dirlist[1].x = dirlist[2].y = dirlist[3].z =  1.0;
    dirlist[4].x = dirlist[5].y = dirlist[6].z = -1.0;
    break;
  case 3:
    if(flag_sobel) break;
  default:
    fprintf(stderr,"[%s] ERROR: unknown #dir %i\n",fcname,ndir);
    return 0;
  }
  if(!flag_sobel) {
    for(n=1; n<ndir; n++) dirlist[n] = niikpt_unit(dirlist[n]);
    if(verbose) {
      fprintf(stderr,"[%s] deformation directions\n",fcname);
      for(n=0; n<ndir; n++) fprintf(stdout,"  %-2i %8.3f %8.3f %8.3f\n",n,dirlist[n].x,dirlist[n].y,dirlist[n].z);
    }
  }

  if((warplist = (nifti_image **)calloc(ndir,sizeof(nifti_image *)))==NULL) {
    fprintf(stderr,"[%s] ERROR: callc warplist\n",fcname);
    return 0;
  }

  switch(flag_sobel) {
  case 1:
    if(verbose) fprintf(stdout,"[%s] sobel gradient images\n",fcname);
    if((sobelimg = niik_image_sobel_filters(movimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_sobel_filters\n");
      return 0;
    }
    for(n=0; n<3; n++) {
      if(!niik_image_type_convert(sobelimg[n],NIFTI_TYPE_FLOAT32)) {
        fprintf(stderr,"ERROR: niik_image_type_convert(sobelimg,NIFTI_TYPE_FLOAT32)\n");
        return 0;
      }
      fs[n]=(float *)sobelimg[n]->data;
    }
  case 2:
    if(verbose) fprintf(stdout,"[%s] sobel gradient images\n",fcname);
    sobelimg = (nifti_image **)calloc(3,sizeof(nifti_image *));
    if((sobelimg[0] = niik_image_sobel_filter(refimg,'m'))==NULL) {
      fprintf(stderr,"ERROR: niik_image_sobel_filters\n");
      return 0;
    }
    if((sobelimg[1] = niik_image_sobel_filter(movimg,'m'))==NULL) {
      fprintf(stderr,"ERROR: niik_image_sobel_filters\n");
      return 0;
    }
    break;
  case 0:
    if((sobelimg = (nifti_image **)calloc(3,sizeof(nifti_image *)))==NULL) {
      fprintf(stderr,"ERROR: sobelimg\n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown flag for gradient %i\n",flag_sobel);
    return 0;
  } /* flag_sobel */

  if(wbsimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: wbsimg is %s 1\n",fcname,nifti_datatype_string(warplist[n]->datatype));
    return 0;
  }

  /* MAIN LOOP */
  for(level=0; level<nlevel; level++) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);

    fprintf(stdout,"  level %-2i img-delta %-5.1f warp-delta %-5.1f FWHM %-5.1f iter %-2i dfm_max %-5.2f  dim %3i  %s\n",
            level+1,idelta[level],wdelta[level],filFWHM[level],maxiter[level],maxdfm[level],
            (int)floor(localsize*wdelta[level] / idelta[level] + 1.5) - (int)floor(-localsize*wdelta[level] / idelta[level]) + 1,
            tmstr);

    /* prepare deformation field */
    if(verbose)   fprintf(stdout,"[niik_image_nregister_bspline] level %i  prepare/update/resample deformation field\n",level+1);
    if(!niik_image_nregister_bspline_update_warp_coeff(wbsimg,refimg,wdelta[level],reg_dist)) {
      fprintf(stderr,"ERROR: niik_image_nregister_update_warp_coeff\n");
      return 0;
    }
    if(wbsimg->datatype!=NIFTI_TYPE_FLOAT32) {
      fprintf(stderr,"[niik_image_nregister_bspline] ERROR: wbsimg is %s\n",nifti_datatype_string(wbsimg->datatype));
      return 0;
    }

    nvox3d1 = wbsimg->nx * wbsimg->ny * wbsimg->nz;

    /* blur image */
    if(verbose) fprintf(stdout,"[niik_image_nregister_bspline] level %i  blur / resample images\n",level+1);
    if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline] level %i    refimg\n",level+1);
    if((crefimg = niik_image_filter_gaussian(refimg,filFWHM[level]*2,filFWHM[level]))==NULL) {
      fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
      return 0;
    }
    if(!niik_image_resample_3d_update(crefimg,idelta[level],idelta[level],idelta[level],-1,-1,-1,NIIK_INTERP_LINEAR)) {
      fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
      return 0;
    }
    if(!niik_image_type_convert(crefimg,NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"ERROR: niik_image_type_convert(crefimg,NIFTI_TYPE_FLOAT32)\n");
      return 0;
    }
    if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline] level %i    refseg\n",level+1);
    if(refseg!=NULL) { /* resample refseg using nearest neighbor interpolation */
      if((crefseg = niik_image_affine_transform_3d(refseg,crefimg,NULL,NIIK_INTERP_NN))==NULL) {
        fprintf(stderr,"ERROR: nifti *niik_image_affine_transform_3d(refseg,crefimg,NULL,NIIK_INTERP_NN)\n");
        return 0;
      }
      if(!niik_image_type_convert(crefseg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"ERROR: niik_image_type_convert\n");
        return 0;
      }
    } else {
      if((crefseg=niik_image_copy_as_type(crefimg,NIFTI_TYPE_UINT8))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
        return 0;
      }
      if(!niik_image_one(crefseg)) {
        fprintf(stderr,"ERROR: niik_image_one\n");
        return 0;
      }
    }
    bimg = (unsigned char *)crefseg->data;
    if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline] level %i    movimg\n",level+1);
    if((cmovimg = niik_image_filter_gaussian(movimg,filFWHM[level]*2,filFWHM[level]))==NULL) {
      fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
      return 0;
    }

    /* testing
    fprintf(stdout,"[niik_image_nregister_bspline] resample 0.3 iso\n");
    if(!niik_image_resample_3d_update(cmovimg,0.3,0.3,0.3,-1,-1,-1,NIIK_INTERP_BSPLINE)){
      fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
      return 0; }
    */

    if(verbose>2) {
      sprintf(fname,"tmp_nregimg_crefimg%i.nii.gz",level);
      fprintf(stdout,"\twriting %s\n",fname);
      niik_image_write(fname,crefimg);
      sprintf(fname,"tmp_nregimg_cmovimg%i.nii.gz",level);
      fprintf(stdout,"\twriting %s\n",fname);
      niik_image_write(fname,cmovimg);
      sprintf(fname,"tmp_nregimg_crefseg%i.nii.gz",level);
      fprintf(stdout,"\twriting %s\n",fname);
      niik_image_write(fname,crefseg);
    }

    if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline] level %i  update warp list %i\n",level+1,ndir);
    if(wbsimg->datatype!=NIFTI_TYPE_FLOAT32) {
      fprintf(stderr,"[niik_image_nregister_bspline] ERROR: wbsimg is %s\n",nifti_datatype_string(wbsimg->datatype));
      return 0;
    }

    for(n=0; n<ndir; n++) {
      /*fprintf(stdout,"\twarp image list %i\n",n);*/
      if((warplist[n] = niik_image_copy(wbsimg))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy %i\n",n);
        return 0;
      }
      if(warplist[n]->datatype!=NIFTI_TYPE_FLOAT32) {
        fprintf(stderr,"[niik_image_nregister_bspline] ERROR: warplist[%i] is %s\n",n,nifti_datatype_string(warplist[n]->datatype));
        return 0;
      }
      wxv[n]=(float *)warplist[n]->data;
      wyv[n]=wxv[n]+nvox3d1;
      wzv[n]=wyv[n]+nvox3d1;
      free(warplist[n]->fname);
      warplist[n]->fname=(char *)calloc(12,sizeof(char));
      sprintf(warplist[n]->fname,"w%i.nii.gz",n+1);
    }
    /*fprintf(stdout,"\twarp image list fin\n");*/

    Ccurr = niik_image_nregister_bspline_obj_func(refimg,refseg,movimg,movseg, wbsimg, sobelimg[0],sobelimg[1],
            -10, /*reg_dist,*/   /* regularization dist */
            0,   /* grad_weight, */  /* gradient term weighting */
            0,0,0,
            refimg->nx-1,refimg->ny-1,refimg->nz-1);
    fprintf(stdout,"[niik_image_nregister_bspline] level %i cc = %9.5f\n",level+1,Ccurr);
    Cprev = Ccurr;

    num_mask_vec=niikvec_init(crefimg->nvox);
    for(k=0; k<crefimg->nvox; k++) {
      num_mask_vec->v[k] = -1;
    }

    /* for each voxel
     * -change the bspline coefficient map
     * -save the best combination
     */
    if(verbose) fprintf(stdout,"[niik_image_nregister_bspline] level %i  main inner loop\n",level+1);

    for(iter=0,dfm=maxdfm[level]; iter<maxiter[level]; dfm*=0.6,iter++) {

      if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline] level %i iter %-2i\n",level+1,iter+1);
      for(k=p=0,ndfm=nskip=0; k<wbsimg->nz; k++) {
        zlo = (k-localsize)*wbsimg->dz / idelta[level];
        zhi = (k+localsize)*wbsimg->dz / idelta[level] + 1.5;
        if(zlo<0) zlo=0;
        if(zhi>=crefimg->nz) zhi=crefimg->nz-1;
        for(j=0; j<wbsimg->ny; j++) {
          ylo = (j-localsize)*wbsimg->dy / idelta[level];
          yhi = (j+localsize)*wbsimg->dy / idelta[level] + 1.5;
          if(ylo<0) ylo=0;
          if(yhi>=crefimg->ny) yhi=crefimg->ny-1;
          for(i=0; i<wbsimg->nx; p++,i++) {
            xlo = (i-localsize)*wbsimg->dx / idelta[level];
            xhi = (i+localsize)*wbsimg->dx / idelta[level] + 1.5;
            if(xlo<0) xlo=0;
            if(xhi>=crefimg->nx) xhi=crefimg->nx-1;

            /* check the mask if data is enough */
            if(verbose>4) fprintf(stdout,"[niik_image_nregister_bspline] level %i iter %-2i [%3i %3i %3i] check mask size %8.3f\n",
                                    level+1,iter+1,i,j,k,num_mask_vec->v[p]);
            if(num_mask_vec->v[p]<0) {
              for(kk=zlo,mm=0; kk<=zhi; kk++) {
                for(jj=ylo; jj<=yhi; jj++) {
                  for(ii=xlo; ii<=xhi; ii++) {
                    nn = ii + jj * crefimg->nx + kk * crefimg->nx * crefimg->ny;
                    if(bimg[nn]) mm++;
                  }
                }
              }
              nn = (xhi-xlo+1) * (yhi-ylo+1) * (zhi-zlo+1);
              num_mask_vec->v[p] = (double)mm/nn;
            }
            if(num_mask_vec->v[p] < 0.5) {
              nskip++;
              continue;
            }

            cw = niikpt_val(wxv[0][p],wyv[0][p],wzv[0][p],0);

            if(flag_sobel==1) {
              if(verbose>4) fprintf(stdout,"[niik_image_nregister_bspline] level %i iter %-2i [%3i %3i %3i] sobel\n",level+1,iter+1,i,j,k);
              dirlist[1]=dirlist[0]=niikpt_zero();
              for(kk=zlo,mm=0; kk<=zhi; kk++) {
                pt.z = kk * crefimg->dz;
                for(jj=ylo; jj<=yhi; jj++) {
                  pt.y = jj * crefimg->dy;
                  for(ii=xlo; ii<=xhi; ii++) {
                    pt.x = ii * crefimg->dx;
                    nn = niik_image_get_index_niikpt(cmovimg,pt);
                    if(nn<0) continue;
                    st = niikpt_val(fs[0][nn],fs[1][nn],fs[2][nn],0);
                    dirlist[1] = niikpt_add(dirlist[1],st);
                  }
                }
              }
              dirlist[1] = niikpt_unit(dirlist[1]);
              dirlist[2] = niikpt_kmul(dirlist[1],-1);
              if(fabs(dirlist[1].x)<1e-5)
                if(fabs(dirlist[1].y)<1e-5)
                  if(fabs(dirlist[1].z)<1e-5)
                    continue;
            } /* using sobel filter */

            if(verbose>4) fprintf(stdout,"[%s] level %i iter %-2i [%3i %3i %3i] %8.3f\n",
                                    fcname,level+1,iter+1,i,j,k,dfm);
            for(n=0; n<ndir; n++) {
              wxv[n][p] = dfm*dirlist[n].x+cw.x;
              wyv[n][p] = dfm*dirlist[n].y+cw.y;
              wzv[n][p] = dfm*dirlist[n].z+cw.z;
            } /* ndir */

            /* CALCULATE THE SIMILARITY */
            if(verbose>4) fprintf(stdout,"[%s] level %i iter %-2i [%3i %3i %3i] similarity\n",
                                    fcname,level+1,iter+1,i,j,k);
            #pragma omp parallel for
            for(n=0; n<ndir; n++) {
              C[n] =
                niik_image_nregister_bspline_obj_func(crefimg,crefseg,cmovimg,crefseg,warplist[n],
                                                      sobelimg[0],sobelimg[1],reg_dist,grad_weight,
                                                      xlo,ylo,zlo,xhi,yhi,zhi) -
                dfm * dfm_fc * (n>0);
            } /* calculate similarity between two images */
            for(n=1,m=0; n<ndir; n++) {
              if(C[n]>C[m]) m=n;
            }
            if(C[m]<0) m=0;
            if(m>0) ndfm++;

            if(verbose>4) {
              for(n=0; n<ndir; n++) {
                if(n==m) {
                  fprintf(stdout,"[%s]   ndir %i %12.8f | %9.4f %9.4f %9.4f  | %9.4f %9.4f %9.4f *\n",
                          fcname,n,C[n],
                          wxv[n][p],wyv[n][p],wzv[n][p],
                          dirlist[n].x,dirlist[n].y,dirlist[n].z);
                } else {
                  fprintf(stdout,"[%s]   ndir %i %12.8f | %9.4f %9.4f %9.4f  | %9.4f %9.4f %9.4f\n",
                          fcname,n,C[n],
                          wxv[n][p],wyv[n][p],wzv[n][p],
                          dirlist[n].x,dirlist[n].y,dirlist[n].z);
                }
              }
            }

            if(verbose>4) {
              fprintf(stdout,"[%s]   ndir %i %12.8f | %9.4f %9.4f %9.4f <- %9.4f %9.4f %9.4f\n",
                      fcname,m,C[m],
                      wxv[m][p],wyv[m][p],wzv[m][p],
                      cw.x,cw.y,cw.z);
            }

            cw = niikpt_move_normal( cw, dirlist[m], 0.5 * dfm);
            for(n=0; n<ndir; n++) {
              wxv[n][p]=cw.x;
              wyv[n][p]=cw.y;
              wzv[n][p]=cw.z;
            } /* each warp list image, add half the warp */

          }
        }
      } /* x,y,z */

      if(verbose>4)
        fprintf(stdout,"[%s] level %i iter %-2i ndfm = %i/%i(%2i%%) skip %i C %6.4f\n",
                fcname,level+1,iter+1,
                ndfm,nvox3d1,(int)floor(100.0*ndfm/nvox3d1+0.5),nskip,Cprev);

      if(verbose) {
        Ccurr = niik_image_nregister_bspline_obj_func(refimg,refseg,movimg,movseg,
                warplist[0],
                sobelimg[0],sobelimg[1],
                -10,
                grad_weight,
                0,0,0,refimg->nx-1,refimg->ny-1,refimg->nz-1);
        fprintf(stdout,"[%s] level %i iter %-2i (%5.2f) ndfm:%i/%i(%2i%%) skip:%i C=%6.4f->%6.4f\n",
                fcname,level+1,iter+1,dfm,
                ndfm,nvox3d1,(int)floor(100.0*ndfm/nvox3d1+0.5),nskip,Cprev,Ccurr);
        Cprev = Ccurr;
      }

      if(0) {
        fprintf(stdout,"[niik_image_nregister_bspline] writing temp files\n");
        tmpimglist[2] = niik_image_copy(refimg);
        fprintf(stdout,"[niik_image_nregister_bspline]   warp image using bspline coefficients\n");
        if(!niik_image_apply_3d_warp_update(movimg,tmpimglist[2],warplist[0],NIIK_WARP_MAP_DISP_BSPLINE,NIIK_INTERP_LINEAR)) {
          fprintf(stderr,"ERROR: niik_image_nregister_bspline_warp_image_update\n");
          return 0;
        }
        tmpimglist[1] = movimg;
        tmpimglist[3] = refimg;
        fprintf(stdout,"[niik_image_nregister_bspline]   combine images\n");
        tmpimglist[0] = niik_image_combine(tmpimglist+1,3,4,140);
        niik_image_type_convert(tmpimglist[0],NIFTI_TYPE_UINT8);
        tmpimglist[0]->cal_min=0;
        tmpimglist[0]->cal_max=200;
        fprintf(stdout,"[niik_image_nregister_bspline]   write tmp_nreg_out.nii.gz\n");
        niik_image_write("tmp_nreg_out.nii.gz",tmpimglist[0]);
        fprintf(stdout,"[niik_image_nregister_bspline]   write tmp_nreg_warp_coeff.nii.gz\n");
        niik_image_write("tmp_nreg_warp_coeff.nii.gz",warplist[0]);
        tmpimglist[2] = niik_image_free(tmpimglist[2]);
        tmpimglist[0] = niik_image_free(tmpimglist[0]);
      }

      /* updating wbsimg (output) */
      if(!niik_image_copy_data(warplist[0],wbsimg)) {
        fprintf(stderr,"ERROR: niik_image_copy_data\n");
        return 0;
      }
    } /* iter */

    /* free memory */
    for(n=0; n<ndir; n++) {
      warplist[n] = niik_image_free(warplist[n]);
    }
    num_mask_vec = niikvec_free(num_mask_vec);
    crefimg = niik_image_free(crefimg);
    crefseg = niik_image_free(crefseg);
    cmovimg = niik_image_free(cmovimg);

    if(verbose>=0) {
      fprintf(stdout,"[niik_image_nregister_bspline] writing temp files\n");
      tmpimglist[2] = niik_image_copy(movimg);
      fprintf(stdout,"[niik_image_nregister_bspline]   warp image using bspline coefficients\n");
      if(!niik_image_apply_3d_warp_update(tmpimglist[2],refimg,wbsimg,NIIK_WARP_MAP_DISP_BSPLINE,NIIK_INTERP_LINEAR)) {
        fprintf(stderr,"ERROR: niik_image_nregister_bspline_warp_image_update\n");
        return 0;
      }
      tmpimglist[1] = movimg;
      tmpimglist[3] = refimg;
      fprintf(stdout,"[niik_image_nregister_bspline]   combine images\n");
      tmpimglist[0] = niik_image_combine(tmpimglist+1,3,4,140);
      niik_image_type_convert(tmpimglist[0],NIFTI_TYPE_UINT8);
      tmpimglist[0]->cal_min=0;
      tmpimglist[0]->cal_max=200;
      fprintf(stdout,"[niik_image_nregister_bspline]   write tmp_nreg_out.nii.gz\n");
      niik_image_write("tmp_nreg_out.nii.gz",tmpimglist[0]);
      fprintf(stdout,"[niik_image_nregister_bspline]   write tmp_nreg_warp_coeff.nii.gz\n");
      niik_image_write("tmp_nreg_warp_coeff.nii.gz",wbsimg);
      tmpimglist[2] = niik_image_free(tmpimglist[2]);
      tmpimglist[0] = niik_image_free(tmpimglist[0]);
    } /* verbose */

  } /* each level */

  /* correct warp */
  if(!niik_image_nregister_bspline_update_warp_coeff(wbsimg,refimg,wdelta[nlevel-1],reg_dist)) {
    fprintf(stderr,"ERROR: niik_image_nregister_update_warp_coeff\n");
    return 0;
  }

  for(n=0; n<3; n++) {
    if(sobelimg[n]==NULL) continue;
    sobelimg[n] = niik_image_free(sobelimg[n]);
  }
  free(sobelimg);
  free(warplist);

  return 1;
} /* niik_image_nregister_bspline */



int niik_image_nregister_bspline_update_warp_coeff(nifti_image *img,nifti_image *refimg,double delta,double regularization_size)
/* update warp coefficient
 * -img is the warp coefficient image
 * -refimg is used for image field of view
 * -delta is the new output image pixel size
 *
 * 1. resamples image to a reference image
 * 2. calculates the actual warp from bspline coefficients for the resampled image
 * 3. iteratively blurs the image and check for negative Jacobian
 *    -> corrects it by blurring the image
 * 4. resample to delta (down-sampling)
 * 5. calculates the new b-spline coefficients
 */
{
  nifti_image
  *tmpimg;
  int
  verbose=1,
  area,size,
  num,
  iter,
  nidx[9],
  n,i,j,k;
  niikpt
  plim,
  p;
  float
  *wx,*wy,*wz,
  *fx,*fy,*fz;
  char fcname[64]="niik_image_nregister_bspline_update_warp_coeff";
  niikmat
  *mk;

  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff] start\n");
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return 0;
  }
  if(img->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"ERROR: img is not double float32\n");
    return 0;
  }
  /* create a warp image */
  if(img->nu!=3 || img->ndim!=5) {
    if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff] new image\n");
    img->ndim=img->dim[0]=5;
    img->dim[5]=img->nu=3;
    img->dim[4]=img->nt=1;
    img->pixdim[5]=img->du=1;
    img->pixdim[4]=img->dt=1;
    free(img->data);
    img->nvox = img->nx*img->ny*img->nz*img->nu;
    img->data=(void *)calloc(img->nvox,sizeof(float));
    img->datatype=NIFTI_TYPE_FLOAT32;
    nifti_datatype_sizes(img->datatype,&(img->nbyper),&(img->swapsize));
    if(refimg==NULL) {
      if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff] resample %9.3f\n",delta);
      if(!niik_image_resample_3d_update(img,delta,delta,delta, -1,-1,-1,NIIK_INTERP_NN)) {
        fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
        return 0;
      }
    }
    return 1;
  } /* finish with new image */

  /* updating */
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff] updateing image\n");
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  if(img->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: img is not float32\n",fcname);
    return 0;
  }
  if(refimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: refimg is not float32\n",fcname);
    return 0;
  }

  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   resample\n");
  if(!niik_image_resample_3d_update(img,
                                    refimg->dx,refimg->dy,refimg->dz,
                                    refimg->nx,refimg->ny,refimg->nz,
                                    NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   %3i %3i %3i \n",img->nx,img->ny,img->nz);
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   warp calculation\n");
  size = img->nx*img->ny*img->nz;
  mk = niikmat_init(img->nz,4);

  #pragma omp parallel for private(n,p,i,j)
  for(k=0; k<img->nz; k++) {
    p.z = k*img->dz / tmpimg->dz;
    for(j=0; j<img->ny; j++) {
      p.y = j*img->dy / tmpimg->dy;
      n=j*img->nx+k*img->nx*img->ny;
      for(i=0; i<img->nx; n++,i++) {
        p.x = i*img->dx / tmpimg->dx;
        niik_image_interpolate_3d_bspline_ijk_update(tmpimg,p,mk->m[k]);
        niik_image_set_voxel(img,n,mk->m[k][0]);
        niik_image_set_voxel(img,n+size,mk->m[k][1]);
        niik_image_set_voxel(img,n+size*2,mk->m[k][2]);
      }
    }
  }
  tmpimg = niik_image_free(tmpimg);
  mk = niikmat_free(mk);
  if(verbose)   fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   iterative correction\n");
  if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]     temp image\n");
  if((tmpimg = niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
    return 0;
  }
  area = img->nx*img->ny;
  size = img->nz*area;
  fx = tmpimg->data;
  fy = fx + size;
  fz = fy + size;
  wx = img->data;
  wy = wx + size;
  wz = wy + size;
  /* limit for Jacobian */
  plim.x = -(1.0-regularization_size)*refimg->dx;
  plim.y = -(1.0-regularization_size)*refimg->dy;
  plim.z = -(1.0-regularization_size)*refimg->dz;
  for(iter=0; iter<2500; iter++) {
    /*if(!niik_image_filter_gaussian_update(img,3,3)){
      fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
      return 0; } */
    for(k=num=0; k<img->nz; k++) {
      n=k*img->nx*img->ny;
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; n++,i++) {
          nidx[0]=(i>0        )?n-1:n;
          nidx[1]=(i<img->nx-1)?n+1:n;
          nidx[2]=(j>0        )?n-img->nx:n;
          nidx[3]=(j<img->ny-1)?n+img->nx:n;
          nidx[4]=(k>0        )?n-area:n;
          nidx[5]=(k<img->nz-1)?n+area:n;
          if(fx[n]-fx[nidx[0]]<=plim.x) {
            num++;
            /*wx[n]=(fx[nidx[0]]+fx[nidx[1]]+fx[nidx[2]]+fx[nidx[3]]+fx[nidx[4]]+fx[nidx[5]])/6.0;*/
            wx[n]=(fx[nidx[0]]+fx[nidx[1]])/2.0;
            if(verbose>1)
              fprintf(stdout,"-x [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n),
                      fx[n-1],fx[n+1],fx[n-img->nx],fx[n+img->nx],fx[n-area],fx[n+area],fx[n]);
          } else if(fy[n]-fy[nidx[2]]<=plim.y) {
            num++;
            /*wy[n]=(fy[nidx[0]]+fy[nidx[1]]+fy[nidx[2]]+fy[nidx[3]]+fy[nidx[4]]+fy[nidx[5]])/6.0;*/
            wy[n]=(fy[nidx[2]]+fy[nidx[3]])/2.0;
            if(verbose>1)
              fprintf(stdout,"-y [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n+size),
                      fy[n-1],fy[n+1],fy[n-img->nx],fy[n+img->nx],fy[n-area],fy[n+area],fy[n]);
            continue;
          } else if(fz[n]-fz[nidx[4]]<=plim.z) {
            num++;
            /*wz[n]=(fz[nidx[0]]+fz[nidx[1]]+fz[nidx[2]]+fz[nidx[3]]+fz[nidx[4]]+fz[nidx[5]])/6.0; */
            wz[n]=(fz[nidx[4]]+fz[nidx[5]])/2.0;
            if(verbose>1)
              fprintf(stdout,"-z [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n+size*2),
                      fz[n-1],fz[n+1],fz[n-img->nx],fz[n+img->nx],fz[n-area],fz[n+area],fz[n]);
            continue;
          } else if(fx[nidx[1]]-fx[n]<=plim.x) {
            num++;
            /*wx[n]=(fx[nidx[0]]+fx[nidx[1]]+fx[nidx[2]]+fx[nidx[3]]+fx[nidx[4]]+fx[nidx[5]])/6.0;*/
            wx[n]=(fx[nidx[0]]+fx[nidx[1]])/2.0;
            if(verbose>1)
              fprintf(stdout,"+x [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n),
                      fx[n-1],fx[n+1],fx[n-img->nx],fx[n+img->nx],fx[n-area],fx[n+area],fx[n]);
            continue;
          } else if(fy[nidx[3]]-fy[n]<=plim.y) {
            num++;
            /*wy[n]=(fy[nidx[0]]+fy[nidx[1]]+fy[nidx[2]]+fy[nidx[3]]+fy[nidx[4]]+fy[nidx[5]])/6.0;*/
            wy[n]=(fy[nidx[2]]+fy[nidx[3]])/2.0;
            if(verbose>1)
              fprintf(stdout,"+y [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n+size),
                      fy[n-1],fy[n+1],fy[n-img->nx],fy[n+img->nx],fy[n-area],fy[n+area],fy[n]);
            continue;
          } else if(fz[nidx[5]]-fz[n]<=plim.z) {
            num++;
            /*wz[n]=(fz[nidx[0]]+fz[nidx[1]]+fz[nidx[2]]+fz[nidx[3]]+fz[nidx[4]]+fz[nidx[5]])/6.0; */
            wz[n]=(fz[nidx[4]]+fz[nidx[5]])/2.0;
            if(verbose>1)
              fprintf(stdout,"+z [%3i %3i %3i] %7.2f <- %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f | %7.2f\n",i,j,k,niik_image_get_voxel(img,n+size*2),
                      fz[n-1],fz[n+1],fz[n-img->nx],fz[n+img->nx],fz[n-area],fz[n+area],fz[n]);
            continue;
          }
        }
      }
    }
    niik_image_copy_data(img,tmpimg);
    if(verbose && !iter) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   iteration %-3i  correct %i / %i\n",iter+1,num,size);
    if(!num) break;
    if(verbose &&  iter) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   iteration %-3i  correct %i / %i\n",iter+1,num,size);
  }
  tmpimg = niik_image_free(tmpimg);
  if(verbose>1) {
    fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   write tmp_nregimg_warp.nii.gz\n");
    niik_image_write("tmp_nregimg_warp.nii.gz",img);
  }
  if(!niik_image_resample_3d_update(img, delta,delta,delta,
                                    refimg->dx*refimg->nx/delta+1,
                                    refimg->dy*refimg->ny/delta+1,
                                    refimg->dz*refimg->nz/delta+1,
                                    NIIK_INTERP_LINEAR)) {
    fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   write tmp_nregimg_warpn.nii.gz\n");
    niik_image_write("tmp_nregimg_warpn.nii.gz",img);
  }
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   bspline coeff\n");
  if(!niik_image_interpolate_convert_3d_bspline_coeff(img)) {
    fprintf(stderr,"ERROR: niik_image_interpoblate_convert_3d_bspline_coeff\n");
    return 0;
  }
  /* convert to float32 */
  if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"[niik_image_nregister_bspline_update_warp_coeff]   write tmp_nregimg_warpnc.nii.gz\n");
    niik_image_write("tmp_nregimg_warpnc.nii.gz",img);
  }
  return 1;
} /* niik_image_nregister_bspline_update_warp_coeff */


#endif /* _FALCON_NREGISTER_BSPLINE_C_ */
