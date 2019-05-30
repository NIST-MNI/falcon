/* Filename:     nifti1_kunio_siena.c
 * Description:  niik implementaiton of SIENA
 * Author:       Kunio Nakamura
 * Date:         October 5, 2012
 */

#ifndef _FALCON_SIENA_C_
#define _FALCON_SIENA_C_

#include "falcon.h"

int niik_image_siena(nifti_image *templateimg,nifti_image *maskimg,nifti_image **imglist,niikmat **matrixlist,int numimglist,
                     double ulim,double uran,int coloroutput,niikmat *out)
/* niik implementation of siena
 * -for now only 2 time-points
 * -templateimg is the subject-specific average template image
 * -maskimg is the binary mask of brain edges (usually inside brain), in the space of templateimg
 * -imglist is the list of intensity-corrected & normalized images
 * -matrixlist is the list of affine matrices that transfrom images in imglist to template space
 * -numimglist is the number of imglist and matrixlist
 * -ulim and uran are the upper intensity limit and range to remove the effect of white matter edge in high-contrast images
 * --probably calculated by bimodal estimation.
 * -coloroutput is a flag for writing color-coded output
 * -out is the output numimglist-by-#mask_voxels-matrix containing the edge shift
 */
{
  int
  i,j,k,m,n,t,ti,nn,
  verbose=2;
  char fcname[64]="niik_image_siena";
  unsigned char *bimg;
  double
  xcthresh=1000,
  xmin=-5,
  xmax=5,
  dx=0.2,
  xi,yi,
  plim=3.0;
  niikpt pt,tt,gr,*norm=NULL,*qt=NULL;
  niikmat
  *A_matrix_for_mode=NULL,
   **invmat=NULL;
  niikvec **y,**d,**xc,*x,*x2,*p;

  niik_fc_display(fcname,1);

  if(numimglist!=2) {
    fprintf(stderr,"[%s] ERROR: asssuming 2 images for now\n",fcname);
    return 0;
  }
  if(templateimg==NULL) {
    fprintf(stderr,"[%s] ERROR: template image is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return 0;
  }
  if(matrixlist==NULL) {
    fprintf(stderr,"[%s] ERROR: matrixlist is null\n",fcname);
    return 0;
  }

  if(niik_image_cmp_dim(templateimg,maskimg)!=0) {
    fprintf(stderr,"[%s] ERROR: templateimg and maskimg\n",fcname);
    return 0;
  }

  if(verbose>1) fprintf(stdout,"[%s] checking datatype\n",fcname);
  if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"[%s] ERROR: maskimg is not uint8, %s\n",fcname,nifti_datatype_string(maskimg->datatype));
    return 0;
  }

  if(verbose>1) fprintf(stdout,"[%s] prepare point\n",fcname);
  if((qt=(niikpt *)calloc(2*numimglist,sizeof(niikpt)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc qt\n",fcname);
    return 0;
  }
  norm=qt+numimglist;

  if(verbose>1) fprintf(stdout,"[%s] inverse matrix\n",fcname);
  invmat=(niikmat **)calloc(numimglist,sizeof(niikmat *));
  for(m=0; m<numimglist; m++) {
    if((invmat[m]=niikmat_inverse(matrixlist[m]))==NULL) {
      fprintf(stderr,"[%s] ERROR: niikmat_inverse %i\n",fcname,m);
      return 0;
    }
  }

  x=niikvec_init_range(xmin,xmax,dx);
  x2=niikvec_init_range(xmin*2.0-dx,xmax*2.0+dx,dx);
  p=niikvec_init(x2->num);
  y=(niikvec **)calloc(numimglist,sizeof(niikvec *));
  d=(niikvec **)calloc(numimglist,sizeof(niikvec *));
  xc=(niikvec **)calloc(numimglist,sizeof(niikvec *));
  A_matrix_for_mode = niikmat_init(x2->num,x2->num);
  niik_bspline_update_A(A_matrix_for_mode);
  for(m=0; m<numimglist; m++) {
    y[m]=niikvec_init(x->num);
    d[m]=niikvec_init(x->num);
    xc[m]=niikvec_init(x->num*2+1);
  }
  for(n=0; n<x2->num; n++) {
    p->v[n]=1.0-NIIK_Heaviside(fabs(x2->v[n])-plim,1.25);
  }

  bimg=maskimg->data;
  if(verbose>1) fprintf(stdout,"[%s] main loop\n",fcname);
  for(n=nn=0; n<templateimg->nvox; n++) {
    if(bimg[n]==0) continue;
    pt=niik_image_get_pt_from_index(templateimg,n,&i,&j,&k);
    gr=niikpt_unit(niik_image_sobel_filter_niikpt(templateimg,pt));
    tt=niikpt_add(pt,gr);
    /*if(nn>=out->row) {
      fprintf(stdout,"[%s] reached the last mask, %3i %3i %3i\n",fcname,i,j,k);
      break; }*/

    /*fprintf(stdout,"\t[%3i %3i %3i] %9.3f %9.3f %9.3f %5.0f\n",i,j,k,gr.x,gr.y,gr.z,gr.w);*/
    /*if(i==54 && j==120 && k==117){ verbose=10; }*/

    for(m=0; m<numimglist; m++) {
      qt[m]=niikpt_affine_transform(invmat[m],pt);
      norm[m]=niikpt_unit(niikpt_sub(niikpt_affine_transform(invmat[m],tt),qt[m]));
      /* fprintf(stdout,"\t  %9.3f %9.3f %9.3f | %9.3f %9.3f %9.3f\n",qt[m].x,qt[m].y,qt[m].z,norm[m].x,norm[m].y,norm[m].z); */

      /* get intensity */
      if(!niik_image_interp_along_normal(imglist[m],NIIK_INTERP_LINEAR,qt[m],norm[m],x,y[m])) {
        fprintf(stderr,"[%s] niik_image_interp_along_normal [%i,%i,%i] %i\n",fcname,i,j,k,m);
        return 0;
      }

      /* central difference */
      if(!niikvec_copy_update(y[m],d[m])) {
        fprintf(stderr,"[%s] niikvec_copy_update for first derivative %i\n",fcname,m);
        return 0;
      }
      niik_central_difference_double_vector(d[m]->v,y[m]->num);

      /* modify by intensity (remove white matter edge) */
      for(t=0; t<x->num; t++) {
        d[m]->v[t] *= (1.0-NIIK_Heaviside(y[m]->v[t]-ulim,uran)) / dx;
      }

      /* modify by first derivative (remove the other side of sulci) */
      for(t=x->num/2; t>0; t--) {
        if(x->v[t]>0) continue;
        if(y[m]->v[t-1]>y[m]->v[t])
          d[m]->v[t-1]=0;
      }

      /* calculate cross-correlation */
      if(!niikvec_cross_correlation(d[0],d[m],xc[m])) {
        fprintf(stderr,"[%s] ERROR: niikvec_cross_correlation\n",fcname);
        return 0;
      }

      /* calculate cross-correlation */
      for(t=ti=0; t<x2->num; t++) {
        if(x2->v[t]<-5.0) continue;
        if(x2->v[t]> 5.0) break;
        if(xc[m]->v[ti]<xc[m]->v[t]) ti=t;
      }

      /* get the realistic maximum */
      if(ti==0)
        out->m[m][nn] = 0;
      else if(xc[m]->v[ti]<xcthresh) {
        out->m[m][nn] = 0;
      } else {
        if(!niik_get_mode_bspline_vector_with_A_matrix(xc[m]->v,xc[m]->num,&xi,&yi,A_matrix_for_mode)) {
          fprintf(stderr,"[%s] ERROR: niik_get_mode_bspline_vector_with_A_matrix\n",fcname);
          return 0;
        }
        out->m[m][nn] = x2->v[ti];
        out->m[m][nn] = (xi*dx)+(x2->v[0]);
      }

      if(verbose>2)
        fprintf(stdout,"[%3i,%3i,%3i] %9i  %9.3f %9.3f %9.3f | %9.3f %9.3f %9.3f | %12.7f %4i\n",i,j,k,nn,
                qt[m].x,qt[m].y,qt[m].z,norm[m].x,norm[m].y,norm[m].z,out->m[m][nn],ti);
      if(verbose==10) {
        fprintf(stdout,"  y%i=[",m);
        niikvec_display(y[m]);
        fprintf(stdout,"];\n");
        fprintf(stdout,"  dy%i=[",m);
        niikvec_display(d[m]);
        fprintf(stdout,"];\n");
        fprintf(stdout,"  xc%i=[",m);
        niikvec_display(xc[m]);
        fprintf(stdout,"];\n");
      }
    } /* each time-point, m */
    nn++;

    if(verbose==10) {
      fprintf(stdout,"  x=[");
      niikvec_display(x);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  x2=[");
      niikvec_display(x2);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  p=[");
      niikvec_display(p);
      fprintf(stdout,"];\n");
      exit(0);
    }

  } /* each voxel */

  /* COLOR CODE ATROPHY RATES */
  if(coloroutput) {

  } /* color-coded output */

  free(qt);
  for(m=0; m<numimglist; m++) {
    invmat[m]=niikmat_free(invmat[m]);
  }
  free(invmat);
  for(m=0; m<numimglist; m++) {
    y[m]=niikvec_free(y[m]);
    d[m]=niikvec_free(d[m]);
    xc[m]=niikvec_free(xc[m]);
  }
  p=niikvec_free(p);
  x=niikvec_free(x);
  x2=niikvec_free(x2);
  free(y);
  free(d);
  free(xc);
  A_matrix_for_mode=niikmat_free(A_matrix_for_mode);

  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_siena */


nifti_image *niik_image_siena_color(nifti_image *img,nifti_image *flow,double lim,double imax) {
  char fcname[32]="niik_image_siena_color";
  int
  i,j,
  num,num2,
  verbose=2;
  niikmat *cm=NULL;
  double
  lim2;  /* assume the range is -2 to 2 */
  nifti_image
  **tmpimg=NULL,
    *outimg=NULL;
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((flow==NULL),fcname,"flow is null");
  NIIK_RET0((img->nvox!=flow->nvox),fcname,"different nvox");
  if(imax<0) {
    imax=niik_image_get_percentile(img,NULL,0.99);
  }
  fprintf(stdout,"[%s] color max: %9.3f\n",fcname,imax);
  num=51;
  NIIK_RET0(((cm=niik_colormap_get(NIIK_COLORMAP_ATROPHY,num))==NULL),fcname,"niik_colormap_get");
  num2=num/2;
  tmpimg=(nifti_image **)calloc(3,sizeof(nifti_image *));
  for(i=0; i<3; i++) {
    NIIK_RET0(((tmpimg[i]=niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  }
  NIIK_RET0(((outimg=niik_image_combine(tmpimg,3,6,0))==NULL),fcname,"niik_image_combine");
  for(i=0; i<3; i++) {
    tmpimg[i]=niik_image_free(tmpimg[i]);
  }
  free(tmpimg);
  lim2=lim*2.0;
  if(verbose>0) {
    fprintf(stdout,"[%s] value range: %9.4f %9.4f\n",fcname,-lim,lim);
  }

  NIIK_RET0((!niikmat_kmul(cm,imax)),fcname,"niikmat_kmul");
  /* niikmat_display(cm); */
  for(i=0; i<img->nvox; i++) {
    j=NIIK_IMINMAX((niik_image_get_voxel(flow,i)+lim)/lim2*num,0,num-1);
    if(j-num2>-2) {
      if(j-num2<2) {
        continue;
      }
    }
    niik_image_set_voxel(outimg,i,cm->m[j][0]);
    niik_image_set_voxel(outimg,i+img->nvox,cm->m[j][1]);
    niik_image_set_voxel(outimg,i+img->nvox*2,cm->m[j][2]);
  } /* each voxel */
  cm=niikmat_free(cm);
  return outimg;
} /* niik_image_siena_color */


#endif /* _FALCON_SIENA_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/