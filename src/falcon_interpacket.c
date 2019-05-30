/* Filename:     nifti1_kunio_interpacket.c
 * Description:  interpacket processing
 * Author:       Kunio Nakamura
 * Date:         November 12, 2012
 *
 * nifti_image **niik_image_split_interpacket(nifti_image *img,int num);
 * nifti_image *niik_image_merge_interpacket(nifti_image **imglist,int num);
 */

#ifndef _FALCON_INTERPACKET_C_
#define _FALCON_INTERPACKET_C_

#include "falcon.h"

nifti_image *niik_image_interpacket_motion_correction(nifti_image *img,int num,int dof) {
  char fcname[64]="niik_image_interpacket_motion_correction";
  nifti_image **imglist=NULL,*outimg=NULL;
  double *affpar=NULL;
  niikvec *iaffpar=NULL;
  niikmat *afmat=NULL;
  niikpt ctr;
  int n;
  /* start */
  niik_fc_display(fcname,1);
  if(num!=2) {
    fprintf(stderr,"[%s] ERROR: not yet implemencted for symmetric registration\n",fcname);
    return NULL;
  }
  /* slit into individual packets */
  if((imglist=niik_image_split_interpacket(img,num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_split_interpacket\n",fcname);
    return NULL;
  }
  affpar=(double *)calloc(25,sizeof(double));
  for(n=0; n<num; n++)
    niik_image_type_convert(imglist[n],NIFTI_TYPE_FLOAT32);
  /* pre-defined parameters */
  ctr=niikpt_image_get_centroid(img,NULL);
  affpar[7]=affpar[8]=affpar[9]=affpar[10]=1;
  affpar[14]=ctr.x;
  affpar[15]=ctr.y;
  affpar[16]=ctr.z;
  niik_aregister_display_affine(affpar);
  /* registration
  if(!niik_image_aregister2_test1(imglist[0],NULL,imglist[1],NULL,affpar,dof,NIIK_REGISTER_CC,NULL)){
    fprintf(stderr,"[%s] ERROR: niik_image_aregister2_test1\n",fcname);
    return NULL; }*/
  iaffpar=niikvec_init(25);
  iaffpar->v[1]=20;
  iaffpar->v[2]=20;
  iaffpar->v[3]=20;
  iaffpar->v[4]=20;
  iaffpar->v[5]=20;
  iaffpar->v[6]=20;
  if(!niik_image_aregister2_test2(imglist[0],NULL,imglist[1],NULL,affpar,NIIK_REGISTER_CC,2,3.2,100,2,2.1,iaffpar,NULL)) {
    fprintf(stderr,"[%s] ERROR: niik_image_aregister2_test2\n",fcname);
    return NULL;
  }
  iaffpar=niikvec_free(iaffpar);
  niik_aregister_display_affine(affpar);
  /*
   * half-way transfomration for symmetric interpolation
   */
  affpar[1]/=-2;
  affpar[2]/=-2;
  affpar[3]/=-2;
  affpar[4]/=-2;
  affpar[5]/=-2;
  affpar[6]=-affpar[6]/2.0+img->dz/2.0; /* to account to slice position */
  affpar[7]=(affpar[7]+1.0)/2.0;
  affpar[8]=(affpar[7]+1.0)/2.0;
  affpar[9]=(affpar[7]+1.0)/2.0;
  affpar[10]=(affpar[7]+1.0)/2.0;
  affpar[11]/=2;
  affpar[12]/=2;
  affpar[13]/=2;
  afmat=niik_aregister_matrix_from_affpar(affpar);
  niikmat_inverse_update(afmat);
  /*niikmat_inverse_update(afmat);*/
  /* image transformation */
  if(!niik_image_affine_transform_3d_update(imglist[1],imglist[0],afmat,NIIK_INTERP_BSPLINE)) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
    return NULL;
  }
  niikmat_inverse_update(afmat);
  if(!niik_image_affine_transform_3d_update(imglist[0],imglist[0],afmat,NIIK_INTERP_BSPLINE)) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
    return NULL;
  }
  /* merge packets */
  if((outimg=niik_image_merge_interpacket(imglist,num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_merge_interpacket\n",fcname);
    return NULL;
  }
  /* free memoery */
  for(n=0; n<num; n++)
    imglist[n]=niik_image_free(imglist[n]);
  free(imglist);
  free(affpar);
  afmat=niikmat_free(afmat);
  niik_fc_display(fcname,0);
  return outimg;
} /* niik_iamge_interpacket_motion_correction */


nifti_image **niik_image_split_interpacket(nifti_image *img,int num)
/* split interpackets
 * -updates dz, and s-matrix
 * -but the s-matrix start position will not be consistent
 */
{
  nifti_image **outlist=NULL;
  /* niikmat *newmat=NULL; */
  char
  fname[128],
        fcname[64]="niik_image_split_interpacket";
  int i,area,j,k,n,N[128],nn;
  niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  if(img->ndim>3) {
    fprintf(stderr,"[%s] ERROR: img is not a 3d image, %i\n",fcname,img->ndim);
    return NULL;
  }
  outlist=(nifti_image **)calloc(num,sizeof(nifti_image *));
  fprintf(stdout,"[%s] image nz = %i\n",fcname,img->nz);
  area=img->nx*img->ny;
  for(n=0; n<num; n++) {
    if((outlist[n]=niik_image_copy(img))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
      return NULL;
    }
    outlist[n]->pixdim[3] = outlist[n]->dz = outlist[n]->dz*num;
    if(img->sform_code>0) {
      if(n==0) {
        /* 2013-03-27, Kunio Nakamura, knakamura@mrs.mcgill.ca
         * Scaling by matrix multiplication results in a huge translation which
         * when combined for interpacket motion, looks very wrong in the start position for minc
         *
         NIIK_RET0n(((newmat=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix");
         NIIK_RET0n((!niikmat_multiply_mat2_free1(niikmat_scale_matrix(1,1,num),
         newmat)),fcname,"niikmat_multiply_mat2_free1");
         outlist[n]->sto_xyz=niikmat_make_mat44(newmat);
         NIIK_RET0n((!niikmat_inverse_update(newmat)),fcname,"niikmat_inverse_update");
         outlist[n]->sto_ijk=niikmat_make_mat44(newmat); */
        mat44_display(outlist[n]->sto_xyz);
        outlist[n]->sto_xyz.m[2][2]*=num;
        mat44_display(outlist[n]->sto_xyz);
        outlist[n]->sto_ijk = nifti_mat44_inverse( outlist[n]->sto_xyz ) ;
      } else {
        outlist[n]->sto_xyz=outlist[0]->sto_xyz;
        outlist[n]->sto_ijk=outlist[0]->sto_ijk;
      }
    }
    nn=img->nz/num;
    if(n<(img->nz%num)) nn++;
    fprintf(stdout,"[%s]   image %2i : nz = %-4i\n",fcname,n,nn);
    outlist[n]->nz=outlist[n]->dim[3]=nn;
    outlist[n]->nvox=outlist[n]->nx*outlist[n]->ny*outlist[n]->nz;
    free(outlist[n]->data);
    outlist[n]->data=(void *)calloc(outlist[n]->nvox,outlist[n]->nbyper);
    N[n]=0;
  }
  for(k=j=0; k<img->nz; k++) {
    n=k%num;
    for(i=0; i<area; i++) {
      niik_image_set_voxel(outlist[n],N[n]++,niik_image_get_voxel(img,j++));
    }
  }
  if(0) {
    for(n=0; n<num; n++) {
      sprintf(fname,"test_img%i.nii.gz",n+1);
      niik_image_write(fname,outlist[n]);
    }
  } /* save images */
  niik_fc_display(fcname,0);
  mat44_display(outlist[n]->sto_xyz);
  return outlist;
} /* niik_image_split_interpacket */


nifti_image *niik_image_merge_interpacket(nifti_image **imglist,int num) {
  nifti_image *outimg=NULL;
  char fcname[128]="niik_image_merge_interpacket";
  int n,i,N[128],j,k,area;
  niik_fc_display(fcname,1);
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return NULL;
  }
  for(n=0; n<num; n++) {
    if(imglist[n]->ndim>3) {
      fprintf(stderr,"[%s] ERROR: imglist[%i] is not 3d, %i\n",fcname,n,imglist[n]->ndim);
      return NULL;
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
    return NULL;
  }
  for(n=outimg->nz=0; n<num; n++) {
    outimg->nz+=imglist[n]->nz;
  }
  outimg->dz/=num;
  free(outimg->data);
  outimg->dim[3]=outimg->nz;
  outimg->nvox=outimg->nx*outimg->ny*outimg->nz;
  outimg->data=(void *)calloc(outimg->nvox,outimg->nbyper);
  area=outimg->nx*outimg->ny;
  for(n=0; n<num; n++) N[n]=0;
  for(k=j=0; k<outimg->nz; k++) {
    n=k%num;
    /*fprintf(stdout,"%i %i \n",k,n);*/
    for(i=0; i<area; i++) {
      niik_image_set_voxel(outimg,j++,niik_image_get_voxel(imglist[n],N[n]++));
    }
    /*fprintf(stdout,"%i %i %i\n",k,n,N[n]);*/
  }
  niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_merge_interpacket */


#endif /* _FALCON_INTERPACKET_C_ */
