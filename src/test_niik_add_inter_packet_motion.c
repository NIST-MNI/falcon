/* Filename:     test_niik_add_inter_packet_motion.c
 * Description:  to study the effect of inter-packet motion
 * Author:       Kunio Nakamura
 * Date:         November 12, 2012
 */

#include "falcon.h"

int main(int argc,char *argv[],char *envp[]) {
  char fcname[256]="test_niik_add_inter_packet_motion";
  nifti_image *img=NULL,*maskimg=NULL,**tmpimglist=NULL,*outimg=NULL;
  niikmat *afmat=NULL;
  double
  maxR=1.0,maxT=1.0,
  *affpar=NULL;
  niikpt ctr;
  int n,num=2;

  if(argc!=4) {
    fprintf(stdout,"  usage: <input.nii> <mask.nii> <out.nii>\n");
    exit(0);
  }

  srand(time(NULL));

  fprintf(stdout,"  read image  %s\n",argv[1]);
  if((img=nifti_image_read(argv[1],1))==NULL) {
    fprintf(stderr,"[%s] ERROR: reading %s\n",fcname,argv[1]);
    exit(0);
  }

  fprintf(stdout,"  read image  %s\n",argv[2]);
  if((maskimg=nifti_image_read(argv[2],1))==NULL) {
    fprintf(stderr,"[%s] ERROR: reading %s\n",fcname,argv[2]);
    exit(0);
  }

  affpar=(double *)calloc(20,sizeof(double));
  affpar[7]=affpar[8]=affpar[9]=affpar[10]=1.0;

  ctr = niikpt_image_get_centroid(img,maskimg);
  fprintf(stdout,"[%s] center  %9.4f %9.4f %9.4f\n",fcname,ctr.x,ctr.y,ctr.z);
  affpar[14]=ctr.x;
  affpar[15]=ctr.y;
  affpar[16]=ctr.z;

  affpar[1]=maxR*(rand()%3600)/3600.0 - maxR/2.0;
  affpar[2]=maxR*(rand()%3600)/3600.0 - maxR/2.0;
  affpar[3]=maxR*(rand()%3600)/3600.0 - maxR/2.0;

  affpar[4]=maxT*(rand()%3600)/3600.0 - maxT/2.0;
  affpar[5]=maxT*(rand()%3600)/3600.0 - maxT/2.0;
  affpar[6]=maxT*(rand()%3600)/3600.0 - maxT/2.0;

  niik_aregister_display_affine(affpar);
  if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
    fprintf(stderr,"[%s] niik_aregister_matrix_from_affpar\n",fcname);
    exit(0);
  }

  /* niikmat_identity_update(afmat); */

  fprintf(stdout,"[%s] split images\n",fcname);
  if((tmpimglist=niik_image_split_interpacket(img,num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_split_interpacket\n",fcname);
    exit(0);
  }

  n=0;
  fprintf(stdout,"[%s] transform image %i\n",fcname,n);
  if(!niik_image_affine_transform_3d_update(tmpimglist[n],tmpimglist[n],afmat,NIIK_INTERP_BSPLINE)) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
    exit(0);
  }

  fprintf(stdout,"[%s] inverse matrix\n",fcname);
  niikmat_inverse_update(afmat);

  n=1;
  fprintf(stdout,"[%s] transform image %i\n",fcname,n);
  if(!niik_image_affine_transform_3d_update(tmpimglist[n],tmpimglist[n],afmat,NIIK_INTERP_BSPLINE)) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
    exit(0);
  }

  fprintf(stdout,"[%s] merge images\n",fcname);
  if((outimg=niik_image_merge_interpacket(tmpimglist,num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_merge_interpacket\n",fcname);
    exit(0);
  }

  fprintf(stdout,"  writing %s\n",argv[3]);
  if(!niik_image_write(argv[3],outimg)) {
    fprintf(stderr,"[%s] ERROR: writing %s\n",fcname,argv[3]);
    exit(0);
  }

  img=niik_image_free(img);
  outimg=niik_image_free(outimg);
  maskimg=niik_image_free(maskimg);
  exit(0);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/