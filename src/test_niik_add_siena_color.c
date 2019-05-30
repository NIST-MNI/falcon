/* Filename:
 * Description:
 * Author:       Kunio Nakamura
 * Date:         November 24, 2012
 */

#include "falcon.h"

int main(int argc,char *argv[],char *envp[]) {
  char fcname[256]="test_niik_add_siena_color";
  nifti_image *img=NULL,*maskimg=NULL;
  niikmat *cmap=NULL;
  double
  dval,cmax,
       dc,dmax=2.5;
  int n,m,num=501;

  if(argc!=5) {
    fprintf(stdout,"  usage: <img.nii> <siena.nii> <Max> <out.nii>\n");
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

  dmax = atof(argv[3]);
  cmap=niik_colormap_get(NIIK_COLORMAP_ATROPHY,num);
  if(!niik_image_convert_to_color_image(img)) {
    fprintf(stderr,"[%s] ERROR: niik_image_convert_to_color_image\n",fcname);
    exit(0);
  }

  cmax = niik_image_get_percentile(img,NULL,0.9);
  fprintf(stdout,"  color max  %.3f\n",cmax);
  for(n=0; n<num; n++) {
    cmap->m[n][0]*=cmax;
    cmap->m[n][1]*=cmax;
    cmap->m[n][2]*=cmax;
  }

  dc=2.0*dmax/(num-1.0);
  fprintf(stdout,"  %.4f %.3f %i\n",dc,dmax,num);
  for(n=0; n<maskimg->nvox; n++) {
    dval = niik_image_get_voxel(maskimg,n);
    if(fabs(dval) < 1e-3) continue;
    m = (dmax + niik_image_get_voxel(maskimg,n)) / dc;
    if(m<0) m=0;
    else if(m>=num) m=num-1;
    /* fprintf(stdout,"%i %i %f -> %f %f %f\n",n,m,dval,cmap->m[m][0],cmap->m[m][1],cmap->m[m][2]); */
    niik_image_set_voxel(img,n,cmap->m[m][0]);
    niik_image_set_voxel(img,n+maskimg->nvox,cmap->m[m][1]);
    niik_image_set_voxel(img,n+2*maskimg->nvox,cmap->m[m][2]);
  }

  niik_image_type_convert(img,NIFTI_TYPE_UINT16);
  fprintf(stdout,"  writing %s\n",argv[4]);
  if(!niik_image_write(argv[4],img)) {
    fprintf(stderr,"[%s] ERROR: writing %s\n",fcname,argv[3]);
    exit(0);
  }

  cmap=niikmat_free(cmap);
  img=niik_image_free(img);
  maskimg=niik_image_free(maskimg);
  exit(0);
}
