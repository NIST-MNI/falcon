/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)


void usage() {
  fprintf(stdout,"niik_stx_registration\n");
  fprintf(stdout,"  usage: img.nii out.nii\n\n");
}


int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *out=NULL,
   *img=NULL;
  int i,j,k,m,n;
  char fcname[32]="test";
  double v;

  if(argc==1) {
    usage();
    exit(0);
  }

  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[1]);
  NIIK_RET0(((img=niik_image_read(argv[1]))==NULL),fcname,"reading pdw, niik_image_read");
  out=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32);

  for(k=n=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(m=0; m<3; n++,m++) {
          v=niik_image_get_voxel(img,n);
          niik_image_set_voxel(out,i+j*img->nx+k*img->nx*img->ny+m*img->nx*img->ny*img->nz,v);
        }
      }
    }
  }

  fprintf(stdout,"[%s] writing image:        %s\n",fcname,argv[2]);
  NIIK_RET0((!niik_image_write(argv[2],out)),fcname,"niik_image_write output");

  img=niik_image_free(img);
  out=niik_image_free(out);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_calc_T2 */
