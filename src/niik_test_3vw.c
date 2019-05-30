/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_3vw.h"

int niikmat_convert_from_xfm_tgtmatrix(niikmat *m,nifti_image *srcimg,niikmat *tgt_sto_ijk) {
  char fcname[64]="niikmat_convert_from_xfm_tgtmatrix";
  niikmat *sm=NULL;
  NIIK_RET0((m==NULL),fcname,"m is null");
  NIIK_RET0((srcimg==NULL),fcname,"srcimg is null");
  NIIK_RET0((tgt_sto_ijk==NULL),fcname,"tgtimg is null");
  NIIK_RET0(((sm = niikmat_mat44_matrix(srcimg->sto_xyz))==NULL),
            fcname,"niikkmat_mat44_matrix");
  NIIK_RET0((!niikmat_multiply_mat1_free2
             (sm,
              niikmat_scale_matrix(1.0/srcimg->dx,1.0/srcimg->dy,1.0/srcimg->dz))),fcname,"niikmat_multiply_mat1_free");
  NIIK_RET0((!niikmat_multiply_mat1_free2(m,sm)),fcname,"niikmat_multiply_mat1_free2");
  NIIK_RET0((!niikmat_multiply_mat2(tgt_sto_ijk,m)),fcname,"niikmat_multiply_mat2");
  return 1;
} /* niikmat_convert_from_xfm_tgtmatrix */

void usage() {
  fprintf(stderr,"  usage: in.mnc stx.xfm\n");
}

int main(int argc,char *argv[],char *envp[]) {
  niik_image_3vw *v;
  int nc=1,sc=1;

  v=niik_image_3vw_init();
  v->verbose=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--help",6)) {
        usage();
        exit(1);
      } else if(!strncmp(argv[nc],"-v",2)) {
        v->verbose=2;
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",__func__,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc!=3) {
    fprintf(stderr,"  usage: in.mnc stx.xfm\n");
    exit(1);
  }

  fprintf(stdout,"[%s] reading image:        %s\n",__func__,argv[1]);
  NIIK_EXIT(((v->img=niik_image_read(argv[1]))==NULL),__func__,"reading mask, niik_image_read",1);
  NIIK_EXIT(((v->stx=niikmat_read_xfm(argv[2]))==NULL),__func__,"reading matrix",1);
  niikmat_convert_from_xfm_tgtmatrix(v->stx,v->img,niikmat_translate_matrix(96,132,78));
  v->ctr=niikpt_val(95,116,75,1);
  niik_image_3vw_update(v);

  fprintf(stdout,"[%s] width  %9.4f\n",__func__,v->width);

  v=niik_image_3vw_free(v);

  exit(0);
} // main
