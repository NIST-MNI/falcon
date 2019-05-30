/* Filename:     xfm2niikmat.c
 * Description:
 * Author:       Kunio Nakamura
 * Date:         January 9, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *xfm2niikmat_version[] = {
  "  xfm2niikmat version history\n"
  "  version 0.0.0, January 9, 2013, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *xfm2niikmat_help[] = {
  "  xfm2niikmat usage:   <mat.xfm> <movimg.nii> <target.nii> <out.matrx> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niikmath general tool\n");
  fprintf(stdout,"  usage: xfm2niikmat <mat.xfm> <movimg.nii> <target.nii> <out.matrx> [options]\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *tgtimg=NULL,
   *movimg=NULL;
  char
  *outcheckname=NULL;
  const char *fcname="xfm2niikmat";
  int
  nc=1,sc=1;

  niikmat
  *tmpmat,
  *afmat=NULL;

  if(argc==1) {
    usage();
    exit(0);
  }

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*xfm2niikmat_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*xfm2niikmat_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*xfm2niikmat_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        usage();
        exit(0);
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc!=5) {
    fprintf(stderr,"[%s] ERROR: wrong usage\n",fcname);
    exit(0);
  }

  fprintf(stdout,"[%s] reading xfm     %s\n",fcname,argv[1]);
  NIIK_EXIT(((afmat = niikmat_read_xfm(argv[1]))==NULL),
            fcname,"niikmat_read_xfm",9);
  fprintf(stdout,"[%s] reading image   %s\n",fcname,argv[2]);
  NIIK_EXIT(((movimg=niik_image_read(argv[2]))==NULL),
            fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading tgtimg  %s\n",fcname,argv[3]);
  NIIK_EXIT(((tgtimg=niik_image_read(argv[3]))==NULL),
            fcname,"niik_image_read",9);

  /* delta*inv(S2)*M*S*inv(delta) */
  tmpmat = niikmat_mat44_matrix(movimg->sto_xyz);
  niikmat_multiply_mat1_free2(tmpmat,
                              niikmat_scale_matrix(1.0/movimg->dx,1.0/movimg->dy,1.0/movimg->dz));
  niikmat_multiply_mat1_free2(afmat,
                              tmpmat);
  niikmat_multiply_mat2_free1(niikmat_mat44_matrix(tgtimg->sto_ijk),
                              afmat);
  niikmat_multiply_mat2_free1(niikmat_scale_matrix(tgtimg->dx,tgtimg->dy,tgtimg->dz),
                              afmat);
  niikmat_display(afmat);

  fprintf(stdout,"[%s] writing output  %s\n",fcname,argv[4]);
  if(!niikmat_write(argv[4],afmat)) {
    fprintf(stderr,"[%s] ERROR: niikmat_write %s\n",fcname,argv[4]);
    exit(0);
  }

  tgtimg=niik_image_free(tgtimg);
  movimg=niik_image_free(movimg);
  afmat=niikmat_free(afmat);
  exit(0);
} /* xfm2niikmat */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/