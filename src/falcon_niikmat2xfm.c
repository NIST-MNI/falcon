/* Filename:     niikmat2xfm.c
 * Description:
 * Author:       Kunio Nakamura
 * Date:         January 9, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (2)
#define MICRO_VERSION (0)

static char *niikmat2xfm_version[] = {
  "  niikmat2xfm version history\n"
  "\n  version 0.0.0, January 9, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
  "\n  version 0.1.0, August 16, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -corrected usage\n"
  "\n  version 0.2.0, October 11, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -uses Kunio's niikmat_inverse to calculate matrix inversion instead of nifti_io.c\n"
  "   to avoid floating point precision problem, (A^-1*A) was not identity\n"
};

static char *niikmat2xfm_help[] = {
  "  niikmat2xfm usage:   <mat.mat> <movimg.nii> <tgtimg.nii> <out.xfm> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niikmath general tool\n");
  fprintf(stdout,"  usage: niikmat2xfm <mat.mat> <movimg.nii> <tgtimg.nii> <out.xfm> [options]\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *tgtimg=NULL,
   *movimg=NULL;
  char
  *outcheckname=NULL,
   fcname[128]="niikmat2xfm";
  int
  nc=1,sc=1;
  niikmat
  *afmat=NULL;
  int flag_inv=0;

  if(argc==1) {
    usage();
    exit(0);
  }

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niikmat2xfm_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*niikmat2xfm_help);
        exit(0);
      } else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niikmat2xfm_help);
        exit(0);
      } else if(!strncmp(argv[nc],"-inv",4)) {
        flag_inv=1;
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

  fprintf(stdout,"[%s] reading mat     %s\n",fcname,argv[1]);
  NIIK_EXIT(((afmat = niikmat_read(argv[1]))==NULL),fcname,"niikmat_read",9);
  fprintf(stdout,"[%s] reading image   %s\n",fcname,argv[2]);
  NIIK_EXIT(((movimg=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading image   %s\n",fcname,argv[3]);
  NIIK_EXIT(((tgtimg=niik_image_read(argv[3]))==NULL),fcname,"niik_image_read",9);

  if(flag_inv) {
    niikmat_inverse_update(afmat);
  }

  NIIK_EXIT((!niikmat_convert_to_xfm(afmat,movimg,tgtimg)),
            fcname,"niikmat_convert_to_xfm",9);

  fprintf(stdout,"[%s] writing output  %s\n",fcname,argv[4]);
  NIIK_EXIT((!niikmat_write_as_linear_xfm(argv[4],afmat)),
            fcname,"niikmat_write_as_linear_xfm",9);

  /*
  MNI Transform File
  %Thu Jan 10 10:29:12 2013>>> /opt/minc/bin/minctracc /var/tmp/mritoself_11057//ASSERT_273-NOH-1_ale_927309_m0_t1g.mnc.gz_float_t1_nlm_nu_cropS.mnc /var/tmp/mritoself_11057//ASSERT_273-NOH-1_ale_927309_m12_t1g.mnc.gz_float_t1_nlm_nu_cropT.mnc img1-to-img2.xfm -est_center -transformation /var/tmp/mritoself_11057//img1-to-img2_tmp1.xfm -lsq6 -mi -groups 256 -threshold 121 0 -step 4.3 4.3 4.3 -simplex 1.5
  %(mni_autoreg 0.99.60)

  Transform_Type = Linear;
  Linear_Transform =
  0.997713267803192 -0.0660243928432465 -0.0144597198814154 -0.773668766021729
  0.0661535039544106 0.997772037982941 0.00864013098180294 1.46318197250366
  0.0138570452108979 -0.00957693345844746 0.99985808134079 -5.23042392730713;
  */

  movimg=niik_image_free(movimg);
  afmat=niikmat_free(afmat);
  exit(0);
} /* niikmat2xfm */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/