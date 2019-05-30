/* Filename:     niikaregistre.c
 * Description:  Kunio's nifti1 affine registration program
 * Author:       Kunio Nakamura
 * Date:         December 3, 2012
 */

#include "falcon.h"

#define NIIKAREGISTER_MAJOR_VERSION (0)
#define NIIKAREGISTER_MINOR_VERSION (0)
#define NIIKAREGISTER_MICRO_VERSION (0)

static char *niikaregister_version[] = {
  "  niikaregister version history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *niikaregister_help[] = {
  "  niikaregister usage:   <img.nii> <ref.nii> <out.matrx> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
  "  -refseg <refseg.nii>              : input reference mask\n"
  /*  "  -refseg-list <mask1>,<mask2>...   : input reference mask\n"*/
  "  -outcheck <check.nii>             : output image for checking registration\n"
  "  -hconf <hconf.txt>                : hierarchical staging\n"
  "      Please create a text file with the following format:\n"
  "      blur_FWHM sampling_distance optimization_tolerance max_iteration  #seeds    delta_rx delta_ry delta_rz    delta_tx delta_ty delta_tz    delta_sg delta_sx delta_sy delta_sz      delta_skew_xy delta_skew_xz delta_skew_yz    delta_center_x delta_center_y delta_center_z \n"
  "         blur_FWHM         : blurring filter's FWHM\n"
  "         sampling_distance : sample distance (larger is faster, smaller may be more accurate?)\n"
  "         toplerance        : optimization tolerance (larger is faster but less accurate)\n"
  "         max_iteration     : iteration for each optimization for each seed\n"
  "         #seeds            : # seeds for each level\n"
  "         delta_xx          : search region for rotation (r), translation (t), scaling (global, and x,y,z), skewing, and center\n"
  "      For multi-level processing, add multiple lines so that the number of levels equals the number of lines\n"
  "  -iconf <iconf.txt>         : initial seeds\n"
  "     Please create a text file with a matrix\n"
  "     0 rx0 ry0 rz0 tx0 ty0 tz0 sg0 sx0 sy0 sz0 kx0 ky0 kz0 cx0 cy0 cz0\n"
  "     0 rx1 ry1 rz1 tx1 ty1 tz1 sg1 sx1 sy1 sz1 kx1 ky1 kz1 cx1 cy1 cz1 ...\n"
  "     -Each line is a seed.\n"
  "     -The number of seeds here must be equal or larger than #seeds in hconf.\n"
};

void usage() {
  fprintf(stdout,"niikmath general tool\n");
  fprintf(stdout,"  usage: niikaregister <img.nii> <ref.nii> <out.matrix> [options]\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *outimg=NULL,
   **imglist=NULL,
     *refimg=NULL,
      *refseg=NULL,
       **refseglist=NULL,
         *crefimg=NULL,
          *cmovimg=NULL,
           *movimg=NULL;
  char
  *outcheckname=NULL,
   fcname[128]="niikaregister";
  int
  verbose=2,
  n,
  nlevel,nl,
  *maxiter=NULL,
   *nseed=NULL,ns,
    areg_cost=NIIK_REGISTER_CC,
    nmi_num[2],
    numrefseglist=0,
    nc=1,sc=1;
  niikvec
  *tol=NULL,
   *delta=NULL,
    *FWHM=NULL;
  double
  imin,imax,omin,omax,
       *affpar=NULL;

  niikmat
  *daffpar=NULL,
   *iconf=NULL,
    *hconf=NULL,
     *afmat=NULL;
  niikpt movctr,refctr;
  nmi_obj *nmiobj=NULL;

  if(argc==1) {
    usage();
    exit(0);
  }

  affpar=(double *)calloc(25,sizeof(double));     /* affine parameters */
  affpar[7]=affpar[8]=affpar[9]=affpar[10]=1;
  nmi_num[0]=nmi_num[1]=0;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niikaregister_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*niikaregister_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-nmi-hist-num",13)) {
        NIIK_RET0((argc>=nc+2),fcname,"missing parameters");
        nmi_num[0] = atoi(argv[++nc]);
        nmi_num[1] = atoi(argv[++nc]);
        fprintf(stdout,"[%s] histogram # bin = %i %i\n",fcname,nmi_num[0],nmi_num[1]);
      } /* histogram bin size */

      else if(!strncmp(argv[nc],"-refseg-list",12)) {
        NIIK_RET0((argc>=nc+1),fcname,"missing parameters");
        nc++;
        if(nc>=argc) {
          usage();
          exit(0);
        }
        if(refseglist!=NULL) {
          fprintf(stderr,"[%s] ERROR: refseglist is already used\n",fcname);
          exit(0);
        }
        if((refseglist = niik_image_read_multiple(argv[nc],&numrefseglist))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_multiple\n",argv[nc]);
          exit(0);
        }
      }

      else if(!strncmp(argv[nc],"-movctr",7)) {
        movctr.x = atof(argv[++nc]);
        movctr.y = atof(argv[++nc]);
        movctr.z = atof(argv[++nc]);
        movctr.w = 1;
      }

      else if(!strncmp(argv[nc],"-refseg",7)) {
        nc++;
        if(nc>=argc) {
          usage();
          exit(0);
        }
        if(refseg!=NULL) {
          fprintf(stderr,"[%s] ERROR: refseg is already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading refimg  %s\n",fcname,argv[nc]);
        if((refseg=nifti_image_read(argv[nc],1))==NULL) {
          fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      }

      else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      }

      else if(!strncmp(argv[nc],"-hconf",5)) {
        hconf=niikmat_read(argv[++nc]);
      } /* hierarchy configuration */

      else if(!strncmp(argv[nc],"-iconf",5)) {
        iconf=niikmat_read(argv[++nc]);
      } /* seed configuration */

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niikaregister_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-nmi",4)) {
        areg_cost=NIIK_REGISTER_NMI;
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

  if(argc!=4) {
    fprintf(stderr,"[%s] ERROR: wrong usage\n",fcname);
    exit(0);
  }

  if(hconf==NULL) {
    fprintf(stderr,"[%s] ERROR: hconf is null\n",fcname);
    exit(0);
  }
  if(iconf==NULL) {
    fprintf(stderr,"[%s] ERROR: iconf is null\n",fcname);
    exit(0);
  }

  fprintf(stdout,"[%s] reading image   %s\n",fcname,argv[1]);
  if((movimg=niik_image_read(argv[1]))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[1]);
    exit(0);
  }

  fprintf(stdout,"[%s] reading refimg  %s\n",fcname,argv[2]);
  if((refimg=niik_image_read(argv[2]))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[2]);
    exit(0);
  }

  if(movctr.w==0) {
    movctr=niikpt_image_get_centroid(movimg,NULL);
  }
  affpar[14]=movctr.x;
  affpar[15]=movctr.y;
  affpar[16]=movctr.z;
  fprintf(stdout,"[%s] using center %9.5f %9.5f %9.5f (movctr)\n",fcname,movctr.x,movctr.y,movctr.z);

  refctr=niikpt_image_get_centroid(refimg,NULL);
  affpar[14]=refctr.x;
  affpar[15]=refctr.y;
  affpar[16]=refctr.z;
  fprintf(stdout,"[%s] using center %9.5f %9.5f %9.5f (refctr)\n",fcname,refctr.x,refctr.y,refctr.z);


  if(areg_cost==NIIK_REGISTER_NMI) {
    imin=niik_image_get_min(movimg,NULL);
    imax=niik_image_get_max(movimg,NULL);
    if(refseglist==NULL) {
      omin=niik_image_get_min(refimg,refseg);
      omax=niik_image_get_max(refimg,refseg);
    } else {
      omin=niik_image_get_min(refimg,refseglist[0]);
      omax=niik_image_get_max(refimg,refseglist[0]);
    }
    if(!(nmi_num[0] * nmi_num[1])) {
      nmi_num[0] = nmi_num[1] = 32;
      nmi_num[0] *= imax / niik_image_get_percentile(movimg,NULL,0.9);
      nmi_num[1] *= omax / niik_image_get_percentile(refimg,NULL,0.9);
    }
    if(!niik_aregister_nmi_obj_update_var(imin,imax,omin,omax,nmi_num[0],nmi_num[1])) {
      fprintf(stderr,"[%s] ERROR: niik_aregister_nmi_obj_update_var\n",fcname);
      exit(0);
    }
    nmiobj=niik_aregister_g_nmi_obj_get();
    niik_aregister_nmi_obj_display(nmiobj);
    niik_image_aregister_set_nmi_obj(nmiobj);
  }


  /* parameters */
  if(verbose>=2) fprintf(stdout,"parameters %i %i \n",hconf->col,hconf->row);
  nlevel=hconf->row;
  if(refseglist!=NULL) {
    if(nlevel != numrefseglist) {
      fprintf(stderr,"[%s] ERROR: refseg list (%i) and #level (%i) did not match\n",fcname,numrefseglist,nlevel);
      exit(0);
    }
  }
  if(verbose>=2) fprintf(stdout,"nlevel %i \n",nlevel);
  nseed=(int *)calloc(nlevel,sizeof(int));
  maxiter=(int *)calloc(nlevel,sizeof(int));
  FWHM=niikvec_init(nlevel);
  delta=niikvec_init(nlevel);
  tol=niikvec_init(nlevel);
  daffpar=niikmat_init(nlevel,20);
  for(nl=0,ns=0; nl<nlevel; nl++) {
    ns=NIIK_IMAX(ns,(int)hconf->m[nl][4]);
  }
  /*fprintf(stdout,"ns %i \n",ns);*/
  if(ns>iconf->row) {
    fprintf(stderr,"ERROR: #seed %i > iconf %i\n",ns,iconf->row);
    exit (0);
  }

  for(nl=0; nl<nlevel; nl++) {
    /*fprintf(stdout,"nlevel %i \n",nl);*/
    FWHM ->v[nl]=hconf->m[nl][0];
    delta->v[nl]=hconf->m[nl][1];
    tol  ->v[nl]=hconf->m[nl][2];
    maxiter [nl]=(int)hconf->m[nl][3];
    nseed   [nl]=(int)hconf->m[nl][4];
    for(n=1; n<17; n++) {
      daffpar->m[nl][n]=hconf->m[nl][n+5];
    }
  }
  for(n=0; n<ns; n++) {
    iconf->m[n][14]+=movctr.x;
    iconf->m[n][15]+=movctr.y;
    iconf->m[n][16]+=movctr.z;
    iconf->m[n][4] -= movctr.x - refctr.x;
    iconf->m[n][5] -= movctr.y - refctr.y;
    iconf->m[n][6] -= movctr.z - refctr.z;
  }

  if(refseglist!=NULL) {
    if(!niik_image_aregister2_multilevel_multimask(refimg,refseglist,movimg,afmat,affpar,areg_cost,nmiobj,nlevel,FWHM,delta,tol,nseed,maxiter,daffpar,iconf)) {
      fprintf(stderr,"[niik_image_aregister2_test1] ERROR: niik_image_aregister2_multilevel\n");
      return 0;
    }
  } else {
    if(!niik_image_aregister2_multilevel(refimg,refseg,movimg,afmat,affpar,areg_cost,nmiobj,nlevel,FWHM,delta,tol,nseed,maxiter,daffpar,iconf)) {
      fprintf(stderr,"[niik_image_aregister2_test1] ERROR: niik_image_aregister2_multilevel\n");
      return 0;
    }
  }

  niik_aregister_display_affine(affpar);
  if(afmat==NULL) afmat=niikmat_init(4,4);
  niik_aregister_matrix_from_affpar_update(afmat,affpar);

  fprintf(stdout,"[%s] writing %s\n",fcname,argv[3]);
  if(!niikmat_write(argv[3],afmat)) {
    fprintf(stderr,"[%s] ERROR: niikmat_write %s\n",fcname,argv[3]);
    exit(0);
  }

  if(outcheckname!=NULL) {
    fprintf(stdout,"    transform image for checking: %s\n",outcheckname);
    if((outimg=niik_image_affine_transform_3d(movimg,refimg,afmat,NIIK_INTERP_BSPLINE))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d\n",fcname);
      exit(0);
    }
    imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
    imglist[0]=outimg;
    imglist[1]=refimg;
    if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_combine\n",fcname);
      exit(0);
    }
    if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
      exit(0);
    }
    outimg->cal_min=0;
    outimg->cal_max=200;
    fprintf(stdout,"[%s] writing output  %s\n",fcname, outcheckname);
    if(!niik_image_write(outcheckname,outimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outcheckname);
      exit(0);
    }
    outimg=niik_image_free(outimg);
    free(imglist);
  }

  crefimg=niik_image_free(crefimg);
  cmovimg=niik_image_free(cmovimg);

  refimg=niik_image_free(refimg);
  movimg=niik_image_free(movimg);
  free(affpar);
  free(daffpar);
  exit(0);
} /* niikaregister */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/