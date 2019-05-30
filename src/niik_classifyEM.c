/* Filename:     niik_classifyEM.c
 * Description:  Kunio's nifti1 classification program using EM
 * Author:       Kunio Nakamura
 * Date:         December 7, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *niik_classify_version[] = {
  "  niik_classify EM version history\n"
  "  0.0.0  December 7, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
};

static char *niik_classify_help[] = {
  "  simple classification program\n"
  "  usage: <in.nii> <anat.nii> <out.nii> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help                   : show this usage\n"
  "  --version                  : show version info\n"
  "\n"
  "  <in.nii> can be 4-dimensional in u-dimension\n"
  "\n"
};


void usage() {
  fprintf(stdout,"%s",*niik_classify_help);
  return;
}

nifti_image *niik_classifyEM(nifti_image *img,nifti_image *anat);

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *tmpimg=NULL,
   *tmpimgs[9],
   *label=NULL,
    *out=NULL,
     *anat=NULL,
      *img=NULL;
  char
  *outcheckname=NULL,
   *outlabelname=NULL,
    fcname[20]="niik_classifyEM";
  int
  maxiter=5,
  num_class=4,
  verbose=2,
  n,
  nc=1,sc=1;
  double
  dlist[12],
        aUniP=0.5,
        aFWHM=5.0;

  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  if(argc==1) {
    usage();
    exit(0);
  }

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);
  fprintf(stdout,"  niik_classifier: version %i.%i.%i\n",MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  fprintf(stdout,"  niik_classifier: executed at %s\n",tmstr);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niik_classify_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*niik_classify_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      } /* outcheck */

      else if(!strncmp(argv[nc],"-outlabel",9)) {
        outlabelname=argv[++nc];
        fprintf(stdout,"[%s] out label name %s\n",fcname,outlabelname);
      } /* outlabel */

      else if(!strncmp(argv[nc],"-nclass",7)) {
        num_class=atoi(argv[++nc]);
      } /* class */

      else if(!strncmp(argv[nc],"-aFWHM",6)) {
        aFWHM=atof(argv[++nc]);
      } /* FWHM */

      else if(!strncmp(argv[nc],"-aUniP",6)) {
        aUniP=atof(argv[++nc]);
      } /* FWHM */

      else if(!strncmp(argv[nc],"-iter",5)) {
        maxiter=atoi(argv[++nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niik_classify_help);
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

  NIIK_EXIT((argc<4),fcname,"too few argments",9);
  NIIK_EXIT((argc>4),fcname,"too many argments",9);

  fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading anat img %s\n",fcname,argv[2]);
  NIIK_EXIT(((anat=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] anat processing   blur %0.3f  UniformMin Prob %0.3f\n",fcname,aFWHM,aUniP);
  NIIK_EXIT((!niik_image_filter_gaussian_update(anat,11,aFWHM)),fcname,"niik_image_filter_gaussian_update",9);
  for(n=0; n<anat->nvox; n++) {
    niik_image_set_voxel(anat,n,niik_image_get_voxel(anat,n)*(1.0-aUniP)+aUniP);
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] writing output image test_anat.nii.gz\n",fcname);
    NIIK_EXIT((!niik_image_write("test_anat.nii.gz",anat)),fcname,"niik_image_write",9);
  }

  /* input must be mask, t1p, pdw, t2w, and flr in u-dimension */
  if(verbose>1) fprintf(stdout,"[%s] calling niik_classifyEM\n",fcname);
  NIIK_EXIT(((out=niik_classifyEM(img,anat))==NULL),fcname,"niik_classifyEM",9);

  if(outlabelname!=NULL) {/* label */
    NIIK_RET0(((label=niik_image_init(img->nx,img->ny,img->nz,1,1,1,1,
                                      img->dx,img->dy,img->dz,1,1,1,1,NIFTI_TYPE_UINT8))==NULL),
              fcname,"niik_image_init");
    for(n=0; n<label->nvox; n++) {
      if(niik_image_get_voxel(img,n)<0.5) continue;
      dlist[0]=niik_image_get_voxel(out,n);
      dlist[1]=niik_image_get_voxel(out,n+label->nvox);
      dlist[2]=niik_image_get_voxel(out,n+label->nvox*2);
      dlist[3]=niik_image_get_voxel(out,n+label->nvox*3);
      if     (dlist[3]>dlist[0] && dlist[3]>dlist[1] && dlist[3]>dlist[2]) niik_image_set_voxel(label,n,8);
      else if(dlist[0]>dlist[1] && dlist[0]>dlist[1] && dlist[0]>dlist[2]) niik_image_set_voxel(label,n,1);
      else if(dlist[1]>dlist[0] && dlist[1]>dlist[2] && dlist[1]>dlist[2]) niik_image_set_voxel(label,n,2);
      else if(dlist[2]>dlist[0] && dlist[2]>dlist[1] && dlist[2]>dlist[3]) niik_image_set_voxel(label,n,4);
    }
    if(verbose>0) fprintf(stdout,"[%s] writing output image %s\n",fcname,outlabelname);
    NIIK_EXIT((!niik_image_write(outlabelname,label)),fcname,"niik_image_write",9);
    label=niik_image_free(label);
  } /* label */


  /* one type of output:
   * includes input images and output probability maps
   */
  for(n=0; n<out->nvox; n++) {
    niik_image_mul_voxel(out,n,100.0);
  }
  tmpimgs[0]=img;
  tmpimgs[1]=out;
  NIIK_EXIT(((tmpimg=niik_image_combine(tmpimgs,2,5,0))==NULL),fcname,"niik_image_combine",9);
  if(verbose>0) fprintf(stdout,"[%s] writing output image %s\n",fcname,argv[3]);
  NIIK_EXIT((!niik_image_write(argv[3],tmpimg)),fcname,"niik_image_write",9);


  img=niik_image_free(img);
  anat=niik_image_free(anat);
  tmpimg=niik_image_free(tmpimg);

  niik_fc_display(fcname,0);
  exit(0);
} /* niik_classify */


void niik_classify_disp_stats(niikmat *mu,niikmat *sd,int nimg) {
  fprintf(stdout,"  mean / stdv matrix\n");
  fprintf(stdout,"          mask    t1p     pdw      t2w      flr\n");
  fprintf(stdout,"  CSF:  ");
  niik_display_double_vector(mu->m[0],nimg);
  fprintf(stdout,"   GM:  ");
  niik_display_double_vector(mu->m[1],nimg);
  fprintf(stdout,"   WM:  ");
  niik_display_double_vector(mu->m[2],nimg);
  fprintf(stdout,"  LES:  ");
  niik_display_double_vector(mu->m[3],nimg);
  fprintf(stdout,"   BG:  ");
  niik_display_double_vector(mu->m[4],nimg);
  fprintf(stdout,"  CSF:  ");
  niik_display_double_vector(sd->m[0],nimg);
  fprintf(stdout,"   GM:  ");
  niik_display_double_vector(sd->m[1],nimg);
  fprintf(stdout,"   WM:  ");
  niik_display_double_vector(sd->m[2],nimg);
  fprintf(stdout,"  LES:  ");
  niik_display_double_vector(sd->m[3],nimg);
  fprintf(stdout,"   BG:  ");
  niik_display_double_vector(sd->m[4],nimg);
  return;
}


nifti_image *niik_classifyEM(nifti_image *img,nifti_image *anat)
/* -img needs dim in u
 *    mask(currently not used),t1p,pdw,t2w,flr
 * -anat has dim in u
 *    csf,gm,wm for normal anatomy
 *
 */
{
  char fcname[32]="niik_classifyEM";
  int verbose=2;
  nifti_image *out=NULL;
  niikmat *mu=NULL,*sd=NULL;
  niikvec *p=NULL,*v=NULL;
  int
  nclass=5,
  nimg=-1;
  int
  check_voxel=(58+76*192+44*192*256),
  i,   /*vox idx*/
  n,m, /*class,img idx*/
  n3,  /*#voxel in 3d*/
  ni,  /*img idx*/
  nmsk=0,nt1p=1,npdw=2,nt2w=3,nflr=4,
  nc=0,
  ng=1,
  nw=2,
  nl=3,
  nb=4;
  unsigned char *mask=NULL;
  double
  sv,st,sx;
  int
  iter,maxiter=1;

  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim<=4),fcname,"img is missing u-dim");

  nimg=img->nu;
  NIIK_RET0((nimg!=5),fcname,"img needs (mask,t1p,pdw,t2w,flr)");

  fprintf(stdout,"[%s] #class = %i\n",fcname,nclass);
  fprintf(stdout,"[%s] #img   = %i\n",fcname,nimg);
  fprintf(stdout,"[%s] initialization\n",fcname);
  NIIK_RET0(((out=niik_image_init(img->nx,img->ny,img->nz,1,nclass,1,1,
                                  img->dx,img->dy,img->dz,1,1,1,1,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_image_init");
  n3=img->nx*img->ny*img->nz;

  mu=niikmat_init(nclass,nimg);
  sd=niikmat_init(nclass,nimg);
  p=niikvec_init(nclass);
  v=niikvec_init(nclass);
  NIIK_RET0(((mask=(unsigned char *)calloc(n3,sizeof(char)))==NULL),fcname,"calloc for mask");
  for(i=0; i<n3; i++) {
    mask[i]=(niik_image_get_voxel(img,i)>0);
  }

  ni=0;
  mu->m[nc][ni]=  1;
  mu->m[ng][ni]= 1;
  mu->m[nw][ni]= 1;
  mu->m[nl][ni]=  1;
  mu->m[nb][ni]= 0; /*mask*/
  ni=1;
  mu->m[nc][ni]= 30;
  mu->m[ng][ni]=70;
  mu->m[nw][ni]=90;
  mu->m[nl][ni]= 75;
  mu->m[nb][ni]= 0; /*t1p*/
  ni=2;
  mu->m[nc][ni]= 70;
  mu->m[ng][ni]=85;
  mu->m[nw][ni]=65;
  mu->m[nl][ni]= 80;
  mu->m[nb][ni]= 0; /*pdw*/
  ni=3;
  mu->m[nc][ni]=120;
  mu->m[ng][ni]=75;
  mu->m[nw][ni]=70;
  mu->m[nl][ni]= 90;
  mu->m[nb][ni]= 0; /*t2w*/
  ni=4;
  mu->m[nc][ni]= 40;
  mu->m[ng][ni]=85;
  mu->m[nw][ni]=75;
  mu->m[nl][ni]=100;
  mu->m[nb][ni]= 0; /*flr*/

  ni=0;
  sd->m[nc][ni]=  0;
  sd->m[ng][ni]= 0;
  sd->m[nw][ni]= 0;
  sd->m[nl][ni]= 0;
  sd->m[nb][ni]=  0; /*mask*/
  ni=1;
  sd->m[nc][ni]= 35;
  sd->m[ng][ni]=10;
  sd->m[nw][ni]=10;
  sd->m[nl][ni]=25;
  sd->m[nb][ni]= 10; /*t1p*/
  ni=2;
  sd->m[nc][ni]= 35;
  sd->m[ng][ni]=10;
  sd->m[nw][ni]=10;
  sd->m[nl][ni]=10;
  sd->m[nb][ni]= 10; /*pdw*/
  ni=3;
  sd->m[nc][ni]= 30;
  sd->m[ng][ni]=15;
  sd->m[nw][ni]=15;
  sd->m[nl][ni]=10;
  sd->m[nb][ni]= 10; /*t2w*/
  ni=4;
  sd->m[nc][ni]= 40;
  sd->m[ng][ni]=20;
  sd->m[nw][ni]=20;
  sd->m[nl][ni]=20;
  sd->m[nb][ni]= 10; /*flr*/
  p->v[nc]=0.05;
  p->v[ng]=0.50;
  p->v[nw]=0.40;
  p->v[nl]=0.05;
  p->v[nb]=0.001;
  niikvec_set_all(p,1);

  niik_classify_disp_stats(mu,sd,nimg);

  for(iter=0; iter<maxiter; iter++) {

    /* update voxel-wise probabilities */
    for(i=0; i<n3; i++) { /*voxel*/
      /*if(mask[i]==0) continue;*/
      NIIK_RET0((niikvec_set_all(v,1)),fcname,"niikvec_set_all");
      for(n=0; n<nclass; n++) { /*class*/
        for(m=1; m<nimg; m++) { /*img*/
          if(i==check_voxel)
            fprintf(stdout,"[%3i %3i %3i   %i %i] %1.9f   %3i\n",
                    i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,
                    n,m,
                    NIIK_GaussPDF(niik_image_get_voxel(img,m*n3+i)-mu->m[n][m],sd->m[n][m]),
                    (int)niik_image_get_voxel(img,m*n3+i));
          if(m==nflr && n==nc)
            v->v[n]*=(1.0-NIIK_Heaviside(niik_image_get_voxel(img,m*n3+i)-
                                         mu->m[nc][m],sd->m[nc][m]));
          else if(m==nflr && n==nl)
            v->v[n]*=NIIK_Heaviside(niik_image_get_voxel(img,m*n3+i)-
                                    mu->m[nl][m],sd->m[nl][m]);
          else
            v->v[n]*=p->v[n]*NIIK_GaussPDF(niik_image_get_voxel(img,m*n3+i)-mu->m[n][m],sd->m[n][m]);

          /*anatomy : csf,gm,wm*/
          if(anat!=NULL) {
            if(n<3)        v->v[n]*=niik_image_get_voxel(anat,n*n3+i);
            else if(n==nl) v->v[n]*=niik_image_get_voxel(anat,nw*n3+i);
          }
        } /*img*/
        if(i==check_voxel)
          fprintf(stdout,"[%3i %3i %3i   %i  ] %1.9f\n",
                  i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,
                  v->v[n]);
      } /*class*/
      sv=niik_get_sum_from_double_vector(v->v,nclass);
      if(i==check_voxel) {
        fprintf(stdout,"[%3i %3i %3i      ] %1.9f\n",
                i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,
                sv);
        niikvec_display(v);
      }
      if(sv>1e-5) {
        for(n=0; n<nclass; n++) { /*class*/
          niik_image_set_voxel(out,n*n3+i,v->v[n]/sv);
        } /*class*/
      } else {
        for(n=0; n<nclass; n++) { /*class*/
          niik_image_set_voxel(out,n*n3+i,0);
        } /*class*/
      } /* can't tell for sure */
    } /*voxel*/

    /* update group stats */
    for(n=0; n<nclass; n++) { /*class*/
      for(m=1; m<nimg; m++) { /*img*/
        mu->m[n][m]=sd->m[n][m]=0;
        st=sx=0;
        for(i=0; i<n3; i++) { /*voxel*/
          /*if(mask[i]==0) continue; */
          sv=niik_image_get_voxel(out,n*n3+i);
          sx+=sv*niik_image_get_voxel(img,m*n3+i);
          st+=sv;
        } /*vox*/
        mu->m[n][m]=sx/st;
        sx=0;
        for(i=0; i<n3; i++) { /*voxel*/
          /*if(mask[i]==0) continue;*/
          sv=niik_image_get_voxel(out,n*n3+i);
          sx+=sv*NIIK_SQ(niik_image_get_voxel(img,m*n3+i)-mu->m[n][m]);
        } /*vox*/
        sd->m[n][m]=sqrt(sx/st);
      } /*img*/

      for(i=0,p->v[n]=0; i<n3; i++) { /*voxel*/
        p->v[n]+=niik_image_get_voxel(out,n*n3+i);
      }
    } /*class*/

    for(n=0,sv=0; n<nclass; n++) {
      sv+=p->v[n];
    }
    for(n=0; n<nclass; n++) {
      p->v[n]/=sv;
    }

    niik_classify_disp_stats(mu,sd,nimg);
  } /* iter */

  mu=niikmat_free(mu);
  sd=niikmat_free(sd);
  p=niikvec_free(p);
  v=niikvec_free(v);
  free(mask);
  mask=NULL;
  if(verbose>0) niik_fc_display(fcname,0);
  return out;
} /* niik_classifyEM */


/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/