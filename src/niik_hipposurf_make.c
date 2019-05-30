/* FILENAME:     niik_hipposurf_make.c
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 * DATE:         February 6, 2015
 */

#include "falcon.h"
#include "falcon_hipposurf.h"

#define MAJOR_VERSION (NIIK_MAJOR_VERSION)
#define MINOR_VERSION (NIIK_MINOR_VERSION)
#define MICRO_VERSION (NIIK_MICRO_VERSION)

static char *niik_hipposurf_make_about_features [] = {
  "niik_hipposurf_make_about_features\n"
  "\nregional features, organized as index=0 for all vertices, then index=1 ...\n"
  "[0-3]    x,y,z coordinates\n"
  "[4]      area\n"
  "[5]      area per unit of angle in spherical space\n"
  "[6-7]    curvature, small and large-scale\n"
  "[8-9]    volume (triangle to the centroid), regular and absolute\n"
  "[10]     ratio of Euclidian edge length to the spherical edge length\n"
  "[11]     dot product of principal axis and vector from center\n"
  "[12-13]  same as 10 with second and third principal axes\n"
  "[14]     magnitude (ssq) of 12-13\n"
  "\nglobal features\n"
  "principal axes(pa):    x,y,z for 3 axes\n"
  "eigenvalues of pa:     3 values\n"
  "ratio of eigenvalues:  y to ssq of x&z\n"
  "area:                  object surface area\n"
  "volume:                object volume\n"
  "edge:                  mean edge length\n"
};

void usage() {
  fprintf(stdout,"niik_hipposurf_make\n");
  fprintf(stdout,"  usage: [options] <hippo_mask.mnc> <out>\n");
  fprintf(stdout," <hippo_mask.mnc>      : input mask; label 2  = left, label 4 = right or\n");
  fprintf(stdout,"                       : input mask; label 16 = left, label 3 = right\n");
  fprintf(stdout," <out>                 : output objects, <out>_lhippo.off and <out>_rhippo.off\n");
  fprintf(stdout,"\n  optional usage:\n");
  fprintf(stdout," -ref <lh> <rh>        : reference left hippo and right hippo\n");
  fprintf(stdout,"                       : extract features and resample\n");
  fprintf(stdout,"                       : additional outputs are <out>_lhippo_text.txt and \n                         <out>_rhippo_text.txt\n"
          fprintf(stdout," -elen <len>           : final edge length [default = 0.5]\n");
          fprintf(stdout,"\n");
}

int niik_hipposurf_make_write_features(kobj *obj,kobj *ref,char *fname) {
  FILE *fp=NULL;
  kvert
  *v=NULL;
  niikmat
  *vin=NULL,
   *vout=NULL;
  int
  i,j,
  ndim=15;
  niikpt
  *pa,
  ctr;

  vin=niikmat_init(ndim,obj->nvert);
  vout=niikmat_init(ndim+3,ref->nvert);
  vout->row=ndim;
  NIIK_RETURN((!niik_hipposurf_extract_feature(obj,vin)),"niik_hipposurf_extract_feature",1);
  NIIK_RETURN((!niik_off_resample_parametric_surface(obj,ref,vin,vout)),
              "niik_off_resample_parametric_surface",1);
  // distance from average hippo
  vout->row=ndim+3;
  for(v=ref->vert,i=0; v!=NULL; v=v->next,i++) {
    vout->m[ndim  ][i]=vout->m[0][i]-v->v.x;
    vout->m[ndim+1][i]=vout->m[1][i]-v->v.y;
    vout->m[ndim+2][i]=vout->m[2][i]-v->v.z;
  }
  /* global measures */
  // principal axes
  for(v=obj->vert,ctr=niikpt_zero(); v!=NULL; v=v->next) {
    ctr=niikpt_add(ctr,v->v);
  }
  ctr=niikpt_kmul(ctr,1.0/obj->nvert);
  NIIK_RETURN(((pa=(niikpt *)calloc(4,sizeof(niikpt)))==NULL),"calloc for pa",0);
  NIIK_RETURN((!niik_hipposurf_principal_axes(obj,pa,ctr)),"niik_hipposurf_principal_axes",0);
  for(i=0; i<3; i++) {
    fprintf(stdout,"[%s] principal axis %i = %9.4f %9.4f %9.4f\n",__func__,i,pa[i].x,pa[i].y,pa[i].z);
    pa[i]=niikpt_unit(pa[i]);
  }
  fprintf(stdout,"[%s] eigval           = %9.4f %9.4f %9.4f\n",__func__,pa[3].x,pa[3].y,pa[3].z);

  fprintf(stdout,"[%s] writing %s\n",__func__,fname);
  NIIK_RETURN(((fp=fopen(fname,"w"))==NULL),"fopen",1);
  for(i=0; i<vout->row; i++) {
    for(j=0; j<vout->col; j++) {
      fprintf(fp,"%9.4f ",vout->m[i][j]);
    }
  }
  for(i=0; i<4; i++) {
    fprintf(fp,"%9.4f %9.4f %9.4f ",pa[i].x,pa[i].y,pa[i].z);
  }
  fprintf(fp,"%9.4f ",pa[3].y/NIIK_SSQ(pa[3].x,pa[3].z));
  fprintf(fp,"%9.4f ",NIIK_SSQ(pa[3].x,pa[3].z));
  fprintf(fp,"%9.4f ",off_get_kobj_area(obj));
  fprintf(fp,"%9.4f ",off_get_kobj_volume(obj));
  fprintf(fp,"%9.4f ",off_get_kobj_mean_edge_length(obj));
  fprintf(fp,"\n");
  free(pa);
  vin=niikmat_free(vin);
  vout=niikmat_free(vout);
  return 1;
} /* niik_hipposurf_make_write_features */

int main(int argc,char *argv[],char *envp[]) {
  niik_hipposurf *h=NULL;
  int
  nc=1,sc=1;
  char fcname[32]="niik_hipposurf_make";
  char fname[4096];
  kobj *lh=NULL;
  kobj *rh=NULL;

  niik_disp_exec_info(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);

  if(argc==1) {
    usage();
    exit(1);
  }

  NIIK_RETURN(((h=niik_hipposurf_init())==NULL),"niik_hipposurf_init",1);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"-u",2)) {
        usage();
        fprintf(stdout,"________________________________________________\n%s\n",*niik_hipposurf_make_about_features);
        exit(1);
      } else if(!strncmp(argv[nc],"-v",2)) {
        NIIK_RETURN((nc+1>=argc),"missing argument(s)",1);
        h->verbose=atoi(argv[++nc]);
        fprintf(stdout,"[%s] verbose level %i\n",fcname,h->verbose);
      } else if(!strncmp(argv[nc],"-elen",5)) {
        NIIK_RETURN((nc+1>=argc),"missing argument(s)",1);
        h->elen=atof(argv[++nc]);
        fprintf(stdout,"[%s] elen %.3f\n",fcname,h->elen);
      } else if(!strncmp(argv[nc],"-ref",4)) {
        NIIK_RETURN((nc+1>=argc),"missing argument(s)",1);
        NIIK_RETURN(((lh=off_kobj_read_offply(argv[++nc]))==NULL),"off_kobj_read",1);
        NIIK_RETURN(((rh=off_kobj_read_offply(argv[++nc]))==NULL),"off_kobj_read",1);
      } else {
        fprintf(stderr,"[%s:%i:%s] ERROR: unknown option %s\n",__FILE__,__LINE__,__func__,argv[nc]);
        exit(1);
      }

      nc++;
    }

    else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  NIIK_RETURN((argc<3),"too few argments",1);
  NIIK_RETURN((argc>3),"too many argments",1);

  fprintf(stdout,"[%s] reading mask:  %s\n",fcname,argv[1]);
  NIIK_RETURN(((h->hippomask=niik_image_read(argv[1]))==NULL),"niik_image_read",1);

  NIIK_RETURN((!niik_hipposurf_make_surf(h)),"niik_hipposurf_make_surf",1);

  sprintf(fname,"%s_lhippo.off",argv[2]);
  fprintf(stdout,"[%s] writing left   %s\n",fcname,fname);
  off_kobj_write_offply(fname,h->lh,0);
  sprintf(fname,"%s_rhippo.off",argv[2]);
  fprintf(stdout,"[%s] writing right  %s\n",fcname,fname);
  off_kobj_write_offply(fname,h->rh,0);


  if(lh!=NULL) {
    sprintf(fname,"%s_lhippo_text.txt",argv[2]);
    niik_hipposurf_make_write_features(h->lh,lh,fname);
  }
  if(rh!=NULL) {
    sprintf(fname,"%s_rhippo_text.txt",argv[2]);
    niik_hipposurf_make_write_features(h->rh,rh,fname);
  }

  h=niik_hipposurf_free(h);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_hipposurf_make */
