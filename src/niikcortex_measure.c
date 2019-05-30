/* FILENAME: niikcortex_calc_thickness.c
 * DESCRIPTION: calculate cortical thickness using simple algorithm
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_point_inline.h"
#include "falcon_cortex.h"

#include <unistd.h>
#include <getopt.h>

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"Various cortical measurements\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  usage: [options] <white.ply> <pial.ply> <out.ply/txt> \n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t--ply        - will output .ply file format with variable thickness\n");
  fprintf(stdout,"\t--white      - use white surface for ply , otherwise - peal\n");
  fprintf(stdout,"\t--pial       - use pieal surface for ply , (default)\n");
  fprintf(stdout,"\t--thickness  - calculate thickness\n");
  fprintf(stdout,"\t--normcos    - Cosine similarities between normals\n");
  fprintf(stdout,"\t--linkcos    - Cosine similarities between normal and link vector\n");
  fprintf(stdout,"\t--all        - All of the above\n");
}


int niikcortex_measure_thickness(kobj *ics, kobj *ocs, double *meas)
/* -calculates vertex-wise cortical thickness
 * -filter_type :  median or mean
 * -filter_size :  size of regional dilation
 * -thk         : ics->nvert-length vector; updated with a list of cortical thickness
 * -ics/ocs :  inner (white) and outer (pial) surfaces
 */
{
  kvert *vi,*vo,**vlist;
  int m,num,vidx;
  double *dlist,*dlist_phi,*dlist_the;
  if(ics==NULL) {
    fprintf(stderr,"ERROR: ics is null\n");
    return 0;
  }
  if(ocs==NULL) {
    fprintf(stderr,"ERROR: ocs is null\n");
    return 0;
  }
  if(meas==NULL) {
    fprintf(stderr,"ERROR: meas is null\n");
    return 0;
  }
  if(ics->nvert!=ocs->nvert) {
    fprintf(stderr,"ERROR: nvert is different %i %i\n",ics->nvert,ocs->nvert);
    return 0;
  }

  for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
    meas[vidx] = niikpt_distance(vi->v,vo->v);
  }
  return 1;
}

int niikcortex_measure_ncos(kobj *ics, kobj *ocs, double *meas)
/* -calculates vertex-wise cortical thickness
 * -filter_type :  median or mean
 * -filter_size :  size of regional dilation
 * -thk         : ics->nvert-length vector; updated with a list of cortical thickness
 * -ics/ocs :  inner (white) and outer (pial) surfaces
 */
{
  kvert *vi,*vo,**vlist;
  int m,num,vidx;
  double *dlist,*dlist_phi,*dlist_the;
  if(ics==NULL) {
    fprintf(stderr,"ERROR: ics is null\n");
    return 0;
  }
  if(ocs==NULL) {
    fprintf(stderr,"ERROR: ocs is null\n");
    return 0;
  }
  if(meas==NULL) {
    fprintf(stderr,"ERROR: meas is null\n");
    return 0;
  }
  if(ics->nvert!=ocs->nvert) {
    fprintf(stderr,"ERROR: nvert is different %i %i\n",ics->nvert,ocs->nvert);
    return 0;
  }

  for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
    meas[vidx] = niikpt_dot(vi->normal, vo->normal);
  }
  return 1;
}

int niikcortex_measure_lcos(kobj *ics, kobj *ocs, double *meas)
/* -calculates vertex-wise cortical thickness
 * -filter_type :  median or mean
 * -filter_size :  size of regional dilation
 * -thk         : ics->nvert-length vector; updated with a list of cortical thickness
 * -ics/ocs :  inner (white) and outer (pial) surfaces
 */
{
  kvert *vi,*vo,**vlist;
  int m,num,vidx;
  double *dlist,*dlist_phi,*dlist_the;
  if(ics==NULL) {
    fprintf(stderr,"ERROR: ics is null\n");
    return 0;
  }
  if(ocs==NULL) {
    fprintf(stderr,"ERROR: ocs is null\n");
    return 0;
  }
  if(meas==NULL) {
    fprintf(stderr,"ERROR: meas is null\n");
    return 0;
  }
  if(ics->nvert!=ocs->nvert) {
    fprintf(stderr,"ERROR: nvert is different %i %i\n",ics->nvert,ocs->nvert);
    return 0;
  }

  for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
    niikpt link = niikpt_unit(niikpt_sub(vo->v, vi->v));
    meas[vidx]  = niikpt_dot(link, vo->normal); 
  }
  return 1;
}


int main(int argc,char *argv[],char *envp[]) {
  kobj *ctx[2];

  double *thk=NULL;
  double *ncos=NULL;
  double *lcos=NULL;

  double *vectors[3];

  const char *fcname="niikcortex_measure.c";
  int clobber=0;

  int output_sph=0;
  int output_ply=0;
  int output_white=0;

  int calc_thickness=0;
  int calc_normcos=0;
  int calc_linkcos=0;
  int calc_all=0;

  const char * meas[3];
  double *val[3];
  int meas_cnt=0;

  const char *in_ics,*in_ocs,*out_file;

  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"sph",              no_argument, &output_sph, 1},
    {"ply",              no_argument, &output_ply, 1},
    {"white",            no_argument, &output_white, 1},
    {"pial",             no_argument, &output_white, 0},
    {"help",             no_argument, 0, 'h'},
    {"version",          no_argument, 0, 'v'},

    {"thickness",        no_argument, &calc_thickness, 1},
    {"normcos",          no_argument, &calc_normcos, 1},
    {"linkcos",          no_argument, &calc_linkcos, 1},
    {"all",              no_argument, &calc_all, 1},

    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "uhv", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
      return 1;
    case '?':
    case 'u':
    case 'h':
    default:
      usage();
      return 1;
    }
  }

  if((argc - optind)<3) {
    usage();
    return 1;
  }
  in_ics=argv[optind];
  in_ocs=argv[optind+1];

  out_file=argv[optind+2];

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading white surface       %s\n",fcname,in_ics);
  NIIK_RETURN(((ctx[0]=off_kobj_read_offply(in_ics))==NULL),"reading white surface",1);

  fprintf(stdout,"[%s] reading pial surface        %s\n",fcname,in_ocs);
  NIIK_RETURN(((ctx[1]=off_kobj_read_offply(in_ocs))==NULL),"reading pial surface",1);

  NIIK_RETURN((ctx[0]->nvert!=ctx[1]->nvert),"wrong number of vertices",1);
  
  if(calc_all||calc_normcos||calc_linkcos)
  {
    /**/
    off_update_kobj_face_normal(ctx[0]);
    off_update_kobj_vert_normal(ctx[0]);
    off_update_kobj_face_normal(ctx[1]);
    off_update_kobj_vert_normal(ctx[1]);
  }

  if(calc_thickness || calc_all) {
        NIIK_RETURN(((thk=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for thickness array",1);
        NIIK_RETURN((! niikcortex_measure_thickness(ctx[0],ctx[1],thk)),"niikcortex_measure_thickness",1);
        meas[meas_cnt]="thickness";
        val[meas_cnt]=thk;
        meas_cnt++;
  }

  if(calc_normcos || calc_all) {
        NIIK_RETURN(((ncos=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for ncos array",1);
        NIIK_RETURN((! niikcortex_measure_ncos(ctx[0],ctx[1],ncos)),"niikcortex_measure_ncos",1);
        meas[meas_cnt]="ncos";
        val[meas_cnt]=ncos;
        meas_cnt++;
  }

  if(calc_linkcos || calc_all) {
        NIIK_RETURN(((lcos=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for lcos array",1);
        NIIK_RETURN((!niikcortex_measure_lcos(ctx[0],ctx[1],lcos)),"niikcortex_measure_lcos",1);
        meas[meas_cnt]="lcos";
        val[meas_cnt]=lcos;
        meas_cnt++;
  }

  fprintf(stdout,"[%s] writing output file        %s\n",fcname,out_file);
  fprintf(stdout,"[%s] total measurements: %d\n",fcname,meas_cnt);

  if(output_ply)
  {
    NIIK_RETURN( (!off_kobj_write_ply_ex(out_file, (output_white?ctx[1]:ctx[0]) ,0, 1,1,0, meas_cnt, meas, val)),"off_kobj_write_ply_ex",1);
  } else {
    NIIK_RETURN( (!niik_write_double_vectors_ex(out_file,val,ctx[0]->nvert,meas_cnt, meas)),"niik_write_double_vectors",1);
  }

  if(thk)  free(thk);
  if(ncos) free(ncos);
  if(lcos) free(lcos);

  ctx[0]=off_kobj_free(ctx[0]);
  ctx[1]=off_kobj_free(ctx[1]);

  niik_fc_display(fcname,0);
  exit(0);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
