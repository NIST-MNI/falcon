/* FILENAME:     niik_hipposurf_resample.c
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 * DATE:         February 8, 2015
 */

#include "falcon.h"
#include "falcon_hipposurf.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"niik_make_hipposurf\n");
  fprintf(stdout,"  usage: [options] <in.off> <ref.off> <out_vals>\n");
  fprintf(stdout," <in.off>          input off (for spherical parameters\n");
  fprintf(stdout," <ref.off>         reference off (for spherical parameters)\n");
  fprintf(stdout," <out_vals.txt>    output values (ref->#vert)\n");
  fprintf(stdout,"\n");
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *in=NULL;
  kobj *ref=NULL;
  kvert *v;
  int
  i,
  nc=1;
  char fcname[32]="niik_hipposurf";
  niikmat *vin,*vout,*vavg;
  int ndim=15;

  niik_disp_exec_info(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);

  if(argc==1) {
    usage();
    exit(1);
  }

  fprintf(stdout,"[%s] reading ref   %s\n",fcname,argv[1]);
  NIIK_RETURN(((ref=off_kobj_read_offply(argv[1]))==NULL),"off_kobj_read",1);
  vout=niikmat_init(ndim,ref->nvert);
  vavg=niikmat_init(ndim,ref->nvert);
  for(nc=3; nc<argc; nc++) {
    fprintf(stdout,"[%s] reading obj   %s\n",fcname,argv[nc]);
    NIIK_RETURN(((in=off_kobj_read_offply(argv[nc]))==NULL),"off_kobj_read",1);
    vin=niikmat_init(ndim,in->nvert);
    NIIK_RETURN((!niik_hipposurf_extract_feature(in,vin)),"niik_hipposurf_extract_feature",1);
    NIIK_RETURN((!niik_off_resample_parametric_surface(in,ref,vin,vout)),
                "niik_off_resample_parametric_surface",1);
    niikmat_add_to_a(vavg,vout);
    off_kobj_free(in);
  }
  niikmat_kmul(vavg,1/(argc-3.0));

  for(v=ref->vert,i=0; v!=NULL; v=v->next,i++) {
    v->v.x = vavg->m[0][i];
    v->v.y = vavg->m[1][i];
    v->v.z = vavg->m[2][i];
  }

  fprintf(stdout,"[%s] writing %s\n",fcname,argv[2]);
  off_kobj_write_offply(argv[2],ref,0);

  niik_fc_display(fcname,0);
  exit(0);
} /* niik_nls_lesion */
