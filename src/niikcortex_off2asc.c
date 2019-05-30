/* FILENAME:     off2asc.c
 * DESCRIPTION:  Conversion from off to text asc file
 * AUTHOR:       Vladimir Fonov
 *
 *
 * reference
 * http://www.danielgm.net/cc/doc/wiki/index.php5?title=FILE_I/O
 */

#include "falcon.h"
#include "falcon_surfaces.h"

int niik_off_write_asc(const char *fname,kobj *obj,double *meas,int output_sph);


static char *prog_version[] = {
  NIIK_VERSION "\n"
  "  program history\n"
  "  0.0.0   : August 16, 2016, Vladimir FONOV <vladimir.fonov@gmail.com>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  off2asc:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"obj2asc\n");
  fprintf(stdout,"  usage: [options] <in.off> [measurements] <out.asc>\n");
  fprintf(stdout,"  options:\n");
  fprintf(stdout,"\t-sph - output spherical coordinates\n");
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int
  nc,sc;
  int output_sph=0;
  const char
  *fcname="off2asc";
  double *meas=NULL;
  int meas_num;
  const char * out_fname;
  niik_fc_display(fcname,1);

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-sph",4)) {
        output_sph=1;
      } else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc>4) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(1);
  } else if(argc<3) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(1);
  }


  NIIK_EXIT(((obj=off_kobj_read_offply(argv[1]))==NULL),fcname,"niik_kobj_read_off",1);
  if(argc>3) {
    meas=niik_read_double_vector(argv[2],&meas_num);

    if(meas_num!=obj->nvert) {
      fprintf(stderr,"[%s] Inconsistent number of measurement: %i , expected %i\n",fcname,meas_num,obj->nvert);
      exit(1);
    }
    out_fname=argv[3];
  } else {
    out_fname=argv[2];
  }

  NIIK_EXIT((!niik_off_write_asc(out_fname,obj,meas,output_sph)),fcname,"niik_off_write_asc",1);

  obj=off_kobj_free(obj);
  if(meas!=NULL)
    free(meas);

  niik_fc_display(fcname,0);
  exit(0);
} /* off2mniasc */

int niik_off_write_asc(const char *fname,kobj *obj,double *meas,int output_sph) {
  const char *fcname="niik_off_write_asc";
  FILE *fp=NULL;
  int i;
  int verbose=niik_verbose();
  kvert *v;
  kface *f;

  niik_fc_display(fcname,1);

  if(verbose>1) {
    fprintf(stdout,"[%s] opening %s\n",fcname,fname);
  }

  /*recalculate normals*/
  /*off_update_kobj_vert_normal(obj);*/

  if((fp=fopen(fname,"w"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",fname);
    return 0;
  }
  fprintf(fp,"//X,Y,Z");
  if(obj->spherecoo && output_sph)
    fprintf(fp,",psi,the");
  if(meas!=NULL)
    fprintf(fp,",thk");
  fprintf(fp,"\n");

  /*output X Y Z coordinates of vertices*/

  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    fprintf(fp,"%15.15f,%15.15f,%15.15f",v->v.x,v->v.y,v->v.z);

    if(obj->spherecoo&& output_sph)
      fprintf(fp,",%15.15f,%15.15f",v->sph.psi,v->sph.the);

    if(meas!=NULL)
      fprintf(fp,",%15.15f",meas[i]);

    fprintf(fp,"\n");
  }

  fclose(fp);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_off_write_asc */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
