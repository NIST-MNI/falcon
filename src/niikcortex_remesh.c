/* FILENAME:     niikcortex_refine.c
 * DESCRIPTION:  Kunio's refine cortical surfaces to fit the data
 * AUTHOR:       Kunio Nakamura
 * DATE:         March 24, 2013
 */
#include "falcon.h"
#include "falcon_cortex.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

int niikcortex_smooth_cortex(niikcortex_deform * dfm);
int niikcortex_remesh_cortex(niikcortex_deform * dfm);

void prog_history() {
  fprintf(stdout,"[falcon_remesh] history\n");
  fprintf(stdout,"  version  0.0.1  February 13, 2019, Vladimir S. FONOV <vladimir.fonov@gmail.com>\n");
}

void usage() {
  fprintf(stdout,"falcon_remesh\n");
  fprintf(stdout,"  usage: [options] <white.off> <pial.off> <out_white.off> <out_pial.off\n\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help                   : show this usage\n");
  fprintf(stdout,"  --version                         : show version info\n");
  fprintf(stdout,"  -debug-keep-tmp                   : keep debug files\n");
  fprintf(stdout,"  Processing options:               \n");
  fprintf(stdout,"  Optimization weights \n");

  fprintf(stdout,"  Additional optimizer parameters \n");
  fprintf(stdout,"  -t1w   <file>                     : T1w scan, for debugging\n");
  fprintf(stdout,"  -depth <n>                        : Quad-tree depth , default 7\n");
  fprintf(stdout,"  -delta <f>                        : time-step (default 0.5)\n");
  fprintf(stdout,"  -apply <f>                        : apply-step (default 0.2)\n");
  fprintf(stdout,"  -iter  <n>                        : maximum number of iterations, default 100\n");
  fprintf(stdout,"  -iter2 <n>                        : maximum number of sub-iterations, default 5\n");
  fprintf(stdout,"\n");
}




int niikcortex_remesh_process(niikcortex_deform *dfm) {
  const char *fcname="niikcortex_remesh_process";
  niik_fc_display(fcname,1);

  NIIK_RET0((dfm==NULL),fcname,"missing self");

  NIIK_RET0((dfm->ctx[0]->nvert!=dfm->ctx[1]->nvert),fcname,"#vert did not match");
  fprintf(stdout,"[%s] parameters\n",fcname);
  fprintf(stdout,"  surface vfe          %i %i %i\n",dfm->ctx[1]->nvert,dfm->ctx[1]->nface,dfm->ctx[1]->nedge);
  fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n",dfm->apply_step);
  fprintf(stdout,"  lambda               %-7.4f %-7.4f\n",dfm->weight->m[0][0],dfm->weight->m[1][0]);
  fprintf(stdout,"  mju                  %-7.4f %-7.4f\n",dfm->weight->m[0][1],dfm->weight->m[1][1]);
  fprintf(stdout,"  max iter             %i\n",dfm->iter);
  fprintf(stdout,"  max iter2            %i\n",dfm->iter2);
  fprintf(stdout,"  quadtree depth       %i\n",dfm->bbox_depth);
  fprintf(stdout,"  tolerance            %-7.4f\n",dfm->tolerance);

  niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_process */


int main(int argc,char *argv[],char *envp[]) {
  niikcortex_deform *dfm=NULL;
  const char *fcname="niikcortex_remesh";
  int n=0,i=0,nc=1,sc=1;
  double *ijk;
  /*niikmat *afmat=NULL;*/
  nifti_image *debugimg;

  if(argc==1) {
    usage();
    exit(0);
  }
  char* timestamp=niik_create_minc_timestamp(argc,argv);
  
#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif


  NIIK_EXIT(((dfm=niikcortex_deform_init())==NULL),fcname,"niikcortex_deform_init",1);
  niik_version_display(fcname, NIIK_MAJOR_VERSION, NIIK_MINOR_VERSION, NIIK_MICRO_VERSION);
  niik_fc_display(fcname,1);

  ijk=(double *)calloc(9,sizeof(double));
  ijk[0]=0;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        prog_history();
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-debug-keep-tmp",15)) {
        dfm->debug_keep_tmp=1;
      }

      else if(!strncmp(argv[nc],"-apply",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->apply_step = atof(argv[nc]);
      } /* delta */

      else if(!strncmp(argv[nc],"-iter2",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->iter2 = atoi(argv[nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-iter",5)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->iter = atoi(argv[nc]);
      } /* iter */

      else if(!strcmp(argv[nc],"-t1w")) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        NIIK_EXIT(((dfm->t1img = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      } /* t1w */

      else if(!strcmp(argv[nc],"-lambda")) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -lambda <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][0] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -lambda <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][0] = atof(argv[nc]);
      }

      else if(!strcmp(argv[nc],"-mju")) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -mju <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][1] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -mju <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][1] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-depth",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->bbox_depth = atoi(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        usage();
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

  if(argc<5) {
    fprintf(stderr,"[%s] ERROR: too few argments (%i != 5)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  } else if(argc>5) {
    fprintf(stderr,"[%s] ERROR: too many argments (%i != 5)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  }
  dfm->debug_pt=ijk;

  fprintf(stdout,"[%s] reading init ics object  %s\n",fcname,argv[1]);
  NIIK_RET1(((dfm->ctx[0]=off_kobj_read_offply(argv[1]))==NULL),fcname,"off_kobj_read_off");
  fprintf(stdout,"[%s] reading init ocs object  %s\n",fcname,argv[2]);
  NIIK_RET1(((dfm->ctx[1]=off_kobj_read_offply(argv[2]))==NULL),fcname,"off_kobj_read_off");

  /* dt */
  dfm->weight->m[CORTEX_ICS][WEIGHT_DT] =
    dfm->weight->m[CORTEX_OCS][WEIGHT_DT] =
      dfm->delta;


  /*Make sure outputs are flushed before long process starts*/
  fflush(stdout);
  fflush(stderr);

  /*display all paramters*/
  niikcortex_remesh_process(dfm);

  /* cortex deformation (niikcortex-deform) */
  NIIK_EXIT((!niikcortex_remesh_cortex(dfm)),
            fcname,"niikcortex_remesh_cortex",9);

  /*append metadata*/
  NIIK_RET1((!off_kobj_add_comment(dfm->ctx[0],timestamp)),fcname,"off_kobj_add_comment");
  NIIK_RET1((!off_kobj_add_comment(dfm->ctx[1],timestamp)),fcname,"off_kobj_add_comment");

  /* writing output for white surface */
  /*fprintf(stdout,"[%s] yellow color for white surface\n",fcname);*/
  /*NIIK_RET1((!off_kobj_add_one_color(dfm->ctx[0],0.8,0.8,0)),fcname,"off_kobj_add_one_color");*/
  fprintf(stdout,"[%s] writing white surface:  %s\n",fcname,argv[3]);
  NIIK_RET1((!off_kobj_write_offply(argv[3],dfm->ctx[0],0)),fcname,"off_kobj_write_off");

  /* writing output for pial surface */
  /*fprintf(stdout,"[%s] red color for pial surface\n",fcname);*/
  /*NIIK_RET1((!off_kobj_add_one_color(dfm->ctx[1],1.0,0.2,0.2)),fcname,"off_kobj_add_one_color");*/
  fprintf(stdout,"[%s] writing pial surface:   %s\n",fcname,argv[4]);
  NIIK_RET1((!off_kobj_write_offply(argv[4],dfm->ctx[1],0)),fcname,"off_kobj_write_off");

  dfm=niikcortex_deform_free(dfm);
  free(ijk);
  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
