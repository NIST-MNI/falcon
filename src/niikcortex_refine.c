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


void prog_history() {
  fprintf(stdout,"[niikcortex_refine] history\n");
  fprintf(stdout,"  version  0.0.0  March 24, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n");
  fprintf(stdout,"  version  0.0.1  September 13, 2018, Vladimir S. FONOV <vladimir.fonov@gmail.com>\n");
}

void usage() {
  fprintf(stdout,"niikcortex_refine\n");
  fprintf(stdout,"  usage: [options] <img> <brain_mask> <csf_mask> <gwi_mask> <avoid mask> <white.off> <pial.off> <out_white.off> <out_pial.off\n\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help                   : show this usage\n");
  fprintf(stdout,"  --version                         : show version info\n");
  fprintf(stdout,"  -debug-keep-tmp                   : keep debug files (obsolete)\n");
  fprintf(stdout,"  -log <file>                       : write convergence log to file\n");
  fprintf(stdout,"  PDE Solver:                       \n");
  fprintf(stdout,"  -backward-euler                   : use backward euler method\n");
  fprintf(stdout,"  -forward-euler                    : use forward euler method (DEFAULT)\n");
  fprintf(stdout,"  -rungekutta                       : use runge-kutta method\n");
  fprintf(stdout,"  -midpoint                         : use midpoint method\n");
  fprintf(stdout,"  Processing options:               \n");
  fprintf(stdout,"  -lesion-mask   <img>              : lesion mask\n");
  fprintf(stdout,"  -white-only                       : deform white surface only\n");
  fprintf(stdout,"  -pial-only                        : deform pial surface only\n");
  fprintf(stdout,"  -both                             : deform both surfaces (default)\n");
  fprintf(stdout,"  -nonctx-mask   <img>              : mask away non-cortex\n");
  fprintf(stdout,"  -priorwm      <img>               : priors for WM\n");
  fprintf(stdout,"  -priorgm      <img>               : priors for GM\n");
  fprintf(stdout,"  -priorcsf     <img>               : priors for CSF\n");
  fprintf(stdout,"  -nomf                             : Don't use multiplication factor (mf)\n");
  fprintf(stdout,"  Optimization weights \n");
  fprintf(stdout,"  -wprox   <white> <pial>           : proximity weights, default 0.2 0.2 \n");
  fprintf(stdout,"  -wgrad   <white> <pial>           : gradient term , default 0.1 0.1 \n");
  fprintf(stdout,"  -wimag   <white> <pial>           : image weights, default 1.0 1.0 \n");
  fprintf(stdout,"  -wprior  <white> <pial>           : prior weights, default 0.0 \n");
  fprintf(stdout,"  -wsmooth <white> <pial>           : tangential surface smoothness weights, default 1.0 1.0 \n");
  fprintf(stdout,"  -wtsmooth <white> <pial>          : thickness smoothness weights, default 0.0 0.0 \n");
  fprintf(stdout,"  -wssmooth <white> <pial>          : simple surface smoothness weights, default 0.0 0.0 \n");
  fprintf(stdout,"  -wbrain  <white> <pial>           : brain mask weights, default 1.0 1.0 \n");
  fprintf(stdout,"  -wvent   <white> <pial>           : ventricles mask weights, default 1.5 1.5 \n");
  fprintf(stdout,"  -wles    <white> <pial>           : lesion mask weights, default 1.5 1.5 \n");
  fprintf(stdout,"  -wabs    <white> <pial>           : absolute thickness weight 0.0 0.0 (not yet implemented) \n");
  fprintf(stdout,"  -wcereb  <white> <pial>           : cerebellum mask weights, default 1.0 1.0 \n");
  fprintf(stdout,"  -wcurv   <white> <pial>           : curvature weights, default 0.0 0.0 \n");
  fprintf(stdout,"  -wflux  <white> <pial>            : flux term , default 0.0 0.0  \n");
  fprintf(stdout,"  Cortical thickness soft constraints \n");
  fprintf(stdout,"  -tmin <white> <pial>              : minimum thickness, default 2.0 2.0 \n");
  fprintf(stdout,"  -tmax <white> <pial>              : maximum thickness, default 5.0 5.0 \n");
  fprintf(stdout,"  -tsigma <white> <pial>            : thickness sigma (bandwidth), default 1.0 1.0 \n");
  fprintf(stdout,"  Update smoothing \n");
  fprintf(stdout,"  -supdate <white> <pial>           : update smoothing sigma, default 0.0 0.0 \n");
  fprintf(stdout,"  Proximity distance constraints \n");
  fprintf(stdout,"  -pmin <val>                       : minimum proximity distance, default 0.6 \n");
  fprintf(stdout,"  Tissue weights \n");
  fprintf(stdout,"  -wmix  <white> <pial>             : tissue mix for ICS and OCS, default 0.5 0.5 \n");

  fprintf(stdout,"  Additional optimizer parameters \n");
  fprintf(stdout,"  -depth <n>                        : Quad-tree depth , default 7\n");
  fprintf(stdout,"  -gradient-FWHM <f> <f>            : blurring kernels for image gradient calculation (voxels), default 1.0 1.0\n");
  fprintf(stdout,"  -divergence-FWHM <f> <f>          : blurring kernels for image divergence calculation (voxels), default 1.0 1.0\n");
  fprintf(stdout,"  -prior-FWHM <f>  <f>              : blurring kernels for prior gradient calculation (voxels), default 1.0 1.0\n");
  fprintf(stdout,"  -delta <f>                        : time-step (default 0.5)\n");
  fprintf(stdout,"  -apply <f>                        : apply-step (default 0.2)\n");
  fprintf(stdout,"  -iter  <n>                        : maximum number of iterations, default 100\n");
  fprintf(stdout,"  -iter2 <n>                        : maximum number of sub-iterations, default 5\n");
  fprintf(stdout,"  -remesh <n>                       : remesh every nth iteration, default 25, 0 to disable\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  Environment variables: \n");
  fprintf(stdout,"   FALCON_TRACE - debug tracing prefix\n");
  fprintf(stdout,"     output will be named  <prefix>_<stage>_<iteration>_<object>.ply - for surface \n");
  fprintf(stdout,"     or  <prefix>_<stage>_<iteration>_<x|y|z>_<slice>.tiff - for volume slice \n");
  fprintf(stdout,"   FALCON_TRACE_SCALE - debug tracing image scale, default 1.0 \n");
  fprintf(stdout,"   FALCON_TRACE_X,FALCON_TRACE_Y,FALCON_TRACE_Z - (n1)[,n2[,n3]] - slices to trace (voxel coordinates) \n");
  fprintf(stdout,"   FALCON_TRACE_SURF - set to any value to dump surfaces \n");
}




int niikcortex_deform_process(niikcortex_deform *dfm) {
  const char *fcname="niikcortex_deform_process";
  niik_fc_display(fcname,1);

  NIIK_RET0((dfm==NULL),fcname,"missing self");

  NIIK_RET0((dfm->ctx[0]->nvert!=dfm->ctx[1]->nvert),fcname,"#vert did not match");
  fprintf(stdout,"[%s] parameters\n",fcname);
  fprintf(stdout,"  t1w image            %s\n",dfm->t1img->fname);
  fprintf(stdout,"  brain mask           %s\n",dfm->brain_mask->fname);
  fprintf(stdout,"  white surface        %s\n",dfm->ctx[0]->fname);
  fprintf(stdout,"  pial surface         %s\n",dfm->ctx[1]->fname);

  if(dfm->prior[0])
    fprintf(stdout,"  WM prior             %s\n", dfm->prior[0]->fname);
  if(dfm->prior[1])
    fprintf(stdout,"  GM prior             %s\n", dfm->prior[1]->fname);
  if(dfm->prior[3])
    fprintf(stdout,"  CSF prior            %s\n", dfm->prior[3]->fname);

  fprintf(stdout,"  surface vfe          %i %i %i\n",dfm->ctx[1]->nvert,dfm->ctx[1]->nface,dfm->ctx[1]->nedge);
  fprintf(stdout,"  deform time step     %-7.4f\n",dfm->delta);
  fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n",dfm->apply_step);
  fprintf(stdout,"  surface sm weights   %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_SURFACE],dfm->weight->m[1][WEIGHT_SURFACE]);
  fprintf(stdout,"  image weights        %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_IMAGE],dfm->weight->m[1][WEIGHT_IMAGE]);
  fprintf(stdout,"  prior weight         %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_PRIOR],dfm->weight->m[1][WEIGHT_PRIOR]);
  fprintf(stdout,"  thickness sm weights %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_THICKNESS_SMOOTHNESS],dfm->weight->m[1][WEIGHT_THICKNESS_SMOOTHNESS]);
  fprintf(stdout,"  brainmask weights    %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_BRAIN_MASK],dfm->weight->m[1][WEIGHT_BRAIN_MASK]);
  fprintf(stdout,"  ventricle weights    %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_VENTRICLE],dfm->weight->m[1][WEIGHT_VENTRICLE]);
  fprintf(stdout,"  cerebellum weights   %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_AVOID],dfm->weight->m[1][WEIGHT_AVOID]);
  fprintf(stdout,"  lesion mask weights  %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_LESION],dfm->weight->m[1][WEIGHT_LESION]);
  fprintf(stdout,"  proximity weights    %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_PROXIMITY],dfm->weight->m[1][WEIGHT_PROXIMITY]);
  fprintf(stdout,"  abs thickness        %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_ABS_THICKINESS],dfm->weight->m[1][WEIGHT_ABS_THICKINESS]);
  fprintf(stdout,"  gradient weights     %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_GRADIENT], dfm->weight->m[1][WEIGHT_GRADIENT]);
  fprintf(stdout,"  curvature weights    %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_CURVATURE],dfm->weight->m[1][WEIGHT_CURVATURE]);
  fprintf(stdout,"  flux weights         %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_DGRADIENT],dfm->weight->m[1][WEIGHT_DGRADIENT]);
  fprintf(stdout,"  thickness min        %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_MIN_THICKNESS],dfm->weight->m[1][WEIGHT_MIN_THICKNESS]);
  fprintf(stdout,"  thickness max        %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_MAX_THICKNESS],dfm->weight->m[1][WEIGHT_MAX_THICKNESS]);
  fprintf(stdout,"  thickness sigma      %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_THICKNESS_SIGMA],dfm->weight->m[1][WEIGHT_THICKNESS_SIGMA]);
  fprintf(stdout,"  update sigma         %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_UPDATE_SIGMA],dfm->weight->m[1][WEIGHT_UPDATE_SIGMA]);
  fprintf(stdout,"  mix weight           %-7.4f %-7.4f\n",dfm->weight->m[0][WEIGHT_MIX],dfm->weight->m[0][WEIGHT_MIX]);

  fprintf(stdout,"  use mf               %s\n",dfm->use_mf?"Yes":"No");
  fprintf(stdout,"  max iter             %i\n",dfm->iter);
  fprintf(stdout,"  max iter2            %i\n",dfm->iter2);
  fprintf(stdout,"  remesh               %i\n",dfm->remesh);
  fprintf(stdout,"  white matter         %9.2f\n",dfm->tissue_val[1]);
  fprintf(stdout,"  white surface value  %9.2f %9.2f\n",dfm->tissue_val[4],dfm->tissue_val[6]);
  fprintf(stdout,"  pial surface value   %9.2f %9.2f\n",dfm->tissue_val[5],dfm->tissue_val[7]);
  fprintf(stdout,"  proximity min dist   %-7.2f\n",dfm->proximity_min_distance);
  fprintf(stdout,"  quadtree depth       %i\n",dfm->bbox_depth);
  fprintf(stdout,"  gradient FWHM        %-7.4f %-7.4f\n", dfm->gradient_FWHM[0],dfm->gradient_FWHM[1]);
  fprintf(stdout,"  divergence FWHM      %-7.4f %-7.4f\n", dfm->divergence_FWHM[0],dfm->divergence_FWHM[1]);
  fprintf(stdout,"  prior FWHM           %-7.4f %-7.4f\n", dfm->prior_FWHM[0],dfm->prior_FWHM[1]);

  if(dfm->gradient_lambda[0]>=0) {
    fprintf(stdout,"  ics gradient lambda      %7.3f\n",dfm->gradient_lambda[0]);
  } else {
    fprintf(stdout,"  ics gradient lambda      TBD\n");
  }
  if(dfm->gradient_lambda[1]>=0) {
    fprintf(stdout,"  ocs gradient lambda      %7.3f\n",dfm->gradient_lambda[1]);
  } else {
    fprintf(stdout,"  ocs gradient lambda      TBD\n");
  }
  fprintf(stdout,"  tolerance            %-7.4f\n",dfm->tolerance);
  if(dfm->ctx_label!=NULL)
    fprintf(stdout,"  non-cortex label     %-i / %i\n",(int)niik_count_zero_from_int_vector(dfm->ctx_label,dfm->ctx[0]->nvert),dfm->ctx[0]->nvert);
  else
    fprintf(stdout,"  non-cortex label     not used\n");
  switch(dfm->cortex_id) {
  case 1:
    fprintf(stdout,"  deform cortex        white surface\n");
    break;
  case 2:
    fprintf(stdout,"  deform cortex        pial surface\n");
    break;
  case 3:
    fprintf(stdout,"  deform cortex        white and pial surfaces\n");
    break;
  default:
    fprintf(stderr,"[%s] ERROR: unknown cortex_id, %i\n",fcname,dfm->cortex_id);
    return 0;
  }

  fprintf(stdout,"  numerical method     %s\n",niik_numerical_method_string(dfm->numerical_method));

  niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_process */


int main(int argc,char *argv[],char *envp[]) {
  niikcortex_deform *dfm=NULL;
  const char *fcname="falcon_cortex_refine";
  int n=0,i=0,nc=1,sc=1;
  double *ijk;
  /*niikmat *afmat=NULL;*/
  nifti_image *debugimg;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  if(argc==1) {
    usage();
    exit(0);
  }

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

      else if(!strncmp(argv[nc],"-backward-euler",15)) {
        dfm->numerical_method=NIIK_NUM_METHOD_BACKWARD_EULER;
      } /* backward-euler method */

      else if(!strncmp(argv[nc],"-forward-euler",14)) {
        dfm->numerical_method=NIIK_NUM_METHOD_FORWARD_EULER;
      } /* forward-euler method */

      else if(!strncmp(argv[nc],"-midpoint",14)) {
        dfm->numerical_method=NIIK_NUM_METHOD_MIDPOINT;
      } /* forward-euler method */
      
      else if(!strncmp(argv[nc],"-rungekutta",11)) {
        dfm->numerical_method=NIIK_NUM_METHOD_RUNGE_KUTTA;
      } /* runge-kutta (4th order) */

      else if(!strncmp(argv[nc],"-debug-keep-tmp",15)) {
        dfm->debug_keep_tmp=1;
      }

      else if(!strncmp(argv[nc],"-nomf",5)) {
        dfm->use_mf=0;
      }
      
      else if(!strncmp(argv[nc],"-gradient-FWHM",13)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -gradient-FWHM\n",fcname);
          exit(1);
        }
        dfm->gradient_FWHM[0] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -gradient-FWHM\n",fcname);
          exit(1);
        }
        dfm->gradient_FWHM[1] = atof(argv[nc]);
      } /* gradient-FWHM */

      else if(!strncmp(argv[nc],"-divergence-FWHM",16)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -divergence-FWHM\n",fcname);
          exit(1);
        }
        dfm->divergence_FWHM[0] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -divergence-FWHM\n",fcname);
          exit(1);
        }
        dfm->divergence_FWHM[1] = atof(argv[nc]);
      } /* divergence-FWHM */

      else if(!strncmp(argv[nc],"-prior-FWHM",11)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -prior-FWHM\n",fcname);
          exit(1);
        }
        dfm->prior_FWHM[0] = atof(argv[nc]);
        
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -prior-FWHM\n",fcname);
          exit(1);
        }
        dfm->prior_FWHM[1] = atof(argv[nc]);
      } /* prior-FWHM */

      else if(!strncmp(argv[nc],"-lesion-mask",12)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading lesion    %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->lesion_mask = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-priorwm",8)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading prior   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->prior[0] = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-priorgm",8)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading prior   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->prior[1] = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-priorcsf",9)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading prior   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->prior[3] = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-white-only",11)) {
        dfm->cortex_id=1;
      } /* white-only */
      else if(!strncmp(argv[nc],"-pial-only",10)) {
        dfm->cortex_id=2;
      } /* pial-only */

      else if(!strncmp(argv[nc],"-delta",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->delta = atof(argv[nc]);
      } /* delta */

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

      else if(!strncmp(argv[nc],"-remesh",7)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->remesh = atoi(argv[nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-nonctx-mask",12)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        NIIK_EXIT(((dfm->nonctx_mask = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      } /* nonctx-mask */

      else if(!strncmp(argv[nc],"-wprox",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprox <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_PROXIMITY] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprox <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_PROXIMITY] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wgrad",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wgrad <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_GRADIENT] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wgrad <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_GRADIENT] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wflux",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wflux <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_DGRADIENT] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wflux <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_DGRADIENT] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wimag",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wimag <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_IMAGE] = atof(argv[nc]);
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wimag <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_IMAGE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wcurv",6)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wcurv <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_CURVATURE] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wcurv <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_CURVATURE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wsmooth",8)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wsmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_SURFACE] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wsmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_SURFACE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wtsmooth",9)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wtsmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_THICKNESS_SMOOTHNESS] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wtsmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_THICKNESS_SMOOTHNESS] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wssmooth",9)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wssmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_SURFACE_SIMPLE] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wssmooth <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_SURFACE_SIMPLE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wabs",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wabs <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_ABS_THICKINESS] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wabs <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_ABS_THICKINESS] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wbrain",7)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wbrain <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_BRAIN_MASK] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wbrain <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_BRAIN_MASK] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wvent",6)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wvent <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_VENTRICLE] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wvent <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_VENTRICLE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wles",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wles <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_LESION] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wles <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_LESION] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wprior",7)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprior <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_PRIOR] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprior <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_PRIOR] = atof(argv[nc]);
      }


      else if(!strncmp(argv[nc],"-wcereb",7)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wcereb <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_AVOID] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wcereb <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_AVOID] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-tmin",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tmin <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_MIN_THICKNESS] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tmin <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_MIN_THICKNESS] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-tmax",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tmax <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_MAX_THICKNESS] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tmax <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_MAX_THICKNESS] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-tsigma",7)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tsigma <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_THICKNESS_SIGMA] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -tsigma <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_THICKNESS_SIGMA] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-supdate",8)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -supdate <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_UPDATE_SIGMA] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -supdate <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_UPDATE_SIGMA] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wmix",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wmix <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_MIX] = atof(argv[nc]);
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wmix <white> <pial>\n",fcname);
          exit(1);
        }
        dfm->weight->m[1][WEIGHT_MIX] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-pmin",5)) {
        if((++nc)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -pmin <dist>\n",fcname);
          exit(1);
        }
        dfm->proximity_min_distance = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-both",5)) {
        dfm->cortex_id=3;
      } /* both surfaces */


      else if(!strncmp(argv[nc],"-depth",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->bbox_depth = atoi(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-log",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -log <fname>\n",fcname);
          exit(1);
        }
        dfm->convergence_log=fopen(argv[nc],"w");
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

  if(argc<10) {
    fprintf(stderr,"[%s] ERROR: too few argments (%i != 11)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  } else if(argc>10) {
    fprintf(stderr,"[%s] ERROR: too many argments (%i != 11)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  }
  dfm->debug_pt=ijk;

  fprintf(stdout,"[%s] reading t1w image        %s\n",fcname,argv[1]);
  NIIK_EXIT(((dfm->t1img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);

  fprintf(stdout,"[%s] reading brainmask image  %s\n",fcname,argv[2]);
  NIIK_EXIT(((dfm->brain_mask=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);

  fprintf(stdout,"[%s] reading CSF image        %s\n",fcname,argv[3]);
  NIIK_EXIT(((dfm->csf_mask=niik_image_read(argv[3]))==NULL),fcname,"niik_image_read",9);

  fprintf(stdout,"[%s] reading GWI image        %s\n",fcname,argv[4]);
  NIIK_EXIT(((dfm->gwi_mask=niik_image_read(argv[4]))==NULL),fcname,"niik_image_read",9);

  fprintf(stdout,"[%s] reading avoid image %s\n",fcname,argv[5]);
  NIIK_EXIT(((dfm->avoid_mask=niik_image_read(argv[5]))==NULL),fcname,"niik_image_read",9);

  fprintf(stdout,"[%s] reading init ics object  %s\n",fcname,argv[6]);
  NIIK_EXIT(((dfm->ctx[0]=off_kobj_read_offply(argv[6]))==NULL),fcname,"off_kobj_read_off",9);
  
  fprintf(stdout,"[%s] reading init ocs object  %s\n",fcname,argv[7]);
  NIIK_EXIT(((dfm->ctx[1]=off_kobj_read_offply(argv[7]))==NULL),fcname,"off_kobj_read_off",9);

  NIIK_EXIT((!niik_image_type_convert(dfm->brain_mask,  NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, brain_mask",1);
  NIIK_EXIT((!niik_image_type_convert(dfm->csf_mask,    NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, csf_mask",1);
  NIIK_EXIT((!niik_image_type_convert(dfm->gwi_mask,    NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, gwi_mask",1);
  NIIK_EXIT((!niik_image_type_convert(dfm->avoid_mask,  NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, avoid_mask",1);

  /*join priors for GM & WM into another one */
  if(dfm->prior[1] && dfm->prior[0]) {
    NIIK_EXIT(((dfm->prior[2] = niik_image_copy_as_type(dfm->prior[0], NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type",1);
    NIIK_EXIT((!niik_image_add_2_images(dfm->prior[2],dfm->prior[1])),fcname,"niik_image_add_2_images",1);

    // DEBUG
    // niik_image_write("debug_prior_0.mnc",dfm->prior[0]);
    // niik_image_write("debug_prior_1.mnc",dfm->prior[1]);
    // niik_image_write("debug_prior_2.mnc",dfm->prior[2]);
    // niik_image_write("debug_prior_3.mnc",dfm->prior[3]);
    // DEBUG
  }



  /* get the tissue values (niikcortex-deform) */
  if( dfm->tissue_val[0]>0 && dfm->tissue_val[1]>0 && dfm->tissue_val[2]>0 && dfm->tissue_val[3]>0 &&
      dfm->tissue_val[4]>0 && dfm->tissue_val[5]>0 && dfm->tissue_val[6]>0 && dfm->tissue_val[7]>0) {
    fprintf(stdout,"[%s] using predefined tissue intensities\n",fcname);
  } else {
    fprintf(stdout,"[%s] niikcortex_estimate_tissue_values\n",fcname);
    NIIK_EXIT((!niikcortex_estimate_tissue_values(dfm->t1img, dfm->brain_mask, dfm->csf_mask, dfm->gwi_mask,
               &dfm->tissue_val[2],  /*CSF*/
               &dfm->tissue_val[1],  /*WM*/
               &dfm->tissue_val[0],  /*GM*/
               &dfm->tissue_val[3],  /*BRAIN*/
               &dfm->tissue_val[4],  /*ICS*/
               &dfm->tissue_val[5],  /*OCS*/
               &dfm->tissue_val[6],  /*rangeICS*/
               &dfm->tissue_val[7],  /*rangeOCS*/
               dfm->weight->m[0][WEIGHT_MIX], dfm->weight->m[1][WEIGHT_MIX])),
              fcname,"niikcortex_estimate_tissue_values",9);
  }

  fprintf(stdout,"[%s] niikcortex_estimate_tissue_values:\n",fcname);
  fprintf(stdout,"           GM WM CSF Brain   %7.3f %7.3f %7.3f %7.3f\n",dfm->tissue_val[0],dfm->tissue_val[1],dfm->tissue_val[2],dfm->tissue_val[3]);
  fprintf(stdout,"           ICS surface       %7.3f (%7.3f)\n",dfm->tissue_val[4],dfm->tissue_val[6]);
  fprintf(stdout,"           OCS surface       %7.3f (%7.3f)\n",dfm->tissue_val[5],dfm->tissue_val[7]);

  /* dt */
  dfm->weight->m[CORTEX_ICS][WEIGHT_DT] =
    dfm->weight->m[CORTEX_OCS][WEIGHT_DT] =
      dfm->delta;


  /*
   DEBUG
   */
#if 0
  debugimg = niik_image_copy(dfm->t1img);
  niikcortex_make_fuzzy(dfm->t1img,debugimg,dfm->tissue_val[4],dfm->tissue_val[6]);
  niik_image_write("debug_ics_membership.mnc",debugimg);
  niikcortex_make_fuzzy(dfm->t1img,debugimg,dfm->tissue_val[5],dfm->tissue_val[7]);
  niik_image_write("debug_ocs_membership.mnc",debugimg);
  niik_image_free(debugimg);
#endif
  /*
   DEBUG
   */


  /*Make sure outputs are flushed before long process starts*/
  fflush(stdout);
  fflush(stderr);

  /*display all paramters*/
  niikcortex_deform_process(dfm);


  /* cortex deformation (niikcortex-deform) */
  NIIK_EXIT((!niikcortex_deform_cortex(dfm)),
            fcname,"niikcortex_deform_cortex",9);

  /*append metadata*/
  NIIK_RET1((!off_kobj_add_comment(dfm->ctx[0], timestamp)),fcname,"off_kobj_add_comment");
  NIIK_RET1((!off_kobj_add_comment(dfm->ctx[1], timestamp)),fcname,"off_kobj_add_comment");

  /* writing output for white surface */
  fprintf(stdout,"[%s] yellow color for white surface\n",fcname);
  NIIK_RET1((!off_kobj_add_one_color(dfm->ctx[0], 0.8, 0.8, 0)),fcname,"off_kobj_add_one_color");
  fprintf(stdout,"[%s] writing white surface:  %s\n",fcname,argv[8]);
  NIIK_RET1((!off_kobj_write_offply(argv[8],dfm->ctx[0], 0)),fcname,"off_kobj_write_off");

  /*FOR DEBUGGING ONLY*/
  /*
  {
    const char *meas_names[]={"mf"};
    const double *meas_vec[]={dfm->mf_list};
    NIIK_RET1((!off_kobj_write_ply_ex(argv[8],dfm->ctx[0], 0, 1, 1, 1, 1, meas_names, meas_vec  )),fcname,"off_kobj_write_off");
  }*/

  /* writing output for pial surface */
  fprintf(stdout,"[%s] red color for pial surface\n",fcname);
  NIIK_RET1((!off_kobj_add_one_color(dfm->ctx[1], 1.0, 0.2, 0.2)),fcname,"off_kobj_add_one_color");
  fprintf(stdout,"[%s] writing pial surface:   %s\n",fcname,argv[9]);
  NIIK_RET1((!off_kobj_write_offply(argv[9], dfm->ctx[1],0)),fcname,"off_kobj_write_off");

  if(dfm->convergence_log) 
    fclose(dfm->convergence_log);
  
  dfm=niikcortex_deform_free(dfm);
  free(ijk);
  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
