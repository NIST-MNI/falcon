/* FILENAME:     niikcortex_shrinkwrap.c
* DESCRIPTION:  Kunio's nifti1 for cortical deformation
* AUTHOR:       Kunio Nakamura
* DATE:         May 8, 2014
*/

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>


#include "falcon.h"
#include "falcon_cortex.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (1)


#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


void prog_version() {
  fprintf(stdout,"  falcon_cortex_shrinkwrap history\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  0.0.1  Sep 14, 2016, Vladimir S. FONOV - moved from niikmath\n");
  fprintf(stdout,"\n");
}

void prog_usage() {
  fprintf(stdout,"falcon_cortex_shrinkwrap\n");
  fprintf(stdout,"  usage: [options] [mask] <out.off> \n");
  fprintf(stdout,"\n");
  fprintf(stdout,"   Shrink-wrap a sphere around a mask, if mask is absent - just produce a sphere of given radius\n");
  fprintf(stdout,"\n  optional usage:\n");
  fprintf(stdout,"  --elen <E>         target edge length [default=2mm]\n");
  fprintf(stdout,"  --start <E>        starting edge length [default=3mm]\n");
  fprintf(stdout,"  --obj <obj.off>    initial object [default=auto]\n");
  fprintf(stdout,"  --radius <R>       initial object's radius [default=estimage from image, if present or 10.0]\n");
  fprintf(stdout,"  --iter <I>         iteration [default=20]\n");
  fprintf(stdout,"  --smooth <d>       Surface smoothing coeff [default=0.01]\n");
  fprintf(stdout,"  --step  <max_step_size> - maximum shrink step length , default 2 \n");
  fprintf(stdout,"  --val   <max_step_size> - same as --step \n");
  fprintf(stdout,"  --label <n> use this label in the mask file, default  any above 0 \n");
  fprintf(stdout,"                     maximum shrink size for each iteration [default=2mm]\n");
  fprintf(stdout,"  --phase=<phase>    change phase (angle) for spherical coordinate [default=0]\n");
  fprintf(stdout,"  --version          show version info\n");
  fprintf(stdout,"  --clobber          overwrite output file\n");
  fprintf(stdout,"  --noremesh         disable remeshing, to keep constant number of nodes\n");
  fprintf(stdout,"  --debug            Added debugging output\n");
  fprintf(stdout,"  --progress         Display progress\n");

}

int main(int argc,char *argv[],char *envp[]) {
  int clobber=0;
  int    iter=20;
  nifti_image *maskimg=NULL;
  kobj *obj=NULL;
  const char *in_mask=NULL;
  const char *in_obj=NULL;
  const char *out_obj=NULL;
  const char *fcname="falcon_cortex_shrinkwrap";
  int i,n;
  double elen=2.0;
  double start_elen=3.0;
  double radius=NIIKMAX;
  double steplen=2.0;
  double phase=0.0;
  double smooth=0.01;
  niikpt  ctr;
  int flag_remesh=1;
  int flag_bbox=1;
  int flag_debug=0;
  int flag_progress=0;
  int label=-1;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument,  &clobber, 1},
    {"debug",   no_argument,  &flag_debug, 1},
    {"noremesh", no_argument, &flag_remesh, 0},
    {"progress", no_argument, &flag_progress, 1},
    {"version", no_argument, 0, 'v'},
    {"help",    no_argument, 0, 'h'},
    {"usage",   no_argument, 0, 'U'},

    {"elen",required_argument,0,     'e'},
    {"start",required_argument,0,    's'},
    {"obj",required_argument,0,      'o'},
    {"radius",required_argument,  0, 'R'},
    {"iter",required_argument,  0,   'I'},
    {"val",required_argument, 0,     'V'},
    {"step",required_argument, 0,    'V'},
    {"phase",required_argument, 0,   'P'},
    {"label",required_argument,  0,  'l'},
    {"smooth",required_argument,0,   'S'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vUhe:o:R:I:V:P:S:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      prog_version();
      return 0;
    case 'e':
      elen=atof(optarg);
      break;
    case 's':
      start_elen=atof(optarg);
      break;
    case 'o':
      in_obj=optarg;
      break;
    case 'R':
      radius=atof(optarg);
      break;
    case 'I':
      iter=atoi(optarg);
      break;
    case 'V':
      steplen=atof(optarg);
      break;
    case 'P':
      phase=atof(optarg);
      break;
    case 'l':
      label=atoi(optarg);
      break;
    case 'S':
      smooth=atof(optarg);
      break;
    case 'h':
    case 'U':
    case '?':
    default:
      prog_usage ();
      return 1;
    }
  }

  if((argc - optind)<1) {
    prog_usage();
    return 1;
  }

  if((argc - optind)>1) {
    in_mask = argv[optind];
    out_obj = argv[optind+1];
  } else {
    in_mask = NULL;
    out_obj = argv[optind];
  }

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);

  if (!clobber && !access (out_obj, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_obj);
    return 1;
  }

  if(in_mask) {
      if( (maskimg=niik_image_read(in_mask))==NULL) {
        fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,in_mask);
        exit(1);
      }
  } else {
    maskimg=NULL;
  }

  if(flag_progress)
  {
    if(!niik_verbose()) niik_set_verbose(1);
  }


  if(label!=-1 && maskimg) {
    double* dimg =NULL;
    int i;

    if((dimg = niik_image_get_voxels_as_double_vector( maskimg ))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_double_vector\n",fcname);
      exit(1);
    }

    for(i=0; i<maskimg->nvox; i++) {
      if( fabs(dimg[i]-label)<0.5 )
        dimg[i]=1.0;
      else dimg[i]=0.0;
    }
    if(!niik_image_set_voxels_from_double_vector(maskimg, dimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_set_voxels_from_double_vector\n",fcname);
      exit(1);
    }
    free(dimg);

    if(flag_debug) {
      fprintf(stdout,"[%s] [DEBUG] Writing itnermediate file to tmp_label.mnc\n",fcname);
      niik_image_write("tmp_label.mnc",maskimg);
    }
  }


#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  /*from niikmath*/

  if(niik_check_double_problem(steplen)) {
    steplen = 2.0;
  }
  if(niik_check_double_problem(phase)) {
    phase = 0.0;
  }

  if(in_obj==NULL) {
    kvert *v;
    fprintf(stdout,"[%s] creating a new spherical object\n",fcname);
    if(maskimg)
      ctr=niikpt_image_get_centroid(maskimg,NULL);
    else
      ctr=niikpt_zero(); /*place it at the origin*/

    if(niik_check_double_problem(radius))
    {
      if(maskimg)
        radius = niikpt_image_get_bounding_sphere(maskimg,1.0,ctr);
      else 
        radius = 10.0;
    }

    fprintf(stdout,"[%s]   centroid = %.3f %.3f %.3f radius = %.3f\n",fcname,ctr.x,ctr.y,ctr.z,radius);
    NIIK_EXIT(((obj=off_make_sphere_from_icosahedron(start_elen,radius,ctr))==NULL),fcname,"off_make_sphere_from_icosahedron",1);

    if(flag_debug) {
      fprintf(stdout,"[%s]  [DEBUG] writing initial object to tmp_init.ply\n",fcname);
      NIIK_EXIT((!off_kobj_write_offply("tmp_init.ply",obj,0)),fcname,"off_kobj_write_off",1);
    }

    /* phase shift */
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->sph.psi += NIIK_DEGREE2RAD(phase);
      v->sph = niiksph_norm(v->sph);
    }
    off_kobj_remove_color(obj);
  } else {
    if((obj=off_kobj_read_offply(in_obj))==NULL) {
      fprintf(stderr,"[%s] ERROR: off_kobj_read_off\n",fcname);
      exit(1);
    }
  }

  if(maskimg && iter>0 )
  {
    if(niik_check_double_problem(elen)) elen=2;

    /*
      * run shrink wrap
      */
    fprintf(stdout,"[%s] shrink wrap\n",fcname);
    NIIK_EXIT((!off_shrinkwrap_kobj_bbox_remesh(maskimg, obj, elen, iter, steplen,
              flag_bbox, flag_remesh, flag_debug, smooth )),fcname,"off_shrinkwrap_kobj_bbox_remesh",1);

  } 
  NIIK_RET1((!off_kobj_add_comment(obj,timestamp)),fcname,"off_kobj_add_comment");

  /*
    * write output
    */
  fprintf(stdout,"[%s] writing output  %s\n", fcname, out_obj);
  NIIK_EXIT((!off_kobj_write_offply(out_obj,obj,0)),fcname,"off_kobj_write_offply",1);

  obj=off_kobj_free(obj);
  if(maskimg)
    maskimg=niik_image_free(maskimg);

  niik_fc_display(fcname,0);
  free(timestamp);
  exit(0);
} /* main */

/*
kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
