#include <float.h>
#include <math.h>

#include "volume_tools.h"
#include "falcon.h"
#include "falcon_surfaces.h"
#include <getopt.h>
#include <unistd.h>


void usage() {
  fprintf(stdout,"Extrasact measurements along surface from a volume\n");
  fprintf(stdout,"  usage: [options] <volume.mnc/nii> <surface.ply/off> <out.ply/txt> \n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t--outside_value <f>- Specify outside value (value to be assigned if surface is outside volume or outside mask)\n");
  fprintf(stdout,"\t--mask <mask> - provide masking volume \n");
  fprintf(stdout,"\t--smooth <fwhm> - apply smoothing kernel (default: disable) - will try to apply it iteratively\n");
  fprintf(stdout,"\t--inside <mask> - value insude mask , default 0.0\n");
  fprintf(stdout,"\t--cubic - Cubic interpolation (default is linear)\n");
  fprintf(stdout,"\t--nearest_neighbor - Nearest neighbor interpolation (default is linear)\n");
  fprintf(stdout,"\t--labels  - Use integer values\n");
  fprintf(stdout,"\t--column <name> - save under this column\n");
  fprintf(stdout,"\t--ply - Save as ply file\n");
  fprintf(stdout,"\t--clobber - overwrite output file\n");
  fprintf(stdout,"\t--verbose - Be more verbose\n");
}

int main(int   argc, char  *argv[] ) {
  VIO_Volume     vol;
  VIO_Volume     vol_mask;
  VIO_Real outside_value=-1.0;

  int cubic=FALSE;
  int nearest_neighbor=FALSE;
  int use_labels=0;
  int clobber=0;
  int verbose=0;
  const char *colname="signal";
  const char *maskfn=NULL;

  kvert *v;
  const char *fcname="falcon_surface_signals";

  char *volume_file, *surface_file,*out_file;
  int i,degrees_continuity=0;
  double values[1];
  int out_ply=0;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  int * measurements_int=NULL;
  double *measurements_double=NULL;

  kobj *obj=NULL;
  int n_meas=0;
  char **meas_name;
  double **meas;

  unsigned char *mask_flag=NULL;
  double *mask_dbl=NULL;
  double smooth=0.0;
  double out_mask_val=0.0;
  double in_mask_val=0.0;

  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"verbose",          no_argument, &verbose, 1},
    {"ply",              no_argument, &out_ply, 1},
    {"labels",           no_argument, &use_labels, 1},
    {"cubic",            no_argument, &cubic, 1},
    {"nearest_neighbor", no_argument, &nearest_neighbor, 1},
    {"outside_value",    required_argument, 0, 'o'},
    {"column",           required_argument, 0, 'c'},
    {"smooth",           required_argument, 0, 's'},
    {"mask",       required_argument, 0, 'm'},
    {"inside",     required_argument, 0, 'I'},
    {"maskval",    required_argument, 0, 'M'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "uhvs:c:o:m:I:M:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,0,0,1);
      return 1;
    case 'o':
      outside_value=atof(optarg);
      break;
    case 'c':
      colname=optarg;
      break;
    case 's':
      smooth=atof(optarg);
      break;
    case 'm':
      maskfn=optarg;
      break;
    case 'M':
      out_mask_val=atof(optarg);
      break;
    case 'I':
      in_mask_val=atof(optarg);
      break;
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

  volume_file  = argv[optind];
  surface_file = argv[optind+1];
  out_file     = argv[optind+2];

  if (cubic && nearest_neighbor) {
    fprintf(stderr,"Specify either -cubic OR -nearest_neighbor.\n");
    return STATUS_ERR;
  }

  if(verbose)
    niik_fc_display(fcname,1);
 
  NIIK_EXIT(((obj=off_kobj_read_offply_ex(surface_file,&n_meas, &meas_name, &meas))==NULL),fcname,"niik_kobj_read_off",1);
  NIIK_RETURN(((mask_flag=(unsigned char *)calloc(obj->nvert,sizeof(unsigned char)))==NULL),"calloc for mask array",1);

  if(maskfn)
  {
      NIIK_RETURN((input_volume(  maskfn, 3, NULL, NC_DOUBLE, FALSE,
                            0.0, FLT_MAX, TRUE, &vol_mask,
                            (minc_input_options *) NULL ) != VIO_OK),"input_volume failed for mask",1);
      
  }

  NIIK_RETURN((input_volume(  volume_file, 3, NULL, NC_DOUBLE, FALSE,
                      0.0, FLT_MAX, TRUE, &vol,
                      (minc_input_options *) NULL ) != VIO_OK),"input_volume failed",1);

  if (cubic)
    degrees_continuity = 2;

  if (nearest_neighbor)
    degrees_continuity = -1;

  if(use_labels) {
    NIIK_RETURN(((measurements_int=(int*)calloc(obj->nvert,sizeof(int)))==NULL),"calloc for measurements array",1);

    for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
      evaluate_volume_in_world(vol,v->v.x,v->v.y,v->v.z,degrees_continuity,TRUE,outside_value,values,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      measurements_int[i] = (int)floor(values[0]+0.5);
    }
  } else {
    NIIK_RETURN(((measurements_double=(double*)calloc(obj->nvert,sizeof(double)))==NULL),"calloc for measurements array",1);

    for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
      evaluate_volume_in_world(vol,v->v.x,v->v.y,v->v.z,degrees_continuity,TRUE,outside_value,values,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      measurements_double[i] = (double)values[0];
    }
  }

  if(maskfn)
  {
    for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
      evaluate_volume_in_world(vol_mask,v->v.x,v->v.y,v->v.z,degrees_continuity,TRUE,in_mask_val+1.0,values,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      mask_flag[i] = fabs((double)values[0] - in_mask_val)<0.5?1:0;

      if(mask_flag[i] == 0) {
        if(use_labels)
          measurements_int[i]=(int)out_mask_val;
        else
          measurements_double[i]=out_mask_val;
      }
    }
  } else {
    /* TODO: evaluate if inside the bounding box of the volume!*/
    memset(mask_flag,1,obj->nvert);
  }

  /* smoothing*/
  if(smooth>0.0)
  {
    if(use_labels)
    {
      fprintf(stderr,"[%s] Can't sample labels and smooth them\n",fcname);
      return 1;
    }

    double it_smooth,it_sigma;
    double mean_elen=off_get_kobj_mean_edge_length(obj);

    /*apply small kernel several times, since each iteration only works with the local neighbourhood*/
    /*HACK: mean_elen*2  - is an arbitrary factor */
    int iterations=ceil((smooth*smooth)/(mean_elen*mean_elen*4));
    int i;
    if(iterations<1) iterations=1;
    it_smooth=sqrt(smooth*smooth/iterations);
    it_sigma=it_smooth/2.355; 

    if(verbose)
      fprintf(stdout,"[%s] Smoothing: %f x %d\n",fcname, it_smooth, iterations);

    if(maskfn)
    {
      if(!off_surface_gauss_smooth_using_vert_with_mask(obj, measurements_double, it_sigma, iterations, mask_flag)) {
        fprintf(stderr,"[%s] ERROR: off_surface_gauss_smooth_using_vert\n",fcname);
        return 1;
      }
    } else {
      if(!off_surface_gauss_smooth_using_vert(obj, measurements_double, it_sigma, iterations)) {
        fprintf(stderr,"[%s] ERROR: off_surface_gauss_smooth_using_vert\n",fcname);
        return 1;
      }
    }
  }

  if(out_ply)
  {
    int i;
    /*preserver measurements if available*/
    /*check if this column is already there*/
    for(i=0;i<n_meas;i++)
    {
      if(!strcmp(colname, meas_name[i])) break;
    }
    if(i == n_meas) {
      n_meas+=2;

      meas_name=(char**)realloc(meas_name,(n_meas)*sizeof(char*));
      meas=(double**)realloc(meas,(n_meas)*sizeof(double*));

      meas_name[i] = strdup(colname);
      meas_name[i+1] = strdup("defined");
      meas[i] = (double*)calloc(obj->nvert,sizeof(double));
      meas[i+1] = (double*)calloc(obj->nvert,sizeof(double));
    }
    if(measurements_double) {
      memmove(meas[i], measurements_double, sizeof(double)*obj->nvert);
      free(measurements_double);
    }
    else 
    {
      int j;
      for(j=0;j<obj->nvert;j++) meas[i][j]=(double)measurements_int[j];
      free(measurements_int);
    }
    int j;
    for(j=0;j<obj->nvert;j++) meas[i+1][j]=(double)mask_flag[j];

    NIIK_RETURN((!off_kobj_write_ply_ex(out_file,obj,0,1,1,1, n_meas, meas_name, meas)),"off_kobj_write_ply_ex",1);

  } else {
    int j;
    /*just save to text file*/
    if(measurements_double)
    {
      const char * headers[]={colname, "defined"};
      double *vectors[2];
      
      vectors[0]=measurements_double;
      NIIK_RETURN(((vectors[1]=(double *)malloc(obj->nvert*sizeof(double)))==NULL),"malloc for mask_dbl array",1);
      for(j=0;j<obj->nvert;j++) vectors[1][j]=(double)mask_flag[j];

      NIIK_RETURN((!niik_write_double_vectors_ex(out_file,vectors,obj->nvert,2, headers)),"niik_write_double_vectors_ex",1);
      free(measurements_double);
      free(vectors[1]);
    }
    else 
    {
      const char * headers[]={colname, "defined"};
      int *vectors[2];
      
      vectors[0]=measurements_int;
      NIIK_RETURN(((vectors[1]=(int *)malloc(obj->nvert*sizeof(int)))==NULL),"malloc for mask_dbl array",1);
      for(j=0;j<obj->nvert;j++) vectors[1][j]=(int)mask_flag[j];

      NIIK_RETURN((!niik_write_int_vectors_ex(out_file,vectors,obj->nvert,2, headers)),"niik_write_int_vectors_ex",1);
      free(measurements_int);
      free(vectors[1]);
    }
  }
  
  if(n_meas>0)
  {
    for(i=0;i<n_meas;i++)
    {
      if(meas[i])      free(meas[i]);
      if(meas_name[i]) free(meas_name[i]);
    }
    free(meas);free(meas_name);
  }
  off_kobj_free(obj);
  delete_volume(vol);

  if(verbose) niik_fc_display(fcname,0);

  free(timestamp);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
