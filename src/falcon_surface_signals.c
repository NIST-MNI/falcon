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
  fprintf(stdout,"\t--outside_value - Specify outside value (value to be assigned if surface is outside volume)\n");
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
  VIO_Real outside_value=-1.0;

  int cubic=FALSE;
  int nearest_neighbor=FALSE;
  int use_labels=0;
  int clobber=0;
  int verbose=0;
  const char *colname="signal";

  kvert *v;
  const char *fcname="falcon_off_signals";

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

  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"verbose",          no_argument, &verbose, 1},
    {"ply",              no_argument, &out_ply, 1},
    {"labels",           no_argument, &use_labels, 1},
    {"cubic",            no_argument, &cubic, 1},
    {"nearest_neighbor", no_argument, &nearest_neighbor, 1},
    {"outside_value",    required_argument, 0, 'o'},
    {"column",           required_argument, 0, 'c'},
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
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,0,0,1);
      return 1;
    case 'o':
      outside_value=atof(optarg);
      break;
    case 'c':
      colname=optarg;
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

  if ( input_volume(  volume_file, 3, NULL, NC_FLOAT, FALSE,
                      0.0, FLT_MAX, TRUE, &vol,
                      (minc_input_options *) NULL ) != VIO_OK )
    return( 1 );

  if (cubic)
    degrees_continuity = 2;

  if (nearest_neighbor)
    degrees_continuity = -1;

  if(use_labels) {
    measurements_int=(int*)calloc(obj->nvert,sizeof(int));
    for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
      evaluate_volume_in_world(vol,v->v.x,v->v.y,v->v.z,degrees_continuity,TRUE,outside_value,values,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      measurements_int[i] = (int)floor(values[0]+0.5);
    }

    
  } else {
    measurements_double=(double*)calloc(obj->nvert,sizeof(double));

    for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
      evaluate_volume_in_world(vol,v->v.x,v->v.y,v->v.z,degrees_continuity,TRUE,outside_value,values,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      measurements_double[i] = (double)values[0];
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
    if(i==n_meas) {
      n_meas++;
      meas_name=(char**)realloc(meas_name,(n_meas)*sizeof(char*));
      meas=(double**)realloc(meas,(n_meas)*sizeof(double*));
      meas_name[i] = strdup(colname);
      meas[i] = (double*)calloc(obj->nvert,sizeof(double));
    }
    if(measurements_double) {
      memmove(meas[i], measurements_double,sizeof(double)*obj->nvert);
      free(measurements_double);
    }
    else 
    {
      int j;
      for(j=0;j<obj->nvert;j++) meas[i][j]=(double)measurements_int[j];
      free(measurements_int);
    }
    NIIK_RETURN((!off_kobj_write_ply_ex(out_file,obj,0,1,1,1, n_meas, meas_name, meas)),"off_kobj_write_ply_ex",1);

  } else {
    /*just save to text file*/
    if(measurements_double)
    {
      NIIK_RETURN((!niik_write_double_vector_ex(out_file,measurements_double,obj->nvert, colname)),"niik_write_double_vector",1);
      free(measurements_double);
    }
    else 
    {
      NIIK_RETURN((!niik_write_int_vector_ex(out_file,measurements_int,obj->nvert,colname)),"niik_write_int_vector",1);
      free(measurements_int);
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
