/* FILENAME:     off_split.c
 * DESCRIPTION:  Split off file into disconnected components
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_cortex.h"
/*#include  <volume_io.h>*/

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input.off> <scalar.txt> <output.off> apply colour map\n"
          "\t--clobber clobber output files\n"
          "\t--min <f>\n"
          "\t--max <f>\n"
          "\t--column <name> - column name to use\n"
          "\t--column_id <name> - column id to use\n"
          "\t--spectral - use spectral colour map (default)\n"
          "\t--atrophy  - use atrophy colour map\n"
          "\t--summer   - use summer colour map\n"
          "\t--jacobian - use jacobian colour map\n"
          "\t--gray     - use gray colour map\n"
          ,name);
}

int main(int argc, char **argv) {
  const char *fcname="off_colour";
  int clobber=0;
  int verbose=0;
  int c;
  int i;
  int cmap=NIIK_COLORMAP_SPECTRAL;

  double omin=0.0;
  double omax=0.0;

  const char *in_off=NULL;
  const char *in_txt=NULL;
  const char *out_off=NULL;

  double *meas_in=NULL;
  int meas_num_in;
  int column_id=0;
  int n;
  kobj *obj=NULL;
  kobj **out_obj;
  const char *column_name=NULL;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"min", required_argument, 0, 'i'},
    {"max", required_argument, 0, 'a'},
    {"column", required_argument,0,'c'},
    {"column_id", required_argument,0,'C'},
    {"spectral", no_argument,&cmap,NIIK_COLORMAP_SPECTRAL},
    {"atrophy", no_argument,&cmap,NIIK_COLORMAP_ATROPHY},
    {"summer", no_argument,&cmap,NIIK_COLORMAP_SUMMER},
    {"jacobian", no_argument,&cmap,NIIK_COLORMAP_JACOBIAN},
    {"gray", no_argument,&cmap,NIIK_COLORMAP_GRAYSCALE},

    {0, 0, 0, 0}
  };
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "i:a:c:C:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'i':
      omin=atof(optarg);
      break;
    case 'a':
      omax=atof(optarg);
      break;
    case 'c':
      column_name=optarg;
      column_id=-1;
      break;
    case 'C':
      column_id=atoi(optarg);
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<3) {
    show_usage(argv[0]);
    return 1;
  }

  in_off =argv[optind];
  in_txt =argv[optind+1];
  out_off=argv[optind+2];


  if (!clobber && !access (out_off, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_off);
    return 1;
  }

  niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"niik_kobj_read_off",1 );
  NIIK_EXIT(((meas_in=(double*)niik_read_vector_ex(in_txt,0,&meas_num_in,column_id,column_name))==NULL),fcname,"niik_read_vector_ex",1);

  if(meas_num_in!=obj->nvert) {
    fprintf(stderr,"[%s] Inconsistent number of measurments, expected %i , got %i\n",fcname,meas_num_in,obj->nvert);
    return 1;
  }

  if(omin==0.0 && omax==0.0) { /*TODO: use different way to set range*/
    omax = niik_get_max_from_double_vector(meas_in,meas_num_in);
    omin = niik_get_min_from_double_vector(meas_in,meas_num_in);
    fprintf(stdout,"[%s] Using range %8.4f %8.4f\n",fcname,omin,omax);
  }

  NIIK_EXIT( (niikcortex_add_color( obj, meas_in, omin, omax, cmap ))==0,      fcname,"niikcortex_add_color",1 );
  NIIK_EXIT((!off_kobj_write_offply(out_off,obj,0)),fcname,"off_kobj_write_off",1);

  off_kobj_free(obj);
  free(meas_in);

  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
