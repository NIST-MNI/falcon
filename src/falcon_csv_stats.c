/* FILENAME:     falcon_csv_stats.c
 * DESCRIPTION:  Simple stats on csv file, can be easily done by R or python
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <scalar.txt> \n"
          "\t--clobber clobber output files\n"
          "\t--column <name> - column name to use\n"
          "\t--column_id <name> - column id to use\n"
          "\t--min - use spectral colour map (default)\n"
          ,name);
}

int main(int argc, char **argv) {
  const char *fcname="falcon_csv_stats";
  int clobber=0;
  int verbose=0;

  const char *in_txt=NULL;
  const char *out_txt=NULL;

  double *meas_in=NULL;
  int meas_num_in;
  int column_id=0;
  int n,c;
  const char *column_name=NULL;
  enum {OP_MIN=0,OP_MAX,OP_RANGE,OP_MEAN,OP_MEDIAN,OP_SD};
  int op=OP_RANGE;
  double omin,omax;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"min", no_argument, &op, OP_MIN},
    {"max", no_argument, &op, OP_MAX},
    {"range", no_argument, &op, OP_RANGE},
    {"mean", no_argument, &op, OP_MEAN},
    {"median", no_argument, &op, OP_MEDIAN},
    {"sd", no_argument, &op, OP_SD},
    {"column", required_argument,0,'c'},
    {"column_id", required_argument,0,'C'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "o:C:C:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
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

  if((argc - optind)<1) {
    show_usage(argv[0]);
    return 1;
  }

  in_txt =argv[optind];

  // if (!clobber && !access (out_off, F_OK)) {
  //   fprintf(stderr,"%s Exists!\n", out_off);
  //   return 1;
  // }

  //niik_fc_display(fcname,1);

  NIIK_EXIT(((meas_in=(double*)niik_read_vector_ex(in_txt,0,&meas_num_in,column_id,column_name))==NULL),fcname,"niik_read_vector_ex",1);
  switch(op)
  { 
    case OP_MIN:
      omin = niik_get_min_from_double_vector(meas_in,meas_num_in);
      fprintf(stdout,"%.10e\n",omin);
      break;
    case OP_MAX:
      omax = niik_get_min_from_double_vector(meas_in,meas_num_in);
      fprintf(stdout,"%.10e\n",omax);
      break;
    case OP_MEAN:
      omin = niik_get_mean_from_double_vector(meas_in,meas_num_in);
      fprintf(stdout,"%.10e\n",omin);
      break;
    case OP_SD:
      omin = niik_get_stdv_from_double_vector(meas_in,meas_num_in);
      fprintf(stdout,"%.10e\n",omin);
      break;
    case OP_MEDIAN:
      omin = niik_median_quicksort_double(meas_in,meas_num_in);
      fprintf(stdout,"%.10e\n",omin);
      break;
    case OP_RANGE:
    default:
      omax = niik_get_max_from_double_vector(meas_in,meas_num_in);
      omin = niik_get_min_from_double_vector(meas_in,meas_num_in);
      fprintf(stdout,"%.10e %.10e\n",omin,omax);
  }
  free(meas_in);
  return 0;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
