/* FILENAME:     ply2off.c
 * DESCRIPTION:  Conversion from ply to off file
 * AUTHOR:       Vladimir Fonov
 *
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include "rply.h"

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

/*needed for ply file format, probably should be used everywhere else too*/
#include <locale.h>

static char *prog_version[] = {
  NIIK_VERSION "\n"
  "  program history\n"
  "  0.0.0   : October 4, 2018, Vladimir FONOV <vladimir.fonov@gmail.com>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  ply2off:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"ply2off\n");
  fprintf(stdout,"  usage: [options] <in.ply> [measurements] <out.off>\n");
  fprintf(stdout,"  options:\n");
  fprintf(stdout,"\t--clobber - clobber output\n");
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int output_sph=0;
  int output_normal=0;
  int output_edge=0;
  int clobber=0;
  const char  *fcname="ply2off";
  double *meas=NULL;
  int meas_num;

  const char *in_ply=NULL;
  const char *out_fname=NULL;
  const char *old_locale;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"help",no_argument, 0, 'h'},
    {"version",no_argument, 0, 'v'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hv", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"%s\n",*prog_version);
      return 0;
    case '?':
    case 'h':
    default:
      usage ();
      return 1;
    }
  }

  if((argc - optind)<2) {
    usage();
    return 1;
  }

  in_ply = argv[optind];
  out_fname=argv[optind+1];

  if (!clobber && !access (out_fname, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_fname);
    return 1;
  }

  niik_fc_display(fcname,1);

  old_locale = setlocale(LC_NUMERIC, NULL);
  setlocale(LC_NUMERIC, "C");

  NIIK_EXIT(((obj=off_kobj_read_ply(in_ply))==NULL),fcname,"off_kobj_read_ply",1);
  NIIK_EXIT((!off_kobj_write_off(out_fname, obj, 0)),fcname,"off_kobj_write_off",1);

  setlocale(LC_NUMERIC, old_locale);
  obj=off_kobj_free(obj);
  niik_fc_display(fcname,0);
  return 0;
} /* ply2off */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
