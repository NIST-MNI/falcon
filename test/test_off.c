/* FILENAME:     test_off.c
 * DESCRIPTION:  test basic object routines
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include  <volume_io.h>

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input.ply> <output.ply> \n",
          name);
}


int main(int argc, char **argv) {
  const char *fcname="off_test";
  int clobber=0;
  int verbose=0;
  int c;
  int i;
  float fwhm=1;
  const char *in_off=NULL;
  const char *out_off=NULL;
  int n;
  kobj *obj=NULL;
  kobj *obj2=NULL;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"verbose", no_argument, &verbose, 1},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "i", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<2) {
    show_usage(argv[0]);
    return 1;
  }

  in_off =argv[optind];
  out_off=argv[optind+1];

  niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"off_kobj_read_offply",1 );
  /*let's see if we are leeking memory*/
  for(i=0;i<1000;i++)
  {
    obj2=off_kobj_copy(obj);
    off_kobj_free(obj);
    obj=obj2;
  }
  NIIK_EXIT((!off_kobj_add_comment(obj2,timestamp)),fcname,"off_kobj_add_comment",1);

  NIIK_EXIT((!off_kobj_write_offply(out_off,obj2,0)),fcname,"off_kobj_write_offply",1);

  off_kobj_free(obj2);

  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
