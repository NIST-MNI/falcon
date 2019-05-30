/* FILENAME:     off_smooth.c
 * DESCRIPTION:  Apply iterativs smoothing to the mesh
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
  fprintf(stdout,"Usage: %s <input.off> <output.off> \n"
          "\t--clobber clobber output file\n"
          "\t--delta <f> smoothing delta (weight of neibourhood) (default 0.5)\n"
          "\t--nei <n> smoothing neighbourhood size (default 4)\n",
          name);
}



int main(int argc, char **argv) {
  const char *fcname="falcon_off_smooth";
  int clobber=0;
  int verbose=0;
  int c;
  int i;
  double delta=0.5;
  int nei=4;

  const char *in_off=NULL;
  const char *out_off=NULL;
  int n;

  kobj *obj=NULL;
  off_curvature_t obj_curv;
  /*kobj **out_obj;*/
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"delta",   required_argument,0, 'd'},
    {"nei",     required_argument,0, 'n'},

    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "d:n:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'd':
      delta=atof(optarg);
      break;
    case 'n':
      nei=atoi(optarg);
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



  if (!clobber && !access (out_off, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_off);
    return 1;
  }

  if(verbose)
    niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"niik_kobj_read_off",1 );
  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);

  if(verbose)
  {
    NIIK_EXIT(!init_curvature(&obj_curv,obj),fcname,   "init_curvature",1);
    NIIK_EXIT(!update_curvature(&obj_curv, obj),fcname,"update_curvature",1);
    printf("Mean curvature before=%f\n",mean_curvature(&obj_curv, obj));
  }

  NIIK_EXIT(!niik_off_apply_surface_smoothing(obj,nei,delta),fcname,"niik_off_apply_surface_smoothing",1);

  if(verbose)
  {
    NIIK_EXIT(!update_curvature(&obj_curv, obj),fcname,"update_curvature",1);
    printf("Mean curvature after=%f\n",mean_curvature(&obj_curv, obj));
  }

  NIIK_EXIT((!off_kobj_add_comment(obj,timestamp)),fcname,"off_kobj_add_comment",1);

  off_kobj_add_one_color(obj,1,0,0);
  NIIK_EXIT((!off_kobj_write_offply(out_off,obj,0)), fcname,"off_kobj_write_off",1);

  off_kobj_free(obj);
  free_curvature(&obj_curv);

  if(verbose)
    niik_fc_display(fcname,0);
    
  free(timestamp);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
