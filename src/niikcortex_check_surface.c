/* FILENAME:     niikcortex_check_off.c
 * DESCRIPTION:  Check surface for self-intersections and mark them in red
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
          "\t--fix  fix self intersections\n",
          name);
}



int main(int argc, char **argv) {
  const char *fcname="niikcortex_check_off.c";
  int clobber=0;
  int verbose=0;
  int c;
  int i;
  int fix_self_intersections=0;

  const char *in_off=NULL;
  const char *out_off=NULL;
  int n;

  bbox *bb=NULL;
  kobj *obj=NULL;
  off_curvature_t obj_curv;
  /*kobj **out_obj;*/
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"fix", no_argument, &fix_self_intersections, 1},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "", long_options, &option_index);

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

  if (!clobber && !access (out_off, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_off);
    return 1;
  }
  niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"niik_kobj_read_off",1 );
  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);
  off_update_kobj_kface_pminmax(obj);

  bb=off_bbox_init(7,320);
  off_create_bbox_from_kobj(bb,obj);

  if((n=off_count_self_intersection_add_color(bb,obj,1))<0) {
    fprintf(stderr,"[%s] ERROR: off_count_self_intersection_add_color(bb,obj,1)\n",fcname);
    exit(1);
  }
  fprintf(stdout,"[%s] surface %i intersection(s)\n",fcname,n);

  if(fix_self_intersections && n>0) {
    fprintf(stdout,"[%s] Fixing intersections : %d\n",fcname,n);
    off_correct_self_intersection_with_elen(bb,obj,-1);
    off_create_bbox_from_kobj(bb,obj);
    if((n=off_count_self_intersection_add_color(bb,obj,1))<0) {
      fprintf(stderr,"[%s] ERROR: off_count_self_intersection_add_color(bb,obj,1)\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] surface %i intersection(s)\n",fcname,n);
  }

  /*append metadata*/
  NIIK_EXIT((!off_kobj_add_comment(obj,timestamp)),fcname,"off_kobj_add_comment",1);

  NIIK_EXIT((!off_kobj_write_offply(out_off,obj,0)), fcname,"off_kobj_write_off",1);

  bb=off_bbox_free(bb);
  off_kobj_free(obj);
  niik_fc_display(fcname,0);

  free(timestamp);
  return 0;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
