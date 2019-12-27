/* FILENAME:     transform_off.c
 * DESCRIPTION:  Apply arbitrary xfm transformation to off file
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"

#include  <volume_io.h>

#include <unistd.h>
#include <getopt.h>


void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input.off/ply> <xfm> <output.off/ply> \n"
          "\t--clobber clobber output files\n"
          "\t--invert_transform \n",
          name);
}

/*TODO: reorient faces if the transform is left-handed?*/
int apply_transform_kobj(kobj *obj,
                         VIO_General_transform* transform ) {
  kvert *v;
  for(v=obj->vert; v!=NULL; v=v->next) {
    VIO_Real x,y,z;
    general_transform_point( transform,
                             (VIO_Real) v->v.x,
                             (VIO_Real) v->v.y,
                             (VIO_Real) v->v.z,
                             &x, &y, &z );
    v->v.x=x;
    v->v.y=y;
    v->v.z=z;
  }
  return 1;
}

int main(int argc, char **argv) {
  const char *fcname="falcon_transform_surface";
  int clobber=0;
  int verbose=0;
  int invert_transform=0;
  int c;
  float fwhm=1;
  const char *in_off=NULL;
  const char *out_off=NULL;
  const char *in_xfm=NULL;

  VIO_General_transform   transform;
  kobj *obj=NULL;

  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"verbose",          no_argument, &verbose, 1},
    {"invert_transform", no_argument, &invert_transform, 1},
    {0, 0, 0, 0}
  };
  char* timestamp = niik_create_minc_timestamp(argc,argv);

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

  if((argc - optind)<3) {
    show_usage(argv[0]);
    return 1;
  }

  in_off =argv[optind];
  in_xfm =argv[optind+1];
  out_off=argv[optind+2];


  if (!clobber && !access (out_off, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_off);
    return 1;
  }
  niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL),              fcname,"niik_kobj_read_off",1);
  NIIK_EXIT( (input_transform_file( in_xfm, &transform ) != VIO_OK),fcname,"input_transform_file",1 );

  if( invert_transform )
    invert_general_transform( &transform );

  NIIK_EXIT( !apply_transform_kobj( obj, &transform ),fcname,"apply_transform_kobj",1 );
  NIIK_EXIT((!off_kobj_write_offply(out_off,obj,0)),fcname,"off_kobj_write_off",1);

  obj=off_kobj_free(obj);
  delete_general_transform( &transform);

  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
