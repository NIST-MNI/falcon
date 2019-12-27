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

#include "NrrdIO.h"

void show_usage (const char *name) {
  fprintf(stdout,
          "Transform surface using antsApplyTransformsToPoints\n"
          " will use system call to antsApplyTransformsToPoints\n"
          " transformations are specified as accepted by antsApplyTransformsToPoints\n"
          "\n"
          "Usage: %s <input.off/ply> [transform spec 1] [transform spec2] <output.off/ply> \n"
          "\t--clobber clobber output files\n"
          "\n"
          "transform spec - are in the format as understood by antsApplyTransformsToPoints,\n"
          "   either simple transformFileName or  \"[transformFileName,useInverse]\" \n"
          ,
          name);
}


int main(int argc, char **argv) {
  const char *fcname="falcon_transform_surface_ants";
  int clobber=0;
  int verbose=0;
  int c,i;
  int n_transforms=0;
  double *raw_coordinates=NULL;

  char tmp_nrrd_in[1024];
  int nrrd_in=-1;
  char tmp_nrrd_out[1024];
  int nrrd_out=-1;
  const char *tmpdir=getenv("TMPDIR");
  
  const char *in_off=NULL;
  const char *out_off=NULL;
  const char **in_transformations=NULL;

  VIO_General_transform   transform;
  kobj *obj=NULL;

  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"verbose",          no_argument, &verbose, 1},
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

  if((argc - optind)<2) {
    show_usage(argv[0]);
    return 1;
  }
  n_transforms=argc-optind-2;

  in_off =argv[optind];

  in_transformations=(const char**)calloc(n_transforms,sizeof(const char*));

  for(i=0;i<n_transforms;++i)
  {
    in_transformations[i]=argv[optind+1+i];
  }
  out_off=argv[argc-1];

  if (!clobber && !access (out_off, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_off);
    return 1;
  }
  niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL),              fcname,"niik_kobj_read_off",1);

  /*
  Convert phisical points to separate array, store in a temporary .mhd or csv file, apply antsApplyTransformsToPoints , then read back
  */
  if(tmpdir==NULL) tmpdir="/tmp";
  sprintf(tmp_nrrd_in,"%s/inXXXXXX.nrrd",tmpdir);
  sprintf(tmp_nrrd_out,"%s/outXXXXXX.nrrd",tmpdir);
  
  if( (nrrd_in=mkstemps(tmp_nrrd_in,5))!=-1 && 
      (nrrd_out=mkstemps(tmp_nrrd_out,5))!=-1 )
  {
    int r;
    char cmdline[3000];
    Nrrd *nval = nrrdNew();
    double *raw_coordinates;
    
    close(nrrd_in);
    close(nrrd_out);

    nrrdAlloc_va(nval, nrrdTypeDouble, 2, 3, obj->nvert );
    raw_coordinates = (double*)nval->data;

    kvert *v;
    for(v=obj->vert,i=0; v!=NULL; v=v->next,++i) {
      raw_coordinates[i*3   ] = v->v.x;
      raw_coordinates[i*3+ 1] = v->v.y;
      raw_coordinates[i*3+ 2] = v->v.z;
    }
    nrrdSave(tmp_nrrd_in, nval, NULL);
    nrrdNuke(nval);

    /*generate command line for antsApplyTransformsToPoints*/
    sprintf(cmdline,"antsApplyTransformsToPoints -d 3 -p 1 -i %s -o %s -f ",tmp_nrrd_in, tmp_nrrd_out);

    for(i=0;i<n_transforms;i++)
    {
      char tmp[1024];
      sprintf(tmp," -t %s ",in_transformations[i] );
      strncat(cmdline,tmp,sizeof(cmdline)-1);
    }
    if(verbose)
      printf("Running:\n%s\n",cmdline);

    if((r=system(cmdline))!=0)
      fprintf(stderr,"%s\nReturn code:%d\n",argv[0], r);
    nval = nrrdNew();

    if (nrrdLoad(nval, tmp_nrrd_out, NULL)) {
        char *err = biffGetDone(NRRD);
        fprintf(stderr, "%s: trouble reading \"%s\":\n%s", argv[0], tmp_nrrd_out, err);
        free(err);
      } else {
         if(verbose) 
         {
            printf("%s: \"%s\" is a %d-dimensional nrrd of type %d (%s)\n", 
              argv[0], tmp_nrrd_out, nval->dim, nval->type,
              airEnumStr(nrrdType, nval->type));
            printf("%s: the array contains %d elements, each %d bytes in size\n",
              argv[0], (int)nrrdElementNumber(nval), (int)nrrdElementSize(nval));
         }
         raw_coordinates = (double*)nval->data;

        for(v=obj->vert,i=0; v!=NULL; v=v->next,++i) {
          v->v.x = raw_coordinates[i*3   ];
          v->v.y = raw_coordinates[i*3+ 1];
          v->v.z = raw_coordinates[i*3+ 2];
        }
        nrrdNuke(nval);

        NIIK_EXIT((!off_kobj_write_offply(out_off,obj,0)),fcname,"off_kobj_write_off",1);

        /*remove temp files*/
        unlink(tmp_nrrd_in);
        unlink(tmp_nrrd_out);
      }
  } else {
    fprintf(stderr,"Failed to create temp files!\n");
    perror("mkstemp");
  }

  obj=off_kobj_free(obj);

  niik_fc_display(fcname,0);
  free(raw_coordinates);
  free(in_transformations);
  free(timestamp);
  return 0;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
