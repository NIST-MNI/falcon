/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 * DATE:
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (8)

static char *prog_help[] = {
  "  niik_inpaint:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
  "  -sd <SD.mnc>                      : output standard deviation SD.mnc image\n"
  "  -float32                          : output 32b float image\n"
  "  -float64                          : output 64b float image\n"
  "  -uint16                           : output 16b unsigned short image\n"
  "  -uint32                           : output 32b unsigned int image\n"
  "  -uint64                           : output 64b unsigned long int image\n"
};

void usage() {
  fprintf(stdout,"  usage: [options] <in.mnc> [...] <output.mnc>\n\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *img=NULL,
   *stdimg=NULL,
    *outimg=NULL;
  double
  d=0;
  int
  i,n,
  datatype=-1,
  nc,sc;
  char
  *stdfname=NULL,
   fcname[32]="niik_avg";

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-float16",8)) {
        datatype=NIFTI_TYPE_FLOAT16;
      } else if(!strncmp(argv[nc],"-float32",8)) {
        datatype=NIFTI_TYPE_FLOAT32;
      } else if(!strncmp(argv[nc],"-float64",8)) {
        datatype=NIFTI_TYPE_FLOAT64;
      } else if(!strncmp(argv[nc],"-double",7)) {
        datatype=NIFTI_TYPE_FLOAT64;
      } else if(!strncmp(argv[nc],"-uint64",7)) {
        datatype=NIFTI_TYPE_UINT64;
      } else if(!strncmp(argv[nc],"-uint32",7)) {
        datatype=NIFTI_TYPE_UINT32;
      } else if(!strncmp(argv[nc],"-uint16",7)) {
        datatype=NIFTI_TYPE_UINT16;
      }

      else if(!strncmp(argv[nc],"-sd",3)) {
        NIIK_EXIT((argc<nc+2),fcname,"missing args",1);
        stdfname=argv[++nc];
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  NIIK_EXIT((argc<3),fcname,"too few arguments",1);

  niik_fc_display(fcname,1);
  if( access( argv[argc-1], F_OK ) != -1 ) {
    fprintf(stderr,"[%s] ERROR: file exists, %s\n",fcname,argv[argc-1]);
    exit(1);
  }
  fprintf(stdout,"[%s] %i images\n",fcname,argc-2);
  d=(double)argc-2.0;

  for(i=1; i<argc-1; i++) {
    fprintf(stdout,"[%s] reading:    %s\n",fcname,argv[i]);
    NIIK_EXIT(((img=niik_image_read(argv[i]))==NULL),
              fcname,
              "reading niik_image_read",1);
    if(outimg==NULL) {
      NIIK_EXIT(((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL),fcname,"niik_image_cop_as_tyep",1);
    } else {
      for(n=0; n<outimg->nvox; n++) {
        niik_image_add_voxel(outimg,n,niik_image_get_voxel(img,n));
      }
    }
    img=niik_image_free(img);
  }
  for(n=0; n<outimg->nvox; n++) {
    niik_image_set_voxel(outimg,n,niik_image_get_voxel(outimg,n)/d);
  }

  if(stdfname!=NULL) { /* SD image */
    NIIK_EXIT(((stdimg=niik_image_copy_as_type(outimg,NIFTI_TYPE_FLOAT64))==NULL),fcname,"niik_image_cop_as_tyep",1);
    niik_image_clear(stdimg);

    for(i=1; i<argc-1; i++) {
      fprintf(stdout,"[%s] reading:    %s\n",fcname,argv[i]);
      NIIK_EXIT(((img=niik_image_read(argv[i]))==NULL),
                fcname,
                "reading niik_image_read",1);
      for(n=0; n<outimg->nvox; n++) {
        niik_image_add_voxel(stdimg,n,
                             pow(niik_image_get_voxel(img,n) - niik_image_get_voxel(outimg,n), 2.0));
      }
      img=niik_image_free(img);
    }
    for(n=0; n<outimg->nvox; n++) {
      niik_image_set_voxel(stdimg,n,sqrt(niik_image_get_voxel(stdimg,n)/d));
    }
    if(datatype>=0) {
      fprintf(stdout,"[%s] type convert\n",fcname);
      NIIK_EXIT((!niik_image_type_convert(stdimg,datatype)),fcname,"niik_image_type_convert",1);
    }
    fprintf(stdout,"[%s] writing:    %s\n",fcname,stdfname);
    NIIK_EXIT((!niik_image_write(stdfname,stdimg)),fcname,"niik_image_write",1);
    stdimg=niik_image_free(stdimg);
  }

  if(datatype>=0) {
    fprintf(stdout,"[%s] type convert\n",fcname);
    NIIK_EXIT((!niik_image_type_convert(outimg,datatype)),fcname,"niik_image_type_convert",1);
  }

  fprintf(stdout,"[%s] writing:    %s\n",fcname,argv[argc-1]);
  NIIK_EXIT((!niik_image_write(argv[argc-1],outimg)),fcname,"niik_image_write",1);
  outimg=niik_image_free(outimg);
  niik_fc_display(fcname,0);
  exit(0);
} /* main */


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/