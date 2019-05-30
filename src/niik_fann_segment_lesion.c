/* FILENAME:     niik_fann_segment_lesion.c
 * DESCRIPTION:  segment lesion using fann
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_fann_segment_lesion.h"

void usage() {
  fprintf(stderr,"  usage: <case_list_file.txt>\n");
  fprintf(stderr,"\n  optional usage:\n");
  fprintf(stderr,"  -tmpfann <filename.dat>      temporary fann train dataset file\n");
  fprintf(stderr,"  -patch <patch_size>          patch size [default=2] means 5x5x5 patch\n");
  fprintf(stderr,"  -trained <out.net>           trained network file\n");
  fprintf(stderr,"  -nvoxles <#lesion_voxel>     number of leison voxels per case [default=100]\n");
  fprintf(stderr,"                               if -1, then all lesion voxels and corresponding number of nabt voxels\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti1_kunio_fann_segment_lesion *s=NULL;
  int nc=1,sc=1;

  niik_fc_display((char *)__func__,1);
  NIIK_EXIT(((s=nifti1_kunio_fann_segment_lesion_init())==NULL),__func__,"nifti1_kunio_fann_segment_lesion_init",1);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--help",6)) {
        usage();
        exit(1);
      } else if(!strncmp(argv[nc],"-help",5)) {
        usage();
        exit(1);
      } else if(!strncmp(argv[nc],"-tmpfann",8)) {
        s->tmpfann_filename=argv[++nc];
        fprintf(stdout,"[%s] tmp fann filename = %s\n",__func__,s->tmpfann_filename);
      } else if(!strncmp(argv[nc],"-patch",6)) {
        s->patch=atoi(argv[++nc]);
        fprintf(stdout,"[%s] patch = %i\n",__func__,s->patch);
      } else if(!strncmp(argv[nc],"-trained",8)) {
        s->trained_data_filename=argv[++nc];
      } else if(!strncmp(argv[nc],"-nvoxles",8)) {
        s->nvox_les=atoi(argv[++nc]);
        fprintf(stdout,"[%s] nvox les = %i\n",__func__,s->nvox_les);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        usage();
        exit(1);
      } else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",__func__,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc!=2) {
    usage();
    exit(1);
  }

  fprintf(stdout,"[%s] case list file: %s\n",__func__,argv[1]);
  s->case_list_file=argv[1];
  NIIK_EXIT((!nifti1_kunio_fann_segment_lesion_read_case_list(s)),__func__,"nifti1_kunio_fann_segment_lesion_read_case_list",1);
  NIIK_EXIT((!nifti1_kunio_fann_segment_lesion_train_create_temp_fann(s)),__func__,"nifti1_kunio_fann_segment_lesion_train_create_temp_fann",1);
  s->num_neurons_hidden=3;
  NIIK_EXIT((!nifti1_kunio_fann_segment_lesion_train_from_temp_fann(s)),__func__,"nifti1_kunio_fann_segment_lesion_train_create_temp_fann",1);

  s=nifti1_kunio_fann_segment_lesion_free(s);
  niik_fc_display((char *)__func__,0);
  exit(0);
} // main

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/