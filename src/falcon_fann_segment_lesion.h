#ifndef _FALCON_FANN_SEGMENT_LESION_H_
#define _FALCON_FANN_SEGMENT_LESION_H_

#include "falcon.h"
#include "floatfann.h"

typedef struct {
  char *case_list_file;
  char **case_list;
  int num_case;
  int nvox_per_case;
  int num_data;
  int patch;
  float les_thresh;
  unsigned int verbose;
  char *tmpfann_filename;
  unsigned int num_layers;
  unsigned int num_input;
  unsigned int num_neurons_hidden;
  unsigned int num_output;
  struct fann *ann;
  float desired_error;
  unsigned int max_epochs;
  unsigned int epochs_between_reports;
  char *trained_data_filename;
  int nvox_les;
} nifti1_kunio_fann_segment_lesion; /* class */

nifti1_kunio_fann_segment_lesion *nifti1_kunio_fann_segment_lesion_init() {
  nifti1_kunio_fann_segment_lesion *s=NULL;
  NIIK_RET0(((s=(nifti1_kunio_fann_segment_lesion *)calloc(1,sizeof(nifti1_kunio_fann_segment_lesion)))==NULL),__func__,"calloc");
  s->case_list_file=NULL;
  s->num_case=0;
  s->case_list=NULL;
  s->nvox_per_case=2000;
  s->patch=2;
  s->les_thresh=0.25;
  s->verbose=1;
  s->tmpfann_filename=NULL;
  s->trained_data_filename=NULL;
  s->num_layers=3;
  s->num_data=-1;
  s->num_input=
    s->num_neurons_hidden=
      s->num_output=0;
  s->desired_error = 0.001;
  s->max_epochs = 5e5;
  s->epochs_between_reports = 1e3;
  s->ann=NULL;
  s->nvox_les=100;
  return s;
}

nifti1_kunio_fann_segment_lesion *nifti1_kunio_fann_segment_lesion_free(nifti1_kunio_fann_segment_lesion *s) {
  if(s==NULL) return s;
  free(s);
  return NULL;
}

int nifti1_kunio_fann_segment_lesion_read_case_list(nifti1_kunio_fann_segment_lesion *s) {
  FILE *fp=NULL;
  char *str;
  int n=0,nstr=0;
  niik_fc_display((char *)__func__,1);
  NIIK_RET0((s==NULL),__func__,"s is null");
  NIIK_RET0((s->case_list_file==NULL),__func__,"case list file is null");
  NIIK_RET0(((fp=fopen(s->case_list_file,"r"))==NULL),__func__,"fopen");
  while(!feof(fp)) {
    fgetc(fp);
    nstr++;
  }
  fclose(fp);
  str=(char *)calloc(nstr+1,sizeof(char));
  NIIK_RET0(((fp=fopen(s->case_list_file,"r"))==NULL),__func__,"fopen");
  while(!feof(fp)) {
    str[n]=fgetc(fp);
    if(str[n]=='\n') str[n]=',';
    n++;
  }
  fclose(fp);
  s->case_list=niik_csstring_to_list(str,&s->num_case);
  for(n=0; n<s->num_case; n++) {
    if(!isalnum(s->case_list[n][0])) break;
    // fprintf(stdout,"%5i %s %i\n",n,s->case_list[n],s->case_list[n][0]);
  }
  s->num_case=n;
  free(str);
  return 1;
}

int nifti1_kunio_fann_segment_lesion_train_create_temp_fann(nifti1_kunio_fann_segment_lesion *s) {
  char fname[4096];
  int n,nl,nm,nskip;
  FILE *fp=NULL;
  nifti_image *mask=NULL;
  nifti_image *t1w=NULL;
  nifti_image *tww=NULL;
  nifti_image *flr=NULL;
  nifti_image *les=NULL;
  nifti_image *iles=NULL;
  nifti_image *pwm=NULL;
  int
  mi,mj,mk,mm,
  idx[3],
  i,j,k;

  niik_fc_display((char *)__func__,1);
  NIIK_RET0((s==NULL),__func__,"s is null");
  NIIK_RET0((s->case_list==NULL),__func__,"case list is null");
  NIIK_RET0((s->num_case<=0),__func__,"case list is empty");
  NIIK_RET0((s->tmpfann_filename==NULL),__func__,"temp fann file name is null");

  fprintf(stdout,"[%s] #cases %i\n",__func__,s->num_case);

  s->num_input=(int)(1+2*pow(s->patch*2+1,3));
  s->num_output=1;
  if(niik_access_file(s->tmpfann_filename)) {
    niik_fc_display((char *)__func__,0);
    return 1;
  }

  fprintf(stdout,"[%s:%i:%s] creating tmp fann file, %s\n",__FILE__,__LINE__,__func__,s->tmpfann_filename);
  NIIK_RET0(((fp=fopen(s->tmpfann_filename,"w"))==NULL),__func__,"fopen");

  if(s->nvox_les<0) {
    for(n=nl=0; n<s->num_case; n++) {
      sprintf(fname,"%s/stx/%s_ct2f.mnc",s->case_list[n],s->case_list[n]);
      NIIK_EXIT(((les=niik_image_read(fname))==NULL),__func__,"read les",1);
      NIIK_EXIT(((iles=niik_image_threshold_new(les,s->les_thresh))==NULL),__func__,"threshold image",1);
      nl+=niik_image_count_mask(iles);
      iles=niik_image_free(iles);
      les=niik_image_free(les);
    }
    fprintf(stdout,"[%s] #data  %i\n",__func__,nl);
    s->nvox_per_case=-1;
    fprintf(fp,"%i %i %i\n",nl*2,(int)(1+2*pow(s->patch*2+1,3)),(int)1);
    s->num_data = nl*2;
  } else {
    fprintf(stdout,"[%s] #data  %i\n",__func__,s->num_case*s->nvox_per_case);
    fprintf(fp,"%i %i %i\n",(int)(s->num_case*s->nvox_per_case),(int)(1+2*pow(s->patch*2+1,3)),(int)1);
    s->num_data = s->num_case * s->nvox_per_case;
  }

  fprintf(stdout,"[%s] les thresh %.3f\n",__func__,s->les_thresh);

  for(n=0; n<s->num_case; n++) {
    sprintf(fname,"%s/stx/%s_t1p_brain.mnc",s->case_list[n],s->case_list[n]);
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] reading image %s\n",__FILE__,__LINE__,__func__,fname);
    NIIK_EXIT(((t1w=niik_image_read(fname))==NULL),__func__,"read t1w",1);
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] reading image %s\n",__FILE__,__LINE__,__func__,fname);
    sprintf(fname,"%s/stx/%s_flr_brain.mnc",s->case_list[n],s->case_list[n]);
    NIIK_EXIT(((flr=niik_image_read(fname))==NULL),__func__,"read flr",1);

    sprintf(fname,"%s/stx/%s_ct2f.mnc",s->case_list[n],s->case_list[n]);
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] reading image %s\n",__FILE__,__LINE__,__func__,fname);
    NIIK_EXIT(((les=niik_image_read(fname))==NULL),__func__,"read les",1);
    NIIK_EXIT(((iles=niik_image_threshold_new(les,s->les_thresh))==NULL),__func__,"threshold image",1);
    nl=niik_image_count_mask(iles);
    sprintf(fname,"%s/%s_t1p_brain_mask_reg2icbm.mnc",s->case_list[n],s->case_list[n]);
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] reading image %s\n",__FILE__,__LINE__,__func__,fname);
    NIIK_EXIT(((mask=niik_image_read(fname))==NULL),__func__,"read mask",1);
    nm=niik_image_count_mask(mask);
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] %s %9i %9i\n",__FILE__,__LINE__,__func__,s->case_list[n],nm,nl);

    // lesion sampling
    if(s->nvox_les<0) {
      nskip=0;
    } else {
      nskip=nl/(s->nvox_les+10);
    }
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] lesion sampling %i / %i\n",__FILE__,__LINE__,__func__,nskip,nl);
    for(i=j=k=0; i<mask->nvox; i++) {
      if(niik_image_get_voxel(iles,i)==0.0) continue;
      if(s->nvox_les>0) {
        if(k>=s->nvox_les) break;
      } else if(k>=nl) break;
      j++;
      if(j>=nskip) {
        j=0;
        k++;
        fprintf(fp,"1 ");
        idx[0]=i%les->nx;
        idx[1]=(i/les->nx)%les->ny;
        idx[2]=i/les->nx/les->ny;
        for(mk=-s->patch; mk<=s->patch; mk++) {
          for(mj=-s->patch; mj<=s->patch; mj++) {
            for(mi=-s->patch; mi<=s->patch; mi++) {
              mm=(idx[0]+mi)+(idx[1]+mj)*les->nx+(idx[2]+mk)*les->nx*les->ny;
              fprintf(fp,"%.5f ",niik_image_get_voxel(t1w,mm));
              fprintf(fp,"%.5f ",niik_image_get_voxel(flr,mm));
            }
          }
        }
        fprintf(fp,"\n%.5f\n",niik_image_get_voxel(iles,i));
        //fprintf(fp,"\n%.5f\n",niik_image_get_voxel(les,i));
        // fprintf(stdout,"%i %i\n",i,(int)niik_image_get_voxel(iles,i));
      }
    }

    // NABT sampling
    if(s->nvox_les<0) {
      nskip=(nm-nl)/(nl+10);
    } else {
      nskip=(nm-nl)/(s->nvox_per_case-s->nvox_les+10);
    }
    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] nabt sampling %i / %i\n",__FILE__,__LINE__,__func__,nskip,nm-nl);
    for(i=j=0; i<mask->nvox; i++) {
      if(niik_image_get_voxel(mask,i)==0.0) continue;
      if(niik_image_get_voxel(iles,i)>0.0) continue;
      if(s->nvox_per_case>0) {
        if(k>=s->nvox_per_case) break;
      } else if(k>=nl*2) break;
      j++;
      if(j>=nskip) {
        j=0;
        k++;
        fprintf(fp,"1 ");
        idx[0]=i%les->nx;
        idx[1]=(i/les->nx)%les->ny;
        idx[2]=i/les->nx/les->ny;
        for(mk=-s->patch; mk<=s->patch; mk++) {
          for(mj=-s->patch; mj<=s->patch; mj++) {
            for(mi=-s->patch; mi<=s->patch; mi++) {
              mm=(idx[0]+mi)+(idx[1]+mj)*les->nx+(idx[2]+mk)*les->nx*les->ny;
              fprintf(fp,"%.5f ",niik_image_get_voxel(t1w,mm));
              fprintf(fp,"%.5f ",niik_image_get_voxel(flr,mm));
            }
          }
        }
        fprintf(fp,"\n%.5f\n",niik_image_get_voxel(iles,i));
        //fprintf(fp,"\n%.5f\n",niik_image_get_voxel(les,i));
      }
    }

    if(s->verbose>1)fprintf(stdout,"[%s:%i:%s] free images\n",__FILE__,__LINE__,__func__);
    les=niik_image_free(les);
    t1w=niik_image_free(t1w);
    flr=niik_image_free(flr);
    iles=niik_image_free(iles);
    mask=niik_image_free(mask);
  }

  fclose(fp);

  niik_fc_display((char *)__func__,0);
  return 1;
}


int nifti1_kunio_fann_segment_lesion_train_from_temp_fann(nifti1_kunio_fann_segment_lesion *s) {

  niik_fc_display((char *)__func__,1);
  NIIK_RET0((s==NULL),__func__,"s is null");
  NIIK_RET0((s->tmpfann_filename==NULL),__func__,"tmpfann filename is null");
  NIIK_RET0((s->trained_data_filename==NULL),__func__,"trained data filename is null");
  NIIK_RET0((niik_access_file(s->trained_data_filename)),__func__,"file exists");
  NIIK_RET0((s->num_input<=0),__func__,"num input is bad");
  NIIK_RET0((s->num_output<=0),__func__,"num output is bad");
  NIIK_RET0((s->num_layers<=0),__func__,"num layers is bad");
  NIIK_RET0((s->num_neurons_hidden<=0),__func__,"num neurons hidden is bad");

  fprintf(stdout,"num_layers = %i\n",s->num_layers);
  fprintf(stdout,"num_input  = %i\n",s->num_input);
  fprintf(stdout,"num_output = %i\n",s->num_output);
  fprintf(stdout,"num_neurons_hidden = %i\n",s->num_neurons_hidden);
  if(s->num_data!=-1) fprintf(stdout,"num_data   = %i\n",s->num_data);
  NIIK_RET0(((s->ann = fann_create_standard(s->num_layers, s->num_input,
                       s->num_neurons_hidden, s->num_output))==NULL),__func__,"fann_create_standard");
  fann_set_activation_function_hidden(s->ann, FANN_SIGMOID_SYMMETRIC);
  fann_set_activation_function_output(s->ann, FANN_SIGMOID_SYMMETRIC);
  fann_print_connections(s->ann);
  fann_train_on_file(s->ann, s->tmpfann_filename, s->max_epochs, s->epochs_between_reports, s->desired_error);

  fprintf(stdout,"[%s] saving %s",__func__,s->trained_data_filename);
  fann_save(s->ann, s->trained_data_filename);
  fann_destroy(s->ann);

  niik_fc_display((char *)__func__,0);
  return 1;
}

#endif // _FALCON_FANN_SEGMENT_LESION_H_
