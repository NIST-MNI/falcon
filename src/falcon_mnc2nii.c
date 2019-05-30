/* Filename:    nifti1_kunio_mnc2nii.c
 * Description: based on from http://minc.sourcearchive.com/documentation/2.0.18-1build1/mnc2nii_8c-source.html
 * Author:      Kunio Nakamura
 * Date:        December 18, 2012
 */
#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <hdf5.h>
#include "minc.h"
#include "minc2.h"
#include "nifti1_io.h"
#include "nifti1_local.h"       /* Our local definitions */
#include "falcon.h"

#ifndef MINC2 /*tehcnically it should be determined by config script*/
#define MINC2 1
#endif

struct midimension;

/* This list is in the order in which dimension lengths and sample
 * widths are stored in the NIfTI-1 structure.
 */
static const char *dimnames[MAX_NII_DIMS] = {
  MIvector_dimension,
  MItime,
  MIzspace,
  MIyspace,
  MIxspace,
  NULL,
  NULL,
  NULL
};

void test_xform(mat44 m, int i, int j, int k) {
  double x, y, z;
  x = m.m[DIM_X][DIM_I] * i + m.m[DIM_X][DIM_J] * j + m.m[DIM_X][DIM_K] * k
      + m.m[DIM_X][3];
  y = m.m[DIM_Y][DIM_I] * i + m.m[DIM_Y][DIM_J] * j + m.m[DIM_Y][DIM_K] * k
      + m.m[DIM_Y][3];
  z = m.m[DIM_Z][DIM_I] * i + m.m[DIM_Z][DIM_J] * j + m.m[DIM_Z][DIM_K] * k
      + m.m[DIM_Z][3];
  printf("%d %d %d => ", i, j, k);
  printf("%f %f %f\n", x, y, z);
}

/* Explicitly set all of the fields of the NIfTI I/O header structure to
 * something reasonable. Right now this is overkill since a simple memset()
 * would do the same job, but I want this function to help me keep track
 * of all of the header fields and to allow me to easily override a default
 * if it becomes useful.
 */
void init_nifti_header(nifti_image *nii_ptr) {
  int i, j;
  nii_ptr->ndim = 0;
  nii_ptr->nx = nii_ptr->ny = nii_ptr->nz = nii_ptr->nt = nii_ptr->nu =
                                nii_ptr->nv = nii_ptr->nw = 1;
  for (i = 0; i < MAX_NII_DIMS; i++) {
    /* Fix suggested by Hyun-Pil Kim (hpkim@ihanyang.ac.kr):
       Use 1 as the default, not zero */
    nii_ptr->dim[i] = 1;
  }
  nii_ptr->nvox = 0;
  nii_ptr->nbyper = 0;
  nii_ptr->datatype = DT_UNKNOWN;
  nii_ptr->dx = nii_ptr->dy = nii_ptr->dz = nii_ptr->dt = nii_ptr->du =
                                nii_ptr->dv = nii_ptr->dw = 0.0;
  for (i = 0; i < MAX_NII_DIMS; i++) {
    nii_ptr->pixdim[i] = 0.0;
  }
  nii_ptr->scl_slope = 0.0;
  nii_ptr->scl_inter = 0.0;
  nii_ptr->cal_min = 0.0;
  nii_ptr->cal_max = 0.0;
  nii_ptr->qform_code = NIFTI_XFORM_UNKNOWN;
  nii_ptr->sform_code = NIFTI_XFORM_UNKNOWN;
  nii_ptr->freq_dim = 0;
  nii_ptr->phase_dim = 0;
  nii_ptr->slice_dim = 0;
  nii_ptr->slice_code = 0;
  nii_ptr->slice_start = 0;
  nii_ptr->slice_end = 0;
  nii_ptr->slice_duration = 0.0;
  nii_ptr->quatern_b = 0.0;
  nii_ptr->quatern_c = 0.0;
  nii_ptr->quatern_d = 0.0;
  nii_ptr->qoffset_x = 0.0;
  nii_ptr->qoffset_y = 0.0;
  nii_ptr->qoffset_z = 0.0;
  nii_ptr->qfac = 0.0;
  nii_ptr->toffset = 0.0;
  nii_ptr->xyz_units = NIFTI_UNITS_MM; /* Default spatial units */
  nii_ptr->time_units = NIFTI_UNITS_SEC; /* Default time units */
  nii_ptr->nifti_type = FT_ANALYZE;
  nii_ptr->intent_code = 0;
  nii_ptr->intent_p1 = 0.0;
  nii_ptr->intent_p2 = 0.0;
  nii_ptr->intent_p3 = 0.0;
  memset(nii_ptr->intent_name, 0, sizeof (nii_ptr->intent_name));
  memset(nii_ptr->descrip, 0, sizeof (nii_ptr->descrip));
  memset(nii_ptr->aux_file, 0, sizeof (nii_ptr->aux_file));
  nii_ptr->fname = NULL;
  nii_ptr->iname = NULL;
  nii_ptr->iname_offset = 0;
  nii_ptr->swapsize = 0;
  nii_ptr->byteorder = 1; /* default order (LSB) */
  nii_ptr->data = NULL;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      nii_ptr->qto_xyz.m[i][j] = 0.0;
      nii_ptr->qto_ijk.m[i][j] = 0.0;
      nii_ptr->sto_xyz.m[i][j] = 0.0;
      nii_ptr->sto_ijk.m[i][j] = 0.0;
    }
  }
  nii_ptr->num_ext = 0;
}

/* private function from from libminc2.  This function is private partially
   because it's parameters are somewhat bizarre.  It would be a good idea to
   rework them into a more rational and easily described form.
*/

extern void restructure_array(int ndims,
                              unsigned char *array,
                              const unsigned long *lengths_perm,
                              int el_size,
                              const int *map,
                              const int *dir);

#define RTN_ERR(code,msg,fcname) if ((code) == MI_ERROR) {fprintf(stderr,"[%s] ERROR %s\n",fcname,msg);return NULL;}


int niik_mat44_from_cosines_start_step(mat44 mat,double dx,double dy,double dz,double *ocosx,double *ocosy,double *ocosz,double *ostart,double *ostep) {
  char fcname[64]="niik_mat44_from_cosines_start_step";
  double step[3],start[3];
  double cosines[3][3];
  int i,j;
  int verbose=0;
  if(verbose>=1) niik_fc_display(fcname,1);
  if(verbose>=1) fprintf(stdout,"[%s] pixel size   : %12.8f %12.8f %12.8f\n",fcname,dx,dy,dz);
  mat44_display(mat);
  step[0]=dx;
  step[1]=dy;
  step[2]=dz;
  if(verbose>=1) fprintf(stdout,"[%s] step         : %12.8f %12.8f %12.8f\n",fcname,step[0],step[1],step[2]);
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      cosines[j][i]=mat.m[i][j]/step[j];
    }
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] cosines0     : %12.8f %12.8f %12.8f\n",fcname,cosines[0][0],cosines[0][1],cosines[0][2]);
    fprintf(stdout,"[%s] cosines1     : %12.8f %12.8f %12.8f\n",fcname,cosines[1][0],cosines[1][1],cosines[1][2]);
    fprintf(stdout,"[%s] cosines2     : %12.8f %12.8f %12.8f\n",fcname,cosines[2][0],cosines[2][1],cosines[2][2]);
  }
  for(i=0; i<3; i++) {
    start[i]=0.0;
    for(j=0; j<3; j++) {
      start[i]+=cosines[i][j]*mat.m[j][3];
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] start        : %12.8f %12.8f %12.8f\n",fcname,start[0],start[1],start[2]);
  for(i=0; i<3; i++) {
    ocosx[i] = cosines[0][i];
    ocosy[i] = cosines[1][i];
    ocosz[i] = cosines[2][i];
    ostep[i] = step[i];
    ostart[i] = start[i];
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}


nifti_image *niik_image_read_minc2(const char *minc_name) {
  mihandle_t
  minc_volume;
  int
  i,n=0,
    num_dim,
    icod,jcod,kcod,
    result,
    verbose=niik_verbose();
  misize_t nvox;
  char fcname[64]="niik_image_read_minc2";
  double cosine[9];
  double *delta=NULL;
  double voxel_position[9];
  double world_position[9];
  /*static char *dimorder[] = { "time", "zspace","yspace","xspace" };*/

  misize_t *start,*size;
  misize_t *count;
  midimhandle_t *dims=NULL;
  nifti_image *img=NULL;
  niikpt pt,qt;
  mat44 mm;
  size_t  minc_history_length=0;

  int do_reorient=1;

  if(minc_name==NULL) {
    fprintf(stderr,"[%s] ERROR: minc_name is null\n",fcname);
    return NULL;
  }
  if(verbose>1) niik_fc_display(fcname,1);

  if(verbose>1) {
    fprintf(stdout,"[%s] miopen_volume %s\n",fcname,minc_name);
  }

  if((result = miopen_volume(minc_name, MI2_OPEN_READ, &minc_volume)) != MI_NOERROR) {
    fprintf(stderr, "[%s] ERROR: opening input file: %d.\n", fcname,result);
    return NULL;
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] miopen_volume: %i  %i\n",fcname,result,MI_NOERROR);
  }

  if((result = miget_volume_voxel_count (minc_volume, &nvox)) != MI_NOERROR) {
    fprintf(stderr, "[%s] ERROR: miget_vol_dimensions: %d.\n", fcname,result);
    return NULL;
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] miget_volume_voxel_count %i\n",fcname,(int)nvox);
  }
  /*if((result = miset_apparent_dimension_order_by_name(minc_volume, 4, dimorder)) != MI_NOERROR) {
    fprintf(stderr, "[%s] ERROR: miset_apparent_dimension_order_by_name: %d.\n", fcname,result);
    return NULL; }
    if(verbose>=1) { fprintf(stdout,"[%s] miset_apparent_dimension_order_by_name %i\n",fcname,result); } */

  if((miget_volume_dimension_count( minc_volume, MI_DIMCLASS_SPATIAL,
                                    MI_DIMATTR_ALL, &num_dim) ) == MI_ERROR ) {
    fprintf(stderr, "[%s] ERROR: miget_volume_dimension_count: %d\n", fcname,result);
    return NULL;
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] miget_volume_dimension_count %i\n",fcname,num_dim);
  }

  dims  = (midimhandle_t *)calloc(num_dim,sizeof(midimhandle_t));
  start = (misize_t *)calloc(num_dim,sizeof(misize_t));
  count = (misize_t *)calloc((num_dim>10)?num_dim:10,sizeof(misize_t));
  size  = (misize_t *)calloc(num_dim,sizeof(misize_t));
  delta = (double  *)calloc((num_dim>10)?num_dim:10,sizeof(double));
  for(i=0; i<10; i++) {
    count[i]=1;
  }
  if(( result = miget_volume_dimensions(minc_volume,
                                        MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
                                        MI_DIMORDER_FILE,
                                        num_dim,
                                        dims) ) == MI_ERROR ) {
    fprintf(stderr, "[%s] ERROR: miget_volume_dimensions: %d\n", fcname,result);
    return NULL;
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] miget_volume_dimensions %i\n",fcname,result);
  }

  for(i=0; i<num_dim; i++) {
    start[i]=0;
    if(verbose>1) {
      fprintf(stdout,"[%s] dim %i\n",fcname,i);
    }
    if(( miget_dimension_size ( dims[i],  &count[i] )) == MI_ERROR ) {
      fprintf(stderr, "[%s] ERROR: miget_dimension_size: %d\n", fcname,result);
      return NULL;
    }
    if(verbose>1) {
      fprintf(stdout,"[%s]   miget_dimension_size %i\n",fcname,result);
    }
    if(verbose>1) {
      fprintf(stdout,"[%s]   count[%i] %i\n",fcname,i,(int)count[i]);
    }
    size[i]=(unsigned long)count[i];
    if(( result = miget_dimension_width ( dims[i],
                                          &delta[i] ) ) == MI_ERROR ) {
      fprintf(stderr, "[%s] ERROR: miget_dimension_width: %d\n", fcname,result);
      return NULL;
    }
    if(verbose>1) {
      fprintf(stdout,"[%s]   delta[%i] %12.7f\n",fcname,i,delta[i]);
    }
  }

  if((img = niik_image_init(count[2],count[1],count[0],count[3],count[4],count[5],count[6],
                            delta[2],delta[1],delta[0],delta[3],delta[5],delta[5],delta[6],
                            NIFTI_TYPE_FLOAT32)) == NULL ) {
    fprintf(stderr, "[%s] ERROR: niik_image_init\n", fcname);
    return NULL;
  }

  if(verbose>1) {
    fprintf(stdout,"  img->nvox  %i\n",(int)img->nvox);
    /*fprintf(stdout,"\tminc_volume = %i\n",minc_volume);*/
  }

  for(i=0; i<num_dim; i++) {
    start[i]=0;
    if(verbose>1) {
      fprintf(stdout,"[%s] dim %i  count = %i\n",fcname,i,(int)count[i]);
    }
    if( (result = miget_dimension_cosines( dims[i], cosine)) == MI_ERROR) {
      fprintf(stderr, "[%s] ERROR: miget_dimension_cosines: %d\n", fcname,result);
      return NULL;
    }
    if(verbose>=1) {
      fprintf(stdout,"[%s]   cosine[%i] %12.7f %12.7f %12.7f\n",fcname,i,cosine[0],cosine[1],cosine[2]);
    }
    img->qto_xyz.m[2-i][0] = cosine[0];
    img->qto_xyz.m[2-i][1] = cosine[1];
    img->qto_xyz.m[2-i][2] = cosine[2];
    if(( result = miget_dimension_start( dims[i],
                                         MI_ORDER_FILE,
                                         &delta[i]) ) == MI_ERROR ) {
      fprintf(stderr, "[%s] ERROR: miget_dimension_start: %d\n", fcname,result);
      return NULL;
    }
    if(verbose>1) {
      fprintf(stdout,"[%s]   start[%i] %12.7f\n",fcname,i,delta[i]);
    }
    img->qto_xyz.m[i][3] = delta[i];
  }
  img->qform_code = NIFTI_XFORM_SCANNER_ANAT;
  if(verbose>1) {
    fprintf(stdout,"\tqto_xyz = \n");
    mat44_display(img->qto_xyz);
  }
  img->sto_xyz=img->qto_xyz;

  if((result = miget_real_value_hyperslab( minc_volume,
               MI_TYPE_FLOAT,
               start,
               size,
               img->data )) == MI_ERROR ) {
    fprintf(stderr, "[%s] ERROR: miget_real_value_hyperslab: %d\n", fcname,result);
    return NULL;
  }

  nifti_mat44_to_orientation(img->qto_xyz,&icod,&jcod,&kcod);
  if(verbose>2) {
    fprintf(stdout,"  x = %s\n",nifti_orientation_string(icod));
    fprintf(stdout,"  y = %s\n",nifti_orientation_string(jcod));
    fprintf(stdout,"  z = %s\n",nifti_orientation_string(kcod));
  }

  if(verbose>1) {
    fprintf(stdout,"\tsto_xyz = \n");
    mat44_display(img->sto_xyz);
  }

  if(do_reorient) {
    if(verbose>1) fprintf(stdout,"  do reorientation\n");
    if     ((fabs(img->qto_xyz.m[0][0])>fabs(img->qto_xyz.m[1][0])) &&
            (fabs(img->qto_xyz.m[0][0])>fabs(img->qto_xyz.m[2][0]))) {
      if(verbose>1) fprintf(stdout,"  x is good\n");
    } else if((fabs(img->qto_xyz.m[1][0])>fabs(img->qto_xyz.m[0][0])) &&
              (fabs(img->qto_xyz.m[1][0])>fabs(img->qto_xyz.m[2][0]))) {
      if(verbose>1) fprintf(stdout,"  x and y swap\n");
      NIIK_RET0((!niik_image_restructure(img,"+y+x+z")),fcname,"niik_image_restructure xy swap");
      if(verbose>1) {
        fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
        fprintf(stdout,"\tqto_xyz = \n");
        mat44_display(img->qto_xyz);
      }
    } else if((fabs(img->qto_xyz.m[2][0])>fabs(img->qto_xyz.m[0][0])) &&
              (fabs(img->qto_xyz.m[2][0])>fabs(img->qto_xyz.m[1][0]))) {
      if(verbose>1) fprintf(stdout,"  x and z swap\n");
      NIIK_RET0((!niik_image_restructure(img,"+z+y+x")),fcname,"niik_image_restructure xz swap");
      if(verbose>1) {
        fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
        fprintf(stdout,"\tqto_xyz = \n");
        mat44_display(img->qto_xyz);
      }
    }

    if     ((fabs(img->qto_xyz.m[1][1])>fabs(img->qto_xyz.m[0][1])) &&
            (fabs(img->qto_xyz.m[1][1])>fabs(img->qto_xyz.m[2][1]))) {
      if(verbose>1) fprintf(stdout,"  y is good\n");
    } else if((fabs(img->qto_xyz.m[0][1])>fabs(img->qto_xyz.m[1][1])) &&
              (fabs(img->qto_xyz.m[0][1])>fabs(img->qto_xyz.m[2][1]))) {
      if(verbose>1) fprintf(stdout,"  x and y swap\n");
      NIIK_RET0((!niik_image_restructure(img,"+y+x+z")),fcname,"niik_image_restructure xy swap");
      if(verbose>1) {
        fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
        fprintf(stdout,"\tqto_xyz = \n");
        mat44_display(img->qto_xyz);
      }
    } else if((fabs(img->qto_xyz.m[2][1])>fabs(img->qto_xyz.m[0][1])) &&
              (fabs(img->qto_xyz.m[2][1])>fabs(img->qto_xyz.m[1][1]))) {
      if(verbose>1) fprintf(stdout,"  y and z swap\n");
      NIIK_RET0((!niik_image_restructure(img,"+x+z+y")),fcname,"niik_image_restructure yz swap");
      if(verbose>1) {
        fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
        fprintf(stdout,"\tqto_xyz = \n");
        mat44_display(img->qto_xyz);
      }
    }
    if(verbose>1) fprintf(stdout,"  done with reorientation\n");
  }

  /* actually calculate the voxel size */
  img->sform_code = NIFTI_XFORM_SCANNER_ANAT;
  if(verbose>1) fprintf(stdout,"[%s] calculating world locations\n",fcname);
  voxel_position[0] = voxel_position[1] = voxel_position[2] = 0;
  miconvert_voxel_to_world ( minc_volume, voxel_position, world_position);
  if(verbose>1) fprintf(stdout,"%12.7f %12.7f %12.7f   ->    %12.7f %12.7f %12.7f\n",
                          voxel_position[0],voxel_position[1],voxel_position[2],
                          world_position[0],world_position[1],world_position[2] );
  for(i=0; i<3; i++) img->sto_xyz.m[i][3]=world_position[i];
  qt = niikpt_val(img->sto_xyz.m[0][3],img->sto_xyz.m[1][3],img->sto_xyz.m[2][3],0);

  voxel_position[0] = 10.0;
  voxel_position[1] = voxel_position[2] = 0;
  miconvert_voxel_to_world ( minc_volume,voxel_position,world_position);
  if(verbose>1) fprintf(stdout,"%12.7f %12.7f %12.7f   ->    %12.7f %12.7f %12.7f  (%7.3f %7.3f %7.3f)\n",
                          voxel_position[0],voxel_position[1],voxel_position[2],
                          world_position[0],world_position[1],world_position[2],
                          world_position[0]-img->sto_xyz.m[0][3],
                          world_position[1]-img->sto_xyz.m[1][3],
                          world_position[2]-img->sto_xyz.m[2][3] );
  pt = niikpt_sub(niikpt_val(world_position[0],world_position[1],world_position[2],0),qt);
  if     ((fabs(pt.x)>=fabs(pt.y)) && (fabs(pt.x)>=fabs(pt.z))) {
    n=0;
  } else if((fabs(pt.y)>=fabs(pt.x)) && (fabs(pt.y)>=fabs(pt.z))) {
    n=1;
  } else if((fabs(pt.z)>=fabs(pt.x)) && (fabs(pt.z)>=fabs(pt.y))) {
    n=2;
  }
  for(i=0; i<3; i++) img->sto_xyz.m[i][n]=(world_position[i]-img->sto_xyz.m[i][3])/10.0;

  voxel_position[1] = 10.0;
  voxel_position[0] = voxel_position[2] = 0;
  miconvert_voxel_to_world ( minc_volume,voxel_position,world_position);
  if(verbose>1) fprintf(stdout,"%12.7f %12.7f %12.7f   ->    %12.7f %12.7f %12.7f  (%7.3f %7.3f %7.3f)\n",
                          voxel_position[0],voxel_position[1],voxel_position[2],
                          world_position[0],world_position[1],world_position[2],
                          world_position[0]-img->sto_xyz.m[0][3],
                          world_position[1]-img->sto_xyz.m[1][3],
                          world_position[2]-img->sto_xyz.m[2][3] );
  pt = niikpt_sub(niikpt_val(world_position[0],world_position[1],world_position[2],0),qt);
  if     ((fabs(pt.x)>=fabs(pt.y)) && (fabs(pt.x)>=fabs(pt.z))) {
    n=0;
  } else if((fabs(pt.y)>=fabs(pt.x)) && (fabs(pt.y)>=fabs(pt.z))) {
    n=1;
  } else if((fabs(pt.z)>=fabs(pt.x)) && (fabs(pt.z)>=fabs(pt.y))) {
    n=2;
  }
  for(i=0; i<3; i++) img->sto_xyz.m[i][n]=(world_position[i]-img->sto_xyz.m[i][3]) / 10.0;

  voxel_position[2] = 10.0;
  voxel_position[0] = voxel_position[1] = 0;
  miconvert_voxel_to_world ( minc_volume,voxel_position,world_position);
  if(verbose>1) fprintf(stdout,"%12.7f %12.7f %12.7f   ->    %12.7f %12.7f %12.7f  (%7.3f %7.3f %7.3f)\n",
                          voxel_position[0],voxel_position[1],voxel_position[2],
                          world_position[0],world_position[1],world_position[2],
                          world_position[0]-img->sto_xyz.m[0][3],
                          world_position[1]-img->sto_xyz.m[1][3],
                          world_position[2]-img->sto_xyz.m[2][3] );
  pt = niikpt_sub(niikpt_val(world_position[0],world_position[1],world_position[2],0),qt);
  if     ((fabs(pt.x)>=fabs(pt.y)) && (fabs(pt.x)>=fabs(pt.z))) {
    n=0;
  } else if((fabs(pt.y)>=fabs(pt.x)) && (fabs(pt.y)>=fabs(pt.z))) {
    n=1;
  } else if((fabs(pt.z)>=fabs(pt.x)) && (fabs(pt.z)>=fabs(pt.y))) {
    n=2;
  }
  for(i=0; i<3; i++) img->sto_xyz.m[i][n]=(world_position[i]-img->sto_xyz.m[i][3])/10.0;
  img->sto_xyz.m[3][3]=1;
  mm=img->sto_xyz;

  if(verbose>=1) {
    fprintf(stdout,"\tsto_xyz = \n");
    mat44_display(img->sto_xyz);
  }
  img->dx = img->pixdim[1] = sqrt(NIIK_SQ(img->sto_xyz.m[0][0]) +
                                  NIIK_SQ(img->sto_xyz.m[1][0]) +
                                  NIIK_SQ(img->sto_xyz.m[2][0])) ;
  img->dy = img->pixdim[2] = sqrt(NIIK_SQ(img->sto_xyz.m[0][1]) +
                                  NIIK_SQ(img->sto_xyz.m[1][1]) +
                                  NIIK_SQ(img->sto_xyz.m[2][1])) ;
  img->dz = img->pixdim[3] = sqrt(NIIK_SQ(img->sto_xyz.m[0][2]) +
                                  NIIK_SQ(img->sto_xyz.m[1][2]) +
                                  NIIK_SQ(img->sto_xyz.m[2][2]));

  for(i=0; i<3; i++) img->qto_xyz.m[i][3]=img->sto_xyz.m[i][3];
  mm=img->sto_xyz;

  if(img->sto_xyz.m[0][0]<0) {
    if(verbose>1) fprintf(stdout,"  flip x\n");
    NIIK_RET0((!niik_image_restructure(img,"-x+y+z")),fcname,"niik_image_restructure x flip");
    pt=niikpt_affine_transform_m44(mm,niikpt_val(img->nx-1,0,0,0));
    /*niikpt_disp(pt);*/
    for(n=0; n<3; n++) mm.m[n][0]*=-1;
    mm.m[0][3]=pt.x;
    mm.m[1][3]=pt.y;
    mm.m[2][3]=pt.z;
    img->sto_xyz=mm;
    if(verbose>=1) {
      fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
      fprintf(stdout,"\tsto_xyz = \n");
      mat44_display(img->sto_xyz);
    }
  }
  if(img->sto_xyz.m[1][1]<0) {
    if(verbose>1) fprintf(stdout,"  flip y\n");
    NIIK_RET0((!niik_image_restructure(img,"+x-y+z")),fcname,"niik_image_restructure y flip");
    pt=niikpt_affine_transform_m44(mm,niikpt_val(0,img->ny-1,0,0));
    /*niikpt_disp(pt);*/
    for(n=0; n<3; n++) mm.m[n][1]*=-1;
    mm.m[0][3]=pt.x;
    mm.m[1][3]=pt.y;
    mm.m[2][3]=pt.z;
    img->sto_xyz=mm;
    if(verbose>1) {
      fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
      fprintf(stdout,"\tsto_xyz = \n");
      mat44_display(img->sto_xyz);
    }
  }
  if(img->sto_xyz.m[2][2]<0) {
    if(verbose>1) fprintf(stdout,"  flip z\n");
    NIIK_RET0((!niik_image_restructure(img,"+x+y-z")),fcname,"niik_image_restructure z flip");
    pt=niikpt_affine_transform_m44(mm,niikpt_val(0,0,img->nz-1,0));
    /*niikpt_disp(pt);*/
    for(n=0; n<3; n++) mm.m[n][2]*=-1;
    mm.m[0][3]=pt.x;
    mm.m[1][3]=pt.y;
    mm.m[2][3]=pt.z;
    img->sto_xyz=mm;
    if(verbose>=1) {
      fprintf(stdout,"\tdim %3i %3i %3i\n",img->nx,img->ny,img->nz);
      fprintf(stdout,"\tsto_xyz = \n");
      mat44_display(img->qto_xyz);
    }
  }

  img->qto_ijk = nifti_mat44_inverse(img->qto_xyz);
  nifti_mat44_to_quatern(img->qto_xyz,
                         &img->quatern_b,&img->quatern_c,&img->quatern_d,
                         &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                         NULL,NULL,NULL,
                         &img->qfac);
  img->sto_ijk = nifti_mat44_inverse(img->sto_xyz);
  img->qform_code=0;


  if(verbose>1) {
    fprintf(stdout,"[%s] final q and s matrices\n",fcname);
    fprintf(stdout,"\tqto_xyz = \n");
    mat44_display(img->qto_xyz);
    fprintf(stdout,"\tsto_xyz = \n");
    mat44_display(img->sto_xyz);
  }

  if(verbose>1) {
    fprintf(stdout,"\tsto_xyz = \n");
    mat44_display(img->sto_xyz);
  }

  /*if(verbose>=1) {fprintf(stdout,"\tminc_volume = %i\n",minc_volume);}*/
  /*img->dx = img->pixdim[1] = img->sto_xyz.m[0][0] / img->qto_xyz.m[0][0];
    img->dy = img->pixdim[2] = img->sto_xyz.m[1][1] / img->qto_xyz.m[1][1];
    img->dz = img->pixdim[3] = img->sto_xyz.m[2][2] / img->qto_xyz.m[2][2];*/
  if(verbose>1)fprintf(stdout,"[%s]   pixel size %12.9f %12.9f %12.9f\n",fcname,img->dx,img->dy,img->dz);

  /*if((result = miget_volume_range(minc_volume,
                                  &delta[0],
                                  &delta[1]) == MI_ERROR)){
    fprintf(stderr,"[%s] ERROR: miget_volume_range\n",fcname);
    return NULL; }*/

  if(verbose>1) {
    /*fprintf(stdout,"\tminc_volume = %i\n",minc_volume);*/
    fprintf(stdout,"[%s] max, min = %12.7f %12.7f\n",
            fcname,delta[0],delta[1]);
  }

  if(0) {
    img->dx = img->pixdim[1] =
                NIIK_DMAX(fabs(img->sto_xyz.m[0][0]),NIIK_DMAX(fabs(img->sto_xyz.m[1][0]),fabs(img->sto_xyz.m[2][0]))) /
                NIIK_DMAX(fabs(img->qto_xyz.m[0][0]),NIIK_DMAX(fabs(img->qto_xyz.m[1][0]),fabs(img->qto_xyz.m[2][0])));
    img->dy = img->pixdim[2] =
                NIIK_DMAX(fabs(img->sto_xyz.m[0][1]),NIIK_DMAX(fabs(img->sto_xyz.m[1][1]),fabs(img->sto_xyz.m[2][1]))) /
                NIIK_DMAX(fabs(img->qto_xyz.m[0][1]),NIIK_DMAX(fabs(img->qto_xyz.m[1][1]),fabs(img->qto_xyz.m[2][1])));
    img->dz = img->pixdim[3] =
                NIIK_DMAX(fabs(img->sto_xyz.m[0][2]),NIIK_DMAX(fabs(img->sto_xyz.m[1][2]),fabs(img->sto_xyz.m[2][2]))) /
                NIIK_DMAX(fabs(img->qto_xyz.m[0][2]),NIIK_DMAX(fabs(img->qto_xyz.m[1][2]),fabs(img->qto_xyz.m[2][2])));
    if(verbose>=1)fprintf(stdout,"[%s]   pixel size %12.5f %12.5f %12.5f\n",fcname,img->dx,img->dy,img->dz);
  }
  if(verbose>=1)fprintf(stdout,"[%s]   pixel size %12.5f %12.5f %12.5f\n",fcname,img->dx,img->dy,img->dz);


  if(verbose>1) {
    for(i=0; i<3; i++)
      fprintf(stdout,"[%s] start size [%i] = %5i %5i\n",
              fcname,i,(int)start[i],(int)size[i]);
  }

  free(dims);
  free(count);
  free(delta);
  free(size);
  free(start);

  img->fname=(char *)calloc(4096,sizeof(char));
  img->iname=(char *)calloc(4096,sizeof(char));
  strncpy(img->fname,minc_name,4096);
  strncpy(img->iname,minc_name,4096);

  /*Read MINC history*/
  if(miget_attr_length(minc_volume,"","history",&minc_history_length) == MI_NOERROR) {
    if(minc_history_length>0) {
      char *minc_history=malloc(minc_history_length+1);
      if( miget_attr_values(minc_volume,MI_TYPE_STRING,"","history",minc_history_length+1,minc_history) == MI_NOERROR) {
        img->minc_history=minc_history;
      } else {
        /*can't get minc history*/
        fprintf(stderr,"[%s] Can't read minc history\n",fcname);
        free(minc_history);
      }
    }
  }

  miclose_volume(minc_volume); /* Oops! KN January 27, 2014 */

  if(verbose>=1) niik_fc_display(fcname,0);
  return img;
} /* niik_image_read_minc2 */


nifti_image *niik_image_read_minc(const char *minc_name) {
  /* NIFTI stuff */
  nifti_image *nii_ptr=NULL;
  int nii_dimids[MAX_NII_DIMS];
  int nii_dir[MAX_NII_DIMS];
  int nii_map[MAX_NII_DIMS];
  unsigned long nii_lens[MAX_NII_DIMS];
  int nii_ndims;
  static int nifti_filetype;
  static int nifti_datatype;
  static int nifti_signed = 1;

  /* MINC stuff */
  /*mihandle_t  minc_volume;     MINC volume */
  int mnc_fd;                 /* MINC file descriptor */
  nc_type mnc_type;           /* MINC data type as read */
  int mnc_ndims;              /* MINC image dimension count */
  int mnc_dimids[MAX_VAR_DIMS]; /* MINC image dimension identifiers */
  long mnc_dlen;              /* MINC dimension length value */
  double mnc_dstep;           /* MINC dimension step value */
  int mnc_icv;                /* MINC image conversion variable */
  int mnc_vid;                /* MINC Image variable ID */
  long mnc_start[MAX_VAR_DIMS]; /* MINC data starts */
  long mnc_count[MAX_VAR_DIMS]; /* MINC data counts */
  int mnc_signed;             /* MINC if output voxels are signed */
  double mnc_rrange[2];       /* MINC real range (min, max) */
  double mnc_vrange[2];       /* MINC valid range (min, max) */

  /* Other stuff */
  char out_str[4096];         /* Big string for filename */
  char att_str[65536];         /* Big string for attribute values */
  nc_type att_datatype;
  int att_length;
  int i;                      /* Generic loop counter the first */
  int j;                      /* Generic loop counter the second */
  char *str_ptr;              /* Generic ASCIZ string pointer */
  int r;                      /* Result code. */
  static int qflag = 0;       /* Quiet flag (default is non-quiet) */

  int id;
  double start;
  double step;
  double dircos[MAX_SPACE_DIMS];
  int tmp;

  int verbose=niik_verbose();
  char fcname[64]="niik_image_mnc2nii";
  int header_only = 1, created_tempfile;
  char *tempfile;
  int Is_MINC2_File = 0;

  ncopts = 0;                 /* Clear global netCDF error reporting flag */

  /* Default NIfTI file type is "NII", single binary file
   */
  nifti_filetype = FT_UNSPECIFIED;
  nifti_datatype = DT_UNKNOWN;

  if(!nifti_signed) {
    switch (nifti_datatype) {
    case DT_INT8:
      nifti_datatype = DT_UINT8;
      if(verbose>0) fprintf(stdout,"[%s] DT_INT8\n",fcname);
      break;
    case DT_INT16:
      nifti_datatype = DT_UINT16;
      if(verbose>0) fprintf(stdout,"[%s] DT_INT16\n",fcname);
      break;
    case DT_INT32:
      nifti_datatype = DT_UINT32;
      if(verbose>0) fprintf(stdout,"[%s] DT_INT32\n",fcname);
      break;
    }
  }

  switch (nifti_datatype) {
  case DT_INT8:
  case DT_UINT8:
    mnc_type = NC_BYTE;
    if(verbose>0) fprintf(stdout,"[%s] minc byte\n",fcname);
    break;
  case DT_INT16:
  case DT_UINT16:
    mnc_type = NC_SHORT;
    if(verbose>0) fprintf(stdout,"[%s] minc short\n",fcname);
    break;
  case DT_INT32:
  case DT_UINT32:
    mnc_type = NC_INT;
    if(verbose>0) fprintf(stdout,"[%s] minc int\n",fcname);
    break;
  case DT_FLOAT32:
    mnc_type = NC_FLOAT;
    if(verbose>0) fprintf(stdout,"[%s] minc float\n",fcname);
    break;
  case DT_FLOAT64:
    mnc_type = NC_DOUBLE;
    if(verbose>0) fprintf(stdout,"[%s] minc double\n",fcname);
    break;
  }

  strncpy(out_str, minc_name,sizeof(out_str)-1);

  if (nifti_filetype == FT_UNSPECIFIED) {
    nifti_filetype = FT_NIFTI_SINGLE;
  }

  /* Open the MINC file.  It needs to exist.
   */

  /* Expand file */
  tempfile = miexpand_file((char *)minc_name, NULL, header_only, &created_tempfile);
  if (tempfile == NULL) {
    (void) fprintf(stderr, "%s: Error expanding file \"%s\"\n",
                   fcname, minc_name);
    return NULL;
  }

  /* Open the file */

  if((mnc_fd = miopen(tempfile, NC_NOWRITE)) == MI_ERROR) {
    fprintf(stderr, "[%s] ERROR: Can't find input file '%s'\n", fcname,tempfile);
    return NULL;
  }

  /* check whether the file is Version 2 */
#ifdef MINC2
  if (MI2_ISH5OBJ(mnc_fd)) {
    Is_MINC2_File = 1;
  }
#endif

  if(verbose>=1) {
    if (Is_MINC2_File) {
      (void) printf("Version: 2 (HDF5)\n");
    } else {
      (void) printf("Version: 1 (netCDF)\n");
    }
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] mnc_fd = %i\n",fcname,mnc_fd);
  }

  mnc_vid = ncvarid(mnc_fd, MIimage);
  if (mnc_vid == MI_ERROR) {
    miclose(mnc_fd);
    if(verbose>=1) fprintf(stdout,"[%s] niik_image_read_minc2\n",fcname);
    if((nii_ptr = niik_image_read_minc2(minc_name))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_read_minc2\n",fcname);
      return NULL;
    }

    if(verbose>=1) niik_fc_display(fcname,0);
    return nii_ptr;
  }

  /* Find the MINC image variable.  If we can't find it, there is no
   * further processing possible...
   */
  mnc_vid = ncvarid(mnc_fd, MIimage);
  if (mnc_vid < 0) {
    fprintf(stderr, "[%s] ERROR: Can't locate the image variable (mnc_vid=%d)\n", fcname,mnc_vid);
    return NULL;
  }

  /* Find out about the MINC image variable - specifically, how many
   * dimensions, and which dimensions.
   */
  r = ncvarinq(mnc_fd, mnc_vid, NULL, NULL, &mnc_ndims, mnc_dimids, NULL);
  if (r < 0) {
    fprintf(stderr, "[%s] ERROR: Can't read information from image variable\n",fcname);
    return NULL;
  }
  if (mnc_ndims > MAX_NII_DIMS) {
    fprintf(stderr, "[%s] ERROR: NIfTI-1 files may contain at most %d dimensions\n", fcname,
            MAX_NII_DIMS);
    return NULL;
  }

  /* Initialize the NIfTI structure
   */
  nii_ptr = (nifti_image *)calloc(1,sizeof(nifti_image));
  init_nifti_header(nii_ptr);

  /* For now we just use the mnc2nii command line as the description
   * field.  Probably we should use something better, perhaps a
   * combination of some other standard MINC fields that might
   * provide more information.
   */
  str_ptr = nii_ptr->descrip;
  str_ptr[0] = '\0';

  nii_ptr->fname = calloc(strlen(out_str) + 4 + 1,sizeof(char));
  nii_ptr->iname = calloc(strlen(out_str) + 4 + 1,sizeof(char));
  strcpy(nii_ptr->fname, out_str);
  strcpy(nii_ptr->iname, out_str);
  strcat(nii_ptr->fname, ".nii");
  strcat(nii_ptr->iname, ".nii");

  miget_image_range(mnc_fd, mnc_rrange); /* Get real range */
  miget_valid_range(mnc_fd, mnc_vid, mnc_vrange); /* Get voxel range */

  if(verbose>=1) {
    fprintf(stdout,"[%s] miget_valid_range mnc_vrange %17.12f %17.12f\n",fcname,mnc_vrange[0],mnc_vrange[1]);
    fprintf(stdout,"[%s] miget_valid_range mnc_rrange %17.12f %17.12f\n",fcname,mnc_rrange[0],mnc_rrange[1]);
  }

  if (mnc_vrange[1] != mnc_vrange[0] && mnc_rrange[1] != mnc_rrange[0]) {
    if(verbose>=1) {
      fprintf(stdout,"[%s] nifti scale    %17.12f %17.12f\n",fcname,nii_ptr->scl_slope,nii_ptr->scl_inter);
    }
  } else {
    nii_ptr->scl_slope = 0.0;
  }

  nii_ptr->nvox = 1;          /* Initial value for voxel count */

  /* Find all of the dimensions of the MINC file, in the order they
   * will be listed in the NIfTI-1/Analyze file.  We use this to build
   * a map for restructuring the data according to the normal rules
   * of NIfTI-1.
   */
  nii_ndims = 0;
  for (i = 0; i < MAX_NII_DIMS; i++) {
    if (dimnames[i] == NULL) {
      nii_dimids[nii_ndims] = -1;
      continue;
    }

    nii_dimids[nii_ndims] = ncdimid(mnc_fd, dimnames[i]);
    if (nii_dimids[nii_ndims] == -1) {
      continue;
    }

    /* Make sure the dimension is actually used to define the image.
     */
    for (j = 0; j < mnc_ndims; j++) {
      if (nii_dimids[nii_ndims] == mnc_dimids[j]) {
        nii_map[nii_ndims] = j;
        break;
      }
    }

    if (j < mnc_ndims) {
      mnc_dlen = 1;
      mnc_dstep = 0;

      ncdiminq(mnc_fd, nii_dimids[nii_ndims], NULL, &mnc_dlen);
      ncattget(mnc_fd, ncvarid(mnc_fd, dimnames[i]), MIstep, &mnc_dstep);

      if (mnc_dstep < 0) {
        nii_dir[nii_ndims] = -1;
        mnc_dstep = -mnc_dstep;
      } else {
        nii_dir[nii_ndims] = 1;
      }

      nii_lens[nii_ndims] = mnc_dlen;
      nii_ndims++;
    }

    nii_ptr->dim[dimmap[i]] = (int) mnc_dlen;
    nii_ptr->nvox *= mnc_dlen;

    nii_ptr->pixdim[dimmap[i]] = (float) mnc_dstep;
  }

  /* Here we do some "post-processing" of the results. Make certain that
   * the nt value is never zero, and make certain that ndim is set to
   * 4 if there is a time dimension and 5 if there is a vector dimension
   */

  if (nii_ptr->dim[3] > 1 && nii_ndims < 4) {
    nii_ndims = 4;
  }

  if (nii_ptr->dim[4] > 1) {
    nii_ptr->intent_code = NIFTI_INTENT_VECTOR;
    nii_ndims = 5;
  }

  nii_ptr->ndim = nii_ndims; /* Total number of dimensions in file */
  nii_ptr->nx = nii_ptr->dim[0];
  nii_ptr->ny = nii_ptr->dim[1];
  nii_ptr->nz = nii_ptr->dim[2];
  nii_ptr->nt = nii_ptr->dim[3];
  nii_ptr->nu = nii_ptr->dim[4];

  nii_ptr->dx = nii_ptr->pixdim[0];
  nii_ptr->dy = nii_ptr->pixdim[1];
  nii_ptr->dz = nii_ptr->pixdim[2];
  nii_ptr->dt = nii_ptr->pixdim[3];
  nii_ptr->du = 1; /* MINC files don't define a sample size for a vector_dimension */

  nii_ptr->nifti_type = nifti_filetype;

  if (nifti_datatype == DT_UNKNOWN) {
    nii_ptr->datatype = DT_FLOAT32; /* Default */
    mnc_type = NC_FLOAT;
    mnc_signed = 1;
  } else {
    nii_ptr->datatype = nifti_datatype;
  }


  /* Load the direction_cosines and start values into the NIfTI-1
   * sform structure.
   *
   */
  for (i = 0; i < MAX_SPACE_DIMS; i++) {
    id = ncvarid(mnc_fd, mnc_spatial_names[i]);

    if (id < 0) {
      continue;
    }

    /* Set default values */
    start = 0.0;
    step = 1.0;
    dircos[DIM_X] = dircos[DIM_Y] = dircos[DIM_Z] = 0.0;
    dircos[i] = 1.0;

    miattget(mnc_fd, id, MIstart, NC_DOUBLE, 1, &start, &tmp);
    miattget(mnc_fd, id, MIstep,  NC_DOUBLE, 1, &step, &tmp);
    miattget(mnc_fd, id, MIdirection_cosines, NC_DOUBLE, MAX_SPACE_DIMS,
             dircos, &tmp);
    ncdiminq(mnc_fd, ncdimid(mnc_fd, mnc_spatial_names[i]), NULL,
             &mnc_dlen);

    nii_ptr->qform_code = NIFTI_XFORM_SCANNER_ANAT;
    nii_ptr->qto_xyz.m[i][0] = dircos[0];
    nii_ptr->qto_xyz.m[i][1] = dircos[1];
    nii_ptr->qto_xyz.m[i][2] = dircos[2];
    nii_ptr->qto_xyz.m[i][3] = start;

    if (step < 0) {
      step = -step;
      start = start - step * (mnc_dlen - 1);
    }

    nii_ptr->sto_xyz.m[0][i] = step * dircos[0];
    nii_ptr->sto_xyz.m[1][i] = step * dircos[1];
    nii_ptr->sto_xyz.m[2][i] = step * dircos[2];

    nii_ptr->sto_xyz.m[0][3] += start * dircos[0];
    nii_ptr->sto_xyz.m[1][3] += start * dircos[1];
    nii_ptr->sto_xyz.m[2][3] += start * dircos[2];

    if( ncattinq(mnc_fd, id, MIspacetype, &att_datatype, &att_length) != MI_ERROR &&
        att_datatype== NC_CHAR &&
        att_length<sizeof(att_str) )
      miattgetstr(mnc_fd, id, MIspacetype, sizeof(att_str), att_str);
    else
      *att_str=0; /*empty */

    /* Try to set the S-transform code correctly.
     */
    if (!strcmp(att_str, MI_TALAIRACH)) {
      nii_ptr->sform_code = NIFTI_XFORM_TALAIRACH;
    } else if (!strcmp(att_str, MI_CALLOSAL)) {
      /* TODO: Not clear what do do here... */
      nii_ptr->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    } else {                /* MI_NATIVE or unknown */
      nii_ptr->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    }
  }

  /* So the last row is right... */
  nii_ptr->sto_xyz.m[3][0] = 0.0;
  nii_ptr->sto_xyz.m[3][1] = 0.0;
  nii_ptr->sto_xyz.m[3][2] = 0.0;
  nii_ptr->sto_xyz.m[3][3] = 1.0;

  nii_ptr->sto_ijk = nifti_mat44_inverse(nii_ptr->sto_xyz);

  nii_ptr->qto_xyz.m[3][0] = 0.0;
  nii_ptr->qto_xyz.m[3][1] = 0.0;
  nii_ptr->qto_xyz.m[3][2] = 0.0;
  nii_ptr->qto_xyz.m[3][3] = 1.0;

  nii_ptr->qto_ijk = nifti_mat44_inverse(nii_ptr->qto_xyz);
  nifti_mat44_to_quatern( nii_ptr->qto_xyz,
                          &nii_ptr->quatern_b, &nii_ptr->quatern_c, &nii_ptr->quatern_d,
                          &nii_ptr->qoffset_x, &nii_ptr->qoffset_y, &nii_ptr->qoffset_z,
                          NULL, NULL, NULL,
                          &nii_ptr->qfac);

  nifti_datatype_sizes(nii_ptr->datatype,
                       &nii_ptr->nbyper, &nii_ptr->swapsize);


  if (!qflag && verbose>=1) {
    nifti_image_infodump(nii_ptr);
  }

  /* Now load the actual MINC data. */
  if(verbose>=1) {
    fprintf(stdout,"  nbyper = %i\n",nii_ptr->nbyper);
    fprintf(stdout,"  nvox     %i\n",nii_ptr->nvox);
  }

  nii_ptr->data = malloc(nii_ptr->nbyper * nii_ptr->nvox);
  if (nii_ptr->data == NULL) {
    fprintf(stderr, "[%s] ERROR: Out of memory.\n",fcname);
    return NULL;
  }

  if (!qflag && verbose>=1) {
    fprintf(stderr, "MINC type %d signed %d\n", mnc_type, mnc_signed);
  }

  mnc_icv = miicv_create();
  miicv_setint(mnc_icv, MI_ICV_TYPE, mnc_type);
  miicv_setstr(mnc_icv, MI_ICV_SIGN, (mnc_signed) ? MI_SIGNED : MI_UNSIGNED);
  miicv_setdbl(mnc_icv, MI_ICV_VALID_MAX, mnc_vrange[1]);
  miicv_setdbl(mnc_icv, MI_ICV_VALID_MIN, mnc_vrange[0]);
  miicv_setint(mnc_icv, MI_ICV_DO_NORM, 1);

  miicv_attach(mnc_icv, mnc_fd, mnc_vid);

  /* Read in the entire hyperslab from the file.
   */
  for (i = 0; i < mnc_ndims; i++) {
    ncdiminq(mnc_fd, mnc_dimids[i], NULL, &mnc_count[i]);
    mnc_start[i] = 0;
  }

  for(i=0; i<MAX_VAR_DIMS && verbose>=1; i++) {
    if(mnc_start[i]>0 || mnc_count[i]>0)
      fprintf(stdout,"[%s]   mnc_start[%i] = %i   |   mnc_count[%i] = %i\n",fcname,i,(int)mnc_start[i],i,(int)mnc_count[i]);
  }

  r = miicv_get(mnc_icv, mnc_start, mnc_count, nii_ptr->data);
  if (r < 0) {
    fprintf(stderr, "[%s] ERROR: Read error\n",fcname);
    return NULL;
  }

  /* Shut down the MINC stuff now that it has done its work.
   */
  miicv_detach(mnc_icv);
  miicv_free(mnc_icv);
  miclose(mnc_fd);

  if (!qflag && verbose>=1) {
    /* Debugging stuff - just to check the contents of these arrays.
     */
    for (i = 0; i < nii_ndims; i++) {
      printf("%d: %ld %d %d\n",
             i, nii_lens[i], nii_map[i], nii_dir[i]);
    }
    printf("bytes per voxel %d\n", nii_ptr->nbyper);
    printf("# of voxels %d\n", nii_ptr->nvox);
  }

  /* Rearrange the data to correspond to the NIfTI dimension ordering.
   */
  restructure_array(nii_ndims,
                    nii_ptr->data,
                    nii_lens,
                    nii_ptr->nbyper,
                    nii_map,
                    nii_dir);

  if (!qflag && verbose>=1) {
    /* More debugging stuff - check coordinate transform.
     */
    test_xform(nii_ptr->sto_xyz, 0, 0, 0);
    test_xform(nii_ptr->sto_xyz, 10, 0, 0);
    test_xform(nii_ptr->sto_xyz, 0, 10, 0);
    test_xform(nii_ptr->sto_xyz, 0, 0, 10);
    test_xform(nii_ptr->sto_xyz, 10, 10, 10);
  }

  if(verbose>=1) fprintf(stdout, "Calling NIFTI-1 Write routine\n");
  /*nifti_image_write(nii_ptr);*/

  if(verbose>=1) {
    fprintf(stdout,"  dim    %i %i    \n",nii_ptr->ndim,nii_ptr->dim[0]);
    fprintf(stdout,"  dim[1] %i %i    \n",nii_ptr->nx,nii_ptr->dim[1]);
    fprintf(stdout,"  dim[2] %i %i    \n",nii_ptr->ny,nii_ptr->dim[2]);
    fprintf(stdout,"  dim[3] %i %i    \n",nii_ptr->nz,nii_ptr->dim[3]);
  }
  nii_ptr->dim[0] = nii_ptr->ndim;
  nii_ptr->dim[1] = nii_ptr->nx;
  nii_ptr->dim[2] = nii_ptr->ny;
  nii_ptr->dim[3] = nii_ptr->nz;
  nii_ptr->dim[4] = nii_ptr->nt;
  nii_ptr->dim[5] = nii_ptr->nu;
  nii_ptr->dim[6] = nii_ptr->nv;
  nii_ptr->dim[7] = nii_ptr->nw;

  nii_ptr->pixdim[1] = nii_ptr->dx;
  nii_ptr->pixdim[2] = nii_ptr->dy;
  nii_ptr->pixdim[3] = nii_ptr->dz;
  nii_ptr->pixdim[4] = nii_ptr->dt;
  nii_ptr->pixdim[5] = nii_ptr->du;
  nii_ptr->pixdim[6] = nii_ptr->dv;
  nii_ptr->pixdim[7] = nii_ptr->dw;

  for(i=1,nii_ptr->nvox=1; i<=nii_ptr->ndim; i++) {
    nii_ptr->nvox *= nii_ptr->dim[i];
  }
  if(verbose>=1) nifti_image_infodump( nii_ptr );

  NIIK_RET0((!niik_image_type_convert_scl(nii_ptr,nii_ptr->datatype,1)),fcname,
            "niik_image_type_convert_scl");

  /* delete temp file */
  if(!strcmp(tempfile,minc_name)) {
    if(verbose>=1) {
      fprintf(stdout,"[%s] tempfile  = %s\n",fcname,tempfile);
      fprintf(stdout,"[%s] minc_name = %s\n",fcname,minc_name);
    }
  } else {
    if(verbose>=1) {
      fprintf(stdout,"[%s] not same?\n",fcname);
      fprintf(stdout,"[%s] tempfile  = %s\n",fcname,tempfile);
      fprintf(stdout,"[%s] minc_name = %s\n",fcname,minc_name);
      fprintf(stdout,"[%s] deleting %s\n",fcname,tempfile);
    }
    remove(tempfile);
  }
  free(tempfile);

  return nii_ptr;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/