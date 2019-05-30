#ifndef _FALCON_NII2MNC_C_
#define _FALCON_NII2MNC_C_

#include "falcon.h"

#include <limits.h>
#include <float.h>
#include <minc.h>
#include <minc2.h>
#include "volume_io.h"

#undef X /* These are used in nifti1_io */
#undef Y
#undef Z

#include <time_stamp.h>
#include "nifti1_io.h"
#include "dbh.h"
#include "nifti1_local.h"

extern void test_xform(mat44 m, int i, int j, int k);

void niik_miinit_default_range(mitype_t mitype, double *valid_max, double *valid_min) {
  switch (mitype) {
  case MI_TYPE_BYTE:
    *valid_min = CHAR_MIN;
    *valid_max = CHAR_MAX;
    break;
  case MI_TYPE_SHORT:
    *valid_min = SHRT_MIN;
    *valid_max = SHRT_MAX;
    break;
  case MI_TYPE_INT:
    *valid_min = INT_MIN;
    *valid_max = INT_MAX;
    break;
  case MI_TYPE_UBYTE:
    *valid_min = 0;
    *valid_max = UCHAR_MAX;
    break;
  case MI_TYPE_USHORT:
    *valid_min = 0;
    *valid_max = USHRT_MAX;
    break;
  case MI_TYPE_UINT:
    *valid_min = 0;
    *valid_max = UINT_MAX;
    break;
  case MI_TYPE_FLOAT:
    *valid_min = -FLT_MAX;
    *valid_max = FLT_MAX;
    break;
  case MI_TYPE_DOUBLE:
    *valid_min = -DBL_MAX;
    *valid_max = DBL_MAX;
    break;
  default:
    *valid_min = 0;
    *valid_max = 1;
    break;
  }
}

static void find_data_range(int datatype,
                            long nvox,
                            void *data,
                            double range[2]) {
  int i;
  range[0] = DBL_MAX;
  range[1] = -DBL_MAX;
  for (i = 0; i < nvox; i++) {
    double tmp;
    switch (datatype) {
    case DT_INT8:
      tmp = (double) ((char *)data)[i];
      break;
    case DT_UINT8:
      tmp = (double) ((unsigned char *)data)[i];
      break;
    case DT_INT16:
      tmp = (double) ((short *)data)[i];
      break;
    case DT_UINT16:
      tmp = (double) ((unsigned short *)data)[i];
      break;
    case DT_INT32:
      tmp = (double) ((int *)data)[i];
      break;
    case DT_UINT32:
      tmp = (double) ((unsigned int *)data)[i];
      break;
    case DT_FLOAT32:
      tmp = (double) ((float *)data)[i];
      break;
    case DT_FLOAT64:
      tmp = (double) ((double *)data)[i];
      break;
    default:
      fprintf(stderr, "Data type %d not handled\n", datatype);
      break;
    }
    if (tmp < range[0]) {
      range[0] = tmp;
    }
    if (tmp > range[1]) {
      range[1] = tmp;
    }
  }
}

int niik_mat44_to_cosines_start_step(mat44 mat,double dx,double dy,double dz,double *ocosx,double *ocosy,double *ocosz,double *ostart,double *ostep) {
  char fcname[]="niik_mat44_to_cosines_start_step";
  double step[3],start[3];
  double cosines[3][3];
  int i,j;
  int verbose=0;
  if(verbose>=1) niik_fc_display(fcname,1);
  if(verbose>=1) fprintf(stdout,"[%s] pixel size   : %12.8f %12.8f %12.8f\n",fcname,dx,dy,dz);
  if(verbose>=1) mat44_display(mat);
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


int niik_image_write_minc2(char *outname,nifti_image *img) {
  int result;
  mihandle_t hvol;
  midimhandle_t dim[3];
  mivolumeprops_t props;
  const char *fcname="niik_image_write_minc2";
  int
  mi_type=MI_TYPE_UNKNOWN,
  i,
  verbose=niik_verbose();
  unsigned long start[3],count[3];
  double
  dstart[3],
         dstep[3],
         imin,imax,
         cosines1[3],
         cosines2[3],
         cosines3[3];

  if(verbose>=1) niik_fc_display(fcname,1);

  if(verbose>=1) fprintf(stdout,"[%s] minew_volume_props\n",fcname);
  result = minew_volume_props(&props);

  if(verbose>=1) fprintf(stdout,"[%s] miset_props_compression_type\n",fcname);
  result = miset_props_compression_type(props,
                                        MI_COMPRESS_ZLIB);

  if(verbose>=1) fprintf(stdout,"[%s] miset_props_zlib_compression\n",fcname);
  result = miset_props_zlib_compression(props, 3);

  if(verbose>=1) fprintf(stdout,"[%s] miset_props_multi_resolution\n",fcname);
  result = miset_props_multi_resolution(props, 1, 3);
  NIIK_RET0((result < 0),fcname,"miset_props_multi_resolution");

  if(verbose>=1) fprintf(stdout,"[%s] micreate_dimension (z)\n",fcname);
  NIIK_RET0(( micreate_dimension("zspace",MI_DIMCLASS_SPATIAL,
                                 MI_DIMATTR_REGULARLY_SAMPLED,
                                 img->nz,&dim[0]) == MI_ERROR),
            fcname,
            "micreate_dimension (z)");

  if(verbose>=1) fprintf(stdout,"[%s] micreate_dimension (y)\n",fcname);
  NIIK_RET0(( micreate_dimension("yspace",MI_DIMCLASS_SPATIAL,
                                 MI_DIMATTR_REGULARLY_SAMPLED,
                                 img->ny,&dim[1]) == MI_ERROR),
            fcname,
            "micreate_dimension (y)");

  if(verbose>=1) fprintf(stdout,"[%s] micreate_dimension (x)\n",fcname);
  NIIK_RET0(( micreate_dimension("xspace",MI_DIMCLASS_SPATIAL,
                                 MI_DIMATTR_REGULARLY_SAMPLED,
                                 img->nx,&dim[2]) == MI_ERROR),
            fcname,
            "micreate_dimension (x)");

  if(img->sform_code!=NIFTI_XFORM_UNKNOWN) {
    NIIK_RET0((!niik_mat44_to_cosines_start_step(img->sto_xyz,(double)img->dx,(double)img->dy,(double)img->dz,cosines1,cosines2,cosines3,dstart,dstep)),fcname,"niik_mat44_to_cosines_start_step");
  } else if(img->qform_code!=NIFTI_XFORM_UNKNOWN) {
    NIIK_RET0((!niik_mat44_to_cosines_start_step(img->qto_xyz,(double)1.0,(double)1.0,(double)1.0,cosines1,cosines2,cosines3,dstart,dstep)),fcname,"niik_mat44_to_cosines_start_step");
    dstep[0]=(double)img->dx;
    dstep[1]=(double)img->dy;
    dstep[2]=(double)img->dz;
  } else {
    fprintf(stdout,"[%s] Warning, step, start, and cosines are estimated!\n",fcname);
    dstep [0]=dstep [1]=dstep [2]=1;
    dstart[0]=dstart[1]=dstart[2]=0;
    cosines1[0]=1;
    cosines1[1]=0;
    cosines1[2]=0;
    cosines2[0]=0;
    cosines2[1]=1;
    cosines2[2]=0;
    cosines3[0]=0;
    cosines3[1]=0;
    cosines3[2]=1;
  }
  /* fprintf(stdout,"[%s] dstep = %f %f %f\n",fcname,dstep[0],dstep[1],dstep[2]); */

  NIIK_RET0((miset_dimension_start(dim[0],dstart[2])==MI_ERROR),fcname,"miset_dimension_start z");
  NIIK_RET0((miset_dimension_start(dim[1],dstart[1])==MI_ERROR),fcname,"miset_dimension_start y");
  NIIK_RET0((miset_dimension_start(dim[2],dstart[0])==MI_ERROR),fcname,"miset_dimension_start x");

  NIIK_RET0((miset_dimension_separation(dim[0],img->dz)==MI_ERROR),fcname,"miset_dimension_separation z");
  NIIK_RET0((miset_dimension_separation(dim[1],img->dy)==MI_ERROR),fcname,"miset_dimension_separation y");
  NIIK_RET0((miset_dimension_separation(dim[2],img->dx)==MI_ERROR),fcname,"miset_dimension_separation x");

  if(verbose>=1) {
    fprintf(stdout,"[%s] cosines from s-matrix\n",fcname);
    mat44_display(img->sto_xyz);
  }

  if(verbose>=1) fprintf(stdout,"[%s] miset_dimension_cosines z:  %.9f %.9f %.9f\n",fcname,cosines1[0],cosines1[1],cosines1[2]);
  NIIK_RET0((miset_dimension_cosines(dim[2],cosines1)==MI_ERROR),fcname,"miset_dimension_cosines z");
  if(verbose>=1) fprintf(stdout,"[%s] miset_dimension_cosines y:  %.9f %.9f %.9f\n",fcname,cosines2[0],cosines2[1],cosines2[2]);
  NIIK_RET0((miset_dimension_cosines(dim[1],cosines2)==MI_ERROR),fcname,"miset_dimension_cosines y");
  if(verbose>=1) fprintf(stdout,"[%s] miset_dimension_cosines x:  %.9f %.9f %.9f\n",fcname,cosines3[0],cosines3[1],cosines3[2]);
  NIIK_RET0((miset_dimension_cosines(dim[0],cosines3)==MI_ERROR),fcname,"miset_dimension_cosines x");

  switch(img->datatype) {
  case NIFTI_TYPE_UINT8:
    mi_type=MI_TYPE_UBYTE;
    break;
  case NIFTI_TYPE_UINT16:
    mi_type=MI_TYPE_USHORT;
    break;
  case NIFTI_TYPE_UINT32:
    mi_type=MI_TYPE_UINT;
    break;
  case NIFTI_TYPE_FLOAT32:
    mi_type=MI_TYPE_FLOAT;
    break;
  case NIFTI_TYPE_FLOAT64:
    mi_type=MI_TYPE_DOUBLE;
    break;
  default:
    fprintf(stdout,"[%s] ERROR: unknown datatype, %s %i\n",fcname,nifti_datatype_string(img->datatype),img->datatype);
    break;
  }

  if(verbose>=1) fprintf(stdout,"[%s] micreate_volume\n",fcname);
  NIIK_RET0(( micreate_volume(outname,
                              3,
                              dim,
                              mi_type,
                              MI_CLASS_REAL,
                              props,&hvol) == MI_ERROR),
            fcname,
            "micreate_volume");

  for(i=0; i<3; i++) {
    start[i]=0;
    count[i]=img->dim[3-i];
    if(verbose>=1) {
      fprintf(stdout,"[%s] start,count [%i]   %i %i\n",fcname,i,(int)start[i],(int)count[i]);
    }
  }

  if(verbose>=1) fprintf(stdout,"[%s] micreate_volume_image\n",fcname);
  NIIK_RET0( ( micreate_volume_image(hvol) == MI_ERROR),
             fcname,
             "micreate_volume_image");

  imin = niik_image_get_min(img,NULL);
  imax = niik_image_get_max(img,NULL);
  niik_miinit_default_range(mi_type,&imax,&imin);
  if((img->datatype==NIFTI_TYPE_FLOAT32) ||
      (img->datatype==NIFTI_TYPE_FLOAT64)) {
    imin=niik_image_get_min(img,NULL);
    imax=niik_image_get_max(img,NULL);
  }

  if(verbose>=1) fprintf(stdout,"[%s] miset_volume_min max  %12.7f %12.7f\n",fcname,imin,imax);
  NIIK_RET0((      miset_volume_range(hvol,imax,imin)==MI_ERROR),fcname,"miset_volume_range");

  if(verbose>=1) fprintf(stdout,"[%s] miset_real_value_hyperslab\n",fcname);
  NIIK_RET0( (miset_real_value_hyperslab(hvol,
                                         mi_type,
                                         (const misize_t *)start,(const misize_t *)count,
                                         (void *)img->data) == MI_ERROR),
             fcname,"miset_real_value_hyperslab");

  if(verbose>=1) niik_fc_display(fcname,0);
  /*free memory*/
  mifree_volume_props( props );

  if(img->minc_history != NULL) {
    if(miset_attr_values(hvol,MI_TYPE_STRING,"","history",strlen(img->minc_history)+1,img->minc_history ) == MI_ERROR)
      fprintf(stderr,"[%s] miset_attr_values failed!\n",fcname);
  }

  miclose_volume( hvol );

  return 1;
} /* minc2 */



int niik_image_write_minc(char *outname,nifti_image *nii_ptr) {
  const char *fcname="niik_image_write_minc";
  /* MINC stuff */
  int mnc_fd; /* MINC file descriptor */
  nc_type mnc_mtype; /* MINC memory data type */
  int mnc_msign; /* MINC !0 if signed data */
  static nc_type mnc_vtype; /* MINC voxel data type */
  static int mnc_vsign; /* MINC !0 if signed data */
  int mnc_ndims; /* MINC image dimension count */
  int mnc_dimids[MAX_VAR_DIMS]; /* MINC image dimension identifiers */
  int mnc_iid; /* MINC Image variable ID */
  long mnc_start[MAX_VAR_DIMS]; /* MINC data starts */
  long mnc_count[MAX_VAR_DIMS]; /* MINC data counts */
  char *mnc_hist=NULL; /* MINC history */
  double mnc_vrange[2]; /* MINC valid min/max */
  double mnc_srange[2]; /* MINC image min/max */
  double mnc_time_step;
  double mnc_time_start;
  int mnc_spatial_axes[MAX_NII_DIMS];
  double mnc_starts[MAX_SPACE_DIMS];
  double mnc_steps[MAX_SPACE_DIMS];
  double mnc_dircos[MAX_SPACE_DIMS][MAX_SPACE_DIMS];
  VIO_Transform mnc_xform;
  VIO_General_transform mnc_linear_xform;
  int mnc_vid; /* Dimension variable id */
  /*struct analyze75_hdr ana_hdr;*/

  /* Other stuff */
  int i; /* Generic loop counter the first */
  int j; /* Generic loop counter the second */
  /* char *str_ptr; */ /* Generic ASCIZ string pointer */
  int r; /* Result code. */
  static int qflag = 0; /* Quiet flag (default is quiet) */
  static int rflag = 1; /* Scan range flag */
  static int oflag = DIMORDER_ZYX;
  static int flip[MAX_SPACE_DIMS];

  static char *mnc_ordered_dim_names[MAX_SPACE_DIMS];
  /*    FILE *fp;
    int ss;
    int must_swap;*/

  int verbose = 2;

  return niik_image_write_minc2(outname,nii_ptr);

  if(verbose>=1) niik_fc_display(fcname,1);

  mnc_hist = time_stamp(1, &outname);

  ncopts = 0; /* Clear global netCDF error reporting flag */
  mnc_vtype = NC_NAT;

  /* Set up the spatial axis correspondence for the call to
   * convert_transform_to_starts_and_steps()
   */
  switch (oflag) {
  default:
  case DIMORDER_ZYX:
    mnc_ordered_dim_names[DIM_X] = MIxspace;
    mnc_ordered_dim_names[DIM_Y] = MIyspace;
    mnc_ordered_dim_names[DIM_Z] = MIzspace;
    mnc_spatial_axes[DIM_X] = DIM_X;
    mnc_spatial_axes[DIM_Y] = DIM_Y;
    mnc_spatial_axes[DIM_Z] = DIM_Z;
    break;
  case DIMORDER_ZXY:
    mnc_ordered_dim_names[DIM_X] = MIyspace;
    mnc_ordered_dim_names[DIM_Y] = MIxspace;
    mnc_ordered_dim_names[DIM_Z] = MIzspace;
    mnc_spatial_axes[DIM_X] = DIM_Y;
    mnc_spatial_axes[DIM_Y] = DIM_X;
    mnc_spatial_axes[DIM_Z] = DIM_Z;
    break;
  case DIMORDER_XYZ:
    mnc_ordered_dim_names[DIM_X] = MIzspace;
    mnc_ordered_dim_names[DIM_Y] = MIyspace;
    mnc_ordered_dim_names[DIM_Z] = MIxspace;
    mnc_spatial_axes[DIM_X] = DIM_Z;
    mnc_spatial_axes[DIM_Y] = DIM_Y;
    mnc_spatial_axes[DIM_Z] = DIM_X;
    break;
  case DIMORDER_XZY:
    mnc_ordered_dim_names[DIM_X] = MIyspace;
    mnc_ordered_dim_names[DIM_Y] = MIzspace;
    mnc_ordered_dim_names[DIM_Z] = MIxspace;
    mnc_spatial_axes[DIM_X] = DIM_Y;
    mnc_spatial_axes[DIM_Y] = DIM_Z;
    mnc_spatial_axes[DIM_Z] = DIM_X;
    break;
  case DIMORDER_YZX:
    mnc_ordered_dim_names[DIM_X] = MIxspace;
    mnc_ordered_dim_names[DIM_Y] = MIzspace;
    mnc_ordered_dim_names[DIM_Z] = MIyspace;
    mnc_spatial_axes[DIM_X] = DIM_X;
    mnc_spatial_axes[DIM_Y] = DIM_Z;
    mnc_spatial_axes[DIM_Z] = DIM_Y;
    break;
  case DIMORDER_YXZ:
    mnc_ordered_dim_names[DIM_X] = MIzspace;
    mnc_ordered_dim_names[DIM_Y] = MIxspace;
    mnc_ordered_dim_names[DIM_Z] = MIyspace;
    mnc_spatial_axes[DIM_X] = DIM_Z;
    mnc_spatial_axes[DIM_Y] = DIM_X;
    mnc_spatial_axes[DIM_Z] = DIM_Y;
    break;
  }

  switch (nii_ptr->datatype) {
  case DT_INT8:
    mnc_msign = 1;
    mnc_mtype = NC_BYTE;
    mnc_vrange[0] = CHAR_MIN;
    mnc_vrange[1] = CHAR_MAX;
    break;
  case DT_UINT8:
    mnc_msign = 0;
    mnc_mtype = NC_BYTE;
    mnc_vrange[0] = 0;
    mnc_vrange[1] = UCHAR_MAX;
    break;
  case DT_INT16:
    mnc_msign = 1;
    mnc_mtype = NC_SHORT;
    mnc_vrange[0] = SHRT_MIN;
    mnc_vrange[1] = SHRT_MAX;
    break;
  case DT_UINT16:
    mnc_msign = 0;
    mnc_mtype = NC_SHORT;
    mnc_vrange[0] = 0;
    mnc_vrange[1] = USHRT_MAX;
    break;
  case DT_INT32:
    mnc_msign = 1;
    mnc_mtype = NC_INT;
    mnc_vrange[0] = INT_MIN;
    mnc_vrange[1] = INT_MAX;
    break;
  case DT_UINT32:
    mnc_msign = 0;
    mnc_mtype = NC_INT;
    mnc_vrange[0] = 0;
    mnc_vrange[1] = UINT_MAX;
    break;
  case DT_FLOAT32:
    mnc_msign = 1;
    mnc_mtype = NC_FLOAT;
    mnc_vrange[0] = -FLT_MAX;
    mnc_vrange[1] = FLT_MAX;
    break;
  case DT_FLOAT64:
    mnc_msign = 1;
    mnc_mtype = NC_DOUBLE;
    mnc_vrange[0] = -DBL_MAX;
    mnc_vrange[1] = DBL_MAX;
    break;
  default:
    fprintf(stderr, "[%s] Data type %d not handled\n", fcname,nii_ptr->datatype);
    break;
  }

  if (mnc_vtype == NC_NAT) {
    mnc_vsign = mnc_msign;
    mnc_vtype = mnc_mtype;
  }

  /* Open the MINC file. It should not already exist.
   */
  mnc_fd = micreate(outname, NC_NOCLOBBER);
  if (mnc_fd < 0) {
    fprintf(stderr, "[%s] ERROR: Can't create output file '%s'\n", fcname,outname);
    return 0;
  }

  if(verbose>=2) fprintf(stdout,"[%s] create dim in minc \n",fcname);
  /* Create the necessary dimensions in the minc file
   */

  mnc_ndims = 0;

  if (nii_ptr->nt > 1) {
    mnc_dimids[mnc_ndims] = ncdimdef(mnc_fd, MItime, nii_ptr->nt);
    mnc_start[mnc_ndims] = 0;
    mnc_count[mnc_ndims] = nii_ptr->nt;
    mnc_ndims++;

    r = micreate_std_variable(mnc_fd, MItime, NC_INT, 0, NULL);
    switch (nii_ptr->time_units) {
    case NIFTI_UNITS_UNKNOWN:
    case NIFTI_UNITS_SEC:
      mnc_time_step = nii_ptr->dt;
      mnc_time_start = nii_ptr->toffset;
      break;
    case NIFTI_UNITS_MSEC:
      mnc_time_step = nii_ptr->dt / 1000;
      mnc_time_start = nii_ptr->toffset / 1000;
      break;
    case NIFTI_UNITS_USEC:
      mnc_time_step = nii_ptr->dt / 1000000;
      mnc_time_start = nii_ptr->toffset / 1000000;
      break;
    default:
      fprintf(stderr, "Unknown time units value %d\n",
              nii_ptr->time_units);
      break;
    }
    miattputdbl(mnc_fd, r, MIstart, mnc_time_start);
    miattputdbl(mnc_fd, r, MIstep, mnc_time_step);
    miattputstr(mnc_fd, r, MIunits, "s");
  }

  if (nii_ptr->nz > 1) {
    mnc_dimids[mnc_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[DIM_Z],
                                     nii_ptr->nz);
    mnc_start[mnc_ndims] = 0;
    mnc_count[mnc_ndims] = nii_ptr->nz;
    mnc_ndims++;
    r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[DIM_Z], NC_INT,
                              0, NULL);
    miattputdbl(mnc_fd, r, MIstep, nii_ptr->dz);
    miattputstr(mnc_fd, r, MIunits, "mm");
  }

  if (nii_ptr->ny > 1) {
    mnc_dimids[mnc_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[DIM_Y],
                                     nii_ptr->ny);
    mnc_start[mnc_ndims] = 0;
    mnc_count[mnc_ndims] = nii_ptr->ny;
    mnc_ndims++;
    r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[DIM_Y], NC_INT,
                              0, NULL);
    miattputdbl(mnc_fd, r, MIstep, nii_ptr->dy);
    miattputstr(mnc_fd, r, MIunits, "mm");
  }

  if (nii_ptr->nx > 1) {
    mnc_dimids[mnc_ndims] = ncdimdef(mnc_fd, mnc_ordered_dim_names[DIM_X],
                                     nii_ptr->nx);
    mnc_start[mnc_ndims] = 0;
    mnc_count[mnc_ndims] = nii_ptr->nx;
    mnc_ndims++;
    r = micreate_std_variable(mnc_fd, mnc_ordered_dim_names[DIM_X], NC_INT,
                              0, NULL);
    miattputdbl(mnc_fd, r, MIstep, nii_ptr->dx);
    miattputstr(mnc_fd, r, MIunits, "mm");
  }

  if (nii_ptr->nu > 1) {
    mnc_dimids[mnc_ndims] = ncdimdef(mnc_fd, MIvector_dimension,
                                     nii_ptr->nu);
    mnc_start[mnc_ndims] = 0;
    mnc_count[mnc_ndims] = nii_ptr->nu;
    mnc_ndims++;
  }

  /* Create scalar image-min and image-max variables.
   */
  if(verbose>=2) fprintf(stdout,"[%s] create image-min and image-max\n",fcname);
  micreate_std_variable(mnc_fd, MIimagemax, NC_DOUBLE, 0, NULL);
  micreate_std_variable(mnc_fd, MIimagemin, NC_DOUBLE, 0, NULL);

  /* Create the group variables.
   */
  if(verbose>=2) fprintf(stdout,"[%s] create group var\n",fcname);
  micreate_std_variable(mnc_fd, MIstudy, NC_INT, 0, NULL);
  if (strlen(nii_ptr->descrip) > 0 && strlen(nii_ptr->descrip) < 79 ) {
    int varid = micreate_std_variable(mnc_fd, MIpatient, NC_INT, 0, NULL);
    (void) miattputstr(mnc_fd, varid, MIfull_name,
                       nii_ptr->descrip);
  } else {
    micreate_std_variable(mnc_fd, MIpatient, NC_INT, 0, NULL);
  }
  micreate_std_variable(mnc_fd, MIacquisition, NC_INT, 0, NULL);

  /* Create the MINC image variable. If we can't, there is no
   * further processing possible...
   */
  if(verbose>=2) fprintf(stdout,"[%s] create minc image variable\n",fcname);
  mnc_iid = micreate_std_variable(mnc_fd, MIimage,
                                  mnc_vtype,
                                  mnc_ndims,
                                  mnc_dimids);
  NIIK_RET0((mnc_iid < 0),fcname,"Can't create the image variable");

  miattputstr(mnc_fd, mnc_iid, MIsigntype,
              (mnc_vsign) ? MI_SIGNED : MI_UNSIGNED);

  /* Calculate the starts, steps, and direction cosines. This only
   * be done properly if the file is NIfTI-1 file. If it is an Analyze
   * file we have to resort to other methods...
   */

  if(verbose>=2) fprintf(stdout,"[%s] calculate starts, steps, and direcion cosines\n",fcname);

  if (nii_ptr->nifti_type != 0 &&
      (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN ||
       nii_ptr->qform_code != NIFTI_XFORM_UNKNOWN)) {
    make_identity_transform(&mnc_xform);
    if (nii_ptr->sform_code != NIFTI_XFORM_UNKNOWN) {
      if (!qflag) {
        printf("Using s-form transform:\n");
      }
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          Transform_elem(mnc_xform, i, j) = nii_ptr->sto_xyz.m[i][j];
          if (!qflag) {
            printf("%8.4f, ", nii_ptr->sto_xyz.m[i][j]);
          }
        }
        if (!qflag) {
          printf("\n");
        }
      }
    } else {
      if (!qflag) {
        printf("Using q-form transform:\n");
      }
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          Transform_elem(mnc_xform, i, j) = nii_ptr->qto_xyz.m[i][j];
          if (!qflag) {
            printf("%8.4f, ", nii_ptr->qto_xyz.m[i][j]);
          }
        }
        if (!qflag) {
          printf("\n");
        }
      }
    }
    create_linear_transform(&mnc_linear_xform, &mnc_xform);
    convert_transform_to_starts_and_steps(&mnc_linear_xform,
                                          MAX_SPACE_DIMS,
                                          NULL,
                                          mnc_spatial_axes,
                                          mnc_starts,
                                          mnc_steps,
                                          mnc_dircos);
  } else {
    /* No official transform was found (possibly this is an Analyze
     * file). Just use some reasonable defaults.
     */
    mnc_steps[mnc_spatial_axes[DIM_X]] = nii_ptr->dx;
    mnc_steps[mnc_spatial_axes[DIM_Y]] = nii_ptr->dy;
    mnc_steps[mnc_spatial_axes[DIM_Z]] = nii_ptr->dz;
    mnc_starts[mnc_spatial_axes[DIM_X]] = -(nii_ptr->dx * nii_ptr->nx) / 2;
    mnc_starts[mnc_spatial_axes[DIM_Y]] = -(nii_ptr->dy * nii_ptr->ny) / 2;
    mnc_starts[mnc_spatial_axes[DIM_Z]] = -(nii_ptr->dz * nii_ptr->nz) / 2;
    /* Unlike the starts and steps, the direction cosines do NOT change
     * based upon the data orientation.
     */
    for (i = 0; i < MAX_SPACE_DIMS; i++) {
      for (j = 0; j < MAX_SPACE_DIMS; j++) {
        mnc_dircos[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
  }
  switch (nii_ptr->xyz_units) {
  case NIFTI_UNITS_METER:
    for (i = 0; i < MAX_SPACE_DIMS; i++) {
      mnc_starts[i] *= 1000;
      mnc_steps[i] *= 1000;
    }
    break;
  case NIFTI_UNITS_MM:
    break;
  case NIFTI_UNITS_MICRON:
    for (i = 0; i < MAX_SPACE_DIMS; i++) {
      mnc_starts[i] /= 1000;
      mnc_steps[i] /= 1000;
    }
    break;
  default:
    fprintf(stderr, "Unknown XYZ units %d\n", nii_ptr->xyz_units);
    break;
  }
  /* Now we write the spatial axis information to the file. The starts,
   * steps, and cosines have to be associated with the correct spatial
   * axes. Also, we perform any orientation flipping that was requested.
   */
  for (i = 0; i < MAX_SPACE_DIMS; i++) {
    if (!qflag) {
      printf("%s start: %8.4f step: %8.4f cosines: %8.4f %8.4f %8.4f\n",
             mnc_spatial_names[i],
             mnc_starts[i],
             mnc_steps[i],
             mnc_dircos[i][DIM_X],
             mnc_dircos[i][DIM_Y],
             mnc_dircos[i][DIM_Z]);
    }
    mnc_vid = ncvarid(mnc_fd, mnc_spatial_names[i]);
    /* If we selected "flipping" of the appropriate axis, do it here
     */
    if(verbose>=2) fprintf(stdout,"[%s] flipping\n",fcname);
    if (flip[i]) {
      miattputdbl(mnc_fd, mnc_vid, MIstart,
                  mnc_starts[i]+((mnc_count[i]-1)*mnc_steps[i]));
      miattputdbl(mnc_fd, mnc_vid, MIstep, -mnc_steps[i]);
    } else {
      miattputdbl(mnc_fd, mnc_vid, MIstart, mnc_starts[i]);
      miattputdbl(mnc_fd, mnc_vid, MIstep, mnc_steps[i]);
    }
    ncattput(mnc_fd, mnc_vid, MIdirection_cosines, NC_DOUBLE,
             MAX_SPACE_DIMS, mnc_dircos[i]);
  }

  /* Find the valid minimum and maximum of the data, in order to set the
   * global image minimum and image maximum properly.
   */
  if(verbose>=2) fprintf(stdout,"[%s] find valid min and max of data\n",fcname);
  if (rflag) {
    find_data_range(nii_ptr->datatype,
                    nii_ptr->nvox,
                    nii_ptr->data,
                    mnc_vrange);
    if(verbose>=2) fprintf(stdout,"[%s]   vrange %12.6f %12.6f\n",fcname,mnc_vrange[0],mnc_vrange[1]);
  }
  if (nii_ptr->scl_slope != 0.0) {
    /* Convert slope/offset to min/max
     */
    mnc_srange[0] = (mnc_vrange[0] * nii_ptr->scl_slope) + nii_ptr->scl_inter;
    mnc_srange[1] = (mnc_vrange[1] * nii_ptr->scl_slope) + nii_ptr->scl_inter;
  } else {
    mnc_srange[0] = mnc_vrange[0];
    mnc_srange[1] = mnc_vrange[1];
  }
  if(verbose>=2) fprintf(stdout,"[%s]   mnc_srange %12.6f %12.6f\n",fcname,mnc_srange[0],mnc_srange[1]);

  ncattput(mnc_fd, mnc_iid, MIvalid_range, NC_DOUBLE, 2, mnc_vrange);

  if(verbose>=2) fprintf(stdout,"[%s]   miattputstr\n",fcname);
  miattputstr(mnc_fd, NC_GLOBAL, MIhistory, mnc_hist);

  /* Switch out of definition mode.
   */
  if(verbose>=2) fprintf(stdout,"[%s] switch out of def\n",fcname);
  ncendef(mnc_fd);

  /* Finally, write the values of the image-min, image-max, and image
   * variables.
   */
  if(verbose>=2) fprintf(stdout,"[%s] write values of image-min, image-max, and image var\n",fcname);
  mivarput1(mnc_fd, ncvarid(mnc_fd, MIimagemin), mnc_start, NC_DOUBLE,
            MI_SIGNED, &mnc_srange[0]);
  mivarput1(mnc_fd, ncvarid(mnc_fd, MIimagemax), mnc_start, NC_DOUBLE,
            MI_SIGNED, &mnc_srange[1]);
  mivarput(mnc_fd, mnc_iid, mnc_start, mnc_count, mnc_mtype,
           (mnc_msign) ? MI_SIGNED : MI_UNSIGNED, nii_ptr->data);
  miclose(mnc_fd);

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}


#endif

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/