/* Filename:     nifti1_kunio_minc.c
 * Description:  minc function(s)
 * Author:       Kunio Nakamura
 * Date:         March 7, 2012
 *
 * Reference = minctoraw.c
 */

#include "falcon.h"
#include <minc.h>
#include <limits.h>
#include <float.h>

#ifndef BOOLEAN_DEFAULT
#define BOOLEAN_DEFAULT -1
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif


/* Variables used for argument parsing */
static int output_datatype = INT_MAX;
static int output_signed = INT_MAX;
static double valid_range[2] = {DBL_MAX, DBL_MAX};
static int normalize_output = BOOLEAN_DEFAULT;

nifti_image *niik_read_minc(char *fname) {
  int mincid, imgid, icvid, ndims, dims[MAX_VAR_DIMS];
  nc_type datatype;
  int is_signed;
  long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS], end[MAX_VAR_DIMS];
  long size;
  int idim;
  void *data;
  double temp;

  nifti_image *outimg=NULL;


  /* Open the file */
  mincid = miopen(fname, NC_NOWRITE);

  /* Inquire about the image variable */
  imgid = ncvarid(mincid, MIimage);
  (void) ncvarinq(mincid, imgid, NULL, NULL, &ndims, dims, NULL);
  (void)miget_datatype(mincid, imgid, &datatype, &is_signed);

  /* Get output data type */
  if (output_datatype == INT_MAX) output_datatype = datatype;

  /* Get output sign */
  if (output_signed == INT_MAX) {
    if (output_datatype == datatype)
      output_signed = is_signed;
    else
      output_signed = (output_datatype != NC_BYTE);
  }

  /* Get output range */
  if (valid_range[0] == DBL_MAX) {
    if ((output_datatype == datatype) && (output_signed == is_signed)) {
      (void) miget_valid_range(mincid, imgid, valid_range);
    } else {
      (void) miget_default_range(output_datatype, output_signed,
                                 valid_range);
    }
  }
  if (valid_range[0] > valid_range[1]) {
    temp = valid_range[0];
    valid_range[0] = valid_range[1];
    valid_range[1] = temp;
  }

  /* Set up image conversion */
  icvid = miicv_create();
  (void) miicv_setint(icvid, MI_ICV_TYPE, output_datatype);
  (void) miicv_setstr(icvid, MI_ICV_SIGN, (output_signed ?
                      MI_SIGNED : MI_UNSIGNED));
  (void) miicv_setdbl(icvid, MI_ICV_VALID_MIN, valid_range[0]);
  (void) miicv_setdbl(icvid, MI_ICV_VALID_MAX, valid_range[1]);
  if ((output_datatype == NC_FLOAT) || (output_datatype == NC_DOUBLE)) {
    (void) miicv_setint(icvid, MI_ICV_DO_NORM, TRUE);
    (void) miicv_setint(icvid, MI_ICV_USER_NORM, TRUE);
  } else if (normalize_output) {
    (void) miicv_setint(icvid, MI_ICV_DO_NORM, TRUE);
  }
  (void) miicv_attach(icvid, mincid, imgid);

  /* Set input file start, count and end vectors for reading a slice
     at a time */
  for (idim=0; idim < ndims; idim++) {
    (void) ncdiminq(mincid, dims[idim], NULL, &end[idim]);
  }
  (void) miset_coords(ndims, (long) 0, start);
  (void) miset_coords(ndims, (long) 1, count);
  size = nctypelen(output_datatype);
  fprintf(stdout,"  output_datatype_byte_len = %i\n",(int)size);

  fprintf(stdout,"  ndims = %i\n",ndims);
  for (idim=ndims-2; idim < ndims; idim++) {
    count[idim] = end[idim];
    size *= count[idim];
    fprintf(stdout,"  count[idim=%i] = %li\n",idim,count[idim]);
  }

  /* Allocate space */
  data = malloc(size*(end[0]-start[0]));

  /* Loop over input slices */

  while (start[0] < end[0]) {
    /* Read in the slice */
    (void) miicv_get(icvid, start, count, data);
    /* Write out the slice */
    if (fwrite(data, sizeof(char), (size_t) size, stdout) != size) {
      (void) fprintf(stderr, "Error writing data.\n");
      exit(EXIT_FAILURE);
    }
    /* Increment start counter */
    idim = ndims-1;
    start[idim] += count[idim];
    while ( (idim>0) && (start[idim] >= end[idim])) {
      start[idim] = 0;
      idim--;
      start[idim] += count[idim];
    }
  }       /* End loop over slices */

  /* Clean up */
  (void) miclose(mincid);
  (void) miicv_free(icvid);
  free(data);


  return NULL;
}


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
