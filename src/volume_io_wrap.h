/*
# Developed by Simon Fristed Eskildsen, eskild@gmail.com
# Part of FALCON
#
# Copyright notice:
# This code is copyright Simon Fristed Eskildsen.
# It may not be copied, altered in any way or transmitted
# to others (unless explicitly stated otherwise) without
# the written permission of the author/developer. 
*/

#ifndef VOLUME_IO_WRAP
#define VOLUME_IO_WRAP

#include <volume_io.h>
#include "surface_tools_basic.h"

typedef struct {
  nc_type type;
  byte type_size;
  BOOLEAN sign;
  int sizes[3];
  void *data;
}Volume_wrap;

/* #ifndef NC_NAT */
/* #define NC_NAT -1 */
/* #endif */
/* #ifndef NC_INT */
/* #define NC_INT NC_LONG */
/* #endif */

#define NC_BYTE_SIZE sizeof(byte)
#define NC_SHORT_SIZE sizeof(short)
#define NC_INT_SIZE sizeof(int)
#define NC_FLOAT_SIZE sizeof(float)
#define NC_DOUBLE_SIZE sizeof(double)
#define NC_NAT_SIZE sizeof(byte)
#define NC_CHAR_SIZE sizeof(char)

int volumeToWrap(VIO_Volume *volume, Volume_wrap *wrap);
int volumeToWrap_real(VIO_Volume *volume, Volume_wrap *wrap);
void free_wrap(Volume_wrap *wrap);
void ***alloc_data3D(int sizes[3],byte size_element);
int wrapDataToVolume(VIO_Volume *volume, Volume_wrap *wrap);
void copy_wrap(Volume_wrap *original, Volume_wrap *copy);
int wrapVoxelDataToVolume(VIO_Volume *volume, Volume_wrap *wrap);
int wrapVoxelDataToVolume_real(VIO_Volume *volume, Volume_wrap *wrap);
int short2byte(Volume_wrap *vol);
VIO_Real get_wrap_min(Volume_wrap *wrp);
VIO_Real get_wrap_max(Volume_wrap *wrp);
void make_binary_wrap(Volume_wrap *wrap);
void output_volume_wrap(Volume_wrap *wrap, VIO_STR filename);
#endif

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
*/