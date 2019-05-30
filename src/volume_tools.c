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

#include <float.h>
#include "surface_tools_basic.h"
#include "volume_tools.h"

//#define DEBUG 1

inline float interpolate_volume(
			 VIO_Volume *volume, 
			 point3D *currentPoint
			 )
{
  int x,y,z;
  float g, x_pol=0, y_pol=0, z_pol=0,decimals;

  x = (int)currentPoint->x;
  y = (int)currentPoint->y;
  z = (int)currentPoint->z;
  
  g = get_volume_real_value(*volume, x, y, z, 0, 0);
  
  /* interpolate the volume */
  if ((decimals=currentPoint->x-x) > 0)
    x_pol = decimals*(get_volume_real_value(*volume, x+1, y, z, 0, 0)-g);
  if ((decimals=currentPoint->y-y) > 0)
    y_pol = decimals*(get_volume_real_value(*volume, x, y+1, z, 0, 0)-g);
  if ((decimals=currentPoint->z-z) > 0)
    z_pol = decimals*(get_volume_real_value(*volume, x, y, z+1, 0, 0)-g);
  return g + (x_pol+y_pol+z_pol)/3.0;
  
}// interpolate_volume()

inline float interpolate_volume_nn(
			 VIO_Volume *volume, 
			 point3D *currentPoint
			 )
{
  int x,y,z;
  float g;

  x = (int)currentPoint->x+0.5;
  y = (int)currentPoint->y+0.5;
  z = (int)currentPoint->z+0.5;
  
  g = get_volume_real_value(*volume, x, y, z, 0, 0);
  
 
  return g;
  
}// interpolate_volume()

float interpolate_volume_cubic(VIO_Volume *volume,point3D *currentPoint, BOOLEAN world){
  VIO_Real voxel[3];
  VIO_Real r_value;
  
  voxel[0]=currentPoint->x;
  voxel[1]=currentPoint->y;
  voxel[2]=currentPoint->z;
 
  if (world)
    evaluate_volume_in_world(*volume,currentPoint->x,currentPoint->y,currentPoint->z,2,TRUE,0.0,&r_value,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  else
    evaluate_volume(*volume,voxel,NULL,2,TRUE,0.0,&r_value,NULL,NULL);

  return((float)r_value);
}

/* float trilinear_interpolant(VIO_Volume volume, point3D coord) */
/* { */
/*    int slcind, rowind, colind, slcmax, rowmax, colmax; */
/*    int slcnext, rownext, colnext,sizes[5]; */
/*    float f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2; */
/*    float v000, v001, v010, v011, v100, v101, v110, v111,result; */

/*    get_volume_sizes(volume,sizes); */

/*    /\* Check that the coordinate is inside the volume *\/ */
/*    slcmax = sizes[0] - 1; */
/*    rowmax = sizes[1] - 1; */
/*    colmax = sizes[2] - 1; */
/*    if ((coord.z  < 0) ||  */
/*        (coord.z  > slcmax) || */
/*        (coord.y  < 0) ||  */
/*        (coord.y  > rowmax) || */
/*        (coord.x  < 0) ||  */
/*        (coord.x  > colmax)) { */
/*      return 0; */
/*    } */

/*    /\* Get the whole part of the coordinate *\/  */
/*    slcind = (int) coord.z; */
/*    rowind = (int) coord.y; */
/*    colind = (int) coord.x; */
/*    if (slcind >= slcmax-1) slcind = slcmax-1; */
/*    if (rowind >= rowmax-1) rowind = rowmax-1; */
/*    if (colind >= colmax-1) colind = colmax-1; */

/*    /\* Get the next voxel up *\/ */
/*    slcnext = slcind+1; */
/*    rownext = rowind+1; */
/*    colnext = colind+1; */

/*    /\* Check for case of dimension of length one *\/ */
/*    if (slcmax == 0) { */
/*       slcind = 0; */
/*       slcnext = 0; */
/*    } */
/*    if (rowmax == 0) { */
/*       rowind = 0; */
/*       rownext = 0; */
/*    } */
/*    if (colmax == 0) { */
/*       colind = 0; */
/*       colnext = 0; */
/*    } */

/*    /\* Get the relevant voxels *\/ */
/*    v000 = volume[slcind*sizes[1]*sizes[2] + rowind*sizes[2] + colind]; */
/*    v001 = volume[slcind*sizes[1]*sizes[2] + rowind*sizes[2] + colnext]; */
/*    v010 = volume[slcind*sizes[1]*sizes[2] + rownext*sizes[2] + colind]; */
/*    v011 = volume[slcind*sizes[1]*sizes[2] + rownext*sizes[2] + colnext]; */

/*    v100 = volume[slcnext*sizes[1]*sizes[2] + rowind*sizes[2] + colind]; */
/*    v101 = volume[slcnext*sizes[1]*sizes[2] + rowind*sizes[2] + colnext]; */
/*    v110 = volume[slcnext*sizes[1]*sizes[2] + rownext*sizes[2] + colind]; */
/*    v111 = volume[slcnext*sizes[1]*sizes[2] + rownext*sizes[2] + colnext]; */

/*    /\* Get the fraction parts *\/ */
/*    f0 = coord.z  - slcind; */
/*    f1 = coord.y  - rowind; */
/*    f2 = coord.x - colind; */
/*    r0 = 1.0 - f0; */
/*    r1 = 1.0 - f1; */
/*    r2 = 1.0 - f2; */

/*    /\* Do the interpolation *\/ */
/*    r1r2 = r1 * r2; */
/*    r1f2 = r1 * f2; */
/*    f1r2 = f1 * r2; */
/*    f1f2 = f1 * f2; */
/*    result = */
/*       r0 *  (r1r2 * v000 + */
/*              r1f2 * v001 + */
/*              f1r2 * v010 + */
/*              f1f2 * v011); */
/*    result += */
/*       f0 *  (r1r2 * v100 + */
/*              r1f2 * v101 + */
/*              f1r2 * v110 + */
/*              f1f2 * v111); */
   
/*    return TRUE; */

/* } */


int normalizeVolume_max_min(VIO_Volume *vol, VIO_Volume *normalized){
  int i,j,k,sizes[VIO_MAX_DIMENSIONS],count=0;
  VIO_Real voxelVal=0.0,max=0.0,sum=0.0,mean=0.0,min=FLT_MAX,minMax=0.0;

  *normalized = copy_volume_definition(*vol,NC_FLOAT,FALSE,0.0,1.0);
  
  get_volume_sizes(*vol,sizes);
  for(i=0;i<sizes[0];i++)
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++){
	voxelVal = get_volume_real_value(*vol,i,j,k,0,0);
	if (voxelVal > max)
	  max=voxelVal;
	if (voxelVal < min)
	  min=voxelVal;
	if (voxelVal != 0.0){
	  sum+=voxelVal;
	  count++;
	}
      }


  mean = sum / count;
  minMax=max-min;
/*   fprintf(stderr,"Max: %f\n",max); */
/*   fprintf(stderr,"Min: %f\n",min); */
  for(i=0;i<sizes[0];i++)
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++){
	voxelVal = get_volume_real_value(*vol,i,j,k,0,0);
	set_volume_real_value(*normalized,i,j,k,0,0,(voxelVal - min)/minMax);
      }
  
  set_volume_real_range(*normalized,0.0,1.0);

  return 0;
}// normalizeVolume()

VIO_Volume normalize_volume(VIO_Volume vol, VIO_Real low, VIO_Real high){
  int i,j,k,sizes[VIO_MAX_DIMENSIONS];
  VIO_Real voxelVal=0.0,max,min,range,newrange;
  VIO_Volume normalized;

  normalized = copy_volume_definition(vol,NC_FLOAT,FALSE,0.0,1.0);
  
  get_volume_sizes(vol,sizes);

  get_volume_real_range(vol,&min,&max);

  range=max-min;
  newrange=high-low;

  for(i=0;i<sizes[0];i++)
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++){
	voxelVal = get_volume_real_value(vol,i,j,k,0,0);
	set_volume_real_value(normalized,i,j,k,0,0,((voxelVal - min)/range)*newrange);
      }
  
  set_volume_real_range(normalized,low,high);

  return normalized;
}// normalize_volume()

int gradient_scale_volume(VIO_Volume *gradient, VIO_Volume *wm_membership, VIO_Volume *gm_membership, VIO_Volume *csf_membership, VIO_Volume *result){
  int sizes[3];
  int i,j,k;
  VIO_Real gradient_value,wm_value,gm_value,csf_value,scale_value,mbr_diff,mbr_scale,gradient_scale;

  *result = copy_volume(*gradient);
  get_volume_sizes(*gradient,sizes);

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++){
      for(k=0;k<sizes[2];k++) {
	gradient_value = get_volume_real_value(*gradient, i, j, k,0,0);
	wm_value = get_volume_real_value(*wm_membership, i, j, k,0,0);
	gm_value = get_volume_real_value(*gm_membership, i, j, k,0,0);
	csf_value = get_volume_real_value(*csf_membership, i, j, k,0,0);
	mbr_diff = gm_value+wm_value-csf_value;
	mbr_scale=mbr_diff*(2.0 - 2.0*cos(mbr_diff));
	gradient_scale = cos(gradient_value*PI/2.0);
	scale_value = (1-ABS(mbr_scale))*gradient_scale + ABS(mbr_scale);
	set_volume_real_value( *result, i, j, k, 0, 0, scale_value);
      }
    }
  }
 
  return STATUS_OK;
}

int add_volumes(VIO_Volume *vol1, VIO_Volume *vol2, VIO_Volume *result){
  int sizes[3];
  int i,j,k;
  VIO_Real value;

  *result = copy_volume(*vol1);
  get_volume_sizes(*vol1,sizes);

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++){
      for(k=0;k<sizes[2];k++) {
	value = get_volume_real_value( *vol1, i, j, k,0,0);
	value += get_volume_real_value( *vol2, i, j, k,0,0);
	set_volume_real_value( *result, i, j, k, 0, 0, value);
      }
    }
  }
 
  return STATUS_OK;
}

VIO_Volume or_volumes(VIO_Volume vol1, VIO_Volume vol2){
  int sizes[3];
  int i,j,k;
  VIO_Real value;
  VIO_Volume result;
  
  result = copy_volume(vol1);
  get_volume_sizes(vol1,sizes);

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++){
      for(k=0;k<sizes[2];k++) {
	if ((get_volume_real_value(vol1, i, j, k,0,0)) || (get_volume_real_value(vol2, i, j, k,0,0)))
	  set_volume_real_value( result, i, j, k, 0, 0, 255.0);
	else
	  set_volume_real_value( result, i, j, k, 0, 0, 0.0);
      }
    }
  }
 
  return result;
}

VIO_Volume not_volume(VIO_Volume vol){
  int sizes[3];
  int i,j,k;
  VIO_Real value;
  VIO_Volume result;
  
  result = copy_volume(vol);
  get_volume_sizes(vol,sizes);

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++){
      for(k=0;k<sizes[2];k++) {
	if ((get_volume_real_value(vol, i, j, k,0,0)) )
	  set_volume_real_value( result, i, j, k, 0, 0, 0.0);
	else
	  set_volume_real_value( result, i, j, k, 0, 0, 1.0);
      }
    }
  }
  return result;
}

VIO_Volume sub_volumes_real(VIO_Volume vol1, VIO_Volume vol2){
  int sizes[3];
  int i,j,k;
  VIO_Real val1,val2,diffval,minval,maxval,*data;
  VIO_Volume diff;

  diff = copy_volume( vol1);
  get_volume_sizes( vol1, sizes );

  minval=FLT_MAX;
  maxval=FLT_MIN;

  data=malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(*data));

  for(i=0;i<sizes[0];i++)
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++){
	val1=get_volume_real_value(vol1,i,j,k,0,0);
	val2=get_volume_real_value(vol2,i,j,k,0,0);
	diffval=val1-val2;
	data[i*sizes[1]*sizes[2]+j*sizes[2]+k]=diffval;
	if (diffval > maxval) maxval = diffval;
	if (diffval < minval) minval = diffval;

      }
  
  set_volume_real_range(diff,minval,maxval);

  for(i=0;i<sizes[0];i++)
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++){
	set_volume_real_value(diff, i, j, k, 0, 0, data[i*sizes[1]*sizes[2]+j*sizes[2]+k]);
      }

  free(data);
  return diff;
}


int thresh_vol(VIO_Volume *orig, VIO_Volume *out, VIO_Real thresh){
  int sizes[3];
  int i,j,k;
  *out = copy_volume( *orig);
  get_volume_sizes( *out, sizes );
  

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++) {
	if (get_volume_real_value( *orig, i, j, k,0,0) > thresh)
	{
	  set_volume_real_value( *out, i, j, k, 0, 0, 1); }
	else {
	  set_volume_real_value( *out, i, j, k, 0, 0, 0); }
      }
  }
  
  return(STATUS_OK);
}

VIO_Volume threshold_volume(VIO_Volume orig, VIO_Real min, VIO_Real max){
  int sizes[3];
  int i,j,k;
  VIO_Real value;
  VIO_Volume out;

  get_volume_sizes( orig, sizes );

/*   out = create_volume(3,get_volume_dimension_names(orig),NC_BYTE,FALSE,0,255); */
/*   set_volume_sizes(out,sizes); */
/*   alloc_volume_data(out); */

  out = copy_volume_definition_no_alloc(orig,NC_BYTE,FALSE,0,255);
  set_volume_real_range(out,0,255);
  alloc_volume_data(out);

/*   *out = copy_volume( *orig); */

#ifdef DEBUG
  fprintf(stderr,"min/max: %f/%f\n",min,max);
#endif

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++)
      for(k=0;k<sizes[2];k++) {
	value = get_volume_real_value( orig, i, j, k,0,0);
	if ((value >= min) && (value <= max)){
	  set_volume_real_value( out, i, j, k, 0, 0, 255);
	}else{
	  set_volume_real_value( out, i, j, k, 0, 0, 0);
	}
      }
  }

#ifdef DEBUG
    output_volume("_thresholded.mnc", NC_BYTE,FALSE,0,255,out,NULL,
		   (minc_output_options *)NULL);
#endif

  
  return out;
}




int create_binary_volume(VIO_Volume *vol, VIO_Real *starts, int *sizes){
  *vol = create_volume(3,XYZ_dimension_names,NC_BYTE,FALSE,0.0,1.0);  

  set_volume_sizes(*vol,sizes);

  set_volume_starts(*vol,starts);

  alloc_volume_data(*vol);

  return STATUS_OK;
}

int expand_volume(VIO_Volume *vol, VIO_Volume *expand, point3D *exp){
  int x,y,z,sizes[3],new_sizes[3],expansion[3];
  VIO_Real starts[3],separations[3];

  expansion[0]=(int)exp->x;
  expansion[1]=(int)exp->y;
  expansion[2]=(int)exp->z;

  fprintf(stderr,"Expanding by %d,%d,%d\n",expansion[0],expansion[1],expansion[2]);

  *expand = copy_volume(*vol);

  get_volume_sizes(*vol,sizes);
  
  get_volume_starts(*vol,starts);

  get_volume_separations(*vol, separations);

  //fprintf(stderr,"Separations: (%f,%f,%f)\n",separations[0],separations[1],separations[2]);

  new_sizes[0] = sizes[0] + expansion[0]*2;
  new_sizes[1] = sizes[1] + expansion[1]*2;
  new_sizes[2] = sizes[2] + expansion[2]*2;
  starts[0] -= expansion[0]*VIO_SIGN(separations[0]);
  starts[1] -= expansion[1]*VIO_SIGN(separations[1]);
  starts[2] -= expansion[2]*VIO_SIGN(separations[2]);

  set_volume_sizes(*expand,new_sizes);

  set_volume_starts(*expand,starts);
  
  alloc_volume_data(*expand);

  for (x=0;x<sizes[0];x++)
    for (y=0;y<sizes[1];y++)
      for (z=0;z<sizes[2];z++)
	set_volume_real_value(*expand,x+expansion[0],y+expansion[1],z+expansion[2],0,0,get_volume_real_value(*vol,x,y,z,0,0));
    
  /* write zeros in the expanded part */
  for (x=0;x<new_sizes[0];x++)
    for (y=0;y<new_sizes[1];y++){
      set_volume_real_value(*expand,x,y,0,0,0,0.0);
      set_volume_real_value(*expand,x,y,new_sizes[2]-1,0,0,0.0);
    }
  for (x=0;x<new_sizes[0];x++)
    for (z=0;z<new_sizes[2];z++){
      set_volume_real_value(*expand,x,0,z,0,0,0.0);
      set_volume_real_value(*expand,x,new_sizes[1]-1,z,0,0,0.0);
    }
  for (y=0;y<new_sizes[1];y++)
    for (z=0;z<new_sizes[2];z++){
      set_volume_real_value(*expand,0,y,z,0,0,0.0);
      set_volume_real_value(*expand,new_sizes[0]-1,y,z,0,0,0.0);
    }

  return STATUS_OK;
}

int decrease_volume(VIO_Volume *vol, VIO_Volume *decreased, point3D *dec){
  int x,y,z,sizes[3],new_sizes[3],decrease[3];
  VIO_Real starts[3],separations[3];

  decrease[0]=(int)dec->x;
  decrease[1]=(int)dec->y;
  decrease[2]=(int)dec->z;

  fprintf(stderr,"Decreasing by %d,%d,%d\n",decrease[0],decrease[1],decrease[2]);

  *decreased = copy_volume(*vol);

  get_volume_sizes(*vol,sizes);
  
  get_volume_starts(*vol,starts);

  get_volume_separations(*vol, separations);

  //fprintf(stderr,"Separations: (%f,%f,%f)\n",separations[0],separations[1],separations[2]);

  new_sizes[0] = sizes[0] - decrease[0]*2;
  new_sizes[1] = sizes[1] - decrease[1]*2;
  new_sizes[2] = sizes[2] - decrease[2]*2;
  starts[0] += decrease[0]*VIO_SIGN(separations[0]);
  starts[1] += decrease[1]*VIO_SIGN(separations[1]);
  starts[2] += decrease[2]*VIO_SIGN(separations[2]);

  set_volume_sizes(*decreased,new_sizes);

  set_volume_starts(*decreased,starts);
  
  alloc_volume_data(*decreased);

  for (x=0;x<new_sizes[0];x++)
    for (y=0;y<new_sizes[1];y++)
      for (z=0;z<new_sizes[2];z++)
	set_volume_real_value(*decreased,x,y,z,0,0,get_volume_real_value(*vol,x+decrease[0],y+decrease[0],z+decrease[0],0,0));
    
  return STATUS_OK;
}


int set_intensity_by_mask(VIO_Volume original, VIO_Volume mask, VIO_Real intensity){
  int i,j,k,sizes[3];


#ifdef DEBUG
  fprintf(stdout,"set_intensity_by_mask(): getting volume sizes\n");
#endif

  get_volume_sizes(original,sizes);
#ifdef DEBUG
  fprintf(stdout,"set_intensity_by_mask(): %d %d %d\n",sizes[0],sizes[1],sizes[2]);
  fprintf(stdout,"set_intensity_by_mask(): going into loop\n");
#endif

  for (i=0;i<sizes[0];i++)
    for (j=0;j<sizes[1];j++)
      for (k=0;k<sizes[2];k++)
	if (get_volume_real_value(mask, i, j, k, 0, 0))
	  set_volume_real_value(original,i,j,k,0,0,intensity);

#ifdef DEBUG
  fprintf(stdout,"set_intensity_by_mask(): returning\n");
#endif
  return STATUS_OK;
}

int volume_histogram(VIO_Volume vol, int *histogram, int bins){
  VIO_Volume normalized;
  int i,j,k,sizes[3];
  VIO_Real min,max,value;
  
  get_volume_sizes(vol,sizes);
  normalized = normalize_volume(vol,0.0,(VIO_Real)(bins-1));

  bzero(histogram,bins*sizeof(int));

  for (i=0;i<sizes[0];i++)
    for (j=0;j<sizes[1];j++)
      for (k=0;k<sizes[2];k++){
	value = get_volume_real_value(normalized,i,j,k,0,0);
	histogram[(int)value]++;
      }


  return STATUS_OK;
}

int get_bounding_box(VIO_Volume *vol, point3D *c1, point3D *c2, VIO_Real thresh){
  int sizes[5],i,j,k;
  
  get_volume_sizes(*vol,sizes);

  c1->x = sizes[0];
  c1->y = sizes[1];
  c1->z = sizes[2];
  c2->x = 0;
  c2->y = 0;
  c2->z = 0;
  

  for(i=0;i<sizes[0];i++){
    for(j=0;j<sizes[1];j++){
      for(k=0;k<sizes[2];k++) {
	if (get_volume_real_value( *vol, i, j, k,0,0) > thresh){
	  c1->x = (i < c1->x ? i : c1->x);
	  c1->y = (j < c1->y ? j : c1->y);
	  c1->z = (k < c1->z ? k : c1->z);
	  c2->x = (i > c2->x ? i : c2->x);
	  c2->y = (j > c2->y ? j : c2->y);
	  c2->z = (k > c2->z ? k : c2->z);
	  
	}
      }
    }
  }
  
  return STATUS_OK;

}

VIO_Real min_intensity( VIO_Volume original)
{
  VIO_Real value;
  VIO_Real min;
  int sizes[3];
  int i,j,k;
  //point3D minpos;

  get_volume_sizes(original, sizes);
 
  min=FLT_MAX;

  for(i=1;i<sizes[0]-1;i++)
    for(j=1;j<sizes[1]-1;j++)
      for(k=1;k<sizes[2]-1;k++)
	{	  
	  value=get_volume_real_value(original, i, j, k, 0, 0);
	  if (min<value)
	    {
	      min=value;
	      //SET_3DPOINT(minpos,i,j,k);
	    }
	}

  //PRINTVEC(minpos);
  return (min);
}


VIO_Real max_intensity( VIO_Volume original)
{
  VIO_Real value;
  VIO_Real max;
  int sizes[3];
  int i,j,k;
  //  point3D minpos;

  get_volume_sizes(original, sizes);

  max=FLT_MIN;

   for(i=1;i<sizes[0]-1;i++)
    for(j=1;j<sizes[1]-1;j++)
      for(k=1;k<sizes[2]-1;k++)
	{
	  {
	    value=get_volume_real_value(original, i, j, k, 0, 0);
	      if (value > max )
		{
		max=value;
		//SET_3DPOINT(minpos,i,j,k);
		}
	    }
	}
   //PRINTVEC(minpos);
   return (max);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
*/