#include "minc2-simple.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const size_t vol_size=10;

int create_file_with_nan(const char *minc2_file_path,int use_short)
{
    int ret=0;
    /*will define my own minc volume*/
    /*fastest varying will go first*/
    struct minc2_dimension my_dims[] = 
    {   {.id=MINC2_DIM_X, .length=vol_size, .irregular=0, .step=1.0, .start=-4.5}, 
        {.id=MINC2_DIM_Y, .length=vol_size, .irregular=0, .step=1.0, .start=-4.5}, 
        {.id=MINC2_DIM_Z, .length=vol_size, .irregular=0, .step=1.0, .start=-4.5},
        {.id=MINC2_DIM_END }
    };
    
    minc2_file_handle minc2_file;
    
    minc2_file=minc2_allocate0();
    
    /*define minc2 volume*/
    minc2_define(minc2_file,my_dims, use_short?MINC2_USHORT:MINC2_FLOAT, MINC2_FLOAT);
    /*we will use global scaling (already setup in previous function actually)*/
    /*minc2_set_scaling(minc2_file,0,0);*/
    fprintf(stdout,"Writing to %s\n",minc2_file_path);
    /*create new minc2 file*/
    if(minc2_create(minc2_file,minc2_file_path)==MINC2_SUCCESS)
    {
        size_t i,j,k;
        double voxel[3];
        double world[3];
        float *volume=(float*)calloc(vol_size*vol_size*vol_size,sizeof(float));
        fprintf(stdout,"initializing whole volume...\n");
        /*fill out the whole volume*/
        for(i=0;i<vol_size;i++)
            for(j=0;j<vol_size;j++)
                for(k=0;k<vol_size;k++)
                {
                    voxel[0]=i;voxel[1]=j;voxel[2]=k;
                    /*convert index to xyz coordinates*/
                    minc2_voxel_to_world(minc2_file,voxel,world);
                    
                    if(i==4 && j==4 && k==3) /*inject +INF*/
                    {
                        volume[i+j*vol_size+k*vol_size*vol_size]=1.0/0;
                        printf("Voxel at %ld,%ld,%ld is %f\n",i,j,k,volume[i+j*vol_size+k*vol_size*vol_size]);
                    } else if(i==4 && j==4 && k==4) /*inject -INF*/
                    {
                        volume[i+j*vol_size+k*vol_size*vol_size]=-1.0/0;
                        printf("Voxel at %ld,%ld,%ld is %f\n",i,j,k,volume[i+j*vol_size+k*vol_size*vol_size]);
                    } else if(i==4 && j==4 && k==5) /*inject +INF*/
                    {
                        volume[i+j*vol_size+k*vol_size*vol_size]=sqrt(-1.0);
                        printf("Voxel at %ld,%ld,%ld is %f\n",i,j,k,volume[i+j*vol_size+k*vol_size*vol_size]);
                    } else {
                        volume[i+j*vol_size+k*vol_size*vol_size]=cos(world[0]*3.14/5.0)*cos(world[1]*3.14/5.0)*cos(world[2]*3.14/5.0)+2.0;
                    }
                }
        fprintf(stdout,"saving whole volume...\n");
        if(minc2_save_complete_volume(minc2_file,volume,MINC2_FLOAT)!=MINC2_SUCCESS)
        {
            fprintf(stderr,"Error writing to %s\n",minc2_file_path);
            ret=1;
        }
        if(minc2_close(minc2_file)!=MINC2_SUCCESS)
        {
            fprintf(stderr,"Error closing %s\n",minc2_file_path);
            ret=1;
        }
        free(volume);
        
    } else {
        fprintf(stderr,"Error creating %s\n",minc2_file_path);
        ret=1;
    }
    minc2_free(minc2_file);
    return ret;
}


int read_file_with_nan(const char *minc2_file_path)
{
    int ret=0;
    /*will define my own minc volume*/
    /*fastest varying will go first*/
    minc2_file_handle minc2_file;
    
    minc2_file=minc2_allocate0();
    
    /*define minc2 volume*/
    
    /*we will use global scaling (already setup in previous function actually)*/
    /*minc2_set_scaling(minc2_file,0,0);*/
    fprintf(stdout,"Reading from %s\n",minc2_file_path);
    /*create new minc2 file*/
    if(minc2_open(minc2_file,minc2_file_path)==MINC2_SUCCESS)
    {
        size_t i,j,k;
        double voxel[3];
        double world[3];
        float *volume=(float*)calloc(vol_size*vol_size*vol_size,sizeof(float));
        fprintf(stdout,"Reading whole volume...\n");
        
        if(minc2_load_complete_volume(minc2_file,volume,MINC2_FLOAT)!=MINC2_SUCCESS)
        {
            fprintf(stderr,"Error writing to %s\n",minc2_file_path);
            ret=1;
        }

        
        /*check three NaN variants*/
        /*fill out the whole volume*/
        if(fpclassify(volume[4+4*vol_size+3*vol_size*vol_size])!=FP_INFINITE)
        {
            fprintf(stderr,"Check failed for voxel %d,%d,%d expected %f got %f\n",4,4,3,1.0/0,volume[4+4*vol_size+3*vol_size*vol_size]);
            ret=1;
        }
        if(fpclassify(volume[4+4*vol_size+4*vol_size*vol_size])!=FP_INFINITE)
        {
            fprintf(stderr,"Check failed for voxel %d,%d,%d expected %f got %f\n",4,4,4,-1.0/0,volume[4+4*vol_size+3*vol_size*vol_size]);
            ret=1;
        }
        if(fpclassify(volume[4+4*vol_size+5*vol_size*vol_size])!=FP_NAN)
        {
            fprintf(stderr,"Check failed for voxel %d,%d,%d expected %f got %f\n",4,4,4,sqrt(-1.0),volume[4+4*vol_size+3*vol_size*vol_size]);
            ret=1;
        }
        if(minc2_close(minc2_file)!=MINC2_SUCCESS)
        {
            fprintf(stderr,"Error closing %s\n",minc2_file_path);
            ret=1;
        }
        free(volume);
        
    } else {
        fprintf(stderr,"Error opening %s\n",minc2_file_path);
        ret=1;
    }
    minc2_free(minc2_file);
    return ret;
}


int main(int argc,char **argv)
{
    int ret;
    int use_short=0;
    const char *fname="test_nan.mnc";
    if(argc==1)
        fprintf(stdout,"Using filename %s, with floating point volume\nUsage: %s [fname] [use_float]\n",fname,argv[0]);
    
    if(argc>1)
        fname=argv[1];
    if(argc>2)
        use_short=atoi(argv[2]);
    
    
    ret=create_file_with_nan(fname,use_short);
    ret=ret||read_file_with_nan(fname);
    
    return ret;
}
