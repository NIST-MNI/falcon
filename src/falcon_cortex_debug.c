/* Filename:     falcon_cortex_debug.c
 * Description:  functions for surface debugging and tracing
 * Author:       Vladimir S. FONOV
 * Date:         Oct 17, 2018
 */

#include "falcon.h"

#include "falcon_cortex.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

int falcon_tracing_init(nifti_image * img, 
    cortex_tracing_info *info) 
    {
    const char *falcon_trace, *trace_x,*trace_y,*trace_z;
    char **CSlist=NULL;
    int n,i;
    double imin,imax;
    double scale=1.0;

    memset(info,0, sizeof(cortex_tracing_info));

    if( (falcon_trace = getenv("FALCON_TRACE"))!=NULL) {
      const char *_scale;
      info->prefix=strdup(falcon_trace);/*TODO: use strdup?*/

      info->slices_x.slice_dir=0;
      info->slices_y.slice_dir=1;
      info->slices_z.slice_dir=2;

      if( (_scale = getenv("FALCON_TRACE_SCALE"))!=NULL) {
         scale = atof(_scale);
      }

      imin=niik_image_get_min(img,NULL);
      imax=niik_image_get_max(img,NULL);

      info->slices_x.compression=info->slices_y.compression=info->slices_z.compression=1;

      if((trace_x = getenv("FALCON_TRACE_X"))!=NULL) {
            NIIK_EXIT(((CSlist=niik_csstring_to_list(trace_x,&n))==NULL),fcname,"niik_csstring_to_list",9);
            info->slices_x.slice=(int *)calloc(n,sizeof(int));
            info->slices_x.slice_num=n;
            for(i=0; i<n; i++) {
                info->slices_x.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
      }

      if((trace_y = getenv("FALCON_TRACE_Y"))!=NULL) {
            NIIK_EXIT(((CSlist=niik_csstring_to_list(trace_y,&n))==NULL),fcname,"niik_csstring_to_list",9);
            info->slices_y.slice=(int *)calloc(n,sizeof(int));
            info->slices_y.slice_num=n;
            for(i=0; i<n; i++) {
                info->slices_y.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
      }

      if((trace_z = getenv("FALCON_TRACE_Z"))!=NULL) {
            NIIK_EXIT(((CSlist=niik_csstring_to_list(trace_z,&n))==NULL),fcname,"niik_csstring_to_list",9);
            info->slices_z.slice=(int *)calloc(n,sizeof(int));
            info->slices_z.slice_num=n;
            for(i=0; i<n; i++) {
                info->slices_z.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
      }
      info->slices_x.dmin  = info->slices_y.dmin  = info->slices_z.dmin  = imin;
      info->slices_x.dmax  = info->slices_y.dmax  = info->slices_z.dmax  = imax;
      info->slices_x.scale = info->slices_y.scale = info->slices_z.scale = scale;

      info->dump_slices = trace_x || trace_y || trace_z;
      if(getenv("FALCON_TRACE_SURF")!=NULL) {
        info->dump_surfaces=1;
      }
      return 1;
    }
    return 0;
}


void falcon_tracing_dump_objects(
    cortex_tracing_info *info,
    int iter, const char *task,
    kobj *obj[],
    int nobj)
{
    char tracing_prefix[4096];
    int i;
    if(! info->dump_surfaces )
        return;

    for(i=0;i<nobj;i++)
    {
        sprintf(tracing_prefix,"%s_%s_%03d_%d.ply",info->prefix,task,iter+1,i);
        off_kobj_write_offply(tracing_prefix,obj[i],0); /*TODO: attach scalars too?*/
    }
}

void falcon_tracing_dump(cortex_tracing_info *info,int iter, const char *task, nifti_image * img,bbox *bb) {
    char tracing_prefix[4096];
    if(! info->dump_slices)
        return;
    
    sprintf(tracing_prefix,"%s_%s_%03d_%%s_%%03d.tiff",info->prefix,task,iter+1);
    if(info->slices_z.slice_num) { /* z-slices */
        info->slices_z.fpattern=tracing_prefix;
        niik_tiff_write_slices_obj(&info->slices_z,img,bb);
    } 

    if(info->slices_y.slice_num) { /* y-slices */
        info->slices_y.fpattern=tracing_prefix;
        niik_tiff_write_slices_obj(&info->slices_y,img,bb);
    } 
    if(info->slices_x.slice_num) { /* x-slices */
        info->slices_x.fpattern=tracing_prefix;
        niik_tiff_write_slices_obj(&info->slices_x,img,bb);
    } 
}

void falcon_tracing_free(cortex_tracing_info *info) {
    if(info->slices_x.slice_num) 
        free(info->slices_x.slice);
    if(info->slices_y.slice_num) 
        free(info->slices_y.slice);
    if(info->slices_z.slice_num) 
        free(info->slices_z.slice);
    if(info->prefix)
        free(info->prefix);
}
