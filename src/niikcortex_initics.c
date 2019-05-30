/* FILENAME:     niikcortex_initics.c
 * DESCRIPTION:  Kunio's nifti1 for cortical deformation
 * AUTHOR:       Kunio Nakamura
 * DATE:         May 8, 2014
 */

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include "falcon.h"
#include "falcon_cortex.h"


#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (1)

void prog_version() {
  fprintf(stdout,"  falcon_cortex_initics history\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  0.0.0  May 8, 2014, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n");
  fprintf(stdout,"  0.0.1  August 26, 2014, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n");
  fprintf(stdout,"  0.0.2  Nov 28, 2017, Vladimir S. FONOV <vladimir.fonov@gmail.com>\n");
  fprintf(stdout,"\n");
}

void prog_usage() {
  fprintf(stdout,"  falcon_cortex_initics:\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help                   : show this usage\n");
  fprintf(stdout,"  --version                         : show version info\n");
  fprintf(stdout,"  --dist <dist.nii>                  : distance map image [default = auto]\n");
  fprintf(stdout,"  --debug-keep-tmp                   : keep tmp files [default = don't keep]\n");
  fprintf(stdout,"  --iter <iter>                      : # iterations [default = 300]\n");
  fprintf(stdout,"  --smooth <f>                       : surface smoothing [default = 0.0]\n");
  fprintf(stdout,"  --elen   <f>                       : target average element length [default = 0.8]\n");
  fprintf(stdout,"  --start-elen   <f>                 : starting average element length [default = 1.5]\n");
  fprintf(stdout,"  --noremesh                         : noremesh\n");

}

void usage() {
  fprintf(stdout,"falcon_cortex_initics\n");
  fprintf(stdout,"  usage: [options] <img> <gwi_mask> <laplacemap> <in.off> <out.off>\n\n");
  fprintf(stdout,"\n");
}


typedef struct {
  nifti_image *gwi_mask;           /* gray-white interface mask */
  nifti_image *t1w;                /* t1w image */
  nifti_image *distancemap;        /* distance map */
  nifti_image *laplacemap;         /* laplace map */
  kobj *obj;                       /* surface object */
  double max_distancemap;          /* max distance on distance map */
  double init_len,fin_len;         /* initial and final edge lengths */
  int iter;                        /* max iteration */
  double step;                     /* step size for each iteration */
  int iter_laplace_step;           /* interval between laplacican normals */
  int iter_remesh_step;            /* interval between remeshing */
  int bbox_depth;                  /* depth for bounding box */
  int verbose;                     /* verbose level */
  int debug_keep_tmp; /* debug flag */
  double smoothing;                /* surface smoothing*/
  int noremesh; 
} niikcortex_initics;              /* initial inner cortical surface deformation using shrink-wrpap */


niikcortex_initics *niikcortex_initics_init() {
  niikcortex_initics *ics=NULL;
  const char *fcname="niikcortex_initics_init";
  NIIK_RET0(((ics=(niikcortex_initics *)calloc(1,sizeof(niikcortex_initics)))==NULL),fcname,
            "calloc for niikcortex_initics");
  ics->gwi_mask=
    ics->t1w=
      ics->distancemap=
        ics->laplacemap=
          NULL;
  ics->obj=NULL;
  ics->max_distancemap=10.0;
  ics->init_len = 1.5;
  ics->fin_len = 0.8;
  ics->iter = 300;
  ics->step = 0.2;
  ics->iter_laplace_step = 2;
  ics->iter_remesh_step  = 5;
  ics->bbox_depth=7;
  ics->verbose=2;
  ics->debug_keep_tmp=0;
  ics->smoothing=0.0;
  ics->noremesh=0;
  return ics;
} /* niikcortex_initics_init */

niikcortex_initics *niikcortex_initics_free(niikcortex_initics *ics) {
  if(ics==NULL) return NULL;
  niik_image_free(ics->gwi_mask);
  niik_image_free(ics->t1w);
  niik_image_free(ics->laplacemap);
  niik_image_free(ics->distancemap);
  off_kobj_free(ics->obj);
  free(ics);
  return NULL;
} /* niikcortex_initics_free */


int niikcortex_initics_process_check_deform_interp(niikcortex_initics *ics, kvert *v, double thresh)
/* -returns zero if the vertex is on the object
 * -otherwise return non-zero
 * -copied from nifti1_kunio_cortex_initics.c
 */
{
  const char *fcname="niikcortex_initics_process_check_deform_interp";
  niikpt
    face_pt_list[64],
    edge_pt_list[64];
  double d;
  int
    n,
    xsc=0;
  
  if( niik_image_interpolate_3d_linear(ics->distancemap,v->v)>thresh ) {
    return 0;
  }
  /*
   * check for mid-edge points
   */
  if(ics->verbose>2) {
    fprintf(stdout,"[%s] check mid-edge points\n",fcname);
  }


  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    edge_pt_list[n] = niikpt_avg(v->neiedge[n]->endpts[0]->v,v->neiedge[n]->endpts[1]->v);
    if(niik_image_interpolate_3d_linear(ics->distancemap,edge_pt_list[n])>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for mid-face points
   */
  if(ics->verbose>2) {
    fprintf(stdout,"[%s] check mid-face points\n",fcname);
  }

  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    face_pt_list[n] = niikpt_avg3(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    if(niik_image_interpolate_3d_linear(ics->distancemap,face_pt_list[n])>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for mid (mid-face)-vertex points
   */
  if(ics->verbose>2) {
    fprintf(stdout,"[%s] check mid-face points\n",fcname);
  }

  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[0]->v))>thresh) {
      xsc++;
      continue;
    }
    if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[1]->v))>thresh) {
      xsc++;
      continue;
    }
    if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[2]->v))>thresh) {
      xsc++;
      continue;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for 1/4 and 3/4-edge points
   */
  if(ics->verbose>2) {
    fprintf(stdout,"[%s] check along 1/4 3/4 edge points\n",fcname);
  }

  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_avg(edge_pt_list[n],v->neiedge[n]->endpts[0]->v))>thresh) {
      xsc++;
    }
    if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_avg(edge_pt_list[n],v->neiedge[n]->endpts[1]->v))>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for 3/4-edge points
   */
  if(ics->verbose>2) {
    fprintf(stdout,"[%s] check along 3/4\n",fcname);
  }
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(v->neiface[n]->vert[0]==v) {
      for(d=0.2; d<0.999; d+=0.2) {
        if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
    if(v->neiface[n]->vert[1]==v) {
      for(d=0.2; d<0.998; d+=0.2) {
        if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[2]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
    if(v->neiface[n]->vert[2]==v) {
      for(d=0.2; d<0.998; d+=0.2) {
        if(niik_image_interpolate_3d_linear(ics->distancemap,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
  }
  if(xsc) {
    return 0;
  }
  return 1;
} /* niikcortex_initics_process_check_deform_interp */


double niikcortex_initics_process_check_deform(niikcortex_initics *ics,bbox *bb,kvert *v,double thresh)
/*
 * -moves the vertex position by v->normal * ics->step
 * -check intensity on the neighboring triangles
 * -check self-intersection
 * -if no problem, then keeps that position
 * -otherwise, vertex does not move
 * -retruns step for the actual distance moved
 *
 * -copied from nifti1_kunio_cortex_initics.c
 * -checks for potential collisions
 */
{
  niikpt origpt;
  niikpt update_smooth=niikpt_zero();

  origpt=v->v;

  if(ics->smoothing>0.0) {
    int nei;
    for(nei=0; nei<v->nei; nei++)
      update_smooth=niikpt_add(update_smooth,v->neivert[nei]->v);
    update_smooth=niikpt_kmul(update_smooth,1.0/(double)v->nei);
    update_smooth=niikpt_sub(update_smooth, v->v);
  }

  /* go outward */
  if(niikcortex_initics_process_check_deform_interp(ics,v,thresh)==0) {
    /* check there is surface intersection */
    if(ics->verbose>2) {
      fprintf(stdout,"\t\tmove outward\n");
    }

    v->v=niikpt_move_normal(v->v, v->normal,    ics->step*(1.0-ics->smoothing));
    v->v=niikpt_move_normal(v->v, update_smooth,ics->step*ics->smoothing);

    off_update_kvert_pminmax(v);
    if(off_check_self_intersection_kvert(bb,v))  {
      v->v=origpt;
      off_update_kvert_pminmax(v);
      return 0;
    }
    return ics->step;
  }

  /* go inward */
  if(ics->verbose>2) {
    fprintf(stdout,"\t\tmove inward %.3f\n",thresh);
  }

  v->v=niikpt_move_normal(v->v,v->normal,-ics->step*(1.0-ics->smoothing));
  v->v=niikpt_move_normal(v->v, update_smooth,ics->step*ics->smoothing);
  if(niikcortex_initics_process_check_deform_interp(ics,v,thresh)==0) {
    v->v=origpt;
    return 0;
  }

  /* check there is surface intersection */
  if(ics->verbose>2) {
    fprintf(stdout,"\t\tno object\n");
  }
  off_update_kvert_pminmax(v);
  if(off_check_self_intersection_kvert(bb,v))  {
    v->v=origpt;
    off_update_kvert_pminmax(v);
    return 0;
  }

  if(ics->verbose>2) {
    fprintf(stdout,"\t\tno triangle intersection\n");
  }
  /* check if there WILL be surface intersection */
  if(ics->verbose>2) {
    fprintf(stdout,"\t\tno intersection\n");
  }

  v->v=niikpt_move_normal(origpt,v->normal,-NIIK_DMINMAX(5*ics->step,1,2)*(1.0-ics->smoothing));
  v->v=niikpt_move_normal(v->v, update_smooth,ics->step*ics->smoothing);

  off_update_kvert_pminmax(v);
  if(off_check_self_intersection_kvert(bb,v))  {
    v->v=origpt;
    off_update_kvert_pminmax(v);
    return 0;
  }
  if(ics->verbose>2) {
    fprintf(stdout,"\t\tno triangle intersection\n");
  }

  v->v=niikpt_move_normal(origpt,v->normal,-ics->step*(1.0-ics->smoothing));
  v->v=niikpt_move_normal(v->v, update_smooth,ics->step*ics->smoothing);

  off_update_kvert_pminmax(v);
  return ics->step;
} /* int niikcortex_initics_shrink_wrap_deform */


int niikcortex_initics_process(niikcortex_initics *ics) {
  const char *fcname="niikcortex_initics_process";
  char tempout[64];
  int
    vpos, /* vertex position (index) */
    iter;
  double
      dval,dsum,
       thresh,
       elen;
  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  kvert *v;
  niikpt sf;   /* sobel filter value */
  bbox *bb=NULL;
  cortex_tracing_info trace;

  int debug_tracing=falcon_tracing_init(ics->gwi_mask, &trace);


  niik_fc_display(fcname,1);
  NIIK_RET0((ics->gwi_mask==NULL),fcname,"missing gwi mask");
  NIIK_RET0((ics->laplacemap==NULL),fcname,"missing laplacemap");
  NIIK_RET0((ics->obj==NULL),fcname,"missing object");

  if(ics->distancemap==NULL) {
    fprintf(stdout,"[%s] distance map %.3f\n",fcname,ics->max_distancemap);
    if((ics->distancemap=niik_image_distance_map(ics->gwi_mask,ics->max_distancemap))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_distance_map\n",fcname);
      return 0;
    }
  }

  if(ics->verbose>0) fprintf(stdout,"[%s] bounding box\n",fcname);
  bb=off_bbox_init(ics->bbox_depth,320);
  if(ics->verbose>0) fprintf(stdout,"[%s]   bbox %8.4f\n",fcname,bb->delta);
  off_update_kobj_kface_pminmax(ics->obj);
  off_create_bbox_from_kobj(bb,ics->obj);


  /* MAIN LOOP */
  for(iter=0; iter<ics->iter; iter++) {
    /* prep iter */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    elen = exp(-6.0*iter/(ics->iter-1.0))*(ics->init_len - ics->fin_len) + ics->fin_len - exp(-6.0);
    thresh=cos(19.0*iter/ics->iter) * exp(-3.0*iter/ics->iter) * NIIK_Heaviside(iter-ics->iter/2.0,ics->iter/4.0);

    if(ics->verbose>0) fprintf(stdout,"[%s] iteration %4i %s  elen %7.4f; thresh %7.4f\n",fcname,iter,tmstr,elen,thresh);

    /* normal vector calculation */
    if(ics->verbose>1) fprintf(stdout,"[%s]     update normals\n",fcname);
    off_update_kobj_face_normal(ics->obj);
    off_update_kobj_vert_normal(ics->obj);
    off_smooth_kobj_vert_normal(ics->obj);

    /* normal from laplace map */
    if((iter%ics->iter_laplace_step)==0 && iter>=5) {
      if(ics->verbose>0) fprintf(stdout,"[%s]   laplace-normal\n",fcname);
      for(v=ics->obj->vert; v!=NULL; v=v->next) {
        vpos = niik_image_get_index_niikpt(ics->laplacemap,v->v);
        if(vpos<0) continue;
        sf = niik_image_sobel_filter_voxel(ics->laplacemap,vpos);
        if(fabs(sf.x) + fabs(sf.y) + fabs(sf.z) > 1e-3) {
          v->normal = niikpt_unit(sf);
        }
      }
    }

    if(ics->verbose>0) fprintf(stdout,"[%s]   deformation\n",fcname);
    for(v=ics->obj->vert,dsum=0; v!=NULL; v=v->next) {
      dval=niikcortex_initics_process_check_deform(ics,bb,v,thresh);
      if(fabs(dval)>(ics->step/2.0)) v->v.w=0;
      dsum += dval;
    }

    for(v=ics->obj->vert,dval=0; v!=NULL; v=v->next) {
      dval+=(v->v.w==0);
    }
    if(ics->verbose>0) fprintf(stdout,"[%s]           %4i %8.2f | #v (no motion) %i / %i\n",fcname,iter,dsum,(int)dval,ics->obj->nvert);

    if((iter%ics->iter_remesh_step)==0) {
      
      if(!ics->noremesh) {
        if(ics->verbose>0) fprintf(stdout,"[%s] remesh %.4f\n",fcname,elen);
        NIIK_RET0((!off_remesh_kobj(ics->obj,elen,10,0)),fcname,"off_remesh_kobj");
      } else {
         NIIK_RET0((!off_relax_kobj(ics->obj,elen,10,0)),fcname,"off_relax_kobj");
      }

      if(ics->debug_keep_tmp>0) {
        fprintf(stdout,"[%s]   write temp off\n",fcname);
        sprintf(tempout,"tmp_initics_%i.off",iter);
        off_kobj_write_offply(tempout,ics->obj,0);
      }

      /* SELF-INTERSECTION CORRECTION */
      if(ics->verbose>1) fprintf(stdout,"[%s]   check intersections\n",fcname);
      NIIK_RET0((!off_correct_self_intersection_with_elen(bb,ics->obj,ics->fin_len)),fcname,"off_correct_self_intersection");

      /* set w-mark to NOT do remesh */
      if(ics->verbose>1) fprintf(stdout,"[%s]   w-markers\n",fcname);
      for(v=ics->obj->vert; v!=NULL; v=v->next) {
        v->v.w=1;
      }

      /* bounding box */
      if(ics->verbose>1) fprintf(stdout,"[%s]   update bbox\n",fcname);
      off_update_kobj_kface_pminmax(ics->obj);
      off_create_bbox_from_kobj(bb,ics->obj);
    } /* remesh */

    if(debug_tracing)
    {
      NIIK_RET0((!off_kobj_add_one_color(ics->obj,1,1,0)),fcname,"off_kobj_add_one_color"); /* yellow for white matter-surface */

      falcon_tracing_dump(&trace,iter,"initics",ics->gwi_mask,bb);
    } 
  } /* iter */

  if(!ics->noremesh) {
    if(ics->verbose>0) fprintf(stdout,"[%s] remesh %.4f\n",fcname,ics->fin_len);
    NIIK_RET0((!off_remesh_kobj(ics->obj,ics->fin_len,10,0)),fcname,"off_remesh_kobj");
  }  else {
    NIIK_RET0((!off_relax_kobj(ics->obj,ics->fin_len,10,0)),fcname,"off_relax_kobj");
  }

  /* SELF-INTERSECTION CORRECTION */
  NIIK_RET0((!off_correct_self_intersection_with_elen(bb,ics->obj,ics->fin_len)),fcname,"off_correct_self_intersection");

  if(ics->verbose>0) fprintf(stdout,"[%s] add color",fcname);
  off_kobj_add_one_color(ics->obj,0.8,0.8,0);

  bb=off_bbox_free(bb);

  falcon_tracing_free(&trace);
  niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_initics_process */



int main(int argc,char *argv[],char *envp[]) {
  const char *fcname="niikcortex_initics";
  niikcortex_initics *initics=niikcortex_initics_init();
  int    i;
  int    clobber=0;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"version", no_argument, 0, 'v'},
    {"help",    no_argument, 0, 'h'},
    {"usage",   no_argument, 0, 'U'},
    {"debug-keep-tmp", no_argument, &initics->debug_keep_tmp,1},
    {"noremesh", no_argument, &initics->noremesh,1},
    {"dist",    required_argument,  0,   'd'},
    {"iter",    required_argument,  0,   'i'},
    {"smooth",  required_argument,  0,   'S'},
    {"elen",    required_argument,  0,   'e'},
    {"start-elen",    required_argument,  0,   's'},
    {0, 0, 0, 0}
  };
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vUhd:i:S:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      prog_version();
      return 0;
    case 'd':
      NIIK_EXIT(((initics->distancemap=niik_image_read(optarg))==NULL),fcname,"niik_image_read",1);
      break;
    case 'i':
      initics->iter=atoi(optarg);
      break;
    case 'S':
      initics->smoothing=atof(optarg);
      break;
    case 'e':
      initics->fin_len=atof(optarg);
      break;
    case 's':
      initics->init_len=atof(optarg);
      break;
    case 'h':
    case 'U':
    case '?':
    default:
      prog_usage ();
      return 1;
    }
  }

  if((argc - optind)<5) {
    fprintf(stderr,"[%s] ERROR: too few argments (%i < 5)\n",fcname,argc - optind);
    usage();
    prog_usage();
    return 1;
  }

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);


#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  fprintf(stdout,"[%s] reading image            %s\n",fcname,argv[optind]);
  NIIK_EXIT(((initics->t1w=niik_image_read(argv[optind]))==NULL),fcname,"niik_image_read",1);
  fprintf(stdout,"[%s] reading GWI mask         %s\n",fcname,argv[optind+1]);
  NIIK_EXIT(((initics->gwi_mask=niik_image_read(argv[optind+1]))==NULL),fcname,"niik_image_read",1);
  fprintf(stdout,"[%s] reading laplacian map    %s\n",fcname,argv[optind+2]);
  NIIK_EXIT(((initics->laplacemap=niik_image_read(argv[optind+2]))==NULL),fcname,"niik_image_read",1);
  fprintf(stdout,"[%s] reading initial obj      %s\n",fcname,argv[optind+3]);
  NIIK_EXIT(((initics->obj=off_kobj_read_offply(argv[optind+3]))==NULL),fcname,"off_kobj_read_offply",1);

  NIIK_EXIT((!niikcortex_initics_process(initics)),fcname,"niikcortex_initics_process",1);

  NIIK_EXIT((!off_kobj_add_one_color(initics->obj,1.0,1.0,0.2)),fcname,"off_kobj_add_one_color",1);

  /*append metadata*/
  NIIK_EXIT((!off_kobj_add_comment(initics->obj,timestamp)),fcname,"off_kobj_add_comment",1);

  fprintf(stdout,"[%s] writing output  %s\n",fcname,argv[optind+4]);
  NIIK_EXIT((!off_kobj_write_offply(argv[optind+4],initics->obj,0)),fcname,"off_kobj_write_off",1);

  initics=niikcortex_initics_free(initics);
  niik_fc_display(fcname,0);
  free(timestamp);
  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/