/* FILENAME:    falcon_midsag.c
* DESCRIPTION:  Reimlimentation of falcon_midsag using Kunio's library
* DATE:         Nov 23, 2017
*/

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>


#include "falcon.h"
#include "falcon_cortex.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (1)


#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


kobj* off_create_plane(const niikpt *p,double elen);

void prog_version() {
  fprintf(stdout,"  falcon_midsag history\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  0.0.1  Nov 23, 2017, Vladimir S. FONOV - reimplementing\n");
  fprintf(stdout,"\n");
}

void prog_usage() {
  fprintf(stdout,"falcon_midsag\n");
  fprintf(stdout,"  usage: [options] <mask_input> \n");
  fprintf(stdout,"\n");
  fprintf(stdout,"\n  optional usage:\n");
  fprintf(stdout,"  --elen <E>         target edge length [default=2mm]\n");
  fprintf(stdout,"  --border <E>       border width [default=1mm]\n");
  fprintf(stdout,"  --out <obj.off>    Output Midsagittal surface \n");
  fprintf(stdout,"  --left <left mask>  Output Left mask \n");
  fprintf(stdout,"  --right <right mask>  Output Right mask \n");
  fprintf(stdout,"  --iter <I>         iteration [default=20]\n");
  fprintf(stdout,"  --step  <max_step_size> - maximum update step, default 0.5 \n");
  fprintf(stdout,"  --label <n> use this label in the mask file, default  any above 0 \n");
  fprintf(stdout,"  --version          show version info\n");
  fprintf(stdout,"  --clobber          overwrite output file\n");
  fprintf(stdout,"  --debug            Added debugging output\n");

}



int optimize_slice_position(nifti_image *maskimg, int slice, int searchspace, int *no_of_voxels) {
  int i,j,k;
  double min_weight=1e10;
  int min_index=-1;

  for(i=(slice-searchspace); i<(slice+searchspace); i++) {
    int offset=i;
    double weight=0.0;
    for(k=0; k<maskimg->nz; k++) {
      for(j=0; j<maskimg->ny; j++) {
        weight += niik_image_get_voxel(maskimg,offset);
        offset+=maskimg->nx;
      }
    }

    if(weight<min_weight || min_index<0) {
      min_weight=weight;
      min_index=i;
    }
  }
  *no_of_voxels=floor(min_weight);
  return min_index;
}


kobj * off_find_minimal_cut_plane(double elen,nifti_image *maskimg)
/**
 * find a y-z plane with minimal intersection surface
 *
 */
{
  int min_index;
  int min_weight;
  int n;
  int slice=maskimg->nx/2; /*use COM ?*/
  int searchspace=5; /*?*/
  kobj *obj;
  niikpt p[4];
  kvert *v;

  min_index=optimize_slice_position(maskimg,slice,searchspace,&min_weight);
  printf("Found middle plane at :%d weight:%d\n",min_index,min_weight);
  /*TODO: check the order*/
  niik_index_to_world(maskimg, min_index, 0, 0,&p[0]);
  niik_index_to_world(maskimg, min_index, maskimg->ny, 0,&p[1]);
  niik_index_to_world(maskimg, min_index, maskimg->ny, maskimg->nz,&p[2]);
  niik_index_to_world(maskimg, min_index, 0, maskimg->nz,&p[3]);

  obj=off_create_plane(p,elen);

  /*remesh to target elen*/
  /*off_remesh_kobj(obj,elen,5,1);*/

  return obj;
}

int off_deform_minimal_cut_plane(nifti_image *img,kobj *obj, double len, int maxiter, double stepsize, int flag_bbox, int flag_remesh, int debug)
/*
* deform minimal cut plane to cut as little as possible (i.e ideally cutting midsagittal junktion only)
* -may/may not use bounding box
*
*/
{
  nifti_image *tmpimg = NULL;
  bbox *bb=NULL;
  kvert *v;
  niikpt origpt;
  double
  tol=0.01,
  stepsum,
  err,
  step;
  int
  iter1,
  iter2,
  midx,
  n,
  ct;
  char fname[512];
  const char *fcname="off_deform_minimal_cut_plane";
  int verbose=2;//niik_verbose();
  double FWHM=2.0;
  double tolerance=0.01;
  cortex_tracing_info trace;
  int debug_trace=falcon_tracing_init(img,&trace);

  nifti_image  **grad_img=NULL;
  nifti_image   *blur_img=NULL;

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }

  if(verbose>1) fprintf(stdout,"[off_shrink_wrap_kobj_bbox_remesh] start off_shrinkwrap_kobj_simple\n");

  if(debug) {
    verbose=2;
    /* make a copy image for uint8 datatype */
    if((tmpimg=niik_image_copy(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"ERROR: niik_image_type_convert\n");
      return 0;
    }
  }

  if(flag_bbox) {
    if(verbose>1) fprintf(stdout,"[%s] bounding box\n",fcname);
    bb=off_bbox_init(6,320);
    if(verbose>1) fprintf(stdout,"[%s] bbox = delta=%f depth=%2i\n",fcname,bb->delta,bb->depth);
  } else {
    if(verbose>1) fprintf(stdout,"[%s] not using bounding box\n",fcname);
  }

  NIIK_RET0(((blur_img=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");

  if(FWHM>0) {
    if(verbose>1) fprintf(stdout,"[%s]   gaussian filter\n",fcname);
    NIIK_RET0((!niik_image_filter_gaussian_update(blur_img,11,FWHM)),fcname,"niik_image_filter_gaussian_update");
  } else if(verbose>1) {
    fprintf(stdout,"[%s]   no gaussian filter\n",fcname);
  }

  NIIK_RET0(((grad_img=niik_image_sobel_filters_with_mag(blur_img))==NULL),fcname,"niik_image_sobel_filters_with_mag");
  blur_img=niik_image_free(blur_img);

  /* determine the step size based on the image voxel size */
  /*step = NIIK_DMIN(img->dx,NIIK_DMIN(img->dy,img->dz)) * 0.8;*/
  step = stepsize / 10.0;
  if(verbose>1) fprintf(stdout,"[%s] step %8.4f\n",fcname,step);

  /* check for self intersection */
  if(flag_bbox) {
    if((n=off_count_self_intersection(bb,obj))>0) {
      fprintf(stderr,"ERROR: self-intersection, %i\n",n);
      return 0;
    } else {
      if(verbose>=1) fprintf(stdout,"        no self-intersection\n");
    }
  }

  /* this is the outer loop */
  for(iter2=1; iter2<=maxiter; iter2++) {

    if(verbose>1) fprintf(stdout,"[%s] Deform %3i\n",fcname,iter2);
    /* update face and vert normal + pmin/pmax */
    if(verbose>1) fprintf(stdout,"[%s] update normal/pmin/pmax\n",fcname);

    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_update_kobj_kface_pminmax(obj);
    off_kobj_update_all_index(obj);

    if(flag_bbox) {
      if(verbose>1) fprintf(stdout,"[%s] bounding box recals\n",fcname);
      off_create_bbox_from_kobj(bb,obj);
      
      if(debug_trace)
        falcon_tracing_dump(&trace,iter2,"midsag",img,bb);
    }


    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v.w=1.0;
    }

    for(stepsum=0,iter1=1; stepsum<stepsize; stepsum+=step,iter1++) {
      int total_inside=0;
      for(v=obj->vert,ct=0; v!=NULL; v=v->next) {
        niikpt normal;
        double normal_mag;
        origpt = v->v;

        /*move point towards higher distance from the object (i.e towards positive gradient)*/
        normal.x = niik_image_interpolate_3d(grad_img[1],v->v,NIIK_INTERP_LINEAR);
        normal.y = niik_image_interpolate_3d(grad_img[2],v->v,NIIK_INTERP_LINEAR);
        normal.z = niik_image_interpolate_3d(grad_img[3],v->v,NIIK_INTERP_LINEAR);

        if(niik_image_interpolate_3d(img,v->v,NIIK_INTERP_LINEAR)>0.0) // inside
          total_inside++;

        /*TODO: update normal to be collinear with surface normal?*/

        //normal = niikpt_unit(normal);
        normal_mag = niikpt_mag(normal);

        if(normal_mag>tolerance) {

          v->v = niikpt_move_normal(origpt, normal, -1.0*step);

          if(flag_bbox) {
            off_update_kvert_pminmax(v);
            if(off_check_self_intersection_kvert(bb,v)) {
              midx=1;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else {
              midx=2;
              v->v.w = 0;
              ct++;
            }
          } /* using bounding box */
        } /* vertex is in the object */
      } /* each vertex */

      err=(double)ct/obj->nvert;

      if(verbose>1) fprintf(stdout,"[%s] Deform %3i  %5.2f %4.2e   %8i / %8i inside:%d  \n",fcname,iter1,stepsum,err,ct,obj->nvert,total_inside);
      /*if(err<=tol) break; */
    } /* inner loop */

    /* remesh */
    if(flag_remesh) {
      if(!off_remesh_kobj(obj,len,5,0)) {
        fprintf(stderr,"ERROR: off_remesh_kobj \n");
        return 0;
      }
    }

    if(flag_bbox) {
      if(!off_correct_self_intersection(bb,obj)) {
        fprintf(stderr,"ERROR: off_correct_self_intersection\n");
        return 0;
      }
      if((n=off_count_self_intersection(bb,obj))>0) {
        fprintf(stderr,"ERROR: self-intersection, %i\n",n);
        return 0;
      } else {
        if(verbose>=1) fprintf(stdout,"[off_shrink_wrap_kobj_bbox_remesh] no self-intersection\n");
      }
    } /* using bounding box */

    if(debug) {
      sprintf(fname,"tmp_iter%03i.off",iter2);
      fprintf(stdout,"  writing %s\n",fname);
      if(!off_kobj_write_offply(fname,obj,0)) {
        fprintf(stderr,"ERROR: off_kobj_write_off \n");
        exit(0);
      }
      off_obj2img(tmpimg,obj,iter2+10);
    }
  } /* larger loop */

  if(debug) {
    fprintf(stdout,"\twriting image:    tmp_shrink_wrap.nii.gza\n");
    niik_image_write("tmp_midsag_deform.nii.gz",tmpimg);
    nifti_image_free( tmpimg );
    tmpimg=NULL;
  }

  if(flag_bbox)
    bb=off_bbox_free(bb);

  for(n=0; n<4; n++) {
    grad_img[n] = niik_image_free(grad_img[n]);
  }

  free(grad_img);

  falcon_tracing_free(&trace);

  return 1;
} /* off_shrinkwrap_kobj */

double face_to_point(kface *cf,niikpt pt) {
  /*dumb algorithm, TODO: update?*/
  niikpt ctr=niikpt_avg3(cf->vert[0]->v, cf->vert[1]->v, cf->vert[2]->v);
  return niikpt_distance2(ctr,pt);
}

int check_left_right_side(bbox *bb, niikpt pt,double border) {
  int i,j,k;

  i = floor( (pt.x - bb->origin.x)/bb->delta );
  j = floor( (pt.y - bb->origin.y)/bb->delta );
  k = floor( (pt.z - bb->origin.z)/bb->delta );

  if( i < 0 ) return -1; //left
  else if(i >= bb->xdim ) return 1; // right
  else { //have to check
    int index,n;
    int closest_index   = -1;
    double closest_dist = 1e10;
    niikpt ctr;

    if( j<0 || j>=bb->ydim ||
        k<0 || k>=bb->zdim ) return 0; // out of the bbox - probably between?

    index = i + j*bb->xdim + k*bb->area;
    //search for the closest face?

    for(n=0; n<bb->ndata[index]; n++) {
      kface *cf=bb->data[index][n];
      ctr=niikpt_avg3( cf->vert[0]->v,cf->vert[1]->v, cf->vert[2]->v);

      double dist=niikpt_distance2(pt,ctr);
      double proj=niikpt_dot(niikpt_sub(pt,ctr),cf->normal);

      if( proj>0.0 && (dist<closest_dist || closest_index<0) ) {
        closest_index = n;
        closest_dist  = dist;
      }
    }

    if(closest_index<0) {
      /*bbox was empty, find the closes non-empty and check*/
      int ii,jj,kk;
      int side=0;
      closest_dist=1000;
      kk=k;
      jj=j;

      for(ii=0; ii<bb->xdim; ii++) {
        if(bb->ndata[ ii + jj*bb->xdim + kk*bb->area ]>0) {
          double dist=fabs(ii-i); /*  *(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k);*/

          if( dist < closest_dist ) {
            closest_dist=dist;
            side=i-ii;
          }
        }
      }
      //printf("found side:%d dist:%f",side,sqrt(closest_dist));
      return side;
    }
    ctr=niikpt_avg3( bb->data[index][closest_index]->vert[0]->v,
                     bb->data[index][closest_index]->vert[1]->v,
                     bb->data[index][closest_index]->vert[2]->v);

    if(pt.x< (ctr.x-border/2)) return -1;
    else if(pt.x>(ctr.x+border/2)) return 1;
  }
  return 0;
}


int off_split_left_right(nifti_image* maskimg, kobj * obj,double border,nifti_image* left_img,nifti_image* right_img  ) {
  int k;
  int offset=0;
  niikpt pmin,pmax; // bounding box
  bbox *bb=off_bbox_init(7,320);
  off_create_bbox_from_kobj(bb,obj);

  pmin=off_calc_kobj_pmin(obj);
  pmax=off_calc_kobj_pmax(obj);
  niikpt_disp(pmin);
  niikpt_disp(pmax);

  printf("delta=%f\n",bb->delta);

  for(k=0; k<maskimg->nz; k++) {
    int j,i;
    for(j=0; j<maskimg->ny; j++) {
      for(i=0; i<maskimg->nx; i++) {
        niikpt p;
        double vox;
        niik_index_to_world(maskimg, i, j, k,&p);
        vox=niik_image_get_voxel(maskimg, offset);

        if(vox>0.0) {
          int side=0;

          if     (p.x < pmin.x-border/2) {
            side=-1;
          } else if(p.x > pmax.x+border/2) {
            side=1;
          } else if(p.y < pmin.y || p.y > pmax.y || p.z < pmin.z || p.z > pmax.z) {
            side=0;
          } else {
            side=check_left_right_side(bb,p,border);  /*update tolerance*/
          }

          if       (side<0) {
            niik_image_set_voxel(left_img,offset, 1.0);
            niik_image_set_voxel(right_img,offset,0.0);
          } else if(side>0) {
            niik_image_set_voxel(left_img,offset, 0.0);
            niik_image_set_voxel(right_img,offset,1.0);
          } else {
            niik_image_set_voxel(left_img, offset,0.0);
            niik_image_set_voxel(right_img,offset,0.0);
          }

        } else {
          niik_image_set_voxel(left_img, offset,0);
          niik_image_set_voxel(right_img,offset,0);
        }

        offset++;
      }
    }
  }
  off_bbox_free(bb);
  return 1;
}


/**
 * Create a rectangle object betwen four edges with minimal thickness
 * HACK! TODO:come up with better way to create thin planes
 */
kobj* off_create_plane(const niikpt *p,double elen) {
  const char *fcname="off_create_plane";
  kobj *obj;
  kvert *v,**vlist;

  kface *f;
  kedge *e;
  int n,m,j,i;

#if 1
  double mlen=niikpt_mag(niikpt_sub(p[0],p[1]));
  niikpt norm=niikpt_unit_normal(p[0],p[1],p[3]);
  niikpt vert[8];

  for(n=0; n<4; n++) {
    niikpt to_center=niikpt_unit( niikpt_sub(p[(n+2)%4],p[n]));

    /*vert[n+4] = niikpt_move_normal(p[n], norm, -0.1*elen);
    vert[n].w = 1;*/

    vert[n]=   niikpt_move_normal(niikpt_move_normal(p[n], to_center, elen),norm,0.01*elen);
    /*vert[n]=niikpt_move_normal(niikpt_move_normal(p[n], to_center, elen),norm,0.01);*/
    vert[n+4]=p[n];
    /*vert[n]   = niikpt_move_normal(p[n], norm,  -elen/2.0);
    vert[n+3] = niikpt_move_normal(p[n], norm, elen/2.0);*/
  }
  obj=off_create_rect(vert);

  /*return obj;*/

  /*now remesh to target elen*/
  while(mlen>elen) {
    if(!off_subdiv_kobj(obj)) {
      fprintf(stderr,"ERROR: off_subdiv_kobj\n");
      return NULL;
    }
    mlen = off_get_kobj_mean_edge_length(obj);
    fprintf(stdout,"[%s] mean elen = %9.5f \n",fcname,mlen);
  }

  return obj;
#else
  double dx=niikpt_distance(p[0],p[1]);
  double dy=niikpt_distance(p[1],p[2]);
  int nx=floor(dx/elen);
  int ny=floor(dy/elen);

  obj = off_obj_init();
  obj->vert = v = off_vert_init();
  vlist=(kvert**)calloc(nx*ny,sizeof(kvert*));
  //create regular mesh of points
  for(j=0,n=0; j<ny; j++)
    for(i=0; i<nx; i++) {
      if(i || j) {
        v->next = off_vert_init();
        v->next->prev=v;
        v = v->next;
      }
      /*bilinear interpolation*/
      niikpt px1 = niikpt_wavg(p[0],p[1],i*elen/dx);
      niikpt px2 = niikpt_wavg(p[3],p[2],i*elen/dx);
      niikpt pt  = niikpt_wavg(px1,px2,j*elen/dy);

      v->v=pt;
      vlist[n++]=v;
    }

  obj->face = f = off_face_init();

  for(j=0; j<(ny-1); j++)
    for(i=0; i<(nx-1); i++) {
      if(i || j) {
        f->next = off_face_init();
        f->next->prev = f;
        f=f->next;
      }

      f->vert[0]=vlist[i+j*nx];
      f->vert[1]=vlist[i+1+j*nx];
      f->vert[2]=vlist[i+1+(j+1)*nx];

      f->next = off_face_init();
      f->next->prev = f;
      f=f->next;

      f->vert[0]=vlist[i+1+j*nx];
      f->vert[1]=vlist[i+1+(j+1)*nx];
      f->vert[2]=vlist[i+(j+1)*nx];
    }



  NIIK_RET0((!off_kobj_read_off_make_edge(obj)),fcname,"off_kobj_read_off_make_edge");

  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[0],e)),fcname,"off_kobj_update_vertex_nei");
    if(!e->endpts[1]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[1],e)),fcname,"off_kobj_update_vertex_nei");
  }
  off_kobj_update_num(obj);

  free(vlist);
#endif
} /* off_create_plane() */




int main(int argc,char *argv[],char *envp[]) {
  int clobber=0;

  int    iter=10;

  nifti_image *maskimg=NULL;
  nifti_image *distimg=NULL;
  double max_dist=5;

  kobj *obj=NULL;

  const char *in_mask=NULL;
  const char *out_obj=NULL;
  const char *out_left=NULL;
  const char *out_right=NULL;

  const char *fcname="falcon_midsag";

  int i,n;
  double elen=2.0;
  double border=1.0;
  double steplen=0.5;
  niikpt  ctr;
  int flag_remesh=0;
  int flag_bbox=1;
  int flag_debug=0;
  int label=-1;

  struct option long_options[] = {
    {"clobber", no_argument,  &clobber, 1},
    {"debug",   no_argument,  &flag_debug, 1},
    {"noremesh", no_argument, &flag_remesh, 0},
    {"version", no_argument, 0, 'v'},
    {"help",    no_argument, 0, 'h'},
    {"usage",   no_argument, 0, 'U'},

    {"elen",required_argument,0,     'e'},
    {"border",required_argument,0,   'b'},
    {"out",required_argument,0,      'o'},
    {"left",required_argument,0,     'L'},
    {"right",required_argument,0,    'R'},
    {"iter",required_argument,  0,   'I'},
    {"step",required_argument, 0,    'V'},
    {"label",required_argument,  0,  'l'},

    {0, 0, 0, 0}
  };
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vhue:o:l:L:R:I:V:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      prog_version();
      return 0;
    case 'b':
      border=atof(optarg);
      break;
    case 'e':
      elen=atof(optarg);
      break;
    case 'o':
      out_obj=optarg;
      break;
    case 'I':
      iter=atoi(optarg);
      break;
    case 'V':
      steplen=atof(optarg);
      break;
    case 'l':
      label=atoi(optarg);
      break;
    case 'R':
      out_right=optarg;
      break;
    case 'L':
      out_left=optarg;
      break;
    case 'h':
    case 'U':
    case '?':
    default:
      prog_usage ();
      return 1;
    }
  }

  if((argc - optind)<1) {
    prog_usage();
    return 1;
  }

  in_mask = argv[optind];

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);

  if (!clobber && !access (out_obj, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_obj);
    return 1;
  }

  if( (maskimg=niik_image_read(in_mask))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,in_mask);
    exit(1);
  }

  if(label!=-1) {
    double* dimg =NULL;
    int i;

    if((dimg=niik_image_get_voxels_as_double_vector(maskimg))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_double_vector\n",fcname);
      exit(1);
    }

    for(i=0; i<maskimg->nvox; i++) {
      if(fabs(dimg[i]-label)<0.5)
        dimg[i]=1.0;
      else dimg[i]=0.0;
    }
    if(!niik_image_set_voxels_from_double_vector(maskimg,dimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_set_voxels_from_double_vector\n",fcname);
      exit(1);
    }
    free(dimg);

    if(flag_debug) {
      fprintf(stdout,"[%s] [DEBUG] Writing itnermediate file to tmp_label.mnc\n",fcname);
      niik_image_write("tmp_label.mnc",maskimg);
    }
  }

#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  if(niik_check_double_problem(steplen)) {
    steplen = 0.5;
  }

  if(iter<0) iter=20;

  if(niik_check_double_problem(elen)) elen=2.0;

  {
    kvert *v;
    fprintf(stdout,"[%s] creating a new planar object\n",fcname);

    NIIK_EXIT( ((obj=off_find_minimal_cut_plane(elen,maskimg))==NULL),fcname, "off_find_minimal_cut_plane",1);

    if(flag_debug) {
      fprintf(stdout,"[%s]  [DEBUG] writing initial object to tmp_init.off\n",fcname);
      NIIK_EXIT((!off_kobj_write_offply("tmp_init.off",obj,0)),fcname,"off_kobj_write_off",1);
    }

    off_kobj_remove_color(obj);
  }

  /**
   * convert to distance map
   */
  NIIK_EXIT( ((distimg=niik_image_distance_map(maskimg,max_dist))==NULL),fcname, "niik_image_distance_map",1);
  /*niik_image_write("tmp_dist.mnc",distimg);*/
  /*
    * run surface deformation
    */
  fprintf(stdout,"[%s] deform\n",fcname);
  NIIK_EXIT((!off_deform_minimal_cut_plane(distimg, obj, elen, iter, steplen, flag_bbox, flag_remesh, flag_debug )),fcname,"off_deform_minimal_cut_plane",1);

  /*
    * write output
    */

  if(out_obj) {
    fprintf(stdout,"[%s] writing output  %s\n",fcname,out_obj);
    NIIK_EXIT((!off_kobj_write_offply(out_obj,obj,0)),fcname,"off_kobj_write_off",1);
  }

  if(out_left || out_right) { /*split initial image into left and right parts*/
    nifti_image *left_img=NULL;
    nifti_image *right_img=NULL;
    left_img=niik_image_copy(maskimg); /* copies image */
    right_img=niik_image_copy(maskimg); /* copies image */
    off_update_kobj_face_normal(obj);
    NIIK_EXIT((!off_split_left_right(maskimg, obj,border, left_img,right_img  )),fcname,"off_deform_minimal_cut_plane",1);

    if(out_left)
      niik_image_write(out_left,left_img);
    if(out_right)
      niik_image_write(out_right,right_img);

    /*niik_image_write("test.mnc",maskimg);*/

    nifti_image_free( left_img );
    nifti_image_free( right_img );

  }

  distimg=niik_image_free(distimg);
  obj=off_kobj_free(obj);
  maskimg=niik_image_free(maskimg);
  niik_fc_display(fcname,0);
  free(timestamp);
  exit(0);
} /* main */

/*
kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/