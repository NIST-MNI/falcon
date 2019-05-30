/* FILENAME:     niik_resample_field_sph.c
 * DESCRIPTION:  resample field values on the surfaces using spherical coordinates
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include  <volume_io.h>

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input_field.txt/csv> <input_surface.off/ply> <reference.off/ply> <output_field.txt/csv> \n"
          "\t--clobber clobber output files\n"
          "\t--primitive - use primitive linear search algorithm O(n*m) , nearest neighbour interpolation\n"
          "\t--linear - use linear interpolation (default)\n"
          "\t--nearest - use nearest neighbor interpolation\n"
          ,
          name);
}

inline static double off_diff_radians(double r1, double r2)
/* add r1 and r2 but considers the limits like PI,
 * lim is half of the limit */
{
  double d;

  d=r1-r2;

  if(d<-NIIK_PI) d+=NIIK_PI2;
  if(d>NIIK_PI) d-=NIIK_PI2;

  return d;
}

/*primitive nearest-neighbour algorithm for resampling the surface*/
int resample_sph_primitive(niiktable *meas_in, kobj *obj, niiktable *meas_out, kobj *ref) {
  kvert *v;
  kvert *r;
  kvert *nn;
  double best_guess;

  for(r=ref->vert; r!=NULL; r=r->next) {
    int j;

    best_guess=1e10;
    /*head-on approach in finding a nearest neighbour*/
    for(v=obj->vert; v!=NULL; v=v->next) {
      double dist_psi=off_diff_radians(v->sph.psi,r->sph.psi);
      double dist_the=off_diff_radians(v->sph.the,r->sph.the);
      double dist2=dist_psi*dist_psi+dist_the*dist_the;

      if(dist2<best_guess) {
        nn=v;
        best_guess=dist2;
      }
    }

    for(j=0; j<meas_in->ncol; j++)
    {
      /*for now just take the best index*/  
      meas_out->col[j]->v[r->index-1] = meas_in->col[j]->v[nn->index-1];
    } 

  }
  return 1;
}

/*use compute barycentric coordinates in spherical space*/
/*ref: http://blackpawn.com/texts/pointinpoly/ */
/*returns TRUE if point is inside of the triangle*/
static inline int sph_compute_baricentric_coordinates(kvert *r,kvert *vert[],double *u,double *v) {
  /* Compute vectors*/
  niiksph v0 = niiksph_norm_sym(niiksph_sub(vert[1]->sph, vert[0]->sph));
  niiksph v1 = niiksph_norm_sym(niiksph_sub(vert[2]->sph, vert[0]->sph));
  niiksph v2 = niiksph_norm_sym(niiksph_sub(r->sph,        vert[0]->sph));

  /* Compute dot products*/
  double dot00 = niiksph_dot(v0, v0);
  double dot01 = niiksph_dot(v0, v1);
  double dot02 = niiksph_dot(v0, v2);
  double dot11 = niiksph_dot(v1, v1);
  double dot12 = niiksph_dot(v1, v2);

  /* Compute barycentric coordinates*/
  double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);

  *u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  *v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  /* Check if point is in triangle*/
  return ((*u) >= 0.0) && ((*v) >= 0.0) && ((*u) + (*v) < 1.0);

}

/*resample using quad-tree*/
int resample_sph_qtree(niiktable *meas_in,kobj *obj,niiktable *meas_out,kobj *ref,int linear_mode) {
  kvert *v;
  kvert *r;
  const char *fcname=__func__;

  bbox_sph *off_bb=off_bbox_sph_init(5,NIIK_PI2);
  if(!off_create_bbox_sph_from_kobj(off_bb,obj)) {
    fprintf(stderr,"[%s] Can't initialize qtree\n",fcname);
    return 1;
  }

  for(r=ref->vert; r!=NULL; r=r->next) {
    kface *best_face=NULL;
    kface *best_face2=NULL;
    int psi_min,psi_max,the_min,the_max;
    int i,j,k;
    double best_guess=1e10;

    psi_min = floor( (r->sph.psi - off_bb->origin.psi) / off_bb->delta - 0.5);
    psi_max = floor( (r->sph.psi - off_bb->origin.psi) / off_bb->delta + 0.5);
    the_min = floor( (r->sph.the - off_bb->origin.the) / off_bb->delta - 0.5);
    the_max = floor( (r->sph.the - off_bb->origin.the) / off_bb->delta + 0.5);

    /*here we can have around 0 wrap-around*/
    psi_min = NIIK_IMINMAX(psi_min-1, -1, off_bb->psi_dim);
    psi_max = NIIK_IMINMAX(psi_max+1, -1, off_bb->psi_dim);
    the_min = NIIK_IMINMAX(the_min-1, -1, off_bb->the_dim);
    the_max = NIIK_IMINMAX(the_max+1, -1, off_bb->the_dim);

    for(j=the_min; j<=the_max&& best_face==NULL ; j++) /**/
      for(i=psi_min; i<=psi_max; i++) {
        int idx=((i+off_bb->psi_dim)%off_bb->psi_dim)+
                ((j+off_bb->the_dim)%off_bb->the_dim)*off_bb->psi_dim;

        for(k=0; k<off_bb->ndata[idx]; k++) {
          kface *f=off_bb->data[idx][k];
          int t;
          double u,v;
          /* Check if point is in triangle*/
          if(sph_compute_baricentric_coordinates(r, f->vert, &u, &v)) {
            best_face=f;
            break;
          }
          /*second best guess - find closest face*/
          for(t=0; t<3; t++) {
            double dist_psi=off_diff_radians(f->vert[t]->sph.psi,r->sph.psi);
            double dist_the=off_diff_radians(f->vert[t]->sph.the,r->sph.the);
            double dist2=dist_psi*dist_psi+dist_the*dist_the;

            if(dist2<best_guess) {
              best_guess=dist2; 
              best_face2=f;
            }
          }
        }
      }

    /*we found the face that contains the point in question.*/
    if(best_face==NULL) {
      best_face=best_face2;
      fprintf(stderr,"[%s] Can't find the good match for vertex %i, using second best guess \n", fcname, r->index-1);
    }

    if(linear_mode) {
      double u,v;
      int c;
      sph_compute_baricentric_coordinates(r,best_face->vert,&u,&v);

      /*compute weighted sum*/
      for(c=0;c<meas_in->ncol;c++)
        meas_out->col[c]->v[r->index-1] = u*         meas_in->col[c]->v[ best_face->vert[1]->index-1]+
                                          v*         meas_in->col[c]->v[ best_face->vert[2]->index-1]+
                                          (1.0-u-v)* meas_in->col[c]->v[ best_face->vert[0]->index-1];
    } else { /*find nearest neighbour*/
      int i,c;
      double best_guess=1e10;
      double u,v;
      kvert *nn=NULL;
      for(i=0; i<3; i++) {
        double dist_psi=off_diff_radians(r->sph.psi,best_face->vert[i]->sph.psi);
        double dist_the=off_diff_radians(r->sph.the,best_face->vert[i]->sph.the);
        double dist2=dist_psi*dist_psi+dist_the*dist_the;

        if(dist2<best_guess) {
          nn=best_face->vert[i];
          best_guess=dist2;
        }
      }

      for(c=0;c<meas_in->ncol;c++)
        meas_out->col[c]->v[r->index-1] = meas_in->col[c]->v[nn->index-1];
    }
  }
  off_bb=off_bbox_sph_free(off_bb);
  return 1;
}


int main(int argc, char **argv) {
  const char *fcname="niik_resample_field_sph";
  int clobber=0;
  int verbose=0;
  int use_primitive=0;
  int use_linear=1;
  int c;
  int i;
  float fwhm=1;
  const char *in_off=NULL;
  const char *in_fld=NULL;
  const char *in_ref=NULL;
  const char *out_fld=NULL;
  int n;
  kobj *obj=NULL;
  kobj *ref=NULL;

  niiktable *meas_in;
  niiktable *meas_out;

  struct option long_options[] = {
    {"clobber",   no_argument, &clobber, 1},
    {"verbose",   no_argument, &verbose, 1},
    {"primitive", no_argument, &use_primitive, 1},
    {"linear",    no_argument, &use_linear, 1},
    {"nearest",    no_argument, &use_linear, 0},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<4) {
    show_usage(argv[0]);
    return 1;
  }

  in_fld     = argv[optind];
  in_off     = argv[optind+1];
  in_ref     = argv[optind+2];
  out_fld    = argv[optind+3];

  if (!clobber && !access (out_fld, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_fld);
    return 1;
  }
  niik_fc_display(fcname,1);

  NIIK_EXIT(((obj=off_kobj_read_offply(in_off))==NULL),fcname,"niik_kobj_read_off",1);
  if(!obj->spherecoo) {
    fprintf(stderr,"[%s] Surface doesn't have spherical coordinates mapping: %s\n",fcname,in_off);
    exit(1);
  }

  NIIK_EXIT(((meas_in=niiktable_read(in_fld))==NULL),fcname,"niiktable_read",1);
  if(meas_in->col[0]->num != obj->nvert) {
    fprintf(stderr,"[%s] Inconsistent number of measurement: %i , expected %i\n",fcname, meas_in->col[0]->num, obj->nvert);
    exit(1);
  }

  NIIK_EXIT(((ref=off_kobj_read_offply(in_ref))==NULL),fcname,"niik_kobj_read_off",1);
  if(!ref->spherecoo) {
    fprintf(stderr,"[%s] Surface doesn't have spherical coordinates mapping: %s\n",fcname,in_ref);
    exit(1);
  }

  meas_out = niiktable_init_ex(meas_in->ncol,ref->nvert, meas_in);

  if(use_primitive) {
    NIIK_EXIT((resample_sph_primitive(meas_in,obj,meas_out,ref)==0),fcname,"resample_sph_primitive",1);
  } else {
    NIIK_EXIT((resample_sph_qtree(meas_in,obj,meas_out,ref,use_linear)==0),fcname,"resample_sph_qtree",1);
  }

  NIIK_EXIT((niiktable_write(out_fld, meas_out)==0),fcname,"niiktable_write",1);

  niiktable_free(meas_out);
  niiktable_free(meas_in);
  obj=off_kobj_free(obj);
  ref=off_kobj_free(ref);
  
  niik_fc_display(fcname,0);
  return 0;
}


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
