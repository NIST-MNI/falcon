/* FILENAME:     falcon_ply.c
 * DESCRIPTION:  Writing ply file format
 * AUTHOR:       Vladimir Fonov
 *
 *
 * reference
 * http://www.danielgm.net/cc/doc/wiki/index.php5?title=FILE_I/O
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include "rply.h"

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stddef.h>

/*needed for ply file format, probably should be used everywhere else too*/
#include <locale.h>


#define PLY_CHECK(test)      if((test)==0) { fprintf(stdout,"[%s:%i:%s] PLY ERROR\n",__FILE__,__LINE__,__func__); return 0; }


int niik_off_write_ply(const char *fname,kobj *obj,const double *meas,int output_sph,int output_normal,int output_edge) {
  /*Legacy function*/
  const char *meas_names[]={"thickness"};
  const double *meas_vec[]={meas};
  int n_meas=meas?1:0;

  return off_kobj_write_ply_ex(fname,obj,output_normal,output_sph, output_edge, 1, n_meas, meas_names, meas_vec);
}

int off_kobj_write_ply(const char *fname, kobj *obj, int write_vertex_normal) {
  /*Legacy function*/
  return off_kobj_write_ply_ex(fname, obj, write_vertex_normal, 1, 1, 1, 0, NULL,NULL);
}

int off_kobj_write_ply_ex(const char *fname, kobj *obj,
       int output_normal,int output_sph,int output_edge, int output_color,
       int n_meas, const char **meas_name, const double **meas) {
  const char *fcname=__func__;
  FILE *fp=NULL;
  int i,n;
  int verbose=niik_verbose();
  kvert *v;
  kface *f;
  kedge *e;
  p_ply iply, oply;

  if(verbose>1)
    niik_fc_display(fcname,1);

  oply = ply_create(fname, PLY_LITTLE_ENDIAN, NULL, 0, NULL); /*PLY_ASCII*/
  if(!oply)
    return 0;
  
  /*define all the properties*/
  ply_add_element(oply,"vertex",obj->nvert); /*vertex*/

  ply_add_scalar_property(oply,"x",PLY_FLOAT64); /*maybe float32 will be enough*/
  ply_add_scalar_property(oply,"y",PLY_FLOAT64); /*maybe float32 will be enough*/
  ply_add_scalar_property(oply,"z",PLY_FLOAT64); /*maybe float32 will be enough*/


  if(obj->spherecoo && output_sph) {
    ply_add_scalar_property(oply,"psi",PLY_FLOAT64);
    ply_add_scalar_property(oply,"the",PLY_FLOAT64);
  }

  if(output_normal) {
    ply_add_scalar_property(oply,"Nx",PLY_FLOAT64);
    ply_add_scalar_property(oply,"Ny",PLY_FLOAT64);
    ply_add_scalar_property(oply,"Nz",PLY_FLOAT64);
  }

  /*dump additional scalar fields, vertex level, if needed*/
  for(i=0;i<n_meas;i++) {
    ply_add_scalar_property(oply, meas_name[i], PLY_FLOAT64);
  }

  ply_add_element(oply,"face",obj->nface); /*face*/
  ply_add_list_property(oply,"vertex_index", PLY_UINT8, PLY_UINT);

  if(obj->color && output_color ) {
    ply_add_scalar_property(oply,"red",  PLY_UINT8);
    ply_add_scalar_property(oply,"green",PLY_UINT8);
    ply_add_scalar_property(oply,"blue", PLY_UINT8);
  }

  if(obj->edge!=NULL && output_edge) {
    ply_add_element(oply, "edge", obj->nedge); /*edge*/
    ply_add_scalar_property(oply,"vertex1",PLY_UINT);
    ply_add_scalar_property(oply,"vertex2",PLY_UINT);

    ply_add_scalar_property(oply,"face1",PLY_UINT);
    ply_add_scalar_property(oply,"face2",PLY_UINT);
  }

  /*write comments (metadata)*/
  /*max comment length in one line is 1024*/
  for(i=0; i<obj->n_comments; i++)
  { 
    const int max_line_len=1000;
    int comment_len = strlen(obj->comment[i]);
    if(comment_len<max_line_len)
    {
      if (!ply_add_comment(oply, obj->comment[i])) 
      {
        fprintf(stderr,"ERROR: [%s] ply_add_comment\n",fcname);
        return 0;
      }
    } else {
      /*kludge to split comment into multiple ones*/
      /*TODO: split at white space*/
      const char *cmt=obj->comment[i];
      /*split into multiple lines*/
      while(comment_len>0)
      {
        char tmp[max_line_len+1];
        int cl;
        strncpy(tmp,cmt,max_line_len);
        tmp[max_line_len]='\0';
        cl=strlen(tmp);
        if (!ply_add_comment(oply, tmp)) {
          fprintf(stderr,"ERROR: [%s] ply_add_comment\n",fcname);
          return 0;
        }
        cmt+=cl;
        comment_len-=cl;
      }
    }
  }

  PLY_CHECK(ply_write_header(oply));


  /*add coordinates*/
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    PLY_CHECK(ply_write(oply,v->v.x));
    PLY_CHECK(ply_write(oply,v->v.y));
    PLY_CHECK(ply_write(oply,v->v.z));
    /*add spherical coordinates*/
    if(obj->spherecoo && output_sph) {
      PLY_CHECK(ply_write(oply,v->sph.psi));
      PLY_CHECK(ply_write(oply,v->sph.the));
    }
    /*add normals*/
    if(output_normal) {
      PLY_CHECK(ply_write(oply,v->normal.x));
      PLY_CHECK(ply_write(oply,v->normal.y));
      PLY_CHECK(ply_write(oply,v->normal.z));
    }
    for(n=0;n<n_meas;n++) {
      PLY_CHECK(ply_write(oply,meas[n][i]));
    }
  }

  /*pass face indexes*/
  for(f=obj->face,i=0; f!=NULL; f=f->next,i++) {
    PLY_CHECK(ply_write(oply,3.0));/*always 3 vertex*/
    PLY_CHECK(ply_write(oply,f->vert[0]->index-1));
    PLY_CHECK(ply_write(oply,f->vert[1]->index-1));
    PLY_CHECK(ply_write(oply,f->vert[2]->index-1));

    /*write face colours*/
    if(obj->color && output_color ) {
      if (f->color)
      {
        PLY_CHECK(ply_write(oply,f->color[0]*255));
        PLY_CHECK(ply_write(oply,f->color[1]*255));
        PLY_CHECK(ply_write(oply,f->color[2]*255));
      } else {
        /*output default color ?*/
        PLY_CHECK(ply_write(oply,0));
        PLY_CHECK(ply_write(oply,0));
        PLY_CHECK(ply_write(oply,0));
      }
    }
  }

  if(obj->edge!=NULL && output_edge ) {
    for(e=obj->edge,i=0; e!=NULL; e=e->next,i++) {
      /*adjacent vertex*/
      PLY_CHECK(ply_write(oply,e->endpts[0]->index-1));
      PLY_CHECK(ply_write(oply,e->endpts[1]->index-1));
      /*adjacent face*/
      PLY_CHECK(ply_write(oply,e->adjface[0]->index-1));
      PLY_CHECK(ply_write(oply,e->adjface[1]->index-1));
    }
  }
  return ply_close(oply);
}

typedef struct {
  double *buf;
  size_t cnt;
  size_t size;
} buf_ptr;

static buf_ptr* init_buf_ptr(size_t cnt) {
  buf_ptr* bbb=(buf_ptr*) calloc(sizeof(buf_ptr),1);
  bbb->buf=(double*) calloc(sizeof(double),cnt);
  bbb->cnt=0;
  bbb->size=cnt;
  return bbb;
}

static buf_ptr* free_buf_ptr(buf_ptr* bbb) {
  if(bbb->buf) free(bbb->buf);
  free(bbb);
  return NULL;
}

static buf_ptr* push_entry_buf(buf_ptr* bbb,double v) {
  /*assert ?*/
  bbb->buf[bbb->cnt]=v;
  bbb->cnt++;
  return bbb;
}

static int vertex_cb(p_ply_argument argument) {
  long idata;
  buf_ptr* pdata;
  ply_get_argument_user_data(argument, (void**)&pdata, &idata);
  push_entry_buf(pdata, ply_get_argument_value(argument));
  return 1;
}

static int face_cb(p_ply_argument argument) {
  long length, value_index;
  long idata;
  buf_ptr *pdata;
  ply_get_argument_property(argument, NULL, &length, &value_index);
  /*validate that it's a triangle*/
  if( length!=3 )
    fprintf(stderr,"ERROR : Found non-triangle\n");

  if(value_index<0) /*advancing to another element*/
    return 1;

  /*Ignore edges above 2*/
  if(value_index>3)
    return 1;

  ply_get_argument_user_data(argument, (void**)&pdata, &idata);
  push_entry_buf(pdata, ply_get_argument_value(argument));

  return 1;
}


static int edge_cb(p_ply_argument argument) {
  long idata;
  buf_ptr* pdata;
  ply_get_argument_user_data(argument, (void**)&pdata, &idata);
  push_entry_buf(pdata, ply_get_argument_value(argument));
  return 1;
}


kobj *off_kobj_read_ply(const char * fname) {
  return off_kobj_read_ply_ex(fname, NULL, NULL,NULL);
}

kobj *off_kobj_read_ply_ex(const char * fname, int *n_meas, char ***meas_name, double ***meas) {
  const char *fcname=__func__;
  int verbose=niik_verbose();
  long nvertices=0, ntriangles=0, nedges=0;
  p_ply ply;
  p_ply_element element = NULL;
  int have_face=0,have_vertex=0,have_edge=0;
  int have_x=0,have_y=0,have_z=0;
  int have_nx=0,have_ny=0,have_nz=0;
  int have_psi=0,have_the=0;
  int have_vertex_index=0;
  int have_vertex_indices=0;
  int have_f_red=0,have_f_green=0,have_f_blue=0;
  int have_edge_vertex1=0,have_edge_vertex2=0;
  int have_edge_face1=0,have_edge_face2=0;
  int have_custom_fields=0;
  char **custom_fields=NULL;
  const char *comment;
  kobj *obj=NULL;

  if(verbose>1)
    niik_fc_display(fcname,1);

  ply = ply_open(fname, NULL, 0, NULL);
  if (!ply) {
    fprintf(stderr,"ERROR %s: ply_open %s\n",fcname, fname);
    return NULL;
  }

  if (!ply_read_header(ply)) {
    fprintf(stderr,"ERROR %s: ply_read_header %s\n",fcname, fname);
    return NULL;
  }

  /*make sure this ply file have all what we need*/
  while ((element = ply_get_next_element(ply, element))) {
    p_ply_property property = NULL;
    long ninst = 0;
    const char *en;
    ply_get_element_info(element, &en, &ninst);
    if(!strcmp(en,"vertex")) {
      nvertices=ninst;
      have_vertex=1;
      while ((property = ply_get_next_property(element, property))) {
        const char *pn;
        e_ply_type pt, length_type, value_type;
        ply_get_property_info(property, &pn, &pt, &length_type, &value_type);
        if      (!strcmp(pn,"x")) have_x=1;
        else if (!strcmp(pn,"y")) have_y=1;
        else if (!strcmp(pn,"z")) have_z=1;
        else if (!strcmp(pn,"psi")) have_psi=1;
        else if (!strcmp(pn,"the")) have_the=1;
        else if (!strcmp(pn,"Nx")) have_nx=1;
        else if (!strcmp(pn,"Ny")) have_ny=1;
        else if (!strcmp(pn,"Nz")) have_ny=1;
        else {
          /*if caller wants custom fields*/
          if(n_meas!=NULL) {
            custom_fields=(char**)realloc(custom_fields,(have_custom_fields+1)*sizeof(char*));
            custom_fields[have_custom_fields]=strdup(pn);
            have_custom_fields++;
          }
        }
      }
    } else if(!strcmp(en,"face")) {
      ntriangles=ninst;
      have_face=1;
      while ((property = ply_get_next_property(element, property))) {
        const char *pn;
        e_ply_type pt, length_type, value_type;
        ply_get_property_info(property, &pn, &pt, &length_type, &value_type);
        if(!strcmp(pn,"vertex_index")) have_vertex_index=1;
        else if(!strcmp(pn,"vertex_indices")) have_vertex_indices=1;
        else if(!strcmp(pn,"red"))   have_f_red=1;
        else if(!strcmp(pn,"green")) have_f_green=1;
        else if(!strcmp(pn,"blue"))  have_f_blue=1;
      }
    } else if(!strcmp(en,"edge")) {
      nedges=ninst;
      have_edge=1;
      while ((property = ply_get_next_property(element, property))) {
        const char *pn;
        e_ply_type pt, length_type, value_type;
        ply_get_property_info(property, &pn, &pt, &length_type, &value_type);
        if(!strcmp(pn,"vertex1")) have_edge_vertex1=1;
        else if(!strcmp(pn,"vertex2")) have_edge_vertex2=1;
        else if(!strcmp(pn,"face1")) have_edge_face1=1;
        else if(!strcmp(pn,"face2")) have_edge_face2=1;
      }
    }
  }
  /*create an empty object*/
  obj=off_obj_init();

  /*read comments*/
  comment = NULL;

  while ((comment = ply_get_next_comment(ply, comment)))
  {
    off_kobj_add_comment(obj,comment);
  }

  /*Verify if this is "good" ply file that we can interpret*/
  if(!have_vertex)
    fprintf(stderr,"ERROR %s: need vertex data %s\n",fcname, fname);
  else if(!have_x)
    fprintf(stderr,"ERROR %s: need x data %s\n",fcname, fname);
  else if(!have_y)
    fprintf(stderr,"ERROR %s: need y data %s\n",fcname, fname);
  else if(!have_z)
    fprintf(stderr,"ERROR %s: need z data %s\n",fcname, fname);
  else if(!have_face)
    fprintf(stderr,"ERROR %s: need face data %s\n",fcname, fname);
  else if(!have_vertex_index && !have_vertex_indices)
    fprintf(stderr,"ERROR %s: need vertex index for face data %s\n",fcname, fname);
  else {
    int i;
    buf_ptr *buf_x=NULL,
            *buf_y=NULL,
            *buf_z=NULL;
    buf_ptr *buf_nx=NULL,
            *buf_ny=NULL,
            *buf_nz=NULL;
    buf_ptr *buf_psi=NULL,
            *buf_the=NULL;

    buf_ptr **buf_custom=NULL;

    buf_ptr *buf_v_idx=NULL;
    buf_ptr *buf_f_colour=NULL;
    buf_ptr *buf_e_v_idx=NULL;
    buf_ptr *buf_e_f_idx=NULL;
    /*to start, read all data into buffers, and then convert to off structure*/

    if(verbose>0)
      fprintf(stdout,"Reading: %ld vertex, %ld faces, %ld edges\n",nvertices,ntriangles,nedges);

    buf_x=init_buf_ptr(nvertices);
    buf_y=init_buf_ptr(nvertices);
    buf_z=init_buf_ptr(nvertices);

    ply_set_read_cb(ply, "vertex", "x", vertex_cb, buf_x, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, buf_y, 0);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, buf_z, 0);

    if(have_nx && have_ny && have_nz) {
      buf_nx=init_buf_ptr(nvertices);
      buf_ny=init_buf_ptr(nvertices);
      buf_nz=init_buf_ptr(nvertices);

      ply_set_read_cb(ply, "vertex", "Nx", vertex_cb, buf_nx, 0);
      ply_set_read_cb(ply, "vertex", "Ny", vertex_cb, buf_ny, 0);
      ply_set_read_cb(ply, "vertex", "Nz", vertex_cb, buf_nz, 0);
    }

    if(have_psi && have_the) {
      buf_psi=init_buf_ptr(nvertices);
      buf_the=init_buf_ptr(nvertices);

      ply_set_read_cb(ply, "vertex", "psi", vertex_cb, buf_psi, 0);
      ply_set_read_cb(ply, "vertex", "the", vertex_cb, buf_the, 0);
    }

    if(have_custom_fields>0 && n_meas!=NULL) {
      buf_custom=(buf_ptr **)calloc(sizeof(buf_ptr *),have_custom_fields);
      for(i=0;i<have_custom_fields;i++) {
        buf_custom[i]=init_buf_ptr(nvertices);
        ply_set_read_cb(ply, "vertex", custom_fields[i], vertex_cb, buf_custom[i], 0);
      }
    }

    buf_v_idx = init_buf_ptr(ntriangles*3);
    if(have_vertex_index)
      ply_set_read_cb(ply, "face", "vertex_index", face_cb, buf_v_idx, 0);
    else /*have_vertex_indices*/
      ply_set_read_cb(ply, "face", "vertex_indices", face_cb, buf_v_idx, 0);

    if(have_f_red && have_f_green && have_f_blue) {
      buf_f_colour = init_buf_ptr(ntriangles*3);
      ply_set_read_cb(ply, "face", "red",   vertex_cb, buf_f_colour, 0);
      ply_set_read_cb(ply, "face", "green", vertex_cb, buf_f_colour, 1);
      ply_set_read_cb(ply, "face", "blue",  vertex_cb, buf_f_colour, 2);
    }

    /*init edges*/
    if(have_edge && have_edge_vertex1 && have_edge_vertex2 && have_edge_face1 && have_edge_face2) {
      buf_e_v_idx = init_buf_ptr(nedges*2);
      buf_e_f_idx = init_buf_ptr(nedges*2);

      ply_set_read_cb(ply, "edge", "vertex1", edge_cb, buf_e_v_idx, 0);
      ply_set_read_cb(ply, "edge", "vertex2", edge_cb, buf_e_v_idx, 1);

      ply_set_read_cb(ply, "edge", "face1", edge_cb, buf_e_f_idx, 0);
      ply_set_read_cb(ply, "edge", "face2", edge_cb, buf_e_f_idx, 1);
    }

    if (!ply_read(ply)) {
      ply_close(ply);
      fprintf(stderr,"ERROR %s: ply_read %s\n", fcname, fname);
      off_kobj_free(obj);
      return NULL;
    }

    /*validate*/
    if(buf_x->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of x coords:%d(%d)\n",fcname,(int)buf_x->cnt,(int)nvertices );
    else if(buf_y->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of y coords:%d(%d)\n",fcname,(int)buf_y->cnt,(int)nvertices );
    else if(buf_z->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of z coords:%d(%d)\n",fcname,(int)buf_z->cnt,(int)nvertices );

    else if(buf_nx && buf_nx->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of Nx coords:%d(%d)\n",fcname,(int)buf_nx->cnt,(int)nvertices );
    else if(buf_ny && buf_ny->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of Ny coords:%d(%d)\n",fcname,(int)buf_ny->cnt,(int)nvertices );
    else if(buf_nz && buf_nz->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of Nz coords:%d(%d)\n",fcname,(int)buf_nz->cnt,(int)nvertices );

    else if(buf_psi && buf_psi->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of psi coords:%d(%d)\n",fcname,(int)buf_psi->cnt,(int)nvertices );
    else if(buf_the && buf_the->cnt!=nvertices)
      fprintf(stderr,"ERROR %s: unexpected number of the coords:%d(%d)\n",fcname,(int)buf_the->cnt,(int)nvertices );

    else if(buf_v_idx->cnt!=(ntriangles*3))
      fprintf(stderr,"ERROR %s: unexpected number corners of faces:%d(%d)\n",fcname,(int)buf_v_idx->cnt,(int)(ntriangles*3) );
    else if(buf_f_colour && buf_f_colour->cnt!=(ntriangles*3))
      fprintf(stderr,"ERROR %s: unexpected number of colour info :%d(%d)\n",fcname,(int)buf_f_colour->cnt,(int)(ntriangles*3) );

    else if(buf_e_v_idx && buf_e_v_idx->cnt!=(nedges*2))
      fprintf(stderr,"ERROR %s: unexpected number corners of edges:%d(%d)\n",fcname,(int)buf_e_v_idx->cnt,(int)(nedges*2) );
    else if(buf_e_f_idx && buf_e_f_idx->cnt!=(nedges*2))
      fprintf(stderr,"ERROR %s: unexpected number sides of edges:%d(%d)\n",fcname,(int)buf_e_f_idx->cnt,(int)(nedges*2) );
    else {
      /*convert to OFF format*/
      kvert *v=NULL,*vc=NULL,**vlist;
      kface *f=NULL,*fc=NULL,**flist;
      kedge *e=NULL,*ce=NULL;
      int n;

      obj->fname = strdup(fname);
      obj->nvert = nvertices;
      obj->nface = ntriangles;
      obj->nedge = nedges;
      niikpt p,norm,color;

      if(buf_psi && buf_the) obj->spherecoo=1;
      if(buf_f_colour) obj->color=1;


      vlist = off_vert_list(obj->nvert);
      flist = off_face_list(obj->nface);
      /*Init vertices*/
      for(n=0; n<obj->nvert; n++) {
        p.x=buf_x->buf[n];
        p.y=buf_y->buf[n];
        p.z=buf_z->buf[n];
        p.w=0.0;

        /* assign value and put in the list */
        v=off_vert_init();
        v->v=p;
        v->index=n+1;

        if(buf_nx && buf_ny && buf_nz) {
          norm.x=buf_nx->buf[n];
          norm.y=buf_ny->buf[n];
          norm.z=buf_nz->buf[n];
          norm.w=0.0;
          v->normal = norm;
        }

        if(buf_psi && buf_the) {
          v->sph.the=buf_the->buf[n];
          v->sph.psi=buf_psi->buf[n];
        }

        /*TODO: v->color*/
        if(obj->vert==NULL) {
          obj->vert=v;
          vc=v;
        } else {
          vc->next=v;
          v->prev=vc;
          vc=v;
        }
        vlist[n]=v;
      }

      /*init faces*/
      for(n=0; n < obj->nface; n++) {
        f = off_face_init();
        /*Debug*/
        if((size_t)buf_v_idx->buf[n*3  ]>=obj->nvert) abort();
        if((size_t)buf_v_idx->buf[n*3+1]>=obj->nvert) abort();
        if((size_t)buf_v_idx->buf[n*3+2]>=obj->nvert) abort();

        f -> vert[0] = vlist[(size_t)buf_v_idx->buf[n*3  ]];
        f -> vert[1] = vlist[(size_t)buf_v_idx->buf[n*3+1]];
        f -> vert[2] = vlist[(size_t)buf_v_idx->buf[n*3+2]];

        f -> index = n+1;

        if(buf_f_colour) {
          f->color = (double *)calloc(4,sizeof(double));
          f->color[0] = buf_f_colour->buf[n*3  ]/255.0;
          f->color[1] = buf_f_colour->buf[n*3+1]/255.0;
          f->color[2] = buf_f_colour->buf[n*3+2]/255.0;
        }
        if(obj->face==NULL) {
          obj->face=f;
          fc = f;
        } else {
          fc->next = f;
          f->prev = fc;
          fc = f;
        }
        flist[n] = f;
      }

      /*init edges*/
      if(buf_e_v_idx && buf_e_f_idx ) {
        for(n=0; n<nedges; n++) {

          double *vidx=&(buf_e_v_idx->buf[n*2]);
          double *fidx=&(buf_e_f_idx->buf[n*2]);

          e = off_edge_init();

          e->endpts[0] = vlist[ (size_t)vidx[0] ];
          e->endpts[1] = vlist[ (size_t)vidx[1] ];

          e->adjface[0] = flist[ (size_t)fidx[0] ];
          e->adjface[1] = flist[ (size_t)fidx[1] ];
          e->index = n+1;

          if(obj->edge==NULL) {
            ce = obj->edge = e;
          } else {
            ce->next = e;
            e->prev = ce;
            ce = e;
          }
        }
      }

      /*updated neighbours*/
      for(e=obj->edge; e!=NULL; e=e->next) {
        int m;
        for(m=0; m<3; m++) {
          if(e->adjface[0]->edge[m] == NULL) {
            e->adjface[0]->edge[m] = e;
            break;
          }
        }
        for(m=0; m<3; m++) {
          if(e->adjface[1]->edge[m] == NULL) {
            e->adjface[1]->edge[m] = e;
            break;
          }
        }
      }

      for(e=obj->edge; e!=NULL; e=e->next) {
        if(!e->endpts[0]->nei)  {
          if(!off_kobj_update_vertex_nei( e->endpts[0], e)) {
            fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
            off_kobj_free(obj);
            return NULL;
          }
        }
        if(!e->endpts[1]->nei)  {
          if(!off_kobj_update_vertex_nei( e->endpts[1], e)) {
            fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
            off_kobj_free(obj);
            return NULL;
          }
        }
      }
      off_kobj_update_num(obj);

      free(vlist);
      free(flist);

      /*import custom measurements*/
      if(n_meas!=NULL )
      {
        *n_meas=have_custom_fields;
        if(have_custom_fields>0) {
          /*outp*/
          *meas_name=custom_fields;
          *meas=(double**)calloc(sizeof(double*),have_custom_fields);
          for(i=0;i<have_custom_fields;i++)
          {
            /*validate that we got enough entries first*/
            if(buf_custom[i]->cnt!=nvertices)
              fprintf(stderr,"ERROR %s: unexpected number of the entries for field %s:%d(%d)\n",fcname, custom_fields[i], (int)buf_custom[i]->cnt,(int)nvertices );
            (*meas)[i] = buf_custom[i]->buf; /*output buffer take memory ownership from here*/
          }
          /*free some memory*/
          free(buf_custom);
        } else{
          *meas_name=NULL;
          *meas=NULL;
        }
      } 
    }

    free_buf_ptr(buf_x);
    free_buf_ptr(buf_y);
    free_buf_ptr(buf_z);

    if(buf_nx) free_buf_ptr(buf_nx);
    if(buf_ny) free_buf_ptr(buf_ny);
    if(buf_nz) free_buf_ptr(buf_nz);
    if(buf_psi) free_buf_ptr(buf_psi);
    if(buf_the) free_buf_ptr(buf_the);
    if(buf_f_colour) free_buf_ptr(buf_f_colour);

    free_buf_ptr(buf_v_idx);
    if(buf_e_v_idx) free_buf_ptr(buf_e_v_idx);
    if(buf_e_f_idx) free_buf_ptr(buf_e_f_idx);

    if(verbose>1)
      niik_fc_display(fcname,0);

    ply_close(ply);
    return obj;
  }
  off_kobj_free(obj);
  ply_close(ply);

  if(verbose>1)
    niik_fc_display(fcname,0);

  return NULL;
}

kobj *off_kobj_read_offply_ex(const char * fname, int *n_meas, char ***meas_name, double ***meas)
{
  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7) ||
      !strncmp(fname+(strlen(fname)-4),".off",4)
    ) {
    *n_meas=0; /*OFF file can't store anything else*/
    return off_kobj_read_off(fname);
  } else {
    return off_kobj_read_ply_ex(fname,n_meas,meas_name,meas);
  }
}

kobj *off_kobj_read_offply(const char * fname) {
  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7) ||
      !strncmp(fname+(strlen(fname)-4),".off",4)
    ) {
    return off_kobj_read_off(fname);
  } else {
    return off_kobj_read_ply(fname);
  }
}

int off_kobj_write_offply(const char * fname,kobj * obj,int write_vertex_normal) {
  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7) ||
      !strncmp(fname+(strlen(fname)-4),".off",4)
    ) {
    return off_kobj_write_off(fname, obj, write_vertex_normal);
  } else {
    return off_kobj_write_ply(fname, obj, write_vertex_normal);
  }
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/


