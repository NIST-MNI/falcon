/* Filename:     nifti1_kunio_off.c
 * Description:  object functions
 * Author:       Kunio Nakamura
 * Date:         March 1, 2012
 *
 * Reference:
 *
 * -O'Rourke, Joseph. Computational Geomtry in C. second edition. Cambridge Academic Press 1998.
 * -Geomview: http://www.geomview.org/
 * -Object file format: http://www.geomview.org/docs/html/OFF.html
 * -Botsch Eurographic Symphosium on Geometry Processing 2004
 *
 *
 * geomview's off file
 * -basic operations (read/write/init/free)
 * -normal calculation
 * -color addition/removal
 *
 * other off functions are in other files
 * -nifti1_kunio_off_remesh.c        remesh functions (probably OK?)
 * -nifti1_kunio_off_shrinkwrap.c    basic shrink-wrap algorithm (incomplete)
 * -nifti1_kunio_off_obj2img.c       previously called off2bounds (-simple version)
 * -nifti1_kunio_off_bbox.c          bounding box system (probably OK?)
 * -nifti1_kunio_off_curberille.c    surface detection from a mask image
 * -nifti1_kunio_off_                ....
 *
 *
 *
 */
#include "falcon.h"
#include "falcon_surfaces.h"

#define USE_ASCII_SPHERICAL_COORDINATE_OFF 0


/* constructors for kobj / kvert / kface / kedge */
kobj *off_obj_init() {
  kobj *obj;
  obj = (kobj *)calloc(1,sizeof(kobj));
  obj->fname=NULL;
  obj->face=NULL;
  obj->edge=NULL;
  obj->vert=NULL;
  obj->color=0;
  obj->spherecoo=0;
  obj->nface=obj->nedge=obj->nvert=0;
  return obj;
}

kvert *off_vert_init() {
  kvert *v;
  v = (kvert *)calloc(1,sizeof(kvert));
  v->v=niikpt_zero();
  v->index=0;
  v->nei=0;
  v->neivert=NULL;
  v->neiface=NULL;
  v->neiedge=NULL;
  v->sph=niiksph_zero();
  v->prev=v->next=NULL;
  return v;
}

kvert *off_vert_init_with_v(double x,double y,double z) {
  kvert *v;
  v = off_vert_init();
  v->v.x = x;
  v->v.y = y;
  v->v.z = z;
  return v;
}

kface *off_face_init() {
  kface *f;
  f = (kface *)calloc(1,sizeof(kface));
  f->vert[0]=f->vert[1]=f->vert[2]=NULL;
  f->edge[0]=f->edge[1]=f->edge[2]=NULL;
  f->index=0;
  f->color=NULL;
  f->prev=f->next=NULL;
  return f;
}

kedge *off_edge_init() {
  kedge *e;
  e = (kedge *)calloc(1,sizeof(kedge));
  e->adjface[0]=e->adjface[1]=NULL;
  e->endpts[0]=e->endpts[1]=NULL;
  e->index=0;
  e->prev=e->next=NULL;
  return e;
}

/* destructors for kobj / kvert / kface / kedge */
void off_kface_free(kface *f) {
  if(f==NULL) return;
  if(f->color!=NULL) free(f->color);
  free(f);
}

void off_kedge_free(kedge *e) {
  if(e!=NULL) free(e);
}

void off_kvert_free(kvert *v) {
  if(v==NULL) return;
  if(v->neiface!=NULL) free(v->neiface);
  if(v->neiedge!=NULL) free(v->neiedge);
  if(v->neivert!=NULL) free(v->neivert);
  free(v);
}


kobj *off_kobj_copy(kobj *src) {
  int n;
  kobj* no;
  kface *f,*f2;
  kedge *e,*e2;
  kvert *v,*v2;

  kface **flist;
  kedge **elist;
  kvert **vlist;

  if(src==NULL) return NULL;

  no = off_obj_init();
  if(no==NULL) return NULL;

  if(src->fname!=NULL) no->fname=strdup(src->fname);
  if(src->n_comments>0) 
  {
    int i;
    no->n_comments=src->n_comments;
    no->comment=(char**)calloc(no->n_comments,sizeof(char*));
    for(i=0;i<src->n_comments;i++)
      no->comment[i]=strdup(src->comment[i]);
  }
  no->color     = src->color;
  no->spherecoo = src->spherecoo;

  flist=(kface**)calloc(src->nface,sizeof(kface*));
  vlist=(kvert**)calloc(src->nvert,sizeof(kvert*));
  elist=(kedge**)calloc(src->nedge,sizeof(kedge*));

  if(src->vert!=NULL) {
    kvert *vc;
    for(v=src->vert; v!=NULL; v=v->next) {
      v2=off_vert_init();
      v2->v      = v->v;
      v2->index  = v->index;
      v2->normal = v->normal;
      v2->sph    = v->sph;
      v2->idata  = v->idata;
      if(v->color)
      {
        v2->color=(double*)calloc(4,sizeof(double));
        memcpy(v2->color,v->color,sizeof(double)*4);
      }
      vlist[v2->index-1]=v2;

      /*TODO: v->color*/
      if(no->vert==NULL) {
        no->vert=v2;
        vc=v2;
      } else {
        vc->next=v2;
        v2->prev=vc;
        vc=v2;
      }
    }
  }

  if(src->face!=NULL) {
    kface *fc;
    for(f=src->face; f!=NULL; f=f->next) {
      f2=off_face_init();
      f2->index=f->index;
      f2->normal=f->normal;
      f2->pmax=f->pmax;
      f2->pmin=f->pmin;
      if(f->color)
      {
        f2->color=(double*)calloc(4,sizeof(double));
        memcpy(f2->color,f->color,sizeof(double)*4);
      }
      flist[f2->index-1]=f2;

      if(no->face==NULL) {
        no->face=f2;
        fc=f2;
      } else {
        fc->next=f2;
        f2->prev=fc;
        fc=f2;
      }
    }
  }
  if(src->edge!=NULL) {
    kedge *ec;
    for(e=src->edge; e!=NULL; e=e->next) {
      e2=off_edge_init();
      e2->index=e->index;
      elist[e2->index-1]=e2;

      if(no->edge==NULL) {
        no->edge=e2;
        ec=e2;
      } else {
        ec->next=e2;
        e2->prev=ec;
        ec=e2;
      }
    }
  }
  no->nface = src->nface;
  no->nedge = src->nedge;
  no->nvert = src->nvert;
  /*regenerate all links and neighbours*/

  for(v=src->vert,v2=no->vert; v!=NULL; v=v->next,v2=v2->next) {
    v2->nei     = v->nei;
    v2->neivert = off_vert_list(v2->nei);
    v2->neiface = off_face_list(v2->nei);
    v2->neiedge = off_edge_list(v2->nei);

    for(n=0; n<v2->nei; n++) {
      v2->neivert[n]=vlist[v->neivert[n]->index-1];
      v2->neiface[n]=flist[v->neiface[n]->index-1];
      v2->neiedge[n]=elist[v->neiedge[n]->index-1];
    }
  }

  for(e=src->edge,e2=no->edge; e!=NULL; e=e->next,e2=e2->next) {
    for(n=0;n<2;n++)
    {
      e2->endpts[n] = vlist[e->endpts[n]->index-1];
      e2->adjface[n] = flist[e->adjface[n]->index-1];
    }
  }

  for(f=src->face,f2=no->face; f!=NULL; f=f->next,f2=f2->next) {
    for(n=0;n<3;n++)
    {
      f2->edge[n] = elist[f->edge[n]->index-1];
      f2->vert[n] = vlist[f->vert[n]->index-1];
    }
  }

  free(flist);free(vlist);free(elist);
  return no;
}

kobj *off_kobj_free(kobj *obj) {
  kface *f,*f2;
  kedge *e,*e2;
  kvert *v,*v2;
  if(obj==NULL) return NULL;
  if(obj->fname!=NULL) free(obj->fname);
  if(obj->face!=NULL) {
    for(f=obj->face; f!=NULL; f=f2) {
      f2=f->next;
      off_kface_free(f);
    }
  }
  if(obj->edge!=NULL) {
    for(e=obj->edge; e!=NULL; e=e2) {
      e2=e->next;
      off_kedge_free(e);
    }
  }
  if(obj->vert!=NULL) {
    for(v=obj->vert; v!=NULL; v=v2) {
      v2=v->next;
      off_kvert_free(v);
    }
  }
  if(obj->n_comments>0) 
  {
    int i;
    for(i=0;i<obj->n_comments;i++)
      free(obj->comment[i]);
  }
  if(obj->comment)
    free(obj->comment);

  free(obj);
  return NULL;
}


/*****************************************
 * check linked list
 ****************************************/

int off_kvert_test_link(kvert *vert) {
  kvert *v,*v2;
  int n;
  if(vert->next==NULL) {
    fprintf(stdout,"vert link test: only 1 item in the list\n");
    return 0;
  }
  if(vert->prev!=NULL) {
    fprintf(stdout,"vert link test: first item has a prev connection\n");
    return 0;
  }
  for(v=vert->next,v2=vert,n=2; v!=NULL; v=v->next,n++) {
    if(v->prev!=v2) {
      fprintf(stdout,"vert link test: vertex %i's previous is broken\n",n);
      fprintf(stdout,"  this vertex = %i\n",v->index);
      fprintf(stdout,"  prev vertex = %i\n",v->prev->index);
      fprintf(stdout,"  prev vertex = %i (actual previous)\n",v2->index);
      return 0;
    }
    v2=v;
  }
  return 1;
}

int off_kface_test_link(kface *face) {
  kface *f,*f2;
  int n;
  if(face->next==NULL) {
    fprintf(stdout,"face link test: only 1 item in the list\n");
    return 0;
  }
  if(face->prev!=NULL) {
    fprintf(stdout,"face link test: first item has a prev connection\n");
    return 0;
  }
  for(f=face->next,f2=face,n=2; f!=NULL; f=f->next,n++) {
    if(f->prev!=f2) {
      fprintf(stdout,"face link test: face %i's previous is broken\n",n);
      fprintf(stdout,"  this face = %i\n",f->index);
      fprintf(stdout,"  prev face = %i\n",f->prev->index);
      fprintf(stdout,"  prev face = %i (actual previous)\n",f2->index);
      return 0;
    }
    f2=f;
  }
  return 1;
}


int off_kedge_test_link(kedge *edge) {
  kedge *e,*e2;
  int n;
  if(edge->next==NULL) {
    fprintf(stdout,"edge link test: only 1 item in the list\n");
    return 0;
  }
  if(edge->prev!=NULL) {
    fprintf(stdout,"edge link test: first item has a prev connection\n");
    return 0;
  }
  for(e=edge->next,e2=edge,n=2; e!=NULL; e=e->next,n++) {
    if(e->prev!=e2) {
      fprintf(stdout,"edge link test: edge %i's previous is broken\n",n);
      fprintf(stdout,"  this edge = %i\n",e->index);
      fprintf(stdout,"  prev edge = %i\n",e->prev->index);
      fprintf(stdout,"  prev edge = %i (actual previous)\n",e2->index);
      return 0;
    }
    e2=e;
  }
  return 1;
}


/*****************************************
 * creates list of (kvert *), (kface *), (kedge *)
 ****************************************/
kvert **off_vert_list(int num) {
  kvert **v;
  v = (kvert **)calloc(num,sizeof(kvert *));
  return v;
}

kface **off_face_list(int num) {
  kface **v;
  v = (kface **)calloc(num,sizeof(kface *));
  return v;
}

kedge **off_edge_list(int num) {
  kedge **v;
  v = (kedge **)calloc(num,sizeof(kedge *));
  return v;
}

kvert ***off_vert_matrix(int num1, int num2) {
  kvert ***out;
  int n;
  if((out = (kvert ***)calloc(num1,sizeof(kvert **)))==NULL) {
    fprintf(stderr,"ERROR: calloc num1, %i\n",num1);
    return NULL;
  }
  for(n=0; n<num1; n++) {
    if((out[n] = (kvert **)calloc(num2,sizeof(kvert *)))==NULL) {
      fprintf(stderr,"ERROR: calloc (%i)\n",n);
      return NULL;
    }
  }
  return out;
}


/*****************************************
 * insert / remove from list
 ****************************************/
void off_kvert_remove(kvert *v,kobj *obj) {
  if(obj->vert == v) { /* should be same as v->prev == NULL */
    obj->vert = obj->vert->next;
    obj->vert->prev = NULL;
    return;
  } else if(v->next==NULL) {
    v->prev->next = NULL;
    return;
  }
  v->prev->next = v->next;
  v->next->prev = v->prev;
  return;
}

void off_kface_remove(kface *f,kobj *obj) {
  if(obj->face == f) { /* should be same as v->prev == NULL */
    obj->face = obj->face->next;
    obj->face->prev = NULL;
    return;
  } else if(f->prev==NULL) {
    fprintf(stderr,"ERROR: link is broken (off_kface_remove %i)\n",f->index);
    exit(0);
  } else if(f->next==NULL) {
    f->prev->next = NULL;
    return;
  }
  /*  fprintf(stdout," *** f       = %i\n",f->index);
  fprintf(stdout," *** f->prev = %i\n",f->prev->index);
  fprintf(stdout," *** f->next = %i\n",f->next->index);
  fprintf(stdout," *** f->prev->next = f->next\n");*/
  f->prev->next = f->next;
  /*  fprintf(stdout," *** f->next->prev = f->prev\n");*/
  f->next->prev = f->prev;
  /*  fprintf(stdout," *** done\n");*/
  return;
}

void off_kedge_remove(kedge *e,kobj *obj) {
  if(obj->edge == e) { /* should be same as v->prev == NULL */
    obj->edge = obj->edge->next;
    obj->edge->prev = NULL;
    return;
  } else if(e->prev==NULL) {
    fprintf(stderr,"ERROR: link is broken (off_kedge_remove %i)\n",e->index);
    exit(0);
  } else if(e->next==NULL) {
    e->prev->next = NULL;
    return;
  }
  e->prev->next = e->next;
  e->next->prev = e->prev;
  return;
}

/*************************************************
 * write off file
 *
 * -very important function
 * -for now, write the basics: coordinates, triangle index and triangle color
 *
 *
 * see also
 * kobj *off_kobj_read_off(char * fname);
 *
 **************************************************/
int off_kobj_write_off(const char *fname,kobj *obj,int write_vertex_normal) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  kvert *v;
  kface *f;
  kedge *e;
  int n,verbose=niik_verbose();
  char sname[4096];
  const char *fcname="off_kobj_write_off";

  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return 0;
  }
  /* re-index */
  if(verbose>0) fprintf(stderr,"[%s] re-index\n",fcname);
  if(obj->vert!=NULL) {
    for(v=obj->vert,n=1; v!=NULL; v=v->next) v->index=n++;
    obj->nvert=n-1;
  }
  if(obj->face!=NULL) {
    for(f=obj->face,n=1; f!=NULL; f=f->next) f->index=n++;
    obj->nface=n-1;
  }
  if(obj->edge!=NULL) {
    for(e=obj->edge,n=1; e!=NULL; e=e->next) e->index=n++;
    obj->nedge=n-1;
  }

  if(verbose>0) fprintf(stderr,"[%s] writing %s\n",fcname,fname);
  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7)) {
    if((gp = gzopen(fname,"w"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return 0;
    }
  } else if((fp=fopen(fname,"w"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s \n",fname);
    return 0;
  }
  if(write_vertex_normal) {
    if(fp!=NULL) fprintf(fp,"N ");
    if(gp!=NULL) gzprintf(gp,"N ");
  }
  if(fp!=NULL) {
    if(verbose>0) fprintf(stderr,"[%s] fp type\n",fcname);
    fprintf(fp,"OFF\n");
    fprintf(fp,"%i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }
  if(gp!=NULL) {
    gzprintf(gp,"OFF\n");
    gzprintf(gp,"%i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }

  if(verbose>0) fprintf(stderr,"[%s] writing vert\n",fcname);
  if(fp!=NULL) {
    for(v=obj->vert; v!=NULL; v=v->next) {
      fprintf(fp,"%15.12lf %15.12lf %15.12lf ",v->v.x,v->v.y,v->v.z);
      if(write_vertex_normal) {
        fprintf(fp,"%8.6f %8.6f %8.6f ",v->normal.x,v->normal.y,v->normal.z);
      }
      fprintf(fp,"\n");
    }
    if(verbose>0) fprintf(stderr,"[%s] writing face\n",fcname);
    for(f=obj->face; f!=NULL; f=f->next) {
      fprintf(fp,"3 %i %i %i ",f->vert[0]->index-1,f->vert[1]->index-1,f->vert[2]->index-1);
      if(obj->color) {
        fprintf(fp,"%8.6f %8.6f %8.6f %8.6f ",f->color[0],f->color[1],f->color[2],f->color[3]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  } /* FILE pointer */
  if(gp!=NULL) {
    for(v=obj->vert; v!=NULL; v=v->next) {
      gzprintf(gp,"%15.8lf %15.8lf %15.8lf ",v->v.x,v->v.y,v->v.z);
      if(write_vertex_normal) {
        gzprintf(gp,"%8.6f %8.6f %8.6f ",v->normal.x,v->normal.y,v->normal.z);
      }
      gzprintf(gp,"\n");
    }
    if(verbose>0) fprintf(stderr,"[%s] writing face\n",fcname);
    for(f=obj->face; f!=NULL; f=f->next) {
      gzprintf(gp,"3 %i %i %i ",f->vert[0]->index-1,f->vert[1]->index-1,f->vert[2]->index-1);
      if(obj->color) {
        gzprintf(gp,"%8.6f %8.6f %8.6f %8.6f ",f->color[0],f->color[1],f->color[2],f->color[3]);
      }
      gzprintf(gp,"\n");
    }
    gzclose(gp);
  } /* gzFile pointer */


  if(verbose>0) fprintf(stderr,"[%s] writing edge\n",fcname);
  if(obj->edge!=NULL) {
    if(!strncmp(fname+(strlen(fname)-7),".off.gz",7)) {
      strncpy(sname,fname,strlen(fname)-3);
      sname[strlen(fname)-3]='\0';
    } else {
      sprintf(sname,"%s",fname);
    }
    if(!off_kobj_write_offe(sname,obj)) {
      fprintf(stderr,"ERROR: off_kobj_write_offe(sname,obj)\n");
      return 0;
    }
  } else {
    fprintf(stdout,"      edge file is not written \n");
  }

  /* spherical coordinates */
  if(verbose>0) fprintf(stderr,"[%s] writing spherecoo\n",fcname);
  if(obj->spherecoo==1) {
    NIIK_EXIT((!off_kobj_write_off_sph(fname,obj)),fcname,"off_kobj_write_off_sph",1);
  }

  return 1;
} /* off_kobj_write_off */



/*************************************************
 * read off file
 *
 * -very important function
 *
 * see also
 *
 * kobj *off_kobj_read_off(char * fname);
 * int off_kobj_read_offe(char *fname,kobj *obj);
 * int off_kobj_write_offe(char *fname,kobj *obj);
 * int off_kobj_update_vertex_nei(kvert *vert,kedge *edge);
 * int off_kobj_write_off(char * fname,kobj *obj);
 *
 **************************************************/

kobj *off_kobj_read_off(const char * fname) {
  FILE *fp;
  gzFile gp;
  kobj *obj;
  kvert *v,*vc=NULL,**vlist;
  kface *f,*fc=NULL;
  kedge *e;
  char
  *cp1,
  *tag,str[4096];
  int
  flag_tag_found=0,
  flag_4d_vertex=0,     /* option in off file */
  flag_vertex_normal=0, /* option in off file */
  flag_vertex_color=0,  /* option in off file */
  iter=0,
  i,n,idx[512];
  niikpt p,norm,color;
  double dlist[512];
  int verbose=0;
  const char *fcname="off_kobj_read_off";

  tag = (char *)calloc(20,sizeof(char));

  obj=off_obj_init();

  fp=NULL;
  gp=NULL;

  /* OPEN FILE -either gz or off */
  if(verbose) fprintf(stdout,"  reading %s\n",fname);
  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7)) {
    if((gp = gzopen(fname,"r"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return 0;
    }
  } else if(!strncmp(fname+(strlen(fname)-4),".off",4)) {
    /* 2012-06-20, Kunio check for gzip format */
    sprintf(str,"%s.gz",fname);
    if((gp=gzopen(str,"r"))==NULL) {
      if((fp=fopen(fname,"r"))==NULL) {
        fprintf(stderr,"ERROR: fopen %s\n",fname);
        return NULL;
      }
    }
  } else {
    if((fp=fopen(fname,"r"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",fname);
      return NULL;
    }
  }
  if(verbose) fprintf(stdout,"    reading tag \n");
  if(fp==NULL) {
    if((gzgets(gp,str,4096))==NULL) {
      fprintf(stderr,"ERROR: no first line (gzFile)\n");
      return NULL;
    }
  } else {
    if((fgets(str,4096,fp))==NULL) {
      fprintf(stderr,"ERROR: no first line\n");
      return NULL;
    }
  }
  if((sscanf(str,"%s",tag))!=1) {
    fprintf(stderr,"ERROR: reading tag \n");
    if(fp!=NULL) fclose (fp);
    if(gp!=NULL) gzclose (gp);
    return NULL;
  }

  iter=0;
  while(!flag_tag_found) {
    if( tag[0] == '4' ) {
      flag_4d_vertex = 1;
      tag = tag+1;
    }

    if( tag[0] == 'N' ) {
      flag_vertex_normal = 1;
      tag = tag+1;
    }

    if( tag[0] == 'C' ) {
      flag_vertex_color = 1;
      tag = tag+1;
    }

    if( (tag[0] == 'S') && (tag[0] == 'T') ) {
      fprintf(stderr,"ERROR: texture coordinate is not supported\n");
      exit(0);
    }

    if( (tag[0] == 'O') &&
        (tag[1] == 'F') &&
        (tag[1] == 'F') ) {
      flag_tag_found=1;
    }
    iter++;
    if(iter==10) break;
  }
  if( (tag[0] != 'O') ||
      (tag[1] != 'F') ||
      (tag[1] != 'F') ) {
    fprintf(stderr,"ERROR: wrong tag \"%s\"\n",tag);
    if(fp!=NULL) fclose (fp);
    if(gp!=NULL) gzclose (gp);
    return NULL;
  }

  if(verbose) fprintf(stdout,"\t  tag = %s\n",tag);

  /* create an object now */
  obj->fname = (char *)calloc(strlen(fname)+1,sizeof(char));
  strncpy(obj->fname,fname,strlen(fname));
  obj->fname[strlen(fname)]=0;

  if(verbose) fprintf(stdout,"    reading #vfe \n");
  if(fp!=NULL) {
    if((fscanf(fp,"%i %i %i",&obj->nvert, &obj->nface, &obj->nedge))!=3) {
      fprintf(stderr,"ERROR: reading #v #f #e \n");
      fclose (fp);
      obj=off_kobj_free(obj);
      return NULL;
    }
  } else {
    if((gzgets(gp,str,4096))==NULL) {
      fprintf(stderr,"ERROR: no first line (gzFile)\n");
      gzclose (gp);
      return NULL;
    }
    if((sscanf(str,"%i %i %i",&obj->nvert, &obj->nface, &obj->nedge))!=3) {
      fprintf(stderr,"ERORR: need #vert, #face, #edge here\n");
      gzclose (gp);
      return NULL;
    }
  }
  if(verbose) fprintf(stdout,"\t  #(v,f,e) = %i %i %i\n",obj->nvert, obj->nface, obj->nedge);


  /* reading vertices */
  if(verbose) fprintf(stdout,"    reading vertices \n");
  vlist=off_vert_list(obj->nvert);
  for(n=0; n<obj->nvert; n++) {
    if(fp!=NULL) {
      if((fscanf(fp,"%lf %lf %lf ",&p.x,&p.y,&p.z))!=3) {
        fprintf(stderr,"ERROR: reading point %i \n",n);
        fclose (fp);
        obj=off_kobj_free(obj);
        return NULL;
      }
      if(flag_4d_vertex) {
        if((fscanf(fp,"%lf ",&p.w))!=1) {
          fprintf(stderr,"ERROR: reading point %i \n",n);
          fclose (fp);
          obj=off_kobj_free(obj);
          return NULL;
        }
      }
    } /* FILE pointer */
    else if(gp!=NULL) {
      if((gzgets(gp,str,4096))==NULL) {
        fprintf(stderr,"ERROR: no first line (gzFile)\n");
        gzclose (gp);
        return NULL;
      }
      if((sscanf(str,"%lf %lf %lf",&p.x,&p.y,&p.z))!=3) {
        fprintf(stderr,"ERORR: 3d coordinate for %i\n",n);
        gzclose (gp);
        return NULL;
      }
      if(flag_4d_vertex) {
        fprintf(stderr,"ERROR: 4d vertex is not supported with off.gz\n");
        gzclose (gp);
        return NULL;
      }
    } /* gzFile pointer */
    if(verbose>2) fprintf(stdout,"\t\t%3i  %10.5f %10.5f %10.5f\n",n,p.x,p.y,p.z);
    p.w=0;

    /* assign value and put in the list */
    v=off_vert_init();
    v->v=p;
    v->index=n+1;
    if(obj->vert==NULL) {
      obj->vert=v;
      vc=v;
    } else {
      if(vc==NULL) {
        fprintf(stderr,"ERROR: vc is still null\n");
        return NULL;
      }
      vc->next=v;
      v->prev=vc;
      vc=v;
    }
    vlist[n]=v;

    if(fp!=NULL) {
      if(flag_vertex_normal) {
        if((fscanf(fp,"%lf %lf %lf ",&norm.x,&norm.y,&norm.z))!=3) {
          fprintf(stderr,"ERROR: reading normal %i \n",n);
          fclose (fp);
          obj=off_kobj_free(obj);
          return NULL;
        }
        norm.w=0;
        v->normal = norm;
      } /* vertex normal */
      if(flag_vertex_color) {
        if((fscanf(fp,"%lf %lf %lf %lf ",&color.x,&color.y,&color.z,&color.w))!=4) {
          fprintf(stderr,"ERROR: reading color %i \n",n);
          fclose (fp);
          obj=off_kobj_free(obj);
          return NULL;
        }
        v->color=(double *)calloc(4,sizeof(double));
        v->color[0]=color.x;
        v->color[1]=color.y;
        v->color[2]=color.z;
        v->color[3]=color.w;
      } /* vertex color */
    }  /* FILE pointer */

    else if(gp!=NULL) {
      if(flag_vertex_normal) {
        fprintf(stderr,"ERROR: vertex_normal is not supported with off.gz\n");
        gzclose (gp);
        return NULL;
      }
      if(flag_vertex_color) {
        fprintf(stderr,"ERROR: vertex_color is not supported with off.gz\n");
        gzclose (gp);
        return NULL;
      }
    } /* gzFile pointer */
  } /* each vertex */


  if(verbose>1)
    for(v=obj->vert; v!=NULL; v=v->next) {
      fprintf(stdout,"\t\t%3i  %10.5f %10.5f %10.5f\n",v->index,v->v.x,v->v.y,v->v.z);
    }


  /* reading triangles */
  if(verbose>=1)
    fprintf(stdout,"[%s] reading triangles\n",fcname);

  for(n=0; n<obj->nface; n++) {
    if(verbose>=2) fprintf(stdout,"[%s] reading triangle %i\n",fcname,n);

    /* get the whole line */
    if(fp!=NULL) {
      if((fgets(str,4096,fp))==NULL) {
        fprintf(stderr,"ERROR: fgets line %i\n",n);
        fclose(fp);
        obj=off_kobj_free(obj);
        return NULL;
      }
    } /* FILE pointer */
    else if(gp!=NULL) {
      if((gzgets(gp,str,4096))==NULL) {
        fprintf(stderr,"ERROR: fgets line %i\n",n);
        gzclose(gp);
        obj=off_kobj_free(obj);
        return NULL;
      }
    } /* gzFile pointer */

    if(verbose>=3) fprintf(stdout,"    %s \n",str);

    idx[0] = atoi(str);
    /* fprintf(stdout,"    N = %i \n",idx[0]); */
    if(idx[0]>512) {
      fprintf(stderr,"ERROR: polygon with too many vertices\n");
      exit(0);
    }

    /*
     * Kunio 2012-03-08
     * -this part had a problem when there was no space at the end
     */
    cp1=strchr(str,' ');
    for(i=0; i<idx[0]; i++) {
      while(cp1[1]==' ') cp1=cp1+1; /* I needed to skip white spaces */
      idx[i+1]=atoi(cp1+1);
      /*fprintf(stderr,"  idx[%i]=%i\n",i+1,idx[i+1]);
      fprintf(stderr,"  cp1  #%s",cp1);*/
      if((strchr(cp1+1,' '))==NULL) {
        for(;;) {
          if     (cp1[1]==' ') cp1++;
          else if(cp1[1]=='.') cp1++;
          else if(cp1[1]=='+') cp1++;
          else if(cp1[1]=='-') cp1++;
          else if(cp1[1]>=48&&cp1[1]<=57) cp1++;
          else break;
        }
        if(cp1[1]!='\n') {
          fprintf(stderr,"ERROR: no additional number found\n");
          fprintf(stderr,"       looking for %i vertices\n",idx[0]);
          fprintf(stderr,"  cp1  #%s",cp1);
          if(fp!=NULL) fclose(fp);
          if(gp!=NULL) gzclose(gp);
          obj=off_kobj_free(obj);
          return NULL;
        }
      } else {
        cp1=strchr(cp1+1,' ');
      }
      /* fprintf(stdout,"      vtx = %i \n",idx[i+1]); */
    }

    /* assign values and put in the list */
    f = off_face_init();
    f -> vert[0] = vlist[idx[1]];
    f -> vert[1] = vlist[idx[2]];
    f -> vert[2] = vlist[idx[3]];
    f -> index = n+1;

    if(idx[1] == idx[2]) {
      fprintf(stderr,"ERROR: same vertex, face %i with %i\n",n+1,idx[1]);
      return NULL;
    }

    if(obj->face==NULL) {
      obj->face=f;
      fc=f;
    } else {
      if(fc==NULL) {
        fprintf(stderr,"ERROR: fc is still null\n");
        return NULL;
      }
      fc->next=f;
      f->prev=fc;
      fc=f;
    }

    /* continue to scan the line for color etc */
    if((strchr(cp1+1,' '))==NULL) {
      continue;
    }

    i=sscanf(cp1,"%lf %lf %lf %lf ",&dlist[0],&dlist[1],&dlist[2],&dlist[3]);
    switch(i) {
    case 1:
      fprintf(stderr,"ERROR: case of colormap [1]   index %i, %f\n",n,dlist[0]);
      fprintf(stderr,"ERROR: it's not supported\n");
      fprintf(stderr,"       cp1 = %s\n",cp1);
      exit(0);
    case 3: /* 3-color */
      f->color=(double *)calloc(4,sizeof(double));
      for(i=0; i<3; i++) f->color[i]=dlist[i];
      if(dlist[0]>1 || dlist[1]>1 || dlist[2]>1)
        for(i=0; i<3; i++) f->color[i]=dlist[i]/255.0;
      obj->color=1;
      break;
    case 4: /* 4-color */
      f->color=(double *)calloc(4,sizeof(double));
      for(i=0; i<4; i++) f->color[i]=dlist[i];
      if(dlist[0]>1 || dlist[1]>1 || dlist[2]>1 || dlist[3]>1)
        for(i=0; i<3; i++) f->color[i]=dlist[i]/255.0;
      for(i=0; i<4; i++) f->color[i]=dlist[i];
      obj->color=1;
      break;
    default:
    case 0:
      break;
    }
    /* fprintf(stdout,"\t\t%8.5f %8.5f %8.5f %8.5f\n",dlist[0],dlist[1],dlist[2],dlist[3]); */

    if(verbose>2) {
      fprintf(stdout,"\t%3i %6i %6i %6i\n",idx[0],idx[1],idx[2],idx[3]);
      for(i=0; i<3; i++) {
        fprintf(stdout,"\t  %6.2f %6.2f %6.2f %i\n",
                f->vert[i]->v.x,f->vert[i]->v.y,f->vert[i]->v.z,f ->vert[i]->index);
      }
    }

  } /* each face */
  free(vlist);


  /* add missing colors */
  if(verbose>=1) fprintf(stdout,"[%s] add missing colors\n",fcname);
  if(obj->color) {
    for(f=obj->face; f!=NULL; f=f->next) {
      if(f->color==NULL) {
        f->color=(double *)calloc(4,sizeof(double));
        for(i=0; i<4; i++) f->color[i]=0;
      }
    }
  }

  if(verbose>=1) {
    for(f=obj->face; f!=NULL; f=f->next) {
      fprintf(stdout,"\t%3i %6i %6i %6i\n",idx[0],idx[1],idx[2],idx[3]);
      for(i=0; i<3; i++)
        fprintf(stdout,"\t  %6.2f %6.2f %6.2f %i\n",f->vert[i]->v.x,f->vert[i]->v.y,f->vert[i]->v.z,f ->vert[i]->index);
    }
  }

  if(fp!=NULL) fclose (fp);
  if(gp!=NULL) gzclose(gp);


  /* read for spherical coordinates */
  NIIK_RET0((!off_kobj_read_off_sph(fname,obj)),fcname,"off_kobj_read_off_sph");


  /* construct edges */
  if(verbose>=1) fprintf(stdout,"[%s] construct edges\n",fcname);

  if(!strncmp(fname+(strlen(fname)-7),".off.gz",7)) {
    strncpy(str,fname,strlen(fname)-3);
    /* this creates offe string
       str[strlen(fname)-3]='e';
       str[strlen(fname)-2]='\0'; */
    str[strlen(fname)-3]='\0';
  } else {
    sprintf(str,"%s",fname);
  }

  if(verbose>=1) fprintf(stdout,"[%s] off_kobj_read_offe\n",fcname);
  if(!off_kobj_read_offe(str,obj)) {
    fprintf(stderr,"ERROR: off_kobj_read_offe %s\n",str);
    return NULL;
  }

  /* make neighbor list */
  if(verbose>=1) fprintf(stdout,"[%s] make neighbor list",fcname);
  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)  {
      if(!off_kobj_update_vertex_nei(e->endpts[0],e)) {
        fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
        return NULL;
      }
    }
    if(!e->endpts[1]->nei)  {
      if(!off_kobj_update_vertex_nei(e->endpts[1],e)) {
        fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
        return NULL;
      }
    }
  }
  free(tag);
  off_kobj_update_num(obj);

  /* -removed this because FreeSurfer's object would fail this test
   *  due to small valence (#nei == 3)
   * -Kunio 2012-03-08
   if(!off_kobj_test(obj)){
   fprintf(stderr,"ERROR: off_kobj_test\n");
   fprintf(stderr,"       %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
   return NULL; }*/

  if(verbose) fprintf(stdout,"  finished reading %s\n",fname);
  return obj;
} /* kobj *off_kobj_read_off(char * fname) */




/*************************************************
 * writes offs file
 *   for spherical coordinates (angles)
 *************************************************/

int off_kobj_write_off_sph(const char *fname,kobj *obj) {
  char sname[4096];
  const char *fcname="off_kobj_write_off_sph";
  kvert *v;
  double *buf=NULL;
  int n,verbose=0;
  gzFile gp;
  FILE *fp=NULL;
  if(verbose>1) niik_fc_display(fcname,1);
  NIIK_RET0((fname==NULL),fcname,"fname is null");
  NIIK_RET0((  obj==NULL),fcname,"obj is null");
  sprintf(sname,"%ssp",fname);

  if(!USE_ASCII_SPHERICAL_COORDINATE_OFF) {
    if((gp = gzopen(sname,"wb"))==NULL) {
      fprintf(stderr,"[%s] ERROR: fopen %s\n",fcname,sname);
      return 0;
    }

    buf=(double *)calloc(2*obj->nvert,sizeof(double));

    for(v=obj->vert,n=0; v!=NULL; v=v->next) {
      buf[n++]=v->sph.the;
      buf[n++]=v->sph.psi;
    }
    gzwrite(gp,buf,sizeof(double)*2*obj->nvert);
    gzclose(gp);
    free(buf);
  } else {
    if(verbose>1) fprintf(stdout,"  writing spherical coordinates\n");
    NIIK_RET0(((fp=fopen(sname,"wb"))==NULL),fcname,"fopen");
    for(v=obj->vert,n=0; v!=NULL; v=v->next) {
      fprintf(fp,"%20.15f %20.15f\n",v->sph.the,v->sph.psi);
    }
    fclose(fp);
  }
  if(verbose>1) niik_fc_display(fcname,0);
  return 1;
} /* off_kobj_write_off_sph */

int off_kobj_read_off_sph(const char *fname,kobj *obj) {
  char sname[4096];
  const char *fcname="off_kobj_read_off_sph";
  kvert *v;
  double *buf=NULL;
  int n;
  gzFile gp;
  FILE *fp=NULL;
  const int verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((fname==NULL),fcname,"fname is null");
  NIIK_RET0((  obj==NULL),fcname,"obj is null");
  sprintf(sname,"%ssp",fname);
  if(!USE_ASCII_SPHERICAL_COORDINATE_OFF) {
    if( access(sname,F_OK)==-1) { /* file does not exist, so don't do anything */
      if(verbose>0) fprintf(stdout,"[%s] file does not exist, %s\n",fcname,sname);
      obj->spherecoo=0;
      return 1;
    }
    if((gp = gzopen(sname,"rb"))==NULL) {
      fprintf(stderr,"[%s] ERROR: fopen %s\n",fcname,sname);
      return 0;
    }
    buf=(double *)calloc(2*obj->nvert,sizeof(double));
    if(verbose>0) fprintf(stdout,"[%s]   reading file, %s\n",fcname,sname);
    if( gzread(gp,buf,sizeof(double)*2*obj->nvert) != (sizeof(double)*2*obj->nvert)) {
      fprintf(stderr,"[%s] spherical coordinates file %s is corrupt!\n",fcname,sname);
      free(buf);
      return 0;
    }

    for(v=obj->vert,n=0; v!=NULL; v=v->next) {
      v->sph.the=buf[n++];
      v->sph.psi=buf[n++];
    }
    gzclose(gp);
  } else {
    if((fp=fopen(sname,"rb"))==NULL) {
      fprintf(stderr,"[%s] ERROR: fopen %s\n",fcname,sname);
      return 0;
    }
    for(v=obj->vert,n=0; v!=NULL; v=v->next) {
      if((fscanf(fp,"%lf %lf\n",&v->sph.the,&v->sph.psi))!=2) {
        fprintf(stderr,"[%s:%i:%s] ERROR: fscanf, n=%i\n",__FILE__,__LINE__,__func__,n);
        return 0;
      }
    }
    fclose(fp);
  }
  obj->spherecoo=1;
  free(buf);
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* off_kobj_write_off_sph */


int off_avg_smooth_spherical_coordinates(kobj *obj)
/* smoothes spherical coordinates
 * -simple average with the surrounding
 * -converts to Euclidian coordinates then average
 *  then converts to spherical coordinates */
{
  const char *fcname="off_smooth_spherical_coordinates";
  kvert *v=NULL,**vlist=NULL;
  int m,n,num,nvert;
  niikmat *nv;
  const int verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((obj==NULL),fcname,"obj is null");
  for(v=obj->vert,nvert=0; v!=NULL; v=v->next) {
    nvert++;
  }
  nv=niikmat_init(2,nvert);
  if(verbose>0) fprintf(stdout,"[%s]   nv = %i x %i\n",fcname,2,nvert);
  vlist=(kvert **)calloc(199,sizeof(kvert *));
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    niikpt pt;
    vlist[0]=v;
    num=1;
    num=off_get_local_kvert_check(vlist,199,num);
    nv->m[0][n]=nv->m[1][n]=0;
    pt=niikpt_zero();

    for(m=0; m<num; m++) {
      pt=niikpt_add(pt,niiksph_xyz(vlist[m]->sph));
    }
    pt=niikpt_unit(pt);

    nv->m[0][n]=acos(pt.z);
    nv->m[1][n]=atan2(pt.y,pt.x);
  }

  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    v->sph.the = nv->m[0][n];
    v->sph.psi = nv->m[1][n];
    if(v->sph.the<0) v->sph.the+=NIIK_PI2;
    if(v->sph.psi<0) v->sph.psi+=NIIK_PI2;
  }
  nv=niikmat_free(nv);
  free(vlist);
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* int off_smooth_spherical_coordinates(kobj *obj) */


/*************************************************
 * writes offe file
 *
 * -offe file contains the number of eges (in integer)
 *  followed by the indices of edge's adjacent triangles
 *  and by the indices of edge's end points for each edge
 * -offe file is compressed by zlib
 * -check the status of COMPRESS_OFF_EDGE_FILE in nifti1_kunio.h
 *
 * see also:
 *   int off_kobj_read_offe(char *fname,kobj *obj);
 *
 *************************************************/

int off_kobj_write_offe(const char *fname,kobj *obj) {
  gzFile gp;
  int *ivec,num,n;
  kedge *e;
  char ename[4096];
  const int verbose=0;
  const char *fcname="off_kobj_write_offe";

  if(verbose>=1) niik_fc_display(fcname,1);

  if(!strncmp(fname+strlen(fname)-4,"offe",4)) {
    sprintf(ename,"%s",fname);
  } else {
    sprintf(ename,"%se",fname);
  }

#ifdef COMPRESS_OFF_EDGE_FILE
  if(verbose>=1) fprintf(stdout,"[%s] opening %s\n",fcname,ename);
  if((gp = gzopen(ename,"wb"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",ename);
    return 0;
  }
  for(e=obj->edge,num=0; e!=NULL; e=e->next) {
    num++;
  }
  obj->nedge=num;
  num = obj->nedge*4;
  ivec=(int *)calloc(num,sizeof(int));
  if(verbose>=1) fprintf(stdout,"[%s] num = %i \n",fcname,num);
  for(e=obj->edge,num=0; e!=NULL; e=e->next) {
    ivec[num++]=e->endpts[0]->index;
    ivec[num++]=e->endpts[1]->index;
    ivec[num++]=e->adjface[0]->index;
    ivec[num++]=e->adjface[1]->index;
  }
  if(verbose>=1) fprintf(stdout,"[%s]  num = %i \n",fcname,num);
  if(verbose>=2) niik_display_int_vector(ivec+num-20,20);
  gzwrite(gp,ivec,sizeof(int)*num);
  gzclose(gp);

  n=niik_get_min_index_from_int_vector(ivec,num);
  if(verbose>2) {
    fprintf(stdout,"[%s] min value = %i at %i\n",fcname,ivec[n],n);
    niik_display_int_vector(ivec-20,20);
  }
  free(ivec);
  return 1;
#endif

  fprintf(stderr,"[%s] ERROR: FILE fp routine here\n",fcname);
  return 0;
} /* int off_kobj_write_offe(char *fname,kobj *obj) */


/*************************************************
 * read offe file
 *
 * -if offe file does not exist, then it will create edge list
 *  by running off_kobj_read_off_make_edge
 *
 *
 * see also
 * int off_kobj_write_offe(char *fname,kobj *obj);
 *
 *************************************************/

int off_kobj_read_offe(const char *fname,kobj *obj) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  char ename[4096];
  int *ivec,n,m,num;
  kedge *e,*ce=NULL;
  kface **flist,*f;
  kvert **vlist,*v;
  const int verbose=0;
  const char *fcname="off_kobj_read_offe";

  if(verbose>=1) niik_fc_display(fcname,1);
  sprintf(ename,"%se",fname);
  if(verbose>=1) fprintf(stdout,"[%s] opening %s\n",fcname,ename);
  if(verbose>=1) fprintf(stdout,"[%s] vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  /* zipped file */
#ifdef COMPRESS_OFF_EDGE_FILE
  if(verbose>=1) fprintf(stdout,"[%s] using gzfile\n",fcname);
  gp=gzopen(ename,"rb");
  if(gp==NULL) {
    fp=fopen(ename,"rb");
    if(fp==NULL) {
      if(verbose>=1) {
        fprintf(stdout,"  make edge list!\n");
      }
      if(!off_kobj_read_off_make_edge(obj)) {
        fprintf(stderr,"ERROR: off_kobj_read_off_make_edge \n");
        return 0;
      }
      return 1;
    }
    return 0;
  }

  /* calculate the #edges in the file */
  num=0;

  while(!gzeof(gp)) {
    if(gzread(gp,&n,sizeof(int))!=sizeof(int)) break;   /*VF: trying to catch a wiered bug on MacOS X*/
    num++;
  }
  /*num--;*/

  /*if(num!=(obj->nedge*4))
          fprintf(stdout,"[%s] unexpected number of indexes: %d , expecting %d\n",fcname,num,obj->nedge*4);*/

  ivec = (int *)calloc(num,sizeof(int));
  if(0) {
    gzrewind(gp);
  } else {
    gzclose(gp);
    gp=gzopen(ename,"rb");
  }
  gzread(gp,ivec,num*sizeof(int));
  if(verbose>=1) fprintf(stdout,"[%s]  read gzfile %i \n",fcname,num);
  if(verbose>=1) fprintf(stdout,"    %i %i %i %i \n",ivec[0],ivec[1],ivec[2],ivec[3]);
  if(verbose>=2) niik_display_int_vector(ivec+num-20,20);
  gzclose(gp);
  gp=NULL;
  obj->nedge=num/4; /*VF: is this so?*/
  if(verbose>=1) fprintf(stdout,"[%s] no. edge = %i\n",fcname,obj->nedge);
  for(n=0; n<num; n++) {
    if(n%4<2) {
      if(ivec[n]<=0) {
        fprintf(stderr,"[%s] ERROR: wrong value vert, %i, val=%i\n",fcname,n,ivec[n]);
        return 0;
      }
      if(ivec[n]>obj->nvert) {
        fprintf(stderr,"[%s] ERROR: wrong value vert, %i, val=%i\n",fcname,n,ivec[n]);
        return 0;
      }
    } else {
      if(ivec[n]<=0) {
        fprintf(stderr,"[%s] ERROR: wrong value face, %i\n",fcname,n);
        return 0;
      }
      if(ivec[n]>obj->nface) {
        fprintf(stderr,"[%s] ERROR: wrong value face, %i\n",fcname,n);
        return 0;
      }
    }
  }


#else
  fp=fopen(ename,"rb");
  if(fp==NULL) {
    if(verbose>=1) {
      fprintf(stdout,"  make edge list!\n");
    }
    if(!off_kobj_read_off_make_edge(obj)) {
      fprintf(stderr,"EROR:R off_kobj_read_off_make_edge \n");
      return 0;
    }
    return 1;
  }
#endif


  /* not zipped file */
  if(fp!=NULL) {
    if(verbose>=1) fprintf(stdout,"pass 1\n");
    /*this counts the number of bytes in the file */
    fseek(fp, 0, SEEK_END);
    num = ftell(fp);
    num /= 4*sizeof(int);

    if(obj->nedge==0) obj->nedge = num;
    else if(obj->nedge!=num) {
      fprintf(stderr,"ERROR: #edge %i did not match offe (%i)\n",obj->nedge,num);
      return 0;
    }
    rewind(fp);

    ivec = (int *)calloc(4*obj->nedge,sizeof(int));
    if(verbose>=1) fprintf(stdout,"nedge from offe = %i\n",num);
    if((fread(ivec,1,4*obj->nedge,fp)!=num)) {
      fprintf(stderr,"ERROR: fread, %i \n",num);
      exit(0);
    }
    fclose(fp);
  } /* read from FILE fp */


  if(verbose>=1) fprintf(stdout,"[%s] create lists\n",fcname);
  vlist = (kvert **)calloc(obj->nvert,sizeof(kvert *));
  flist = (kface **)calloc(obj->nface,sizeof(kface *));
  for(f=obj->face,n=0; f!=NULL; f=f->next) flist[n++]=f;
  for(v=obj->vert,n=0; v!=NULL; v=v->next) vlist[n++]=v;

  /* nifti_k_spmat_display_int_vector(ivec,num); */

  if(verbose>=1) fprintf(stdout,"[%s] edge import, %i\n",fcname,num);
  for(n=0,m=1; n<num; m++) {
    if(verbose>=2) fprintf(stdout,"  edge %i | %i %i %i %i\n",m,ivec[n],ivec[n+1],ivec[n+2],ivec[n+3]);
    e = off_edge_init();
    e->endpts[0]=vlist[ivec[n++]-1];
    e->endpts[1]=vlist[ivec[n++]-1];
    e->adjface[0]=flist[ivec[n++]-1];
    e->adjface[1]=flist[ivec[n++]-1];
    e->index=m;
    if(obj->edge==NULL) {
      ce = obj->edge = e;
    } else {
      //if(verbose>1) fprintf(stdout,"  edge %i\n",n);
      if(ce==NULL) {
        fprintf(stderr,"ERROR: ce is still null\n");
        return 0;
      }
      ce->next = e;
      e->prev = ce;
      ce = e;
      /* fprintf(stdout,"  edge %i *\n",n); */
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] edge import, %i, done\n",fcname,num);

  free(vlist);
  free(flist);
  free(ivec);


  /* put edge into face->edge */
  if(verbose>=1) fprintf(stdout,"[%s] put edge into face->edge\n",fcname);
  for(e=obj->edge; e!=NULL; e=e->next) {
    if(verbose>=2) {
      fprintf(stdout,"[%s]   edge %i\n",fcname,e->index);
      if(e->adjface[0]==NULL) {
        fprintf(stderr,"[%s] ERROR: e[%i]->adjface[0] is null\n",fcname,e->index);
        return 0;
      }
      fprintf(stdout,"[%s]   edge %i adj 0 %i\n",fcname,e->index,e->adjface[0]->index);
      if(e->adjface[1]==NULL) {
        fprintf(stderr,"[%s] ERROR: e[%i]->adjface[1] is null\n",fcname,e->index);
        return 0;
      }
      fprintf(stdout,"[%s]   edge %i adj 1 %i\n",fcname,e->index,e->adjface[1]->index);
    }
    for(m=0; m<3; m++) {
      if(e->adjface[0]->edge[m]==NULL) {
        e->adjface[0]->edge[m]=e;
        break;
      }
    }
    for(m=0; m<3; m++) {
      if(e->adjface[1]->edge[m]==NULL) {
        e->adjface[1]->edge[m]=e;
        break;
      }
    }
  }

  if(verbose>1) {
    fprintf(stdout,"[%s] display info\n",fcname);
    for(f=obj->face; f!=NULL; f=f->next) {
      fprintf(stdout,"%5i\n",f->index);
      fprintf(stdout,"%5i %5i %5i\n",f->edge[0]->index,f->edge[1]->index,f->edge[2]->index);
    }
    for(e=obj->edge; e!=NULL; e=e->next) {
      fprintf(stdout,"%5i\n",e->index);
      fprintf(stdout,"  e->endpts[0]  = %5i\n",e->endpts[0]->index);
      fprintf(stdout,"  e->endpts[1]  = %5i\n",e->endpts[1]->index);
      fprintf(stdout,"  e->adjface[0] = %5i\n",e->adjface[0]->index);
      fprintf(stdout,"  e->adjface[1] = %5i\n",e->adjface[1]->index);
    }
  }

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* int off_kobj_read_offe(char *fname,kobj *obj)  */


int off_kobj_display_valency_stats(kobj *obj) {
  niikvec *vec;
  kvert *v;
  int n,m,i;
  if(obj==NULL) {
    fprintf(stderr,"[off_kobj_display_valency_stats] ERROR: obj is null\n");
    return 0;
  }
  vec=niikvec_init(obj->nvert);
  for(v=obj->vert,m=i=n=0; v!=NULL; v=v->next) {
    vec->v[n++]=(double)v->nei;
    if(v->nei==6) {
      m++;
    }
    if(v->nei>=5 && v->nei<=7) {
      i++;
    }
  }
  niik_display_stats_for_double_vector(vec->v,vec->num);
  fprintf(stderr,"[off_kobj_display_valency_stats] valence=6     %9i  %9.3f%%\n",m,(double)m/obj->nvert*100.0);
  fprintf(stderr,"[off_kobj_display_valency_stats] valence=5:7   %9i  %9.3f%%\n",i,(double)i/obj->nvert*100.0);
  vec=niikvec_free(vec);
  return 1;
} /* off_kobj_display_valency_stats */

int off_kobj_show_edge_info(kobj *obj) {
  kedge *e;
  for(e=obj->edge; e!=NULL; e=e->next) {
    fprintf(stdout,"%5i %5i %5i %5i \n",
            e->endpts[0]->index,e->endpts[1]->index,
            e->adjface[0]->index,e->adjface[1]->index);
  }
  return 1;
}











/*************************************************
 * make edge list
 *
 * -very slow for large objects (like cortex)
 * -should use offe files
 *
 *
 * see also
 *
 * kobj *off_kobj_read_off(char * fname);
 * int off_kobj_read_offe(char *fname,kobj *obj);
 * int off_kobj_write_offe(char *fname,kobj *obj);
 *
 **************************************************/

int off_kobj_read_off_make_edge(kobj *obj) {
  kface *f,*face;
  kedge *e,*ec=NULL,*et;
  kvert *vend[2],*vert;
  niikpt pt;
  double dmin=10;
  int i,n;
  int verbose=0;
  const char *fcname="off_kobj_read_off_make_edge";

  NIIK_RET0((obj==NULL),fcname,"obj is a null poointer");
  NIIK_RET0((obj->edge!=NULL),fcname,"obj->edge exists");
  obj->nedge=0;

  if(verbose) {
    fprintf(stdout,"[%s] off_kobj_read_off_make_edge\n",fcname);
  }

  /* go thru face and create an edge */
  for(f=obj->face,n=1; f!=NULL; f=f->next) {
    /*if(f->index==6538) verbose=2;
      else verbose=0;*/
    for(i=0; i<3; i++) {
      vend[0]=f->vert[i];
      vend[1]=f->vert[(i+1)%3];

      if(verbose>=2) fprintf(stdout,"  my end points = %i %i \n",vend[0]->index,vend[1]->index);
      /* put it on the list */
      if(obj->edge==NULL) {
        /* make the edge */
        e=off_edge_init();
        obj->nedge++;
        f->edge[i]=e;
        e->adjface[0]=f;
        e->endpts[0]=vend[0];
        e->endpts[1]=vend[1];
        obj->edge=e;
        e->index=n++;
        ec=e;
      } /* first edge */

      else {
        /* check if it is already in the list */
        for(et=obj->edge; et!=NULL; et=et->next) {
          if(et->adjface[1]!=NULL) {
            et=NULL;
            break;
          }
          if(et->endpts[0]==vend[0]) {
            if(et->endpts[1]==vend[1]) {
              fprintf(stderr,"[%s] ERROR: flipped triangle ? redundancy?\n",fcname);
              fprintf(stderr,"  edge: endpts[0]->index = %i \n",et->endpts[0]->index);
              fprintf(stderr,"  edge: endpts[1]->index = %i \n",et->endpts[1]->index);
              fprintf(stderr,"  vertex 1: xyz = %9.4f %9.4f %9.4f\n",et->endpts[0]->v.x,et->endpts[0]->v.y,et->endpts[0]->v.z);
              fprintf(stderr,"  vertex 2: xyz = %9.4f %9.4f %9.4f\n",et->endpts[1]->v.x,et->endpts[1]->v.y,et->endpts[1]->v.z);
              fprintf(stderr,"\n\n[%s] writing error object, tmp_err.off\n",fcname);
              fprintf(stderr,"[%s]   adding white/red color\n",fcname);
              NIIK_RET0((!off_kobj_add_white_color(obj)),fcname,"off_kobj_add_color");
              f->color[1]=f->color[2]=0.0;
              obj->edge=NULL;
              fprintf(stderr,"[%s]   writing data\n",fcname);
              NIIK_RET0((!off_kobj_write_off("tmp_err.off",obj,0)),fcname,"off_kobj_write_off");

              fprintf(stderr,"[%s]   making smaller data\n",fcname);
              pt=et->endpts[0]->v;
              for(f=obj->face; f!=NULL; f=f->next) {
                if(niikpt_distance(f->vert[0]->v,pt)>dmin)
                  if(niikpt_distance(f->vert[1]->v,pt)>dmin)
                    if(niikpt_distance(f->vert[2]->v,pt)>dmin) {
                      off_kface_remove(f,obj);
                      f=obj->face;
                    }
              }
              f=obj->face;
              if(niikpt_distance(f->vert[0]->v,pt)>dmin)
                if(niikpt_distance(f->vert[1]->v,pt)>dmin)
                  if(niikpt_distance(f->vert[2]->v,pt)>dmin) {
                    obj->face=obj->face->next;
                    obj->face->next->prev=NULL;
                  }
              off_kobj_update_all_index(obj);
              fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
              for(vert=obj->vert; vert!=NULL; vert=vert->next) vert->v.w=0;
              for(face=obj->face; face!=NULL; face=face->next) {
                face->vert[0]->v.w=face->vert[1]->v.w=face->vert[2]->v.w=1;
              }
              for(vert=obj->vert; vert!=NULL; vert=vert->next)
                if(vert->v.w==0)
                  off_kvert_remove(vert,obj);
              off_kobj_update_all_index(obj);
              fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
              obj->edge=NULL;
              off_kobj_update_all_index(obj);
              fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);

              fprintf(stderr,"[%s]   writing data\n",fcname);
              NIIK_RET0((!off_kobj_write_off("tmp_err_small.off",obj,0)),fcname,"off_kobj_write_off");

              return 0;
            }
          } /* check for redundancy */
          if(et->endpts[0]==vend[1])
            if(et->endpts[1]==vend[0])
              break;
        }

        if(et==NULL) { /* not in the list, so add it */
          if(verbose>=2)
            fprintf(stdout,"  not in the list, so add it to top %i %i\n",vend[0]->index,vend[1]->index);
          /* make the edge */
          e=off_edge_init();
          obj->nedge++;
          f->edge[i]=e;
          e->adjface[0]=f;
          e->endpts[0]=vend[0];
          e->endpts[1]=vend[1];
          e->index=n++;
          e->next=obj->edge;
          e->next->prev=e;
          obj->edge=e;

        } else { /* it's already in the list so add the second adjface */
          if(verbose>=2)
            fprintf(stdout,"  already in the list %i, so add the second adjface: et %i %i \n",
                    et->index,et->endpts[0]->index,et->endpts[1]->index);
          et->adjface[1]=f;
          f->edge[i]=et;

          /* move this edge to the end so that search is faster */
          if(et==obj->edge) {
            if(verbose>2) fprintf(stdout,"  move it (first item) to the end of the list\n");
            obj->edge=et->next;
            et->next->prev=NULL;
            if(ec==NULL) {
              fprintf(stderr,"[%s] ERROR: ec is still null\n",fcname);
              return 0;
            }
            ec->next=et;
            et->prev=ec;
            ec=et;
            et->next=NULL;
          } else if(et->next==NULL) {
            if(verbose>2) fprintf(stdout,"  it's already at the end\n");
          } else if(et->prev==NULL) {
            if(verbose>2) fprintf(stdout,"  et->prev is null?\n");
          } else {
            if(verbose>2) fprintf(stdout,"  move it to the end of the list\n");
            et->prev->next=et->next;
            et->next->prev=et->prev;
            ec->next=et;
            et->prev=ec;
            ec=et;
            et->next=NULL;
            if(verbose>2) fprintf(stdout,"  move it to the end of the list *\n");
          }
        } /* already in the list */
      } /* not the first edge */

      if(verbose>2) {
        fprintf(stdout,"    intermediate resutls:\n");
        for(e=obj->edge; e!=NULL; e=e->next) {
          if(e->adjface[1]==NULL)  {
            if(e->prev==NULL) {
              fprintf(stdout,"  e[%5i]: %5i %5i   %5i  undef\n",e->index,
                      e->endpts[0]->index,e->endpts[1]->index,
                      e->adjface[0]->index);
            } else {
              fprintf(stdout,"  e[%5i]: %5i %5i   %5i  undef  prev=%i\n",e->index,
                      e->endpts[0]->index,e->endpts[1]->index,
                      e->adjface[0]->index,e->prev->index);
            }
          } else {
            if(e->prev==NULL) {
              fprintf(stdout,"  e[%5i]: %5i %5i   %5i %5i\n",e->index,
                      e->endpts[0]->index,e->endpts[1]->index,
                      e->adjface[0]->index,e->adjface[1]->index);
            } else {
              fprintf(stdout,"  e[%5i]: %5i %5i   %5i %5i   prev=%i\n",e->index,
                      e->endpts[0]->index,e->endpts[1]->index,
                      e->adjface[0]->index,e->adjface[1]->index,e->prev->index);
            }
          }
          /*if(e->prev!=NULL) fprintf(stdout,"      e->prev = %i\n",e->prev->index);
            if(e->next!=NULL) fprintf(stdout,"      e->next = %i\n",e->next->index);*/
        }
      } /* checking */
    } /* each edge and end vertex */
  } /* each face */

  for(e=obj->edge,obj->nedge=0; e!=NULL; e=e->next) {
    obj->nedge++;
  }
  if(verbose) {
    fprintf(stdout,"  edge content: %i \n",obj->nedge);
  }
  for(e=obj->edge; e!=NULL; e=e->next) {
    if(e->endpts[0]==NULL) {
      fprintf(stderr,"[%s] ERROR: edge %i missing endpts[0] | total #e = %i\n",fcname,e->index,obj->nedge);
      return 0;
    } else if(e->endpts[1]==NULL) {
      fprintf(stderr,"[%s] ERROR: edge %i missing endpts[1] | total #e = %i\n",fcname,e->index,obj->nedge);
      return 0;
    } else if(e->adjface[0]==NULL) {
      fprintf(stderr,"[%s] ERROR: edge %i missing adjface[0] | total #e = %i\n",fcname,e->index,obj->nedge);
      fprintf(stderr,"  e->endpts[0]->v->index     = %i\n",e->endpts[0]->index);
      fprintf(stderr,"  e->endpts[0]->v->v.[x,y,z] = %9.5f %9.5f %9.5f\n",e->endpts[0]->v.x,e->endpts[0]->v.y,e->endpts[0]->v.z);
      fprintf(stderr,"  e->endpts[1]->v->index     = %i\n",e->endpts[1]->index);
      fprintf(stderr,"  e->endpts[1]->v->v.[x,y,z] = %9.5f %9.5f %9.5f\n",e->endpts[1]->v.x,e->endpts[1]->v.y,e->endpts[1]->v.z);
      return 0;
    } else if(e->adjface[1]==NULL) {
      fprintf(stderr,"[%s] ERROR: edge %i missing adjface[1] | total #e = %i\n",fcname,e->index,obj->nedge);
      fprintf(stderr,"  e->endpts[0]->v->index     = %i\n",e->endpts[0]->index);
      fprintf(stderr,"  e->endpts[0]->v->v.[x,y,z] = %9.5f %9.5f %9.5f\n",e->endpts[0]->v.x,e->endpts[0]->v.y,e->endpts[0]->v.z);
      fprintf(stderr,"  e->endpts[1]->v->index     = %i\n",e->endpts[1]->index);
      fprintf(stderr,"  e->endpts[1]->v->v.[x,y,z] = %9.5f %9.5f %9.5f\n",e->endpts[1]->v.x,e->endpts[1]->v.y,e->endpts[1]->v.z);
      return 0;
    }
  }

  if(verbose>1) fprintf(stdout,"  forward MAKE EDGE \n");
  for(e=obj->edge,n=1; e!=NULL; e=e->next,n++) {
    /* fprintf(stdout,"e = %7i\n",e->index); */
    if(e->next==NULL) break;
  }
  if(verbose>1) fprintf(stdout,"  backward MAKE EDGE \n");
  for(e=e; e!=NULL; e=e->prev,n--) {
    /* fprintf(stdout,"e = %7i\n",e->index); */
  }
  if(n) fprintf(stdout,"  index = %i\n",n);
  /*exit(0); */

  return 1;
} /* int off_kobj_read_off_make_edge(kobj *obj) */






/**************************************************
 * make neighbor list
 *
 * -update verte's neilist
 * -called when reading an off file as well as
 *  when mesh configuration changes
 *
 **************************************************/

int off_kobj_update_vertex_nei(kvert *vert,kedge *edge) {
  kface *cf,*flist[4096];
  kedge *ce,*elist[4096];
  kvert *vlist[4096];
  int n;
  const int verbose=0;
  /* Testing, Kunio, 2012-08-23,
     if(edge->index==9532) verbose=1;
     if(vert->index==2661) verbose=1;*/

  if(vert==NULL) {
    fprintf(stderr,"ERROR: vert is a null pointer\n");
    return 0;
  }
  if(edge==NULL) {
    fprintf(stderr,"ERROR: edge is a null pointer\n");
    return 0;
  }

  if(verbose>=1) fprintf(stdout,"\nstarting off_kobj_update_vertex_nei  (v=%i, e=%i)\n",vert->index,edge->index);

  vert->nei=0;
  ce = edge;

  elist[vert->nei]=ce;
  if(ce->endpts[0]==vert) vlist[vert->nei]=ce->endpts[1];
  else                    vlist[vert->nei]=ce->endpts[0];

  if(vert->neivert!=NULL) free(vert->neivert);
  if(vert->neiface!=NULL) free(vert->neiface);
  if(vert->neiedge!=NULL) free(vert->neiedge);

  if(verbose>=1) fprintf(stdout," edge = %i [%i %i] (start)  \n",ce->index,ce->adjface[0]->index,ce->adjface[1]->index);
  cf = ce->adjface[0];
  if(verbose>=1) fprintf(stdout," face = %i  \n",cf->index);
  flist[vert->nei]=cf;

  for(;;) {
    if(verbose>=1) fprintf(stdout," face->edge = %i %i %i  \n",cf->edge[0]->index,cf->edge[1]->index,cf->edge[2]->index);
    for(n=0; n<3; n++) {
      if(verbose>=1) fprintf(stdout,"   face->edge[%i,%i]->endpts = %i %i\n",cf->edge[n]->index,n,
                               cf->edge[n]->endpts[0]->index,
                               cf->edge[n]->endpts[1]->index);
      if(cf->edge[n]==ce) continue;
      if(cf->edge[n]->endpts[0]==vert) break;
      if(cf->edge[n]->endpts[1]==vert) break;
    }
    if(n==3) {
      fprintf(stderr,"ERROR: edge not found\n");
      fprintf(stderr," face=%i -> edge = %i %i %i\n",cf->index,cf->edge[0]->index,cf->edge[1]->index,cf->edge[2]->index);
      fprintf(stderr," looking for %i\n",ce->index);
      fprintf(stderr,"ERROR: edge not found; tri[%i]-> v=%i %i %i | edge %i(pt=%i,%i) %i(pt=%i,%i) %i(pt=%i,%i)\n",
              cf->index,
              cf->vert[0]->index,cf->vert[1]->index,cf->vert[2]->index,
              cf->edge[0]->index,cf->edge[0]->endpts[0]->index,cf->edge[0]->endpts[1]->index,
              cf->edge[1]->index,cf->edge[1]->endpts[0]->index,cf->edge[1]->endpts[1]->index,
              cf->edge[2]->index,cf->edge[2]->endpts[0]->index,cf->edge[2]->endpts[1]->index);
      fprintf(stderr,"       off_kobj_update_vertex_nei v=%i e=%i\n",vert->index,edge->index);
      exit(0);
    }
    if(verbose>=1) fprintf(stdout,"   edge = %i\n",cf->edge[n]->index);
    ce = cf->edge[n];
    vert->nei++;
    if(vert->nei>=4096) {
      fprintf(stderr,"ERROR: too many neighbors %i\n",vert->nei);
      return 0;
    }
    if(ce==edge) break;

    if(verbose>=1) fprintf(stdout," edge = %i [%i %i]\n",ce->index,ce->adjface[0]->index,ce->adjface[1]->index);
    if(ce->adjface[0]==cf) cf = ce->adjface[1];
    else                   cf = ce->adjface[0];

    elist[vert->nei]=ce;
    if(ce->endpts[0]==vert) vlist[vert->nei]=ce->endpts[1];
    else                    vlist[vert->nei]=ce->endpts[0];
    flist[vert->nei]=cf;

    if(verbose>=1) fprintf(stdout," face = %i  \n",cf->index);
  }

  vert->neivert=off_vert_list(vert->nei);
  vert->neiface=off_face_list(vert->nei);
  vert->neiedge=off_edge_list(vert->nei);

  for(n=0; n<vert->nei; n++) {
    vert->neivert[n]=vlist[n];
    vert->neiface[n]=flist[n];
    vert->neiedge[n]=elist[n];
  }

  if(verbose>=1)
    for(n=0; n<vert->nei; n++)
      fprintf(stdout,"\t%i %i %i\n",vert->neivert[n]->index,vert->neiface[n]->index,vert->neiedge[n]->index);

  return 1;
} /* off_kobj_update_vertex_nei */









/***********************************************
 *
 * update object's nvert, nface, nedge
 *
 *
 ***********************************************/

void off_kobj_update_num(kobj *obj) {
  off_kobj_update_nvert(obj);
  off_kobj_update_nface(obj);
  off_kobj_update_nedge(obj);
}

void off_kobj_update_nvert(kobj *obj) {
  kvert *v;
  for(v=obj->vert,obj->nvert=0; v!=NULL; v=v->next) obj->nvert++;
}

void off_kobj_update_nface(kobj *obj) {
  kface *f;
  for(f=obj->face,obj->nface=0; f!=NULL; f=f->next) obj->nface++;
}

void off_kobj_update_nedge(kobj *obj) {
  kedge *e;
  for(e=obj->edge,obj->nedge=0; e!=NULL; e=e->next) obj->nedge++;
}

void off_kobj_update_all_index(kobj *obj) {
  off_kobj_update_vert_index(obj);
  off_kobj_update_face_index(obj);
  off_kobj_update_edge_index(obj);
}

void off_kobj_update_vert_index(kobj *obj) {
  kvert *v;
  int n=0;
  for(v=obj->vert,n=1; v!=NULL; v=v->next) v->index=n++;
  obj->nvert=n-1;
}

void off_kobj_update_face_index(kobj *obj) {
  kface *v;
  int n=0;
  for(v=obj->face,n=1; v!=NULL; v=v->next) v->index=n++;
  obj->nface=n-1;
}

void off_kobj_update_edge_index(kobj *obj) {
  kedge *v;
  int n=0;
  for(v=obj->edge,n=1; v!=NULL; v=v->next) v->index=n++;
  obj->nedge=n-1;
}


/***********************************************
 *
 * update object's color (face)
 *
 * -vert color is supported but never used in my library
 *
 *
 *  off_kobj_gray_kface
 *  off_kobj_add_color
 *  off_kobj_remove_color
 *
 *
 ***********************************************/

int off_kobj_gray_kface(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=f->color[1]=f->color[2]=f->color[3]=0.66;
  }
  return 1;
}

int off_kobj_add_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
      f->color[0]=f->color[1]=f->color[2]=f->color[3]=0.66;
    }
  }
  return 1;
}

int off_kobj_remove_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=0;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color!=NULL) {
      free(f->color);
      f->color=NULL;
    }
  }
  return 1;
}

int off_kobj_add_white_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=f->color[1]=f->color[2]=f->color[3]=1.0;
  }
  return 1;
}

int off_kobj_add_green_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=0;
    f->color[1]=1;
    f->color[2]=0;
    f->color[3]=0;
  }
  return 1;
}

int off_kobj_add_blue_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=0;
    f->color[1]=0;
    f->color[2]=1;
    f->color[3]=0;
  }
  return 1;
}

int off_kobj_add_red_color(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=1;
    f->color[1]=0;
    f->color[2]=0;
    f->color[3]=0;
  }
  return 1;
}

int off_kobj_add_one_color(kobj *obj,double R, double G, double B) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(R>1 || R<0) {
    fprintf(stderr,"ERROR: R is out of range %9.4f\n",R);
    return 0;
  }
  if(G>1 || G<0) {
    fprintf(stderr,"ERROR: G is out of range %9.4f\n",G);
    return 0;
  }
  if(B>1 || B<0) {
    fprintf(stderr,"ERROR: B is out of range %9.4f\n",B);
    return 0;
  }
  obj->color=1;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color==NULL) {
      f->color=(double *)calloc(4,sizeof(double));
    }
    f->color[0]=R;
    f->color[1]=G;
    f->color[2]=B;
    f->color[3]=0;
  }
  return 1;
}

int off_kobj_apply_color_map(kobj *obj,double *var, double dmin,double dmax,int colormap) {
  kface *f;
  double dran;
  int n,num=512;
  niikmat *cm;
  cm = niik_colormap_get(colormap,num);
  /*niikmat_display(cm);*/
  dran = (dmax-dmin)/(num-1.0);
  off_kobj_add_color(obj);
  for(f=obj->face; f!=NULL; f=f->next) {
    n = ((var[f->vert[0]->index-1] + var[f->vert[1]->index-1] + var[f->vert[2]->index-1]) / 3.0 - dmin) / dran;
    if(n<0) n=0;
    else if(n>=num) n=num-1;
    f->color[0] = cm->m[n][0];
    f->color[1] = cm->m[n][1];
    f->color[2] = cm->m[n][2];
    f->color[3] = 0;
  }
  cm = niikmat_free(cm);
  return 1;
}


/***********************************************
 *
 * create a simple icosahedron
 *
 *
 ***********************************************/

kobj *off_create_icosahedron() {
  kobj *obj;
  kvert *v,*vlist[12];
  kface *f;
  kedge *e;
  int n,m;
  int fvlist[60] = { 0,2,10,0,8,2,0,10,6,0,6,4,2,7,10,2,5,7,2,8,5,10,7,11,10,11,6,7,5,3,7,3,11,8,9,5,8,0,4,8,4,9,11,3,1,11,1,6,5,9,3,9,4,1,9,1,3,4,6,1 };
  const char *fcname="off_create_icosahedron";

  if(0) niik_fc_display(fcname,1);

  obj = off_obj_init();
  obj->vert = v = off_vert_init();
  v->v.x = 0.0000000000;
  v->v.y = -0.5257310867;
  v->v.z = 0.8506507874;  /* index=0 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.0000000000;
  v->v.y = -0.5257310867;
  v->v.z = -0.8506507874;   /* index=1 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.0000000000;
  v->v.y = 0.5257310867;
  v->v.z = 0.8506507874;   /* index=2 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.0000000000;
  v->v.y = 0.5257310867;
  v->v.z = -0.8506507874;   /* index=3 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.5257310867;
  v->v.y = -0.8506507874;
  v->v.z = 0;   /* index=4 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.5257310867;
  v->v.y = 0.8506507874;
  v->v.z = 0;   /* index=5 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = -0.5257310867;
  v->v.y = -0.8506507874;
  v->v.z = 0;   /* index=6 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = -0.5257310867;
  v->v.y = 0.8506507874;
  v->v.z = 0;   /* index=7 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.8506507874;
  v->v.y = 0.0;
  v->v.z = 0.5257310867;   /* index=8 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = 0.8506507874;
  v->v.y = 0.0;
  v->v.z = -0.5257310867;   /* index=9 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = -0.8506507874;
  v->v.y = 0.0;
  v->v.z = 0.5257310867;   /* index=10 */
  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v.x = -0.8506507874;
  v->v.y = 0.0;
  v->v.z = -0.5257310867;   /* index=11 */

  for(v=obj->vert,n=0; v!=NULL; v=v->next) vlist[n++]=v;

  /* added spherical coordinates */
  for(n=0; n<12; n++) {
    vlist[n]->sph=niikxyz_sph(vlist[n]->v);
  }
  obj->spherecoo=1;

  m=0;
  obj->face = f = off_face_init();
  for(n=0; n<3; n++) f->vert[n]=vlist[fvlist[m++]];

  while(m<60) {
    f->next = off_face_init();
    f->next->prev = f;
    f=f->next;
    for(n=0; n<3; n++) f->vert[n]=vlist[fvlist[m++]];
  }

  NIIK_RET0((!off_kobj_read_off_make_edge(obj)),fcname,"off_kobj_read_off_make_edge");

  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[0],e)),fcname,"off_kobj_update_vertex_nei");
    if(!e->endpts[1]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[1],e)),fcname,"off_kobj_update_vertex_nei");
  }
  off_kobj_update_num(obj);

  /*fprintf(stdout,"  writing tmp.off\n");
    if(!off_kobj_write_off("tmp.off",obj,0)){
    fprintf(stderr,"ERROR: off_kobj_write_off \n");
    exit(0); }
    exit(0);*/

  return obj;
} /* off_create_icosahedron() */


int off_subdiv_kobj(kobj *obj) {
  kface *f,*nf,*f2,**flist;
  kvert *v,**vlist;
  kedge *e,*e2;
  int i,nvert,nface;
  int verbose=0;

  if(verbose) fprintf(stdout,"  off_subdiv_kobj \n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null poointer\n");
    return 0;
  }

  off_kobj_update_num(obj);
  off_kobj_update_all_index(obj);


  vlist = off_vert_list(obj->nvert + obj->nedge);
  flist = off_face_list(obj->nface * 4);

  /* create a list of new vertices */
  nvert = 0;
  for(e=obj->edge; e!=NULL; e=e->next) {
    v = off_vert_init();
    v->v = off_avg_point_on_edge(e);
    if(obj->spherecoo==1) {
      v->sph=niiksph_avg_lim(e->endpts[0]->sph,e->endpts[1]->sph,NIIK_PI);
    }
    vlist[nvert++]=v;
  }
  for(v=obj->vert; v!=NULL; v=v->next) {
    vlist[nvert++]=v;
  }

  /* create a list of new faces */
  for(f=obj->face,nface=0; f!=NULL; f=f->next) {
    if(!off_match_face_vert_edge(f)) {
      fprintf(stderr,"ERROR: off_match_face_vert_edge \n");
      return 0;
    }

    nf = off_face_init();
    nf->vert[0] = f->vert[0];
    nf->vert[1] = vlist[f->edge[0]->index-1];
    nf->vert[2] = vlist[f->edge[2]->index-1];
    flist[nface++]=nf;

    nf = off_face_init();
    nf->vert[0] = f->vert[1];
    nf->vert[1] = vlist[f->edge[1]->index-1];
    nf->vert[2] = vlist[f->edge[0]->index-1];
    flist[nface++]=nf;

    nf = off_face_init();
    nf->vert[0] = f->vert[2];
    nf->vert[1] = vlist[f->edge[2]->index-1];
    nf->vert[2] = vlist[f->edge[1]->index-1];
    flist[nface++]=nf;

    nf = off_face_init();
    nf->vert[0] = vlist[f->edge[0]->index-1];
    nf->vert[1] = vlist[f->edge[1]->index-1];
    nf->vert[2] = vlist[f->edge[2]->index-1];
    flist[nface++]=nf;
  }

  /* re-create vert/face/edge */
  for(f=obj->face; f!=NULL; f=f2) {
    f2=f->next;
    off_kface_free(f);
  }
  for(e=obj->edge; e!=NULL; e=e2) {
    e2=e->next;
    off_kedge_free(e);
  }
  obj->edge=NULL;
  obj->face=NULL;

  obj->vert = vlist[0];
  for(i=1; i<nvert; i++)   vlist[i]->prev = vlist[i-1];
  for(i=0; i<nvert-1; i++) vlist[i]->next = vlist[i+1];

  obj->face = flist[0];
  for(i=1; i<nface; i++)   flist[i]->prev = flist[i-1];
  for(i=0; i<nface-1; i++) flist[i]->next = flist[i+1];

  /* free nei */
  for(v=obj->vert; v!=NULL; v=v->next) {
    if(v->nei>0) {
      free(v->neivert);
      free(v->neiedge);
      free(v->neiface);
    }
    v->nei=0;
    v->neivert=NULL;
    v->neiface=NULL;
    v->neiedge=NULL;
  }

  if(!off_kobj_read_off_make_edge(obj)) {
    fprintf(stderr,"ERROR: off_kobj_read_off_make_edge\n");
    return 0;
  }

  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)
      if(!off_kobj_update_vertex_nei(e->endpts[0],e)) {
        fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
        return 0;
      }
    if(!e->endpts[1]->nei)
      if(!off_kobj_update_vertex_nei(e->endpts[1],e)) {
        fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(e->endpts[0],e)\n");
        return 0;
      }
  }

  off_kobj_update_num(obj);
  off_kobj_update_all_index(obj);

  free(flist);
  free(vlist);

  return 1;
} /* int off_subdiv_kobj(kobj *obj) */







/***********************************************
 *
 * simple operations
 *
 *
 ***********************************************/

/* may not be working... */
void off_swap_kvert(kvert *e1,kvert *e2) {
  kvert *e;
  e=e1;
  e1=e2;
  e2=e;
}
void off_swap_kface(kface *f1,kface *f2) {
  kface *f;
  f=f1;
  f1=f2;
  f2=f;
}
void off_swap_kedge(kedge *e1,kedge *e2) {
  kedge *e;
  e=e1;
  e1=e2;
  e2=e;
}


/*
 *  match the vert and edge for face (f)
 * so that edge[0] is between vert[0] and vert[1]
 *         edge[1] is between vert[1] and vert[2]
 *         edge[2] is between vert[2] and vert[0]
 */
int off_match_face_vert_edge(kface *f) {
  int i;
  kedge *te;
  if(f==NULL) {
    fprintf(stderr,"ERROR: f is a null pointer \n");
    return 0;
  }

  /*fprintf(stdout,"f[%i]  v = %i %i %i  e = (%i %i) (%i %i) (%i %i)\n",
    f->index,
    f->vert[0]->index,f->vert[1]->index,f->vert[2]->index,
    f->edge[0]->endpts[0]->index,f->edge[0]->endpts[1]->index,
    f->edge[1]->endpts[0]->index,f->edge[1]->endpts[1]->index,
    f->edge[2]->endpts[0]->index,f->edge[2]->endpts[1]->index);*/
  for(i=0; i<3; i++) {
    if     (f->edge[i]->endpts[0]==f->vert[0] && f->edge[i]->endpts[1]==f->vert[1]) break;
    else if(f->edge[i]->endpts[1]==f->vert[0] && f->edge[i]->endpts[0]==f->vert[1]) break;
  }
  if(i==3) return 0;
  if(i!=0) {
    te = f->edge[0];
    f->edge[0] = f->edge[i];
    f->edge[i] = te;
    /*    fprintf(stdout,"SWAPP?\n");
    fprintf(stdout,"  edge 0,1 = %i,%i \n",f->edge[0]->index,f->edge[1]->index);
    off_swap_kedge(f->edge[0],f->edge[i]);
    fprintf(stdout,"  edge 0,1 = %i,%i \n",f->edge[0]->index,f->edge[1]->index);
    exit(0);*/
  }
  for(i=1; i<3; i++) {
    if     (f->edge[i]->endpts[0]==f->vert[1] && f->edge[i]->endpts[1]==f->vert[2]) break;
    else if(f->edge[i]->endpts[1]==f->vert[1] && f->edge[i]->endpts[0]==f->vert[2]) break;
  }
  if(i==3) return 0;
  if(i!=1) {
    /* off_swap_kedge(f->edge[1],f->edge[i]); */
    te = f->edge[1];
    f->edge[1] = f->edge[i];
    f->edge[i] = te;
  }
  /*  fprintf(stdout,"f[%i]  v = %i %i %i  e = (%i %i) (%i %i) (%i %i)\n",
    f->index,
    f->vert[0]->index,f->vert[1]->index,f->vert[2]->index,
    f->edge[0]->endpts[0]->index,f->edge[0]->endpts[1]->index,
    f->edge[1]->endpts[0]->index,f->edge[1]->endpts[1]->index,
    f->edge[2]->endpts[0]->index,f->edge[2]->endpts[1]->index);*/
  return 1;
}

void off_display_vert_info(kvert *v,int type) {
  int n;
  switch(type) {
  default:
  case 0:
    fprintf(stdout,"vert [%i]\n",v->index);
    return;
  case 20:
    fprintf(stdout,"vert[%i] %8.4f %8.4f %8.4f vert\n",v->index,v->v.x,v->v.y,v->v.z);
    fprintf(stdout,"vert[%i] %8.4f %8.4f %8.4f norm\n",v->index,v->normal.x,v->normal.y,v->normal.z);
    fprintf(stdout,"vert[%i] %8.4f %8.4f         the,psi\n",v->index,v->sph.the,v->sph.psi);
  case 10:
    fprintf(stdout,"vert[%i] neiV=[",v->index);
    for(n=0; n<v->nei; n++) {
      fprintf(stdout,"%i ",v->neivert[n]->index);
    }
    fprintf(stdout,"] neiF=[");
    for(n=0; n<v->nei; n++) {
      fprintf(stdout,"%i ",v->neiface[n]->index);
    }
    fprintf(stdout,"] neiE=[");
    for(n=0; n<v->nei; n++) {
      fprintf(stdout,"%i ",v->neiedge[n]->index);
    }
    fprintf(stdout,"]\n");
    return;
  }
  return;
}

void off_display_edge_info(kedge *e,int type) {
  switch(type) {
  case 0:
    fprintf(stdout,"edge [%i]\n",e->index);
    return;
  case 10:
  default:
    fprintf(stdout,"edge[%i] endpts=[%i %i] adjface=[%i %i]\n",e->index,
            e->endpts[0]->index,e->endpts[1]->index,
            e->adjface[0]->index,e->adjface[1]->index);
    return;
  }
  return;
}

void off_display_face_info(kface *f,int type) {
  switch(type) {
  default:
  case 0:
    fprintf(stdout,"face [%i]\n",f->index);
    return;
  case 20:
    fprintf(stdout,"face[%i] v[0]=(%6.2f %6.2f %6.2f) v[1]=(%6.2f %6.2f %6.2f) v[2]=(%6.2f %6.2f %6.2f)\n",f->index,
            f->vert[0]->v.x,f->vert[0]->v.y,f->vert[0]->v.z,
            f->vert[1]->v.x,f->vert[1]->v.y,f->vert[1]->v.z,
            f->vert[2]->v.x,f->vert[2]->v.y,f->vert[2]->v.z);
  case 10:
    fprintf(stdout,"face[%i] v=[%i %i %i] e=[%i (%i %i), %i (%i %i), %i (%i %i)]\n",f->index,
            f->vert[0]->index,f->vert[1]->index,f->vert[2]->index,
            f->edge[0]->index,f->edge[0]->endpts[0]->index,f->edge[0]->endpts[1]->index,
            f->edge[1]->index,f->edge[1]->endpts[0]->index,f->edge[1]->endpts[1]->index,
            f->edge[2]->index,f->edge[2]->endpts[0]->index,f->edge[2]->endpts[1]->index);
    return;
  }
  return;
}

/* returns average point on an edge */
niikpt off_avg_point_on_edge(kedge *e) {
  return niikpt_avg(e->endpts[0]->v,e->endpts[1]->v);
}

/* updates pmin and pmax for a triangle */
void off_update_kface_pminmax(kface *f) {
  f->pmin=f->pmax=f->vert[0]->v;
  f->pmin=niikpt_min(f->pmin,f->vert[1]->v);
  f->pmin=niikpt_min(f->pmin,f->vert[2]->v);
  f->pmax=niikpt_max(f->pmax,f->vert[1]->v);
  f->pmax=niikpt_max(f->pmax,f->vert[2]->v);
}

/* updates pmin and pmax for the neighboring triangles */
void off_update_kvert_pminmax(kvert *v) {
  int n;
  for(n=0; n<v->nei; n++) {
    off_update_kface_pminmax(v->neiface[n]);
  }
}

/* updates pmin and pmax for all triangles in an object */
void off_update_kobj_kface_pminmax(kobj *obj) {
  kface *f;
  for(f=obj->face; f!=NULL; f=f->next) {
    off_update_kface_pminmax(f);
  }
}

/* calculates the min point for an object */
niikpt off_calc_kobj_pmin(kobj *obj) {
  niikpt p;
  kvert *v;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return niikpt_nan();
  }
  p=obj->vert->v;
  for(v=obj->vert; v!=NULL; v=v->next) {
    p=niikpt_min(p,v->v);
  }
  return p;
}

/* calculates the max point for an object */
niikpt off_calc_kobj_pmax(kobj *obj) {
  niikpt p;
  kvert *v;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return niikpt_nan();
  }
  p=obj->vert->v;
  for(v=obj->vert; v!=NULL; v=v->next) {
    p=niikpt_max(p,v->v);
  }
  return p;
}

void off_calc_facelist_bounds(kface **facelist,int num,niikpt *  pmin,niikpt *  pmax) {
  int i;
  if(num<1) return;
  *pmin=facelist[0]->pmin;
  *pmax=facelist[0]->pmax;

  for(i=1; i<num; i++) {
    *pmin=niikpt_min(*pmin,facelist[i]->pmin);
    *pmax=niikpt_max(*pmax,facelist[i]->pmax);
  }
}



/* updates the triangle normal */
void off_update_kface_normal(kface *f) {
  f->normal = niikpt_unit_normal(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v);
}

void off_update_kvert_normal(kvert *v) {
  int n;
  v->normal = v->normal;
  for(n=0; n<v->nei; n++) {
    v->normal = niikpt_add(v->normal,v->neiface[n]->normal);
  }
  v->normal = niikpt_unit(v->normal);
}

int off_update_kobj_face_normal(kobj *obj) {
  kface *f;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return 0;
  }
  for(f=obj->face; f!=NULL; f=f->next) {
    off_update_kface_normal(f);
  }
  return 1;
}

int off_update_kobj_vert_normal(kobj *obj)
/* assuming that face normal is updated */
{
  kvert *v;
  niikpt *nlist;
  int m,n;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return 0;
  }
  if((nlist = niikpt_alloc(obj->nvert))==NULL) {
    fprintf(stderr,"ERROR: niikpt_alloc for nlist\n");
    return 0;
  }
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    nlist[n]=v->normal;
    for(m=0; m<v->nei; m++)
      nlist[n] = niikpt_add(v->neiface[m]->normal,nlist[n]);
  }
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    v->normal=niikpt_unit(nlist[n]);
  }
  free(nlist);
  return 1;
}

int off_smooth_kobj_vert_normal(kobj *obj)
/* assuming that vertex normal is already calcualted once */
{
  char fcname[64]="off_smooth_kobj_vert_normal";
  kvert *v;
  niikpt *nlist;
  int m,n,verbose=0;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return 0;
  }
  nlist = niikpt_alloc(obj->nvert);
  if(nlist==NULL) {
    fprintf(stderr,"ERROR: niikpt_alloc\n");
    return 0;
  }
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    nlist[n]=v->normal;
    if(verbose>=1) fprintf(stdout,"[%s]   vertex %9i / %9i\n",fcname,n,obj->nvert);
    for(m=0; m<v->nei; m++)
      nlist[n] = niikpt_add(v->neivert[m]->normal,nlist[n]);
  }
  if(verbose>=1) fprintf(stdout,"[%s] update\n",fcname);
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    v->normal=niikpt_unit(nlist[n]);
  }
  if(verbose>=1) fprintf(stdout,"[%s] updated\n",fcname);
  free(nlist);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
}

double off_get_kobj_mean_edge_length(kobj *obj) {
  double dsum=0;
  int n=0;
  kedge *e;
  for(e=obj->edge; e!=NULL; e=e->next,n++) {
    dsum += niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
  }
  return dsum/n;
}

int off_display_kobj_edge_stats(kobj *obj) {
  kedge *e;
  double dval,dmin,dsum,dssq,dmax,dmean;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  e=obj->edge;
  dmin=dmax=niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
  for(e=obj->edge,dsum=dssq=0; e!=NULL; e=e->next) {
    dval=niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
    if(dval>dmax) dmax=dval;
    else if(dval<dmin) dmin=dval;
    dsum+=dval;
    dssq+=dval*dval;
  }
  dmean=dsum/obj->nedge;
  fprintf(stdout,"\tedge stats: %9.5f +/- %9.5f (min,max) = %9.5f, %9.5f\n",
          dmean,sqrt(dssq/(obj->nedge)-dmean*dmean),dmin,dmax);
  return 1;
}

double off_get_kobj_volume(kobj *obj) {
  niikpt md;
  kvert *v;
  kface *f;
  double d;
  md=niikpt_zero();
  for(v=obj->vert; v!=NULL; v=v->next) {
    md.x+=v->v.x;
    md.y+=v->v.y;
    md.z+=v->v.z;
  }
  md.x/=obj->nvert;
  md.y/=obj->nvert;
  md.z/=obj->nvert;
  for(d=0,f=obj->face; f!=NULL; f=f->next) {
    d += niikpt_det4(f->vert[0]->v.x,f->vert[1]->v.x,f->vert[2]->v.x,md.x,
                     f->vert[0]->v.y,f->vert[1]->v.y,f->vert[2]->v.y,md.y,
                     f->vert[0]->v.z,f->vert[1]->v.z,f->vert[2]->v.z,md.z,
                     1,1,1,1);
  }
  return d / 6.0;
} /* off_get_kobj_volume */

double off_get_kobj_area(kobj *obj) {
  double dsum=0;
  int n=0;
  kface *f;
  for(f=obj->face; f!=NULL; f=f->next,n++) {
    dsum += niikpt_area2(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v);
  }
  return dsum/2.0;
}

double off_get_kobj_mean_tri_area(kobj *obj) {
  return off_get_kobj_area(obj)/obj->nface;
}

niikpt off_get_kobj_global_min_coord(kobj *obj) {
  niikpt p;
  kvert *v;
  p.x=p.y=p.z=1e30;
  p.w=0;
  for(v=obj->vert; v!=NULL; v=v->next) {
    p=niikpt_min(p,v->v);
  }
  return p;
}

niikpt off_get_kobj_global_max_coord(kobj *obj) {
  niikpt p;
  kvert *v;
  p.x=p.y=p.z=-1e30;
  p.w=0;
  for(v=obj->vert; v!=NULL; v=v->next) {
    p=niikpt_max(p,v->v);
  }
  return p;
}



/************************************************************
 *
 * functions to find local regions
 *
 *
 *
 ************************************************************/
int off_get_local_kvert(kvert **vlist,int num) {
  int nei,n,m,N;
  kvert *v;
  N=num;
  for(n=0; n<N; n++) {
    v=vlist[n];
    /* add vertex
       fprintf(stdout,"  %i index -> %i nei -- num %i \n",n,v->nei,num); */
    for(nei=0; nei<v->nei; nei++) {
      vlist[num++]=v->neivert[nei];
    }
  }
  /* check if it's in the list
     fprintf(stdout,"  check if it's in the list\n");*/
  for(n=0; n<num; n++) {
    for(m=n+1; m<num; m++) {
      if(vlist[n]==vlist[m]) {
        vlist[m]=vlist[--num];
        m--;
      }
    }
  }
  return num;
}

int off_get_local_kvert_check(kvert **vlist,int vnum,int num)
/* -similar to off_get_local_kvert but checks for the number of vlist =vnum
 * -returns negative number for error
 * -on successful return, number of vertices in the list is returned
 */
{
  int nei,n,m,N,verbose=0;
  kvert *v;
  N=num;
  if(verbose>1) fprintf(stdout,"  NUM = %i\n",num);
  for(n=0; n<N; n++) {
    v=vlist[n];
    for(nei=0; nei<v->nei; nei++) {
      if(num>=vnum) {
        fprintf(stdout,"[%s] ERROR: overflow\n",__func__);
        return -1;
      }
      if(verbose)   {
        fprintf(stdout,"  n = %i   num = %i :  nei = %i\n",n,num,nei);
      }
      vlist[num++]=v->neivert[nei];
    }
  }
  /* check if it's in the list
     fprintf(stdout,"  check if it's in the list\n");*/
  for(n=0; n<num; n++) {
    for(m=n+1; m<num; m++) {
      if(vlist[n]==vlist[m]) {
        vlist[m]=vlist[--num];
        m--;
      }
    }
  }
  return num;
}

kvert **off_get_local_kvert_list(kvert *vert,int *num,int local_num)
/* -finds local vertices
 * -uses the local connectivity distance of local_num (in vertex connection, not milimeter)
 * -returns the list of vertices
 * -on output, num will be the number of vertices */
{
  kvert **vlist,**outlist;
  int
  n,m,nv;
  if(vert==NULL) {
    fprintf(stderr,"ERROR: vert is null for off_get_local_kvert_list\n");
    return NULL;
  }
  if(local_num<=0) {
    fprintf(stderr,"ERROR: local_num is undefined %i\n",local_num);
    return NULL;
  }
  for(n=0,nv=1; n<local_num; n++) {
    nv *= 8;
  }
  do {
    vlist = (kvert **)calloc(nv,sizeof(kvert *));
    for(n=0,m=1,vlist[0]=vert; n<local_num; n++) {
      if(m>0)
        m = off_get_local_kvert_check(vlist,nv,m);
    }
    if(m<0) {
      free(vlist);
      nv*=2;
    }
  } while (m<0);
  *num = m;
  outlist=(kvert **)calloc(m,sizeof(kvert *));
  for(n=0; n<m; n++) {
    outlist[n]=vlist[n];
  }
  free(vlist);
  return outlist;
} /* off_get_local_kvert_list */

int off_get_local_kface(kface **flist,int num) {
  int n,m;
  kedge *e;
  for(n=0; n<num; n++) {
    for(m=0; m<3; m++) {
      e=flist[n]->edge[m];
      if(e->adjface[0]==flist[n])
        flist[num++] = e->adjface[1];
      else
        flist[num++] = e->adjface[0];
    }
  }
  /* check if it's in the list */
  for(n=0; n<num; n++) {
    for(m=n+1; m<num; m++) {
      if(flist[n]==flist[m]) {
        flist[m]=flist[--num];
        m--;
      }
    }
  }
  return num;
}


niikpt niikpt_kvert_local_average(kvert *v,int n)
/* simple local average (no area-weighting)
 *
 * see also
 *   niikpt niikpt_kvert_local_average(kvert *v,int n);
 *   -with area-weighting
 */
{
  int m,num;
  kvert *vlist[4096];
  niikpt p;
  vlist[0]=v;
  num=1;
  for(m=0; m<n; m++) {
    if((num=off_get_local_kvert_check(vlist,4096,num))<0) {
      fprintf(stderr,"ERROR: off_get_local_kvert_check(vlist,4096,num)\n");
      return niikpt_problem();
    }
  }
  p = niikpt_zero();
  for(m=0; m<num; m++) {
    p = niikpt_add(p,vlist[m]->v);
  }
  p = niikpt_kmul(p,1.0/num);
  return p;
} /* niikpt_kvert_local_average */

niikpt niikpt_kvert_simple_local_avg(kvert *v) {
  int n;
  niikpt g=niikpt_zero();
  for(n=0; n<v->nei; n++) {
    g=niikpt_add(g,v->neivert[n]->v);
  }
  return niikpt_kmul(g,1.0/(double)v->nei);
}


niikpt niikpt_kvert_local_average2(kvert *v) {
  return niikpt_kvert_local_average_with_tangential_relaxation(v,1);
}

niikpt niikpt_kvert_local_average_with_tangential_relaxation(kvert *v,int use_tangential_relaxation)
/* -based on niikpt_kvert_local_average2
 * -area-weighted average as in remesh function
 * -choice of tangential relaxation */
{
  const double lambda = 0.95;
  niikpt g,p;
  double a;
  int n;
  g=niikpt_zero();
  for(n=0; n<v->nei; n++) {
    a = niikpt_area2(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    p = niikpt_avg3 (v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    g.x += a * p.x;
    g.y += a * p.y;
    g.z += a * p.z;
    g.w += a;
  }
  g.x /= g.w;
  g.y /= g.w;
  g.z /= g.w;
  if(!use_tangential_relaxation) return g;
  g = niikpt_move_normal(g,v->normal,lambda*niikpt_dot(v->normal,niikpt_sub(v->v,g)));
  return g;
} /* niikpt_kvert_local_average_with_tangential_relaxation */


niikpt niikpt_kvert_local_average_with_tangential_relaxation2(kvert *v,double lambda)
/* -based on niikpt_kvert_local_average2
 * -area-weighted average as in remesh function
 * -choice of tangential relaxation */
{
  niikpt g,p;
  double a;
  int n;
  g=niikpt_zero();
  for(n=0; n<v->nei; n++) {
    a = niikpt_area2(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    p = niikpt_avg3 (v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    g.x += a * p.x;
    g.y += a * p.y;
    g.z += a * p.z;
    g.w += a;
  }
  g.x /= g.w;
  g.y /= g.w;
  g.z /= g.w;
  if(lambda<0.0) return g;
  g = niikpt_move_normal(g,v->normal,lambda*niikpt_dot(v->normal,niikpt_sub(v->v,g)));
  return g;
} /* niikpt_kvert_local_average_with_tangential_relaxation2 */


/****************************************************************
 *
 * surface smoothing
 *
 * -various surface smoothing for variable var
 * -simplex averaging, median filtering, and trimmed averaging
 *
 * int off_surface_smooth_using_vert(kobj *obj,double *var,int vnei,double wself);
 *
 ****************************************************************/
int off_surface_smooth_using_vert(kobj *obj, double *var, int vnei, double wself)
/* -surface smoothing for var
 * -var is a obj->nvert-length vector
 * -vnei is the number of neighborhood usually around 1-4
 * -wself is the weighting for the vertex itself
 */
{
  kvert *v,*vlist[65536];
  double *tmpvar;
  int 
    verbose=0,
    vindex,m,n;

  if(verbose) fprintf(stdout,"[off_surface_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (double *)calloc(obj->nvert,sizeof(double));
  if(verbose) fprintf(stdout,"[off_surface_smooth_using_vert] vertex loop\n");
  for(v=obj->vert,vindex=1; v!=NULL; v=v->next,vindex++) {
    v->index=vindex;
  }
  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    tmpvar[vindex] = wself * var[vindex];
    vlist[0]=v;
    n=1;
    for(m=0; m<vnei; m++) {
      if((n=off_get_local_kvert_check(vlist,65536,n))<0) {
        fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
        return 0;
      }
    }
    for(m=0; m<n; m++) {
      tmpvar[vindex] += var[vlist[m]->index-1];
    }
    tmpvar[vindex] /= n+wself;
  }
  if(verbose) fprintf(stdout,"[off_surface_smooth_using_vert] update var\n");
  for(n=0; n<obj->nvert; n++) {
    var[n]=tmpvar[n];
  }
  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_smooth_using_vert] complete\n");
  return 1;
}

int off_surface_median_smooth_using_vert(kobj *obj,double *var,int vnei)
/* -surface median smoothing for var
 * -var is a obj->nvert-length vector
 * -vnei is the number of neighborhood usually around 1-4
 */
{
  kvert *v,*vlist[65536];
  double *tmpvar,dvlist[65536];
  int
  verbose=0,
  vindex,m,n;
  if(verbose) fprintf(stdout,"[off_surface_median_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (double *)calloc(obj->nvert,sizeof(double));
  if(verbose) fprintf(stdout,"[off_surface_median_smooth_using_vert] vertex loop\n");
  for(v=obj->vert,vindex=1; v!=NULL; v=v->next,vindex++) {
    v->index=vindex;
  }
  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    tmpvar[vindex] = var[vindex];
    for(m=0,vlist[0]=v,n=1; m<vnei; m++) {
      if((n=off_get_local_kvert_check(vlist,65536,n))<0) {
        fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
        return 0;
      }
    }
    for(m=0; m<n; m++) {
      dvlist[m]=var[vlist[m]->index-1];
    }
    tmpvar[vindex] = niik_median_quicksort_double(dvlist,n); /* median calculation */
  }
  if(verbose) fprintf(stdout,"[off_surface_median_smooth_using_vert] update var\n");
  for(n=0; n<obj->nvert; n++) {
    var[n]=tmpvar[n];
  }
  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_median_smooth_using_vert] complete\n");
  return 1;
}

int off_surface_trimmed_average_smooth_using_vert(kobj *obj,double *var,int vnei,double trim)
/* -surface smoothing (trimmed average) for var
 * -var is a obj->nvert-length vector
 * -vnei is the number of neighborhood usually around 1-4
 */
{
  kvert *v,*vlist[65536];
  double *tmpvar,dvlist[65536];
  int
    verbose=0,
    vindex,m,n;
  if(verbose) fprintf(stdout,"[off_surface_trimmed_average_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (double *)calloc(obj->nvert,sizeof(double));
  if(verbose) fprintf(stdout,"[off_surface_trimmed_average_smooth_using_vert] vertex loop\n");
  for(v=obj->vert,vindex=1; v!=NULL; v=v->next,vindex++) {
    v->index=vindex;
  }
  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    tmpvar[vindex] = var[vindex];
    for(m=0,n=1,vlist[0]=v; m<vnei; m++) {
      if((n=off_get_local_kvert_check(vlist,65536,n))<0) {
        fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
        return 0;
      }
    }
    for(m=0; m<n; m++) {
      dvlist[m]=var[vlist[m]->index-1];
    }
    if(!niik_get_trimmed_average_from_double_vector(dvlist,n,trim,&tmpvar[vindex])) {
      fprintf(stderr,"ERROR: niik_get_trimmed_average_from_double_vector\n");
      return 0;
    }
  }
  if(verbose) fprintf(stdout,"[off_surface_trimmed_average_smooth_using_vert] update var\n");
  for(n=0; n<obj->nvert; n++) {
    var[n]=tmpvar[n];
  }
  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_trimmed_average_smooth_using_vert] complete\n");
  return 1;
}

/****************************************************************
 *
 * surface 3D variable smoothing, itaratively apply gaussian kernel using only nearest neighbours
 *
 * -simple gaussian filtering
 *
 * int off_surface_smooth_using_vert(kobj *obj,niikpt *var,double sigma,int iter);
 *
 ****************************************************************/
int off_surface_gauss_smooth_using_vert(kobj *obj,  double *  var, double sigma,int iter)
/* -surface smoothing for var
 * -var is a obj->nvert-length vector
 * -sigma - smoothing kernel 
 * -iter - number of iterations
 */
{
  kvert *v;
  double *  tmpvar;
  int 
    verbose=0,
    vindex,m,n,it;

  if(verbose) fprintf(stdout,"[off_surface_gauss_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (double *)calloc(obj->nvert,sizeof(double));
  if(verbose) fprintf(stdout,"[off_surface_gauss_smooth_using_vert] vertex loop\n");

  for(it=0;it<iter;it++) {

    for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
      double total_w = 1.0;
      /*DEBUG*/
      if( (v->index-1) != vindex)
        abort();
      /**/
      tmpvar[vindex] = var[vindex];

#if _OPENMP>=201307
    #pragma omp simd
#endif
      for(n=0; n<v->nei; n++) {
        double dist2 = niikpt_distance2(v->v, v->neivert[n]->v);
        double w = exp( -dist2/(2*sigma*sigma) );
        total_w += w;

        tmpvar[vindex] += w * var[ v->neivert[n]->index-1 ];
      }
      tmpvar[vindex] /= total_w;
    }

#if _OPENMP>=201307
    #pragma omp simd
#endif
    for(n=0; n<obj->nvert; n++) {
      var[n] = tmpvar[n];
    }
  }

  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] complete\n");
  return 1;
}


/****************************************************************
 *
 * surface 3D variable smoothing, itaratively apply gaussian kernel using only nearest neighbours
 *
 * -simple gaussian filtering
 *
 * int off_surface_smooth_using_vert(kobj *obj,niikpt *var,double sigma,int iter);
 *
 ****************************************************************/
int off_surface_gauss_smooth_using_vert_with_mask(kobj *obj,  double *  var, double sigma,int iter,unsigned char *mask)
/* -surface smoothing for var
 * -var is a obj->nvert-length vector
 * -sigma - smoothing kernel 
 * -iter - number of iterations
 */
{
  kvert *v;
  double *  tmpvar;
  int 
    verbose=0,
    vindex,m,n,it;

  if(verbose) fprintf(stdout,"[off_surface_gauss_smooth_using_vert_with_mask] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (double *)calloc(obj->nvert,sizeof(double));
  if(verbose) fprintf(stdout,"[off_surface_gauss_smooth_using_vert_with_mask] vertex loop\n");

  for(it=0;it<iter;it++) {

    for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
      double total_w = 1.0;
      /*DEBUG*/
      if( (v->index-1) != vindex)
        abort();
      /**/
      tmpvar[vindex] = var[vindex];

#if _OPENMP>=201307
    #pragma omp simd
#endif
      for(n=0; n<v->nei; n++) {
        double dist2 = niikpt_distance2(v->v, v->neivert[n]->v);
        double w = exp( -dist2/(2*sigma*sigma) );

        if(mask[v->neivert[n]->index-1]>0)
        {
          total_w += w;

          tmpvar[vindex] += w * var[ v->neivert[n]->index-1 ];
        }
      }
      tmpvar[vindex] /= total_w;
    }

#if _OPENMP>=201307
    #pragma omp simd
#endif
    for(n=0; n<obj->nvert; n++) {
      var[n] = tmpvar[n];
    }
  }

  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_gauss_smooth_using_vert_with_mask] complete\n");
  return 1;
}



/****************************************************************
 *
 * surface 3D field smoothing
 *
 * -simple gaussian filtering
 *
 * int off_surface_smooth_using_vert(kobj *obj,niikpt *var,double sigma);
 *
 ****************************************************************/
int off_surface_field_smooth_using_vert(kobj *obj, niikpt *var, double sigma)
/* -surface smoothing for var
 * -var is a obj->nvert-length vector
 * -vnei is the number of neighborhood usually around 1-4
 * -wself is the weighting for the vertex itself
 */
{
  kvert *v;
  niikpt *tmpvar;
  int 
    verbose=0,
    vindex,m,n;

  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (niikpt *)calloc(obj->nvert,sizeof(niikpt));
  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] vertex loop\n");

  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    double total_w = 1.0;
    /*DEBUG*/
    if( (v->index-1) != vindex)
      abort();
    /**/
    tmpvar[vindex] = var[vindex];
  
    for(n=0; n<v->nei; n++) {
      double dist2 = niikpt_distance2(v->v, v->neivert[n]->v);
      double w = exp(-dist2/(2*sigma*sigma));
      total_w += w;

      tmpvar[vindex].x += w * var[v->neivert[n]->index-1].x;
      tmpvar[vindex].y += w * var[v->neivert[n]->index-1].y;
      tmpvar[vindex].z += w * var[v->neivert[n]->index-1].z;
    }
    tmpvar[vindex].x /= total_w;
    tmpvar[vindex].y /= total_w;
    tmpvar[vindex].z /= total_w;
  }

  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] update var\n");
  for(n=0; n<obj->nvert; n++) {
    var[n]=tmpvar[n];
  }
  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] complete\n");
  return 1;
}


int off_surface_field_smooth_using_vert_thresholded(kobj *obj, niikpt *var, double sigma,double threshold)
/* -surface smoothing for var
 * -var is a obj->nvert-length vector
 * -vnei is the number of neighborhood usually around 1-4
 * -wself is the weighting for the vertex itself
 */
{
  kvert *v;
  niikpt *tmpvar;
  int 
    verbose=0,
    vindex,m,n;

  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] start\n");
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  tmpvar = (niikpt *)calloc(obj->nvert,sizeof(niikpt));
  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] vertex loop\n");

  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    double total_w = 1.0;
    /*DEBUG*/
    if( (v->index-1) != vindex)
      abort();
    /**/
    tmpvar[vindex] = var[vindex];
  
    for(n=0; n<v->nei; n++) {

      if(
        fabs(var[v->neivert[n]->index-1].x)<threshold &&
        fabs(var[v->neivert[n]->index-1].y)<threshold &&
        fabs(var[v->neivert[n]->index-1].z)<threshold
       )
        continue;

      double dist2 = niikpt_distance2(v->v, v->neivert[n]->v);
      double w = exp(-dist2/(2*sigma*sigma));
      total_w += w;

      tmpvar[vindex].x += w * var[v->neivert[n]->index-1].x;
      tmpvar[vindex].y += w * var[v->neivert[n]->index-1].y;
      tmpvar[vindex].z += w * var[v->neivert[n]->index-1].z;
    }
    tmpvar[vindex].x /= total_w;
    tmpvar[vindex].y /= total_w;
    tmpvar[vindex].z /= total_w;
  }

  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] update var\n");
  for(n=0; n<obj->nvert; n++) {
    var[n]=tmpvar[n];
  }
  free(tmpvar);
  if(verbose) fprintf(stdout,"[off_surface_field_smooth_using_vert] complete\n");
  return 1;
}


/***************************************
 *
 * icosahedron
 *
 ***************************************/

kobj *off_make_sphere_from_icosahedron(double elen, double radius, niikpt ctr)
/* makes a sphere from icosahedron
 * --added spherical coordinates (checking the code), knakamura 2014-06-16 */
{
  kobj *obj;
  kvert *v;
  double  dfrac, meanlen;
  const int verbose=1;
  if(verbose>0) niik_fc_display(__func__,1);
  if((obj=off_create_icosahedron())==NULL) {
    fprintf(stderr,"ERROR: off_create_icosahedron\n");
    return NULL;
  }
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v = niikpt_kmul(niikpt_unit(v->v),radius);
  }
  meanlen = off_get_kobj_mean_edge_length(obj);
  if(verbose>0) fprintf(stdout,"[%s] mean elen = %9.5f \n",__func__,meanlen);
  dfrac = (off_get_kobj_mean_edge_length(obj) - elen)/elen;
  while(dfrac>1.0) {
    if(!off_subdiv_kobj(obj)) {
      fprintf(stderr,"ERROR: off_subdiv_kobj\n");
      return NULL;
    }
    if(verbose) fprintf(stdout,"[%s] subdiv vfe = %i %i %i\n",__func__,obj->nvert,obj->nface,obj->nedge);
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v = niikpt_kmul(niikpt_unit(v->v),radius);
    }
    meanlen = off_get_kobj_mean_edge_length(obj);
    dfrac = (meanlen - elen)/elen;
    if(verbose) fprintf(stdout,"[%s] mean elen = %9.5f \n",__func__,meanlen);
  }
  /* image center and radius */
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v = niikpt_move_normal(ctr,niikpt_unit(v->v),radius);
  }
  /* remesh */
  if(verbose>0) fprintf(stdout,"[%s] remeshed %7.3f\n",__func__,elen);
  if(!off_remesh_kobj(obj,elen,10,0)) {
    fprintf(stderr,"ERROR: off_remesh_kobj \n");
    return NULL;
  }
  meanlen = off_get_kobj_mean_edge_length(obj);
  if(verbose>0) {
    fprintf(stdout,"[%s] sphere from icosahedron: vfe = %i %i %i   elen_avg %9.4f\n",__func__,
            obj->nvert,obj->nface,obj->nedge,meanlen);
  }
  if(verbose>0)niik_fc_display(__func__,0);
  return obj;
}


int niik_affine_transform_off(kobj *obj,niikmat *afmat)
/* -applies affine transform to obj
 * -updates vertex only (not normal or anything else) */
{
  kvert *v;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(afmat==NULL) {
    fprintf(stderr,"ERROR: afmat is null\n");
    return 0;
  }
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v = niikpt_affine_transform(afmat,v->v);
  }
  return 1;
}



/****************************************************************
 *
 * curvature map
 *
 ****************************************************************/

#define niik_off_curvature_vert_vmax 65535

int niik_off_curvature_vert_sph(kvert *v,int nei,double *out) {
  kvert *vlist[niik_off_curvature_vert_vmax];
  int m,n;
  niikpt ctr;
  double dsum;
  for(m=0,n=1,vlist[0]=v; m<nei; m++) {
    if((n=off_get_local_kvert_check(vlist,niik_off_curvature_vert_vmax,n))<0) {
      fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
      return 0;
    }
  }
  for(m=0,dsum=0,ctr=niikpt_zero(); m<n; m++) {
    ctr=niikpt_add(ctr,vlist[m]->v);
    dsum+=fabs(vlist[m]->sph.the-v->sph.the);
    dsum+=fabs(vlist[m]->sph.psi-v->sph.psi);
  }
  ctr=niikpt_kmul(ctr,1.0/(double)n);
  *out=niikpt_dot(niikpt_sub(v->v,ctr),v->normal)/dsum*NIIK_PI;
  return 1;
}

int niik_off_curvature_vert(kvert *v,int nei,double *out) {
  kvert *vlist[niik_off_curvature_vert_vmax];
  int m,n;
  niikpt ctr;
  for(m=0,n=1,vlist[0]=v; m<nei; m++) {
    if((n=off_get_local_kvert_check(vlist,niik_off_curvature_vert_vmax,n))<0) {
      fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
      return 0;
    }
  }
  for(m=0,ctr=niikpt_zero(); m<n; m++) {
    ctr=niikpt_add(ctr,vlist[m]->v);
  }
  ctr=niikpt_kmul(ctr,1.0/(double)n);
  *out=niikpt_dot(niikpt_sub(v->v,ctr),v->normal);
  return 1;
}

int niik_off_curvature_vert_meek2000(kvert *v,double *out) {
  int nei;
  double T=0,S=0;
  for(nei=0; nei<v->nei; nei++) {
    T+=niikpt_angle_between_vectors(niikpt_sub(v->neivert[nei]->v,v->v),niikpt_sub(v->neivert[(nei+1)%(v->nei)]->v,v->v));
    S+=niikpt_area2(v->neiface[nei]->vert[0]->v,v->neiface[nei]->vert[1]->v,v->neiface[nei]->vert[2]->v);
  }
  *out=(NIIK_PI2 - NIIK_DEGREE2RAD(T)) / (S/6.0);
  //fprintf(stdout,"%9i  %.5f %.5f\n",v->index,T,S);
  return 1;
}

#undef niik_off_curvature_vert_vmax

int niik_off_curvature_map_update(kobj *obj,niikvec *cm,int nei)
/* calculates curvature map
 * -obj needs updated normal for vertices
 */
{
  kvert *v,*vlist[65536];
  int m,n,vindex;
  niikpt ctr;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(cm==NULL) {
    fprintf(stderr,"ERROR: cm is null\n");
    return 0;
  }
  if(cm->num!=obj->nvert) {
    fprintf(stderr,"ERROR: num is different (obj %i, cm %i)\n",obj->nvert,cm->num);
    return 0;
  }
  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {

    for(m=0,n=1,vlist[0]=v; m<nei; m++) {
      if((n=off_get_local_kvert_check(vlist,65536,n))<0) {
        fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
        return 0;
      }
    }

    for(m=0,ctr=niikpt_zero(); m<n; m++) {
      ctr=niikpt_add(ctr,vlist[m]->v);
    }

    ctr=niikpt_kmul(ctr,1.0/(double)n);
    cm->v[vindex]=niikpt_dot(niikpt_sub(v->v,ctr), v->normal);
  }
  return 1;
}

int niik_off_apply_surface_smoothing(kobj *obj,int nei,double delta)
/* smoothes the surface
 * -obj needs updated normal for vertices
 */
{
  kvert *v,*vlist[65536];
  int m,n,vindex;
  niikpt ctr,*plist;
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  plist=(niikpt *)calloc(obj->nvert,sizeof(niikpt));

  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    for(m=0,n=1,vlist[0]=v; m<nei; m++) {
      if((n=off_get_local_kvert_check(vlist,65536,n))<0) {
        fprintf(stderr,"ERROR: off_get_local_kvert_check\n");
        return 0;
      }
    }
    for(m=0,ctr=niikpt_zero(); m<n; m++) {
      ctr=niikpt_add(ctr,vlist[m]->v);
    }

    ctr=niikpt_kmul(ctr,1.0/(double)n);
    plist[vindex]=niikpt_wavg(v->v, ctr, delta);
  }
  for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
    v->v=plist[vindex];
  }

  free(plist);
  return 1;
}


kvert *niik_off_find_closest_vertex(kobj *obj,niikpt p) {
  kvert *v,*vc;
  double dval,dmin=1e9;
  for(v=vc=obj->vert; v!=NULL; v=v->next) {
    dval=niikpt_distance(v->v,p);
    if(dmin>dval) {
      dmin=dval;
      vc=v;
    }
  }
  return vc;
}



/**
 * Create a rectangle object based on eight points
 *
 */
kobj *off_create_rect(const niikpt *p) {
  kobj *obj;
  kvert *v,*vlist[8];

  int fvlist[36] = { 0,1,2,  2,3,0,  0,3,4, 7,4,3, 3,2,7, 6,7,2, 1,5,2, 6,2,5, 4,7,5, 6,5,7, 0,4,1, 5,1,4  };

  kface *f;
  kedge *e;
  int n,m;
  const char *fcname="off_create_rect";

  obj = off_obj_init();
  obj->vert = v = off_vert_init();
  v->v=p[0];  /* index=0 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[1];   /* index=1 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[2];   /* index=2 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[3];   /* index=3 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[4];  /* index=0 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[5];   /* index=1 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[6];   /* index=2 */

  v->next = off_vert_init();
  v->next->prev=v;
  v = v->next;
  v->v=p[7];   /* index=3 */

  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) vlist[n]=v;

  /* added spherical coordinates */
  for(n=0; n<8; n++) {
    vlist[n]->sph=niikxyz_sph(vlist[n]->v);
  }
  obj->spherecoo=1;

  m=0;
  obj->face = f = off_face_init();
  for(n=0; n<3; n++,m++) f->vert[n]=vlist[ fvlist[m] ];

  while(m<36) {
    f->next = off_face_init();
    f->next->prev = f;
    f=f->next;
    for(n=0; n<3; n++,m++) f->vert[n]=vlist[fvlist[m]];
  }
  NIIK_RET0((!off_kobj_read_off_make_edge(obj)),fcname,"off_kobj_read_off_make_edge");

  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[0],e)),fcname,"off_kobj_update_vertex_nei");
    if(!e->endpts[1]->nei)
      NIIK_RET0((!off_kobj_update_vertex_nei(e->endpts[1],e)),fcname,"off_kobj_update_vertex_nei");
  }
  off_kobj_update_num(obj);

  return obj;
}

int off_kobj_add_comment(kobj *obj,const char *comment)
{
  const char *fcname="off_kobj_add_comment";
  char *pos;
  obj->n_comments += 1;
  NIIK_RET0(((obj->comment = (char**)realloc(obj->comment,sizeof(char*)*obj->n_comments)) ==NULL),fcname,"realloc");  
  NIIK_RET0(((obj->comment[obj->n_comments-1] = strdup(comment)) ==NULL),fcname,"strdup");

  /*strip newlines*/
  if ((pos=strchr(obj->comment[obj->n_comments-1], '\n')) != NULL)
    *pos = '\0';

  return 1;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
