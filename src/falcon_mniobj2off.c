/* FILENAME:     mniobj2off.c
 * DESCRIPTION:  Kunio's conversion from obj to off
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_surfaces.h"

#if 0
#ifdef HAVE_BICPL
#include  <bicpl.h>
#endif
#endif 

static char *prog_version[] = {
  "  program history\n"
  "  0.0.0   : August 25, 2014, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  falcon_obj2off:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"falcon_obj2off\n");
  fprintf(stdout,"  usage: [options] <in.obj> <out.off/.ply>\n");
}

int niik_off_mniobj2off_files(char *iname,char *oname);
kobj *niik_off_mniobj2off(char *iname);

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int
    nc,sc,NC=3;
  const char
    *fcname="falcon_obj2off";
  char* timestamp=niik_create_minc_timestamp(argc,argv);
  int verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc>NC) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(1);
  } else if(argc<NC) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(1);
  }

  /*NIIK_EXIT((!niik_off_mniobj2off_files(argv[1],argv[2])),fcname,"niik_off_mniobj2off_files",1);*/
  NIIK_EXIT(((obj=niik_off_mniobj2off(argv[1]))==NULL),fcname,"niik_off_mniobj2off",1);
  off_kobj_add_comment(obj,timestamp);
  off_kobj_write_offply(argv[2],obj,0);
  
  if(verbose>0) niik_fc_display(fcname,0);
  exit(0);
} /* mniobj2off */


kobj *niik_off_mniobj2off(char *iname) {
  const char *fcname=__func__;
  kobj *obj=NULL;
  kvert **vlist=NULL;
  kface **flist=NULL;
  int verbose=0;
  FILE *ifp=NULL;
  char line1[1024];
  double
    amb,diff,spec1,spec2,opacity;
  int
    color_code,
    n_triangles=0,
    n,tri[3],
    i,n_points;
  char desc;
    niikpt
    objcolor,
    *colorlist=NULL;

  niik_fc_display(fcname,1);
  if(verbose>1) {
    fprintf(stdout,"[%s] opening %s\n",fcname,iname);
  }

  if((ifp=fopen(iname,"r"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",iname);
    return 0;
  }

  NIIK_RET0(((fgets(line1,1024,ifp))==NULL),fcname,"fgets problem");
  if(verbose>5) fprintf(stdout,"line 1 = %s\n",line1);
  sscanf(line1,"%c %lf %lf %lf %lf %lf %i\n",&desc,&amb,&diff,&spec1,&spec2,&opacity,&n_points);

  if(verbose>4) {
    fprintf(stdout,"  desc    %c\n",desc);
    fprintf(stdout,"  amb     %f\n",amb);
    fprintf(stdout,"  diff    %f\n",diff);
    fprintf(stdout,"  spec    %f\n",spec1);
    fprintf(stdout,"  spec    %f\n",spec2);
    fprintf(stdout,"  opac    %f\n",opacity);
    fprintf(stdout,"  #pts    %i\n",n_points);
  }

  obj=off_obj_init();

  if(verbose>4) {
    fprintf(stdout,"[%s] reading points, %i\n",fcname,n_points);
  }
  NIIK_RET0(((vlist=off_vert_list(n_points))==NULL),fcname,"off_vert_list");
  for(i=0; i<n_points; i++) {
    vlist[i]=off_vert_init();
    if((fscanf(ifp,"%lf %lf %lf\n",&vlist[i]->v.x,&vlist[i]->v.y,&vlist[i]->v.z))!=3) {
      fprintf(stderr,"[%s] ERROR: point %i\n",fcname,i);
      return 0;
    }
    vlist[i]->index=i+1;
  }

  if(verbose>4) {
    fprintf(stdout,"[%s] reading normals, %i\n",fcname,n_points);
  }
  for(i=0; i<n_points; i++) {
    if((fscanf(ifp,"%lf %lf %lf\n",&vlist[i]->normal.x,&vlist[i]->normal.y,&vlist[i]->normal.z))!=3) {
      fprintf(stderr,"[%s] ERROR: normal %i\n",fcname,i);
      return 0;
    }
  }
  vlist[0]->next=vlist[1];
  for(i=1; i<n_points-1; i++) {
    vlist[i]->prev=vlist[i-1];
    vlist[i]->next=vlist[i+1];
  }
  vlist[n_points-1]->prev=vlist[n_points-2];

  if((fscanf(ifp,"%i",&n_triangles))!=1) {
    fprintf(stderr,"[%s] ERROR: n_triangles\n",fcname);
    return 0;
  }
  if(verbose>4) {
    fprintf(stdout,"[%s] #triangles = %i\n",fcname,n_triangles);
  }

  if(verbose>4) {
    fprintf(stdout,"[%s] reading colors\n",fcname);
  }
  if((fscanf(ifp,"%i",&color_code))!=1) {
    fprintf(stderr,"[%s] ERROR: color_code %i\n",fcname,color_code);
    return 0;
  }
  if(color_code==0) {
    if(verbose>4) {
      fprintf(stdout,"[%s]   color code is %i -> object color\n",fcname,color_code);
    }
    if((fscanf(ifp,"%lf %lf %lf %lf",&objcolor.x,&objcolor.y,&objcolor.z,&objcolor.w))!=4) {
      fprintf(stderr,"[%s] ERROR: color\n",fcname);
      return 0;
    }
  } else if(color_code==2) {
    if(verbose>4) {
      fprintf(stdout,"[%s]   color code is %i -> point color\n",fcname,color_code);
    }
    NIIK_RET0(((colorlist=niikpt_alloc(n_points))==NULL),fcname,"niikpt_alloc");
    for(i=0; i<n_points; i++) {
      colorlist[i].x=colorlist[i].y=colorlist[i].z=colorlist[i].w=0;
      if((fscanf(ifp,"%lf %lf %lf %lf\n",&colorlist[i].x,&colorlist[i].y,&colorlist[i].z,&colorlist[i].w))!=4) {
        fprintf(stderr,"[%s] ERROR: color %i\n",fcname,i);
        return 0;
      }
    }
  } else {
    fprintf(stdout,"[%s] ERROR: unknown color code %i\n",fcname,color_code);
  } /* color codes */

  /* I don't understand, but ignore */
  for(i=0; i<n_triangles; i++) {
    if((fscanf(ifp,"%i ",&n))!=1) {
      fprintf(stderr,"[%s] ERROR: triangle %i\n",fcname,i);
      return 0;
    }
  }

  if(verbose>4) {
    fprintf(stdout,"[%s] face list %i\n",fcname,n_triangles);
  }
  NIIK_RET0(((flist=off_face_list(n_triangles))==NULL),fcname,"off_vert_list");

  /* read triangles */
  for(i=0; i<n_triangles; i++) {
    if((fscanf(ifp,"%i %i %i ",&tri[0],&tri[1],&tri[2]))!=3) {
      fprintf(stderr,"[%s] ERROR: triangle %i\n",fcname,i);
      return 0;
    }
    //fprintf(stdout,"\t%9i %9i %9i\n",tri[0],tri[1],tri[2]);
    flist[i]=off_face_init();
    flist[i]->vert[0]=vlist[tri[0]];
    flist[i]->vert[1]=vlist[tri[1]];
    flist[i]->vert[2]=vlist[tri[2]];
    flist[i]->index=i+1;
  }
  if(verbose>4) {
    fprintf(stdout,"[%s] double linked list %i\n",fcname,n_triangles);
  }
  flist[0]->next=flist[1];
  for(i=1; i<n_triangles-1; i++) {
    flist[i]->prev=flist[i-1];
    flist[i]->next=flist[i+1];
  }
  flist[n_triangles-1]->prev=flist[n_triangles-2];
  fclose(ifp);

  obj->vert=vlist[0];
  obj->face=flist[0];
  obj->nface=n_triangles;
  obj->nvert=n_points;
  obj->color=0;
  obj->spherecoo=0;

  off_kobj_read_off_make_edge(obj);

  niik_fc_display(fcname,0);
  return obj;
} /* niik_off_mniobj2off */


/*Direct conversion code, without importing into internal data format*/
/*using BICPL library*/
#if 0
int niik_off_mniobj2off_files(char *iname,char *oname) 
{
  const char *fcname=__func__;
  int verbose=0;

  FILE *ofp=NULL;
  VIO_File_formats         format;
  object_struct        **object_list;
  int n_objects=0;

  double
    amb,diff,spec1,spec2,opacity;
  int
    color_code,
    n_triangles=0,
    n,tri[3],
    i,n_points;
  int idx;
  char desc;
    niikpt
    objcolor,
    *colorlist=NULL,
    *normlist=NULL,
    *plist=NULL;

  polygons_struct      *polygons;


  niik_fc_display(fcname,1);
  if(verbose>1) {
    fprintf(stdout,"[%s] opening %s\n",fcname,iname);
  }

  if( input_graphics_file( iname, &format, &n_objects,
                           &object_list ) != VIO_OK ) {
    fprintf(stderr,"ERROR: fopen %s\n",iname);
    return 0;
  }

  if( n_objects != 1) {
    fprintf(stderr,"File %s contains %d objects, only one is supported presently\n",iname,n_objects);
    return 0;
  }

  if(get_object_type( object_list[0] ) != POLYGONS ) {
    fprintf(stderr,"File %s contains unsupported object, only polygons are supported presently\n",iname);
    return 0;
  }

  polygons = get_polygons_ptr(object_list[0]);

  n_points=polygons->n_points;

  if(verbose>4) {
    fprintf(stdout,"[%s] reading points, %i\n",fcname,n_points);
  }
  NIIK_RET0(((plist=niikpt_alloc(n_points))==NULL),fcname,"niikpt_alloc");
  for(i=0; i<n_points; i++) {
    plist[i].x=polygons->points[i].coords[0];
    plist[i].y=polygons->points[i].coords[1];
    plist[i].z=polygons->points[i].coords[2];
    plist[i].w=0;
  }

  if(verbose>4) {
    fprintf(stdout,"[%s] reading normals, %i\n",fcname,n_points);
  }

  NIIK_RET0(((normlist=niikpt_alloc(n_points))==NULL),fcname,"niikpt_alloc");
  for(i=0; i<n_points; i++) {
    normlist[i].x=polygons->normals[i].coords[0];
    normlist[i].y=polygons->normals[i].coords[1];
    normlist[i].z=polygons->normals[i].coords[2];
    normlist[i].w=0;
  }
  n_triangles = polygons->n_items;

  if(verbose>4) {
    fprintf(stdout,"[%s] #triangles = %i\n",fcname,n_triangles);
  }

  if(polygons->colour_flag != ONE_COLOUR) {
    NIIK_RET0(((colorlist=niikpt_alloc(n_triangles))==NULL),fcname,"niikpt_alloc");
    for(i=0; i<n_triangles; i++) {
      colorlist[i].x=get_Colour_r_0_1(polygons->colours[i]);
      colorlist[i].y=get_Colour_g_0_1(polygons->colours[i]);
      colorlist[i].z=get_Colour_b_0_1(polygons->colours[i]);
      colorlist[i].w=0;
    }
  }

  /* write output file */
  if((ofp=fopen(oname,"w"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",oname);
    return 0;
  }
  fprintf(ofp,"OFF\n%i %i %i\n",n_points,n_triangles,0);
  for(i=0; i<n_points; i++) {
    fprintf(ofp,"%15.20f %15.20f %15.20f\n",plist[i].x,plist[i].y,plist[i].z);
  }

  /* read/write triangles */
  for(i=0,idx=0; i<n_triangles; i++) {
    int j,k;
    /*convert at most 3 points per face.... hack?*/
    for(j=idx,k=0; j<polygons->end_indices[i] && k<3; j++,k++) {
      tri[k]=polygons->indices[j];
    }
    idx=polygons->end_indices[i];

    fprintf(ofp,"3 %i %i %i\n",tri[0],tri[1],tri[2]);
  }

  fclose(ofp);
  delete_object_list( n_objects, object_list );

  niik_fc_display(fcname,0);
  return 1;

} /* niik_off_obj2off */
#endif
/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
 */
