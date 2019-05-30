/* FILENAME:     off2mniobj.c
 * DESCRIPTION:  Kunio's conversion from off to mni's obj
 * AUTHOR:       Kunio Nakamura
 *
 *
 * reference
 * http://www.bic.mni.mcgill.ca/users/david/FAQ/polygons_format.txt
 */

#include "falcon.h"
#include "falcon_surfaces.h"

#ifdef HAVE_BICPL
#include  <bicpl.h>
#endif

int niik_off_write_mniobj(char *fname, kobj *obj);


static char *prog_version[] = {
  NIIK_VERSION "\n"
  "  program history\n"
  "  0.0.0   : August 25, 2014, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  off2mniobj:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"obj2mniobj\n");
  fprintf(stdout,"  usage: [options] <in.off> <out.obj>\n");
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int
  nc,sc,NC=3;
  const char
  *fcname="off2mniobj";

  niik_fc_display(fcname,1);

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

  NIIK_EXIT(((obj=off_kobj_read_offply(argv[1]))==NULL),fcname,"niik_kobj_read_off",1);
  NIIK_EXIT((!niik_off_write_mniobj(argv[2],obj)),fcname,"niik_off_wirte_mniobj",1);

  obj=off_kobj_free(obj);

  niik_fc_display(fcname,0);
  exit(0);
} /* off2mniobj */

int niik_off_write_mniobj(char *fname,kobj *obj) {
  const char *fcname=__func__;
#ifdef HAVE_BICPL
  FILE *fp=NULL;
  int
  i,
  verbose=niik_verbose();
  kvert *v;
  kface *f;
  object_struct   *object_list;
  polygons_struct      *polygons;

  niik_fc_display(fcname,1);

  /*recalculate normals ?*/
  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);


  object_list=create_object(POLYGONS);
  polygons = get_polygons_ptr(object_list);

  initialize_polygons_with_size(polygons,
                                make_rgba_Colour(200,100,100,255),  /*reddish?*/
                                NULL,
                                obj->nvert,obj->nface,3);

  i=0;
  for(v=obj->vert; v!=NULL; v=v->next) {
    /*assign XYZ coords*/
    fill_Point(polygons->points[i],v->v.x,v->v.y,v->v.z);
    fill_Point(polygons->normals[i],v->normal.x,v->normal.y,v->normal.z);
    /*assign normals?*/
    i++;
  }

  /*assign indeces*/
  i=0;
  for(f=obj->face; f!=NULL; f=f->next) {
    /*f->vert[0]->index-1, f->vert[1]->index-1, f->vert[2]->index-1);*/
    polygons->indices[i*3]=f->vert[0]->index-1;
    polygons->indices[i*3+1]=f->vert[1]->index-1;
    polygons->indices[i*3+2]=f->vert[2]->index-1;
    i++;
  }
  /*assume that end_indeces are properly set already*/

  if(obj->color) {
    /*assign colours*/
    if(verbose>1) fprintf(stdout,"Outputting per item colour\n");
    set_polygon_per_item_colours(polygons);
    i=0;
    for(f=obj->face; f!=NULL; f=f->next) {
      if(f->color)
        polygons->colours[i]=make_rgba_Colour(f->color[0]*255,f->color[1]*255,f->color[2]*255,255);
      i++;
    }
  }

  output_graphics_file( fname, ASCII_FORMAT, 1, &object_list );

  niik_fc_display(fcname,0);
  return 1;
#else /*HAVE_BICPL*/
  FILE *fp=NULL;
  int
  i,
  verbose=niik_verbose();
  kvert *v;
  kface *f;

  niik_fc_display(fcname,1);

  fprintf(stdout,"WARNING: using implementation without BICPL, broken obj file might be created!\n");

  if(verbose>1) {
    fprintf(stdout,"[%s] opening %s\n",fcname,fname);
  }

  /*recalculate normals*/
  off_update_kobj_vert_normal(obj);

  if((fp=fopen(fname,"w"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",fname);
    return 0;
  }
  fprintf(fp,"P 0.3 0.7 0.5 100 1 %i\n",obj->nvert);

  /*output X Y Z coordinates of vertices*/
  for(v=obj->vert; v!=NULL; v=v->next) {
    fprintf(fp,"%15.15f %15.15f %15.15f\n",
            v->v.x,v->v.y,v->v.z);
  }

  /*output X Y Z values of normals, one per vertex*/
  for(v=obj->vert; v!=NULL; v=v->next) {
    fprintf(fp,"%15.15f %15.15f %15.15f\n",
            v->normal.x,v->normal.y,v->normal.z);
  }

  fprintf(fp,"%i\n",obj->nface);
  /*if(obj->color==0) */
  {
    /* COLOR IS NOT SUPPORTED AT THE MOMENT
     * -this is because Kunio's off object has color associated with triangles not vertices
     * -to do, assign average color to points
     */
    fprintf(fp,"2 ");
    /*0.7 0.7 0.7 1\n"); Assign same colour*/
    for(v=obj->vert; v!=NULL; v=v->next) {
      fprintf(fp,"%15.15f %15.15f %15.15f %15.5f\n",
              0.7,0.7,0.7,1.0);
    }
  }

  /*write indeces of faces*/
  for(i=0; i<obj->nface; i++)
    fprintf(fp,"%i ", i*3);

  fprintf(fp,"\n");
  /*write faces*/
  for(f=obj->face; f!=NULL; f=f->next)
    fprintf(fp,"%i %i %i ",
            f->vert[0]->index-1, f->vert[1]->index-1, f->vert[2]->index-1);
  fprintf(fp,"\n");

  fclose(fp);

  niik_fc_display(fcname,0);
  return 1;
#endif
} /* niik_off_off2mniobj */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
