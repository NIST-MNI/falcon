/* FILENAME:     off2pov.c
 * DESCRIPTION:  Conversion from off to POV RAY Mesh2 format
 * AUTHOR:       Vladimir Fonov
 *
 *
 * reference
 * http://www.povray.org/documentation/3.7.0/r3_4.html#r3_4_5_2_4
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include "falcon_cortex.h"

#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>

int niik_off_write_pov(const char *fname,kobj *obj,double diffuse,double phong,const char* obj_name,int metallic,int sphere,double sphere_radius);


static char *prog_version[] = {
  NIIK_VERSION "\n"
  "  program history\n"
  "  0.0.0   : August 29, 2016, Vladimir FONOV <vladimir.fonov@gmail.com>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  off2pov:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"off2pov\n");
  fprintf(stdout,"  usage: [options] <in.off> <out.pov>\n");
  fprintf(stdout,"  options:\n"
          "\t--object <name> - specify output object name, default brain\n"
          "\t--diffuse <f>   - object diffuse parameter, default 0.6\n"
          "\t--phong   <f>   - object phong   parameter, default 0.7\n"
          "\t--clobber       - clobber output file\n"
          "\t--metallic      - add metallic keyword\n"
          "\t--sphere        - use spherical coordinates instead of x,y,z\n"
          "\t--radius <n>    - sphere radius, default 50mm\n"
          "\t--colorize <map> - colorize surface using field\n"
          "\t--min <f>\n"
          "\t--max <f>\n"
          "\t--column <name> - column name to use\n"
          "\t--column_id <name> - column id to use\n"
           );
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int  nc,sc;
  const char *fcname="off2pov";
  double *meas=NULL;
  int meas_num;
  const char *in_off;
  const char *in_txt=NULL;

  const char * out_pov;
  int clobber=0;
  int sphere=0;
  double sphere_radius=50.0;

  const char * obj_name="brain";
  double diffuse=0.6;
  double phong=0.7;
  int metallic=0;

  int cmap=NIIK_COLORMAP_SPECTRAL;

  double omin=0.0;
  double omax=0.0;

  const char *column_name=NULL;
  int column_id=0;
  int print_range=0;

  niiktable *meas_in;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"metallic", no_argument, &metallic, 1},
    {"diffuse", required_argument, 0, 'd'},
    {"phong",  required_argument, 0,  'p'},
    {"object", required_argument, 0, 'o'},
    {"radius", required_argument, 0, 'r'},
    {"sphere", no_argument, &sphere, 1},

    {"min", required_argument, 0, 'i'},
    {"max", required_argument, 0, 'a'},
    {"colorize", required_argument, 0, 'c'},
    {"column", required_argument,0,'C'},
    {"column_id", required_argument,0,'N'},

    {"spectral", no_argument,&cmap,NIIK_COLORMAP_SPECTRAL},
    {"atrophy", no_argument,&cmap,NIIK_COLORMAP_ATROPHY},
    {"summer", no_argument,&cmap,NIIK_COLORMAP_SUMMER},
    {"jacobian", no_argument,&cmap,NIIK_COLORMAP_JACOBIAN},
    {"gray", no_argument,&cmap,NIIK_COLORMAP_GRAYSCALE},

    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c;

    c = getopt_long (argc, argv, "d:p:o:i:a:c:C:N:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'c':
      in_txt=optarg;
      break;
    case 'i':
      omin=atof(optarg);
      break;
    case 'a':
      omax=atof(optarg);
      break;
    case 'd':
      diffuse=atof(optarg);
      break;
    case 'p':
      phong=atof(optarg);
      break;
    case 'o':
      obj_name=optarg;
      break;
    case 'r':
      sphere_radius=atof(optarg);
      break;
    case 'C':
      column_name=optarg;
      column_id=-1;
      break;
    case 'N':
      column_id=atoi(optarg);
      break;
    case '?':
    default:
      usage();
      return 1;
    }
  }

  if((argc - optind)<2) {
    usage();
    return 1;
  }
  in_off  = argv[optind];
  out_pov = argv[optind+1];

  if (!clobber && !access (out_pov, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_pov);
    return 1;
  }

  niik_fc_display(fcname,1);
  NIIK_EXIT(((obj=off_kobj_read_offply(in_off))==NULL),fcname,"niik_kobj_read_off",1);

  if(in_txt)
  {
   
    NIIK_EXIT(((meas_in=niiktable_read(in_txt))==NULL),fcname,"niik_read_vector_ex",1);

    if(meas_in->col[0]->num!=obj->nvert) {
      fprintf(stderr,"[%s] Inconsistent number of measurments in txt, expected %i, got %i\n",fcname,meas_in->col[0]->num,obj->nvert);
      return 1;
    }
    /*find the right column*/
    if(column_id<0)
    {
      int i;
      for(i=0;i<meas_in->ncol;i++)
        if(!strcmp(column_name, meas_in->name[i])) break;
      if(i==meas_in->ncol)
      {
        fprintf(stderr,"[%s] Can't find column %s\n",fcname,column_name);
        return 1;
      }
      column_id=i;
    }
    if(omin==0.0 && omax==0.0) { /*TODO: use different way to set range*/
      omax = niik_get_max_from_double_vector(meas_in->col[column_id]->v, meas_in->col[column_id]->num);
      omin = niik_get_min_from_double_vector(meas_in->col[column_id]->v, meas_in->col[column_id]->num);
      fprintf(stdout,"[%s] Using range %8.4f %8.4f\n",fcname, omin, omax);
    }

    NIIK_EXIT( (niikcortex_add_color( obj, meas_in->col[column_id]->v, omin, omax, cmap ))==0,      fcname,"niikcortex_add_color",1 );
  }

  if(sphere && !obj->spherecoo)
  {
    fprintf(stderr,"%s missing spherical coordinates!\n", in_off);
    return 1;
  }

  NIIK_EXIT((!niik_off_write_pov(out_pov, obj, diffuse, phong, obj_name, metallic, sphere, sphere_radius)),fcname,"niik_off_write_pov",1);
  obj=off_kobj_free(obj);
  if(meas!=NULL)
    free(meas);

  niik_fc_display(fcname,0);
  exit(0);
} /* off2mniasc */

int niik_off_write_pov(const char *fname,kobj *obj, double diffuse, double phong, const char* obj_name, int metallic, int sphere,double sphere_radius) {
  const char *fcname="niik_off_write_pov";
  FILE *fp=NULL;
  int i;
  int verbose=niik_verbose();
  kvert *v;
  kface *f;
  kedge *e;
  int output_normal=1;

  niik_fc_display(fcname,1);

  if(output_normal) {
    /*update all normals?*/
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
  }


  fp=fopen(fname,"w");
  fprintf(fp,"#declare %s = \n",obj_name);
  fprintf(fp,"mesh2{\n");

  fprintf(fp,"vertex_vectors {\n");
  fprintf(fp,"%d,\n",obj->nvert);

  /*add coordinates*/
  for(v=obj->vert,i=0; v!=NULL; v=v->next, i++) {
    niikpt pt;
    if(sphere)
      pt=niikpt_kmul(niiksph_xyz(v->sph),sphere_radius); 
    else
      pt=v->v;
    fprintf(fp,"<%f,%f,%f>",pt.x,pt.y,pt.z);

    if(v->next)
      fprintf(fp,",");
    fprintf(fp,"\n");
  }
  fprintf(fp,"}\n");

  /*add normals*/
  if(output_normal) {
    fprintf(fp,"normal_vectors {\n");
    fprintf(fp,"%d,\n",obj->nvert);

    /*add coordinates*/
    for(v=obj->vert,i=0; v!=NULL; v=v->next, i++) {
      niikpt pt;
      if(sphere)
        pt=niiksph_xyz(v->sph);
      else
        pt=v->normal;
      
      fprintf(fp,"<%f,%f,%f>",pt.x,pt.y,pt.z);
      if(v->next)
        fprintf(fp,",");
      fprintf(fp,"\n");
    }
    fprintf(fp,"}\n");
  }

  /*write face colours*/
  if(obj->color) {
    fprintf(fp,"texture_list {\n");
    fprintf(fp,"%d,\n",obj->nface);

    /*add colours*/
    for(f=obj->face,i=0; f!=NULL; f=f->next, i++) {
      fprintf(fp,"texture{finish{ diffuse %f phong %f %s} pigment{rgb <%f,%f,%f>}}\n", 
              diffuse, phong, metallic?"metallic":"", f->color[0], f->color[1], f->color[2]);
    }
    fprintf(fp,"}\n");
  } else {
    /*uniform colour*/
    fprintf(fp,"texture_list {\n");
    fprintf(fp,"1,\n");
    fprintf(fp,"texture{finish{ diffuse %f phong %f %s } pigment{rgb <%f,%f,%f>}}\n",
              diffuse, phong, metallic?"metallic":"", 0.8, 0.8, 0.8 );
  }

  fprintf(fp,"face_indices {\n");
  fprintf(fp,"%d,\n", obj->nface);

  /*pass face indexes*/
  for(f=obj->face,i=0; f!=NULL; f=f->next, i++) {
    fprintf(fp,"<%d,%d,%d>", f->vert[0]->index-1, f->vert[1]->index-1, f->vert[2]->index-1);
    if(obj->color) {
      fprintf(fp,",%d", i);
    } else {
      fprintf(fp,",%d", 0);
    }
    if(f->next)
      fprintf(fp,",");
    fprintf(fp,"\n");
  }
  fprintf(fp,"}\n");
  fprintf(fp,"inside_vector <0,0,-1>\n");
  fprintf(fp,"}\n");

  fclose(fp);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_off_write_ply */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
