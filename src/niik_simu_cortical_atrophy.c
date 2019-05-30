/* FILENAME:     niik_simu_cortical_atrophy.c
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 * DATE:         February 3, 2014
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include "falcon_morph.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  niik_simu_cortical_atrophy history\n"
  "  0.0.0 February 3, 2014; knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_describe[] = {
  "  [niik_simu_cortical_atrophy] description\n"
};

static char *prog_help[] = {
  "  niik_simu_cortical_atrophy:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niik_simu_cortical_atrophy\n");
  fprintf(stdout,"  usage: [options] <GM.nii> <white.off> <pial.off> <out.nii>\n\n");
}

int niik_simu_cortical_atrophy_func(kobj **obj,nifti_image *refimg,double radius);

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *img=NULL;
  kobj *obj[2];
  int
  verbose=1,
  nc,sc;
  char
  fcname[32]="niik_simu_cortical_atrophy";
  double radius=3;

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(1);
      } else if(!strncmp(argv[nc],"--describe",10)) {
        fprintf(stdout,"%s",*prog_describe);
        exit(1);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(1);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(1);
      }

      else if(!strncmp(argv[nc],"-verbose",8)) {
        verbose=1;
      } else if(!strncmp(argv[nc],"-quiet",6)) {
        verbose=0;
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

  NIIK_EXIT((argc>5),fcname,"too many arguments",1);
  NIIK_EXIT((argc<5),fcname,"too few arguments",1);

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading reference image:  %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"reading niik_image_read",1);
  fprintf(stdout,"[%s] reading white surf:       %s\n",fcname,argv[2]);
  NIIK_EXIT(((obj[0]=off_kobj_read_offply(argv[2]))==NULL),fcname,"off_kobj_read_off",1);
  fprintf(stdout,"[%s] reading pial surf:        %s\n",fcname,argv[3]);
  NIIK_EXIT(((obj[1]=off_kobj_read_offply(argv[3]))==NULL),fcname,"off_kobj_read_off",1);


  NIIK_EXIT((niik_simu_cortical_atrophy_func((kobj **)obj,img,radius)),fcname,"niik_simu_cortical_atrophy_func",1);

  fprintf(stdout,"[%s] writing output image: %s\n",fcname,argv[4]);
  niik_image_write(argv[4],img);


  img=niik_image_free(img);
  niik_fc_display(fcname,0);

  exit(0);
} /* niik_simu_cortical_atrophy */

int niik_simu_cortical_atrophy_func(kobj **obj,nifti_image *refimg,double radius)
/* return zero on success
 *
 */
{
  nifti_image
  *tmpimg=NULL,
   *maskimg=NULL;
  char fcname[32]="niik_simu_cortical_atrophy_func";
  kvert *v=NULL,*vi=NULL,*vo=NULL;
  kface *f=NULL;
  int
  iter,maxiter=5,
       n,n3;
  unsigned char *bimg=NULL;
  double
  pthresh=0.2,
  R;

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] define ROI\n",fcname);
  NIIK_RET1(((maskimg=niik_image_copy_as_type(refimg,NIFTI_TYPE_UINT8))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET1((!niik_image_clear(maskimg)),fcname,"niik_image_clear");
  NIIK_RET1(((tmpimg=niik_image_copy_as_type(refimg,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET1((!off_obj2img(maskimg,obj[1],255)),fcname,"off_obj2img");
  fprintf(stdout,"[%s]   dilation %.3f\n",fcname,radius);
  NIIK_RET1((!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_DILATE,radius)),fcname,"niik_image_morph_3d_radius");

  fprintf(stdout,"[%s] prep deformation image\n",fcname);
  refimg->ndim=6;
  refimg->nt=refimg->dim[4]=1;
  refimg->nu=refimg->dim[5]=1;
  refimg->nv=refimg->dim[6]=3;
  refimg->nvox=refimg->nx*refimg->ny*refimg->nz*refimg->nv;
  refimg->datatype=NIFTI_TYPE_FLOAT32;
  refimg->nbyper=sizeof(float);
  free(refimg->data);
  NIIK_RET1(((refimg->data=(void *)calloc(refimg->nvox,sizeof(float)))==NULL),fcname,"calloc for refimg->data");
  n3=maskimg->nx*maskimg->ny*maskimg->nz;
  fprintf(stdout,"[%s]     %i %i %i\n",fcname,n3,refimg->nvox,maskimg->nvox);

  fprintf(stdout,"[%s]   normal calculation\n",fcname);
  for(vo=obj[1]->vert,vi=obj[0]->vert; vo!=NULL; vo=vo->next,vi=vi->next) {
    vi->normal=vo->normal=niikpt_unit(niikpt_sub(vo->v,vi->v));
  }
  off_kobj_add_color(obj[1]);
  off_update_kobj_face_normal(obj[1]);
  for(f=obj[1]->face; f!=NULL; f=f->next) {
    f->color[0]=f->normal.x;
    f->color[1]=f->normal.y;
    f->color[2]=f->normal.z;
    f->color[3]=1;
  }

  fprintf(stdout,"[%s]   maskout non-<gray matter>\n",fcname);
  bimg=(unsigned char *)maskimg->data;
  for(n=0; n<maskimg->nvox; n++) {
    if(niik_image_get_voxel(tmpimg,n)<pthresh)  {
      bimg[n]=0;
    }
  }
  if(1) {
    fprintf(stdout,"[%s] writing test_mask.nii.gz\n",fcname);
    niik_image_write("test_mask.nii.gz",maskimg);
  }


  fprintf(stdout,"[%s] deformation image\n",fcname);
  for(iter=0; iter<maxiter; iter++) {
    R = (2.0+radius)-radius/NIIK_IMAX(1,(maxiter-1))*iter;
    fprintf(stdout,"[%s]   iteration %3i  %9.4f\n",fcname,iter+1,R);

    refimg->ndim=6;
    refimg->nt=refimg->dim[4]=1;
    refimg->nu=refimg->dim[5]=1;
    refimg->nv=refimg->dim[6]=3;
    NIIK_RET1((!off_obj2img_color(refimg,obj[1],1.0)),fcname,"off_obj2img");

    /* for(v=obj[1]->vert;v!=NULL;v=v->next){
      n=niik_image_get_index_niikpt(maskimg,v->v);
      niik_image_set_voxel(refimg,n,v->normal.x); n+=n3;
      niik_image_set_voxel(refimg,n,v->normal.y); n+=n3;
      niik_image_set_voxel(refimg,n,v->normal.z);
      }*/

    fprintf(stdout,"[%s]     mask\n",fcname);
    for(n=0; n<maskimg->nvox; n++) {
      if(bimg[n]==0) {
        continue;
      } else {
        niik_image_set_voxel(refimg,n,0);
        n+=n3;
        niik_image_set_voxel(refimg,n,0);
        n+=n3;
        niik_image_set_voxel(refimg,n,0);
      }
    }

    fprintf(stdout,"[%s]     blur %9.3f\n",fcname,R);
    refimg->ndim=5;
    refimg->nt=refimg->dim[4]=1;
    refimg->nu=refimg->dim[5]=3;
    refimg->nv=refimg->dim[6]=1;
    NIIK_RET1((!niik_image_filter_gaussian_update(refimg,5,R)),fcname,"niik_image_filter_gaussian_update");
  }

  maskimg=niik_image_free(maskimg);
  niik_fc_display(fcname,0);
  return 0;
} /* niik_simu_cortical_atrophy_func */


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/