/* FILENAME: niikcortex_calc_thickness.c
 * DESCRIPTION: calculate cortical thickness using simple algorithm
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_cortex.h"

#include <unistd.h>
#include <getopt.h>

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"Simple cortical thickness calculation\n");
  fprintf(stdout,"Will calculate thickness per vortex assuming node-to-node correspondense\n");
  fprintf(stdout,"  usage: [options] <white.ply> <pial.ply> <out.ply/txt> \n\n");
  fprintf(stdout,"Options:\n");

  fprintf(stdout,"\t--sph - will include spherical coordinates in the output\n");
  fprintf(stdout,"\t--ply - will output .ply file format with variable thickness\n");
  fprintf(stdout,"\t--white - use white surface for ply , otherwise - peal\n");
  fprintf(stdout,"\t--pial  - use pieal surface for ply , (default)\n");
  fprintf(stdout,"\t--smooth <fwhm> - apply smoothing kernel (default: disable) - will try to apply it iteratively\n");
  fprintf(stdout,"\t--mask <mask.mnc> - Use mask to mark vertices with unreliable thickness (non-cortex)\n");
  fprintf(stdout,"\t--maskval <mask value> - Use this value inside mask, default 0.0\n");
}


int main(int argc,char *argv[],char *envp[]) {
  kobj *ctx[2];
  double *thk=NULL,*psi=NULL,*the=NULL;
  unsigned char *mask_flag=NULL;
  double *vectors[3];
  double smooth=0.0;
  double mask_val=0.0;
  const char *maskfn=NULL;
  const char *fcname="niikcortex_calc_thickness.c";
  int clobber=0;
  int output_sph=0;
  int output_ply=0;
  int output_white=0;
  nifti_image *maskimg=NULL;

  const char *in_ics,*in_ocs,*out_thickness;
  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"sph",              no_argument, &output_sph, 1},
    {"ply",              no_argument, &output_ply, 1},
    {"white",            no_argument, &output_white, 1},
    {"pial",             no_argument, &output_white, 0},
    {"help",             no_argument, 0, 'h'},
    {"version",          no_argument, 0, 'v'},
    {"smooth",     required_argument, 0, 's'},
    {"mask",       required_argument, 0, 'm'},
    {"maskval",    required_argument, 0, 'M'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "uhvs:m:M:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
      return 1;
    case 's':
      smooth=atof(optarg);
      break;
    case 'm':
      maskfn=optarg;
      break;
    case 'M':
      mask_val=atof(optarg);
      break;
    case '?':
    case 'u':
    case 'h':
    default:
      usage();
      return 1;
    }
  }

  if((argc - optind)<3) {
    usage();
    return 1;
  }
  in_ics=argv[optind];
  in_ocs=argv[optind+1];
  out_thickness=argv[optind+2];

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading white surface       %s\n",fcname,in_ics);
  NIIK_RETURN(((ctx[0]=off_kobj_read_offply(in_ics))==NULL),"reading white surface",1);

  fprintf(stdout,"[%s] reading pial surface        %s\n",fcname,in_ocs);
  NIIK_RETURN(((ctx[1]=off_kobj_read_offply(in_ocs))==NULL),"reading pial surface",1);

  NIIK_RETURN((ctx[0]->nvert!=ctx[1]->nvert),"wrong number of vertices",1);
  NIIK_RETURN(((thk=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for thickness array",1);

  if(maskfn)
  {
     if( (maskimg=niik_image_read(maskfn))==NULL) {
        fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,maskfn);
        exit(1);
     }
     // binarize ?
     //niik_image_type_convert(maskfn,NIFTI_TYPE_UINT8);
 
  }


  if(output_sph && !output_ply) {
    NIIK_RETURN(((psi=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for psi array",1);
    NIIK_RETURN(((the=(double *)calloc(ctx[0]->nvert,sizeof(double)))==NULL),"calloc for the array",1);
  }

  NIIK_RETURN((!niikcortex_calc_thickness(ctx[0], ctx[1], thk,psi,the,1,0)),"niikcortex_calc_thickness",1);

  if(maskimg)
  {
    int i;
    kvert *vi,*vo;
    NIIK_RETURN(((mask_flag=(unsigned char *)calloc(ctx[0]->nvert,sizeof(unsigned char)))==NULL),"calloc for mask array",1);

    for(vi=ctx[0]->vert, vo=ctx[1]->vert,i=0; vi!=NULL && vo!=NULL; vi=vi->next,vo=vo->next,++i)
    { 
      int mask_i = (niik_image_interpolate_3d_linear_xyz(maskimg, vi->v) >= 0.5);
      int mask_o = (niik_image_interpolate_3d_linear_xyz(maskimg, vo->v) >= 0.5);

      if(mask_i||mask_o)
      {
        thk[i]=mask_val;
        mask_flag[i]=0;
      } else {
        mask_flag[i]=1;
      }
    }
  }

  if(smooth>0.0)
  {
    double it_smooth,it_sigma;
    double mean_elen=off_get_kobj_mean_edge_length(ctx[1]);

    /*apply small kernel several times, since each iteration only works with the local neighbourhood*/
    /*HACK: mean_elen*2  - is an arbitrary factor */
    int iterations=ceil((smooth*smooth)/(mean_elen*mean_elen*4));
    int i;
    if(iterations<1) iterations=1;
    it_smooth=sqrt(smooth*smooth/iterations);
    it_sigma=it_smooth/2.355; 

    fprintf(stdout,"[%s] Smoothing: %f x %d\n",fcname, it_smooth, iterations);

    if(mask_flag) 
    {
      if(!off_surface_gauss_smooth_using_vert_with_mask(ctx[1], thk, it_sigma, iterations, mask_flag)) {
        fprintf(stderr,"[%s] ERROR: off_surface_gauss_smooth_using_vert\n",fcname);
        return 1;
      }
    } else {
      if(!off_surface_gauss_smooth_using_vert(ctx[1], thk, it_sigma, iterations)) {
        fprintf(stderr,"[%s] ERROR: off_surface_gauss_smooth_using_vert\n",fcname);
        return 1;
      }
    }
  }

  fprintf(stdout,"[%s] writing output file        %s\n",fcname,out_thickness);

  if(output_ply)
  {
    const char * meas[]={"thickness"};
    double *val[]={thk};
    NIIK_RETURN((!off_kobj_write_ply_ex(out_thickness, (output_white?ctx[1]:ctx[0]) ,0, output_sph,1,0,1, meas, val)),"off_kobj_write_ply_ex",1);
  } else {
    if(output_sph) {
      const char * headers[]={"thickness","psi","the"};
      vectors[0]=thk;
      vectors[1]=psi;
      vectors[2]=the;
      NIIK_RETURN((!niik_write_double_vectors_ex(out_thickness,vectors,ctx[0]->nvert,3, headers)),"niik_write_double_vectors",1);
    } else {
      NIIK_RETURN((!niik_write_double_vector_ex(out_thickness,thk,ctx[0]->nvert,"thickness")),"niik_write_double_vector",1);
    }
    if(output_sph) {
      free(the);
      free(psi);
    }
  }

  free(thk);
  ctx[0]=off_kobj_free(ctx[0]);
  ctx[1]=off_kobj_free(ctx[1]);

  niik_fc_display(fcname,0);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/