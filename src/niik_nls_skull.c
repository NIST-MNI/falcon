/* FILENAME:     niik_nls_skull.c
 * DESCRIPTION:  Kunio's nifti1 implementation of skull extraction
 *               using nonlocal segmentation (NLS) technique
 * AUTHOR:       Kunio Nakamura
 * DATE:         January 10, 2014
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

#define NIIK_NLS_VOX_DIFF    0
#define NIIK_NLS_ABS_DIFF    1
#define NIIK_NLS_CORRELATION 2

typedef struct {          /* class for nonlocal segmentation of skull */
  nifti_image *img;       /* input image */
  nifti_image *skull;     /* output skull */
  nifti_image **template; /* templates (# = ntemplate) */
  nifti_image **mask;     /* masks of skull matching to template */
  int ntemplate;          /* # of templates */
  int nlocal;
  double postprocess_median_radius;
} niik_nls_skull;         /* class for nonlocal segmentation of skull */

niik_nls_skull *niik_nls_skull_init(char *conf,int nlocal,double ppradius) {
  niik_nls_skull *nls=NULL;
  FILE *fp=NULL;
  int
  n,nimg=0,
    verbose=1;
  char
  iname[4096],mname[4096],
        c,fcname[32]="niik_nls_skull_init";
  NIIK_RET0(((nls=(niik_nls_skull *)calloc(1,sizeof(niik_nls_skull)))==NULL),fcname,"calloc failure");
  nls->ntemplate=0;
  nls->img=nls->skull=NULL;
  nls->template=nls->mask=NULL;
  nls->nlocal=nlocal;
  fprintf(stdout,"[%s] window size = %i   ->  %i\n",fcname,nls->nlocal,nlocal*2+1);
  nls->postprocess_median_radius=ppradius;
  fprintf(stdout,"[%s] median filter R = %f  (post-processing)\n",fcname,nls->postprocess_median_radius);
  if(verbose>0) fprintf(stdout,"[%s] opening %s\n",fcname,conf);
  NIIK_RET0(((fp=fopen(conf,"r"))==NULL),fcname,"fopen failure");
  while((c=getc(fp))!=EOF) if(c=='\n') nls->ntemplate++;
  fclose(fp);
  fprintf(stdout,"[%s]      # template files %-i\n",fcname,nls->ntemplate);
  NIIK_RET0(((fp=fopen(conf,"r"))==NULL),fcname,"fopen failure");
  NIIK_RET0(((nls->template=(nifti_image **)calloc(nls->ntemplate,sizeof(nifti_image *)))==NULL),fcname,"calloc for templates");
  NIIK_RET0(((nls->mask=(nifti_image **)calloc(nls->ntemplate,sizeof(nifti_image *)))==NULL),fcname,"calloc for template masks");
  for(n=0; n<nls->ntemplate; n++) {
    nls->template[n]=nls->mask[n]=NULL;
  }
  while(!feof(fp)) {
    if(fscanf(fp,"%s %s",iname,mname)!=2 || feof(fp)) break;

    fprintf(stdout,"[%s]  %3i %s %s\n",fcname,nimg,iname,mname);
    NIIK_RET0(((nls->template[nimg]=niik_image_read(iname))==NULL),fcname,"niik_image_read img");
    // NIIK_RET0((!niik_image_type_convert(nls->template[nimg],NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert");
    NIIK_RET0(((nls->mask[nimg]=niik_image_read(mname))==NULL),fcname,"niik_image_read mask");
    NIIK_RET0((!niik_image_type_convert(nls->mask[nimg],NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert");
    nimg++;
  } /* while */
  fclose(fp);
  return nls;
} /* niik_nls_skull_init  constructor */

niik_nls_skull *niik_nls_skull_free(niik_nls_skull *nls) {
  int n;
  if(nls==NULL) return NULL;
  // fprintf(stdout,"[niik_nls_skull_free] start\n");
  for(n=0; n<nls->ntemplate; n++) {
    nls->template[n]=niik_image_free(nls->template[n]);
    nls->mask[n]=niik_image_free(nls->mask[n]);
  }
  free(nls->template);
  free(nls->mask);
  nls->img=niik_image_free(nls->img);
  free(nls);
  return NULL;
} /* niik_nls_skull_free   destructor */


int niik_image_patch_list(nifti_image *img,nifti_image **lib,int nlib,int i,int j,int k,int len,double *out,int method)
/* patch and feature extraction
 * img is the reference image; for multi-contrast image, use u-dimension
 * lib is the list of libraries
 * i,j,k is the 3d coordinate
 * nlib is the # of images in the library (length of vector)
 * len is the local size (for example, 1 for 3x3x3 patch)
 * out is the output vector
 * method is the choice of metric
 *   0 = intensity absolute difference (no patch)
 *   1 = absolute difference
 *   2 = correlation coefficient
 */
{
  int
  verbose=0,
  area,size=0;
  int
  ilo,jlo,klo,ihi,jhi,khi,
      m,n,u,
      ii,jj,kk,nn;
  char fcname[32]="niik_image_patch_list";
  double x,y,sxy=0,sx=0,sy=0,sxx=0,syy=0;
  area=img->nx*img->ny;
  size=area*img->nz;
  for(n=0; n<nlib; n++) {
    switch(method) {
    default:
      fprintf(stderr,"[%s] ERROR: unknown method, %i\n",fcname,method);
      return 0;
    case NIIK_NLS_VOX_DIFF: /* intensity difference */
      out[n]+=fabs(niik_image_get_voxel(img,i+j*img->nx+k*area) -
                   niik_image_get_voxel(lib[n],i+j*img->nx+k*area) );
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]);
      break;
    case NIIK_NLS_ABS_DIFF: /* absolute difference */
      out[n]=0;
      for(kk=k-len; kk<=k+len; kk++) {
        if(kk<0) continue;
        if(kk>=img->nz) break;
        for(jj=j-len; jj<=j+len; jj++) {
          if(jj<0) continue;
          if(jj>=img->ny) break;
          for(ii=i-len,m=kk*area+jj*img->nx+ii; ii<=i+len; m++,ii++) {
            /* fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]); */
            if(ii<0) continue;
            if(ii>=img->nx) break;
            for(u=0; u<img->nu; u++) {
              out[n]+=fabs(niik_image_get_voxel(img,   m+u*size) -
                           niik_image_get_voxel(lib[n],m+u*size) );
            } /* u-dim */
          }
        }
      }
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]);
      break;
    case NIIK_NLS_CORRELATION: /* correlation */
      sxy=sx=sy=sxx=syy=0.0;
      nn=0;
      ilo=i-len;
      ilo=(ilo<0)?0:ilo;
      jlo=j-len;
      jlo=(jlo<0)?0:jlo;
      klo=k-len;
      klo=(klo<0)?0:klo;
      ihi=i+len;
      ihi=(ihi>=img->nx)?img->nx-1:ihi;
      jhi=j+len;
      jhi=(jhi>=img->ny)?img->ny-1:jhi;
      khi=k+len;
      khi=(khi>=img->nz)?img->nz-1:khi;
      out[n]=0;
      for(u=0; u<img->nu; u++) {
        for(kk=klo; kk<=khi; kk++) {
          if(kk<0) continue;
          if(kk>=img->nz) break;
          for(jj=jlo; jj<=jhi; jj++) {
            if(jj<0) continue;
            if(jj>=img->ny) break;
            for(ii=ilo,m=kk*area+jj*img->nx+ii; ii<=ihi; m++,ii++) {
              if(ii<0) continue;
              if(ii>=img->nx) break;
              x=niik_image_get_voxel(img,m);
              y=niik_image_get_voxel(lib[n],m);
              sx+=x;
              sy+=y;
              sxx+=x*x;
              syy+=y*y;
              sxy+=x*y;
              nn++;
            }
          }
        }
        out[n] += (nn * sxy - sx * sy) / sqrt(nn*sxx-sx*sx) / sqrt(nn*syy-sy*sy);
      } /* each contrast */
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f   | %6.0f %6.0f %6.0f %6.0f\n",i,j,k,out[n],sx,sy,sxy,(double)nn);
      break;
    } /* method */
  } /* lib img */
  return 1;
} /* niik_image_patch_list */


int niik_nls_skull_segment(niik_nls_skull *nls)
/* segment skull */
{
  char fcname[32]="niik_nls_skull_segment";
  int
  debug=0,
  nmax=5, /* number of top candidates used */
  m,n,nn,n3,i,j,k;
  double d;
  niikvec *v;
  int *vidx,checkidx[4];

  niik_fc_display(fcname,1);
  NIIK_RET0((nls==NULL),fcname,"missing nls");
  NIIK_RET0((nls->ntemplate<=0),fcname,"no templates");
  for(n=0; n<nls->ntemplate; n++) {
    for(j=1; j<=3; j++) {
      NIIK_RET0((nls->template[n]->dim[j]!=nls->img->dim[j]),fcname,"#voxel mismatch (templates)");
      NIIK_RET0((nls->mask[n]->    dim[j]!=nls->img->dim[j]),fcname,"#voxel mismatch (mask)");
    }
  }
  NIIK_RET0((nls->skull!=NULL),fcname,"skull already exists");

  // fprintf(stdout,"[%s]   checked\n",fcname);
  n3=nls->img->nx*nls->img->ny*nls->img->nz;
  NIIK_RET0(((nls->skull=niik_image_copy(nls->mask[0]))==NULL),fcname,"niik_image_copy");
  NIIK_RET0((!niik_image_clear(nls->skull)),fcname,"niik_image_clear");
  NIIK_RET0(((v=niikvec_init(nls->ntemplate))==NULL),fcname,"niikvec_init");

  vidx=(int *)calloc(nls->ntemplate,sizeof(int));

  checkidx[0]=1;
  checkidx[1]=93;
  checkidx[2]=206;
  checkidx[3]=130;

  for(k=m=0; k<nls->skull->nz; k++) {
    for(j=0; j<nls->skull->ny; j++) {
      for(i=0; i<nls->skull->nx; m++,i++) {

        if(checkidx[0]>0) {
          debug=0;
          /*if(i!=checkidx[1]) continue;
            if(j!=checkidx[2]) continue;
            if(k!=checkidx[3]) continue;
            debug=1; */
          if(i==checkidx[1])
            if(j==checkidx[2])
              if(k==checkidx[3])
                debug=1;
        }

        // if(niik_image_get_voxel(nls->img,m)>0.5) continue;
        for(n=0,d=0; n<nls->ntemplate; n++) {
          d+=niik_image_get_voxel(nls->mask[n],m);
        }
        if(d<=0.00001) continue;
        if(fabs(d-nls->ntemplate)<=0.00001)  {
          niik_image_set_voxel(nls->skull,m,100.0);
          continue;
        }

        /* patch */
        niik_image_patch_list(nls->img,nls->template,nls->ntemplate,i,j,k,nls->nlocal,v->v,NIIK_NLS_ABS_DIFF);
        if(debug>0) {
          fprintf(stdout,"[%3i %3i %3i] \n",i,j,k);
          fprintf(stdout,"  im          %3i",
                  (int)niik_image_get_voxel(nls->img,m));
          for(n=1; n<nls->img->nu; n++) {
            fprintf(stdout," %3i",(int)niik_image_get_voxel(nls->img,m+n3*n));
          }
          fprintf(stdout,"\n");
          for(n=0; n<nls->ntemplate; n++) {
            fprintf(stdout,"  %2i %4i %i   %3i",n,(int)v->v[n],(int)niik_image_get_voxel(nls->mask[n],m),(int)niik_image_get_voxel(nls->template[n],m));
            for(nn=1; nn<nls->img->nu; nn++) {
              fprintf(stdout," %3i",(int)niik_image_get_voxel(nls->template[n],m+n3*nn));
            }
            fprintf(stdout,"\n");
          }
          niik_sort_double_vector_index(v->v,vidx,v->num);
          for(n=0; n<nls->ntemplate; n++) {
            fprintf(stdout,"\t%2i %3i %6i | %i\n",n,(int)vidx[n],(int)v->v[vidx[n]],(int)niik_image_get_voxel(nls->mask[vidx[n]],m));
          }
        }

        niik_sort_double_vector_index(v->v,vidx,v->num);
        for(n=0,d=0; n<nmax; n++) d += v->v[vidx[n]];
        for(n=0; n<nmax; n++) v->v[vidx[n]] = d / (v->v[vidx[n]]+1e-6); /* inv normalize */
        for(n=0,d=0; n<nmax; n++) d += v->v[vidx[n]];
        for(n=0; n<nmax; n++) v->v[vidx[n]] /= d; /* calculate factors */
        for(n=0,d=0; n<nmax; n++) /* weighted sum with factors */
          d += niik_image_get_voxel(nls->mask[vidx[n]],m) * v->v[vidx[n]];

        if(debug>0) {
          for(n=0; n<nmax; n++) {
            fprintf(stdout,"    %2i [%2i]   %3i | %i\n",vidx[n],n,(int)(100.0*v->v[vidx[n]]),(int)niik_image_get_voxel(nls->mask[vidx[n]],m));
          }
          fprintf(stdout,"  output = %-3i\n",(int)floor(d*100.0+0.5));
        }

        niik_image_set_voxel(nls->skull,m,d*100.0); /* set voxel (probability) */

      }
    }
  }

  if(nls->postprocess_median_radius>0) {
    fprintf(stdout,"[%s] post-process median filter %.3f\n",fcname,nls->postprocess_median_radius);
    NIIK_RET0((!niik_image_filter_median_radius(nls->skull,NULL,nls->postprocess_median_radius)),fcname,"niik_image_filter_median_radius");
  }


  free(vidx);
  v=niikvec_free(v);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_nls_skull_segment   segment1 */



static char *prog_skull_segment_version[] = {
  "  niik_skull_segment history\n"
  "  0.0.0  January 10, 2014, Kunio Nakamura, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_skull_segment_help[] = {
  "  niik_skull_registration:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niik_skull_segment_nls\n");
  fprintf(stdout,"  usage: [options] <img.nii> <conf> <out.nii>\n\n");
  fprintf(stdout,"\n");
}


int main(int argc,char *argv[],char *envp[]) {
  niik_nls_skull *nls=NULL;
  nifti_image *tmpimg,**imglist;
  char
  fcname[20]="niik_nls_skull";
  int
  nlocal=2,
  verbose=2,
  nc=1,sc=1;
  double ppradius=1;

  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  if(argc==1) {
    usage();
    exit(0);
  }

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);
  fprintf(stdout,"  niik_nls_skull: version %i.%i.%i\n",MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  fprintf(stdout,"  niik_nls_skull: executed at %s\n",tmstr);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*prog_skull_segment_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*prog_skull_segment_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*prog_skull_segment_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"--window",8)) {
        nlocal=atoi(argv[++nc]);
      }

      else if(!strncmp(argv[nc],"--median",8)) {
        ppradius=atof(argv[++nc]);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
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

  NIIK_EXIT((argc<4),fcname,"too few argments",9);
  NIIK_EXIT((argc>4),fcname,"too many argments",9);

  fprintf(stdout,"[%s] reading conf file  %s\n",fcname,argv[2]);
  NIIK_EXIT(((nls=niik_nls_skull_init(argv[2],nlocal,ppradius))==NULL),fcname,"niik_nls_skull_init",9);
  fprintf(stdout,"[%s] reading image      %s\n",fcname,argv[1]);
  NIIK_EXIT(((nls->img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);

  /* run segmentation */
  NIIK_EXIT((!niik_nls_skull_segment(nls)),fcname,"niik_nls_skull_segment",9);

  /* write output */
  fprintf(stdout,"[%s] writing output     %s\n",fcname,argv[3]);
  NIIK_EXIT((!niik_image_write(argv[3],nls->skull)),fcname,"niik_image_write",9);

  /*  imglist=niik_image_unmerge3d(nls->img);
  NIIK_EXIT(((tmpimg=niik_image_mask_add_red_color_uint8(imglist[2],nls->skull))==NULL),fcname,"niik_image_mask_add_red_color_uint8",9);
  fprintf(stdout,"[%s] writing output     check.nii.gz\n",fcname);
  NIIK_EXIT((!niik_image_write("check.nii.gz",tmpimg)),fcname,"niik_image_write",9); */

  nls=niik_nls_skull_free(nls);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_nls_skull */
