/* Filename:     niik_bpf.c
 * Description:  Kunio's BPF calculation program
 * Author:       Kunio Nakamura
 * Date:         April 10, 2012
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (3)
#define MICRO_VERSION (0)

#define THRESHOLD_METHOD_RIDLER 1
#define THRESHOLD_METHOD_OTSU   2

static char *niik_bpf_version[] = {
  "  niik_bpf version history\n"
  "  0.0.0  April 10, 2013, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
  "\n  0.1.0  April 24, 2013, knakamura@mrs.bic.mcgill.ca\n"
  "  -force zero partial volume outside brainmask\n"
  "\n  0.2.0  April 26, 2013, knakamura@mrs.bic.mcgill.ca\n"
  "  -non-brain is mean, not mode\n"
  "\n  0.3.0  November 14, 2013, knakamura@mrs.bic.mcgill.ca\n"
  "  -added option to read/handle lesion mask\n"
};

static char *niik_bpf_help[] = {
  "  simple classification program\n"
  "  usage: <in_T1W.nii> <brainmask.nii> <out_parenchyma.nii> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help                   : show this usage\n"
  "  --version                  : show version info\n"
  "\n"
  "  -otsu                      : use Otsu threshold [default -ridler]\n"
  "  -ridler                    : use Ridler threshold [default -ridler]\n"
  "\n"
  "  -lm <lesionmask.nii>       : use lesion mask [default none]\n"
  "\n"
};

double niik_bpf_image_get_non_brain_mode(nifti_image *img,nifti_image *maskimg,double dmin,double dmax,int num,int avgnum);

void usage() {
  fprintf(stdout,"%s",*niik_bpf_help);
  return;
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *img=NULL,
   *tmpimg=NULL,
    *outimg=NULL,
     *nonbrain=NULL,
      *lesionmask=NULL,
       *brainmask=NULL;
  char
  fcname[24]="niik_bpf";
  int
  method_thresh=THRESHOLD_METHOD_RIDLER,
  verbose=0,
  nx,xy,
  n,i,j,k,nn,
  nc=1,sc=1;
  double
  dmin,dmax,dsum,
       bmean,
       nmean,
       thresh=2;
  unsigned char *bimg=NULL,*bmask=NULL;
  double *dimg=NULL;

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
  fprintf(stdout,"[%s]: version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  fprintf(stdout,"[%s]: executed at %s\n",fcname,tmstr);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niik_bpf_version);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niik_bpf_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-verbose",8)) {
        verbose=2;
      } else if(!strncmp(argv[nc],"-ridler",7)) {
        method_thresh=THRESHOLD_METHOD_RIDLER;
      } else if(!strncmp(argv[nc],"-otsu",5)) {
        method_thresh=THRESHOLD_METHOD_OTSU;
      }


      else if(!strncmp(argv[nc],"-lm",3)) {
        NIIK_EXIT(((lesionmask=niik_image_read(argv[++nc]))==NULL),fcname,"niik_image_read lesion mask",9);
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


  /*********************************************
   * READING IMAGES AND BRAIN MASK
   *********************************************/
  fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read t1w",9);

  fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[2]);
  NIIK_EXIT(((brainmask=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read mask",9);
  NIIK_EXIT((!niik_image_threshold(brainmask,0.5)),fcname,"niik_image_threshold",9);
  bmask=brainmask->data;
  fprintf(stdout,"[%s] count = %i\n",fcname,niik_image_count_mask(brainmask));

  if(lesionmask!=NULL) {
    NIIK_EXIT((lesionmask->nvox!=brainmask->nvox),fcname,"mismatch of mask #vox",9);
  }


  /* SEGMENTATION */

  /* NIIK_RET0((!niik_image_display_stats_mask(img,brainmask)),fcname,"niik_image_display_stats_mask"); */
  switch(method_thresh) {
  case THRESHOLD_METHOD_RIDLER:
    NIIK_EXIT((!niik_image_thresh_ridler(img,brainmask,&thresh)),fcname,"niik_image_thresh_ridler",9);
    break;
  case THRESHOLD_METHOD_OTSU:
    NIIK_EXIT((!niik_image_thresh_otsu  (img,brainmask,&thresh)),fcname,"niik_image_thresh_otsu",9);
    break;
  default:
    fprintf(stderr,"[%s] ERROR: unkown method for threshold, %i\n",fcname,method_thresh);
    exit(0);
  }

  fprintf(stdout,"[%s] threshold = %-.4f\n",fcname,thresh);
  NIIK_RET0(((outimg=niik_image_threshold_new(img,thresh))==NULL),fcname,"niik_image_threshold_new");
  NIIK_RET0((!niik_image_mask(outimg,brainmask)),fcname,"niik_image_mask");
  if(lesionmask!=NULL) {
    fprintf(stdout,"[%s] adding lesions\n",fcname);
    NIIK_RET0(((tmpimg=niik_image_copy(outimg))==NULL),fcname,"niik_image_copy");
    for(i=0; i<img->nvox; i++) {
      if(niik_image_get_voxel(lesionmask,i)>0)
        niik_image_set_voxel(tmpimg,i,1);
    }
    fprintf(stdout,"[%s] non-pve volume                :  %-12.6f\n",fcname,niik_image_get_mask_vol(tmpimg));
    tmpimg=niik_image_free(tmpimg);
  } else {
    fprintf(stdout,"[%s] non-pve volume                :  %-12.6f\n",fcname,niik_image_get_mask_vol(outimg));
  }

  // brain intensity
  bmean=niik_image_get_mean(img,outimg);
  fprintf(stdout,"[%s] brain mean        %12.6f\n",fcname,bmean);

  // non-brain intensity
  NIIK_RET0(((nonbrain=niik_image_copy(brainmask))==NULL),fcname,"niik_image_copy");
  NIIK_RET0((!niik_image_maskout(nonbrain,outimg)),fcname,"niik_image_maskout");
  dmin=niik_image_get_min(img,brainmask);
  dmax=niik_image_get_max(img,brainmask);
  if(verbose>1) fprintf(stdout,"[%s] min / max = %9.4f %9.4f\n",fcname,dmin,dmax);
  /*nmean=niik_bpf_image_get_non_brain_mode(img,nonbrain,
                            dmin,dmax,
                            250,5);*/
  nmean=niik_image_get_mean(img,nonbrain);
  fprintf(stdout,"[%s] non-brain mean    %12.6f\n",fcname,nmean);
  nonbrain=niik_image_free(nonbrain);


  /* PARTIAL VOLUME CORRECTION */
  bimg = niik_image_get_voxels_as_uint8_vector(outimg);
  dimg = niik_image_get_voxels_as_double_vector(outimg);
  nx=xy=img->nx;
  xy  *=img->ny;

  for(k=n=0,dsum=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; n++,i++) {
        /* if(bmask[n]==0) continue; /* version 0.1.0 */
        if(bimg[n]>0) {
          dimg[n]=1;
          dsum+=1;
          continue;
        }
        if(lesionmask!=NULL) if(niik_image_get_voxel(lesionmask,i)>0.5) {
            dimg[n]=1;
            dsum+=1;
            continue;
          }
        nn=0;
        if(i>0) if(bimg[n- 1]!=bimg[n]) nn++;
        if(j>0) if(bimg[n-nx]!=bimg[n]) nn++;
        if(k>0) if(bimg[n-xy]!=bimg[n]) nn++;
        if(i<img->nx-1) if(bimg[n+ 1]!=bimg[n]) nn++;
        if(j<img->ny-1) if(bimg[n+nx]!=bimg[n]) nn++;
        if(k<img->nz-1) if(bimg[n+xy]!=bimg[n]) nn++;
        if(!nn) dimg[n]=0;
        else {
          dimg[n]=1.0-niik_pvc(niik_image_get_voxel(img,n),nmean,bmean);
        }
        dsum+=dimg[n];
      }
    }
  }
  free(bimg);
  niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32);
  niik_image_set_voxels_from_double_vector(outimg,dimg);
  free(dimg);

  dsum*=niik_image_get_voxel_size(outimg)*0.001;
  fprintf(stdout,"[%s] brain parenchymal volume      :  %-12.6f\n",fcname,dsum);
  fprintf(stdout,"[%s] beast volume                  :  %-12.6f\n",fcname,niik_image_get_mask_vol(brainmask));
  fprintf(stdout,"[%s] beast parenchymal fraction    :  %-12.8f    %s %s\n",fcname,
          dsum/niik_image_get_mask_vol(brainmask),
          argv[1],argv[2]);

  /* OUTPUT VOLUME */
  fprintf(stdout,"[%s] output image  %s\n",fcname,argv[3]);
  NIIK_RET0((!niik_image_write(argv[3],outimg)),fcname,"niik_image_write");

  /*********************************************
   * clean up
   *********************************************/
  img=niik_image_free(img);
  outimg=niik_image_free(outimg);
  brainmask=niik_image_free(brainmask);

  niik_fc_display(fcname,0);
  return 0;
} /* niik_bpf */


double niik_bpf_image_get_non_brain_mode(nifti_image *img,nifti_image *maskimg,double dmin,double dmax,int num,int avgnum) {
  double
  dx,d,
  *hx,*hy;
  int
  n,
  verbose=2;
  const char *fcname="niik_bpf_image_get_non_brain_mode";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(num<=0) {
    fprintf(stderr,"[%s] ERROR: num is invalid %i\n",fcname,num);
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (%s) memory alloc + init %i\n",fcname,num);
  if((hx=niik_calloc_double_vector(num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_calloc_double_vector\n",fcname);
    return NIIKMAX;
  }
  if((hy=niik_calloc_double_vector(num))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_calloc_double_vector\n",fcname);
    return NIIKMAX;
  }
  dx=(dmax-dmin)/(num-1.0);
  for(n=0; n<num; n++) {
    hx[n]=n*dx+dmin;
  }
  if(verbose) fprintf(stdout,"-d (%s) histogram %8.3f %8.3f %8.3f %i\n",fcname,dmin,dx,dmax,num);
  if(!niik_image_histogram_limits(img,maskimg,hx,hy,num)) {
    fprintf(stderr,"[%s] ERROR: niik_image_histogram\n",fcname);
    return NIIKMAX;
  }
  /* 2012-03-14 Kunio
   * -check if histogram exists */
  for(n=0,d=0; n<num; n++) {
    d+=hy[n];
  }
  if(d<1e-3) {
    fprintf(stderr,"[%s] ERROR: no histogram\n",fcname);
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"-d (%s) writing tmp_x.txt\n",fcname);
    niik_write_double_vector("tmp_x.txt",hx,num);
    niik_write_double_vector("tmp_y0.txt",hy,num);
  }
  if(verbose) fprintf(stdout,"-d (%s) average: vec %i   avg %i\n",fcname,num,avgnum);
  if(!niik_runavg_double_vector(hy,num,avgnum)) {
    fprintf(stderr,"[%s] ERROR: niik_runavg_double_vector\n",fcname);
    return NIIKMAX;
  }
  if(verbose>1) {
    fprintf(stdout,"-d (%s) writing tmp_y1.txt\n",fcname);
    niik_write_double_vector("tmp_y1.txt",hy,num);
  }
  for(n=1; n<num; n++) {
    if(hy[n-1]<=hy[n]) continue;
    break;
  }
  n--;
  /*if((n=niik_get_max_index_double_vector(hy,num*0.7))<0) {
    fprintf(stderr,"[%s] ERROR: niik_get_min_index_double_vector\n",fcname);
    return NIIKMAX; }*/
  d=hx[n];
  free(hx);
  free(hy);
  return d;
} /* niik_bpf_image_get_mode */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/