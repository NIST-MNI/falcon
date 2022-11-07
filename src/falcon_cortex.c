/* Filename:     nifti1_kunio_cortex.c
 * Description:  functions for cortical surface model
 * Author:       Kunio Nakamura
 * Date:         March 26, 2012
 */

#include "falcon_cortex.h"
#include "falcon_morph.h"


/* smallest deformation possibly applied */

int g_niik_cortex_debug = 0;
int  niikcortex_get_debug() {
  return g_niik_cortex_debug;
}
void niikcortex_set_debug(int debug) {
  g_niik_cortex_debug = debug;
}



char *niik_numerical_method_string( int code ) {
  switch( code ) {
  case NIIK_NUM_METHOD_FORWARD_EULER:
    return "Forward Euler";
  case NIIK_NUM_METHOD_BACKWARD_EULER:
    return "Backward Euler";
  case NIIK_NUM_METHOD_MIDPOINT:
    return "Mid-point";
  case NIIK_NUM_METHOD_RUNGE_KUTTA:
    return "Runge-Kutta";
  case NIIK_NUM_METHOD_UNKNOWN:
    return "unknown";
  }
  return "unknown";
}


/**********************************************************************
 *
 * niikcortex_estimate_tissue_values
 *
 * -estimates tissue and surface intensiites
 * -CSF from ventricle mask (mode)
 * -WM from GWI_mask and brain_mask
 *   -mode within GWI_mask and above threshold (mean brain_mask)
 * -GM also from GWI_mask and brain_mask
 *   -mode within GWI_mask and below threshold (mean brain_mask)
 * -ICS above 50% on sobel-gradient image above mean GM
 *
 **********************************************************************/
int niikcortex_estimate_tissue_values(nifti_image *img,nifti_image *brain_mask,
                                      nifti_image *csf_mask,nifti_image *GWI_mask,
                                      double *iCSF, double *iWM, double *iGM,double *iBrain,
                                      double *intICS,double *intOCS,
                                      double *ranICS,double *ranOCS,
                                      double mixICS,double mixOCS)

{
  nifti_image
      *gwi_dil,
      *gradimg,
      *brain_dil,
      *tmpimg;

  double
      gthresh,
      thresh,
      dmax,
      dmin;

  unsigned char
      *bseg,
      *bimg;
  
  int
      i,nvox,
      verbose=niik_verbose();
  
  const char *fcname=__func__;

  if(verbose>0) {
    fprintf(stdout,"[%s] initial values: \n",fcname);
    fprintf(stdout,"  GM,WM,CSF      %5f %5f %5f\n",*iGM,*iWM,*iCSF);
    fprintf(stdout,"  ICS            %5f %5f\n",*intICS,*ranICS);
    fprintf(stdout,"  OCS            %5f %5f\n",*intOCS,*ranOCS);
  }

  if(*iCSF>0.0 && *iWM>0.0 && *iGM>0.0 && *iBrain>0.0 && *intICS>0.0 && *intOCS>0.0 && *ranICS>0.0 && *ranOCS>0.0)
    return 1;

  if(verbose>=1) niik_fc_display(fcname,1);

  NIIK_RET0((img->nx!=brain_mask->nx),fcname,"brain mask width");
  NIIK_RET0((img->nx!=csf_mask->nx),fcname,"ventricle mask width");
  NIIK_RET0((img->nx!=GWI_mask->nx),fcname,"GWI mask width");
  NIIK_RET0((img->ny!=brain_mask->ny),fcname,"brain mask height");
  NIIK_RET0((img->ny!=csf_mask->ny),fcname,"ventricle mask height");
  NIIK_RET0((img->ny!=GWI_mask->ny),fcname,"GWI mask height");
  NIIK_RET0((img->nz!=brain_mask->nz),fcname,"brain mask depth");
  NIIK_RET0((img->nz!=csf_mask->nz),fcname,"ventricle mask depth");
  NIIK_RET0((img->nz!=GWI_mask->nz),fcname,"GWI mask depth");
  NIIK_RET0((brain_mask->datatype!=NIFTI_TYPE_UINT8),fcname,"brain mask is not uint8");

  if(*iBrain<0)
    *iBrain = niik_image_get_mean(img, brain_mask);

  if(verbose>1) fprintf(stdout,"    Brain = %6.2f\n",*iBrain);

  if(*iCSF>0.0 && *iWM>0.0 && *iGM>0.0 && *iBrain>0.0 && *intICS>0.0 && *intOCS>0.0 && *ranICS>0.0 && *ranOCS>0.0) {
    return 1;
  }

  bseg = (unsigned char *)brain_mask->data;
  /* dmax = niik_image_get_percentile(img,brain_mask,0.98); */
  dmax = niik_image_get_max(img, brain_mask);
  dmin = niik_image_get_min(img, brain_mask);
  /* dmax = niik_image_get_max(img,NULL);  */
  if(verbose>1) {
    fprintf(stdout,"[%s] min = %7.2f max = %7.2f\n",fcname,dmin,dmax);
  }
  if(*iCSF<=0.0) /*VF niik_image_get_mode call is fishy*/
    *iCSF = niik_image_get_mode(img,csf_mask,dmin,dmax,128,3);/*(int)(dmax-dmin)/3.0,NIIK_IMIN(10,dmax/40.0)*/

  if(verbose>1) fprintf(stdout,"    CSF = %6.2f\n",*iCSF);
  if(verbose>1) niik_image_write("tmp_csf_mask.mnc",csf_mask);

  if(*iCSF>0.0 && *iWM>0.0 && *iGM>0.0 && *iBrain>0.0 && *intICS>0.0 && *intOCS>0.0 && *ranICS>0.0 && *ranOCS>0.0) {
    return 1;
  }

  /* make a temp image for WM intensity */
  NIIK_RET0(((tmpimg=niik_image_copy_as_type(GWI_mask, NIFTI_TYPE_UINT8))==NULL),fcname,"niik_image_copy_as_type for GWI_mask -> tmpimg");
  bimg = (unsigned char *)tmpimg->data;

  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(img,i)<*iBrain) {
      bimg[i]=0;
    }
  }

  if(*iWM<=0.0)
    *iWM = niik_image_get_mode(img,tmpimg,dmin,dmax,128,3);
  if(verbose>1) fprintf(stdout,"    WM = %6.2f\n",*iWM);
  if(verbose>1) niik_image_write("tmp_wm_mask.mnc",tmpimg);

  /* reset tmp as brain mask */
  NIIK_RET0((!niik_image_copy_data(brain_mask,tmpimg)),fcname,"niik_image_copy_data");
  bimg = (unsigned char *)tmpimg->data;

  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(img,i)>*iBrain) {
      bimg[i]=0;
    }
  }
  NIIK_RET0((!niik_image_maskout(tmpimg,csf_mask)),fcname,"niik_image_maskout tmpimg,venmask");

  if(*iGM<=0.0)
    *iGM = niik_image_get_mode(img,tmpimg,dmin,dmax,128,3);
  if(verbose>1) fprintf(stdout,"    GM = %6.2f\n",*iGM);
  if(verbose>1) niik_image_write("tmp_gm_mask.mnc",tmpimg);

  if(verbose) fprintf(stdout,"    [CSF,GM,WM,Brain] = [%.2f,%.2f,%.2f,%.2f]\n",*iCSF,*iGM,*iWM,*iBrain);

  if(*intICS >= 0.0 && *intOCS >= 0.0 &&
      *ranICS >= 0.0 && *ranOCS) {
    if(verbose>0) {
      fprintf(stdout,"    WM surface      %8.2f  +/-  %8.2f  user defined\n",*intICS,*ranICS);
      fprintf(stdout,"    Pial surface    %8.2f  +/-  %8.2f  user defined\n",*intOCS,*ranOCS);
    }
    return 1;
  }

  /* SURFACE MEAN INTENSITIES */
  NIIK_RET0(((gradimg=niik_image_copy_as_type(img, NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type gradimg");
  NIIK_RET0((!niik_image_filter_gaussian_update(gradimg,9,1.0)),fcname,"niik_image_filter_gaussian_update");
  NIIK_RET0((!niik_image_sobel_filter_update(gradimg,'m')),fcname,"niik_image_sobel_filter_update");
  gthresh = niik_image_get_percentile(gradimg, brain_mask, 0.6);
  if(verbose>2) {
    fprintf(stdout,"\t  gradient threshold = %8.2f\n",gthresh);
    niik_image_write("tmp_gradimg.mnc",gradimg);
  }

  /* ESTIMATED WM SURFACE INTENSITY */
  thresh = *iGM;
  bimg = (unsigned char *) tmpimg->data;
  for(i=nvox=0; i<img->nvox; i++) {
    bimg[i] = bseg[i];
    if( niik_image_get_voxel(img,i)<thresh ) bimg[i]=0;
    else if( niik_image_get_voxel(gradimg,i)<gthresh ) bimg[i]=0;
    nvox+=(bimg[i]>0);
  }
  if(*intICS<0) {
    /* *intICS = niik_image_get_mode(img,tmpimg,0,dmax,dmax/3.0,NIIK_IMIN(10,dmax/40));*/
    /* *intICS = niik_image_get_median(img,tmpimg); */
    /* *intICS = 0.8 * (*iWM) + 0.2 * (*iGM); VF: Kunio set 0.8*WM+0.2*GM  ???*/
    *intICS = mixICS * (*iWM) + (1.0-mixICS) * (*iGM);
  }
  if(*ranICS<0) {
    /**ranICS = fabs(*iWM - *intICS) * 0.5;*/ /*VF: which one ?*/
    *ranICS = 0.5 * (*iWM - *iGM); /* VF: half of the difference? */
  }
  if(verbose) fprintf(stdout,"    WM surface      %8.2f  +/-  %8.2f   %9i\n", *intICS, *ranICS, nvox);
  if(verbose>1) niik_image_write("tmp_wm_ics.mnc", tmpimg);

  /* ESTIMATE PIAL SURFACE INTENSITY */
  NIIK_RET0(((gwi_dil=niik_image_copy_as_type(GWI_mask, NIFTI_TYPE_UINT8))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET0((!niik_image_morph_3d_radius_mask(gwi_dil, NULL, NIIK_MORPH_DILATE,1.24)),fcname,"niik_image_morph_3d_radius_mask");

  NIIK_RET0(((brain_dil=niik_image_copy_as_type(brain_mask, NIFTI_TYPE_UINT8))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET0((!niik_image_morph_3d_radius_mask(brain_dil, NULL, NIIK_MORPH_DILATE,1.7)),fcname,"niik_image_morph_3d_radius_mask");

  thresh = *iCSF;
  bimg = (unsigned char *) tmpimg->data;
  if(verbose) {
    fprintf(stdout,"\tthresh  %9.4f\n",thresh);
    fprintf(stdout,"\tintICS  %9.4f\n",*intICS);
    fprintf(stdout,"\tgthresh %9.4f\n",gthresh);
  }
  for(i=nvox=0; i<img->nvox; i++) {
    bimg[i]=(niik_image_get_voxel(brain_dil,i)>0);
    if     (niik_image_get_voxel(img,i)< thresh) bimg[i]=0;
    else if(niik_image_get_voxel(img,i)>*intICS) bimg[i]=0;
    else if(niik_image_get_voxel(gradimg,i)<gthresh) bimg[i]=0;
    else if(niik_image_get_voxel(gwi_dil,i)>      0) bimg[i]=0;
    nvox+=(bimg[i]>0);
  }
  if(nvox==0) {
    fprintf(stderr,"ERROR: missing mask\n");
    exit(0);
  }

  gwi_dil = niik_image_free( gwi_dil );
  brain_dil = niik_image_free( brain_dil );
  gradimg = niik_image_free( gradimg );

  if(verbose>1) niik_image_write("tmp_ocs_mask.mnc",tmpimg);
  if( *intOCS ) {
    NIIK_RET0((niik_image_count_mask(tmpimg)==0),fcname,"missing mask for intOCS");
    /*niik_image_write("tmp_tmpimg.nii.gz",tmpimg);*/
    *intOCS = niik_image_get_mode(img,tmpimg,0,dmax,dmax/3.0,NIIK_IMIN(10,dmax/40.0));
    /* *intOCS = niik_image_get_median(img,tmpimg);*/
    /**intOCS = 0.7 * (*iGM) + 0.3 * (*iCSF); VF: Kunio's original split*/
    *intOCS = mixOCS * (*iGM) + (1.0-mixOCS) * (*iCSF);
  }
  if(*ranOCS) {
    /* *ranOCS = fabs(*iGM - *intOCS) * 0.5; */
    *ranOCS = 0.5 * (*iGM) - 0.5 * (*iCSF);
  }
  if(verbose>0) fprintf(stdout,"    Pial surface    %8.2f  +/-  %8.2f   %9i\n",*intOCS, *ranOCS, nvox);

  tmpimg = niik_image_free(tmpimg);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_estimate_tissue_values */


/* -bb is the bounding box, need to be initialized but not needed to be updated
 * -ics and ocs are the cortical surfaces
 * -returns the number of intersections
 * -if failure, returns a negative integer */
int niikcortex_off_count_intersection(bbox *bb,kobj *ics, kobj *ocs)
{
  kobj *ctx[2];
  kface *f[2];
  int n,num[2],
      verbose=niik_verbose();
  ctx[0]=ics;
  ctx[1]=ocs;

  if(!off_create_bbox_from_multiple_kobj(bb,ctx,2)) {
    fprintf(stderr,"ERROR: off_create_bbox_from_multiple_kobj(bb,ctx,2)\n");
    return -1;
  }

  if(!off_check_self_intersection_from_bbox(bb)) {
    fprintf(stderr,"ERROR: off_check_self_intersection_from_bbox(bb)\n");
    return -1;
  }

  for(n=0; n<2; n++) {
    for(num[n]=0,f[n]=ctx[n]->face; f[n]!=NULL; f[n]=f[n]->next) {
      if(f[n]->color[0] > f[n]->color[1]) num[n]++;
    }
  }

  n=num[0]+num[1];
  if(n>0 && verbose>0) {
    fprintf(stdout,"[niikcortex_off_count_intersection] intersection %i = %i %i \n",n,num[0],num[1]);
  }
  return n;
}


/**********************************************
* correct self intersection after remeshing
*
**********************************************/
int niikcortex_off_correct_self_intersections(bbox *bb,kobj *obj[],int maxiter)
/* -corrects (removes) surface intersections
 * -smoothes (relocates vertices)
 * -remeshes if necessary
 */
{
  kface **flist,*f;
  kvert *v;
  int iter = 0,
      cnt,n,m,i;
  int verbose = 0;
  int nei=0; /*include neighbours*/
  NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,obj,2)),fcname,"off_create_bbox_from_multiple_kobj");

  cnt =   niikcortex_off_count_intersection(bb, obj[0], obj[1]);

  while(cnt>0 && iter<maxiter) {
    if(verbose>1) fprintf(stdout,"[off_correct_self_intersection] iter %3i   intersections %-i\n",iter,cnt);
    iter++;
    flist = off_face_list(cnt);

    /* find the surface with self-intersections */
    for(i=0,n=0;i<2;i++) {
      for(f=obj[i]->face; f!=NULL; f=f->next) {
        if(f->color[0]>f->color[1]) {
          flist[n++]=f;
        }
      }
    }

    /* smoothes with relocate function */
    for(n=0; n<cnt; n++) {
      for(m=0; m<3; m++) {
        off_remesh_kvert_relocate(flist[n]->vert[m]);
      }
    }

    free(flist);

    cnt = niikcortex_off_count_intersection(bb, obj[0], obj[1]);
    if(iter%20==19) {
      NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,obj,2)),fcname,"off_create_bbox_from_multiple_kobj");
     }
  } /* while (cnt>0) */
  if(cnt>0)
    fprintf(stderr,"[%s] correcting self intersections failed\n",__func__);
  return 1;
}



/*********************************************************
 *
 * nonlinearly transform CLADA image to refimg
 *
 *
 ********************************************************/

nifti_image *niikcortex_gwi_mask_get_CLADA_mask(char *filename,nifti_image *refimg,nifti_image *warpimg,int interp)
/* -returns nonlinearly transformed image
 * -the file should be in $NIIKDIR/data/CLADA/
 * -refimg is the reference image
 * -warpimg is a location map (not a displacement map)
 * -interp is the interpolation */
{
  nifti_image *outimg;
  char fname[4096];
  const char *NIIKDIR;
  int verbose=niik_verbose();
  if(verbose>=1) niik_fc_display("niikcortex_gwi_mask_get_CLADA_mask",1);
  if(( NIIKDIR = get_NIIKDIR() )==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    fprintf(stderr,"       for example,\n");
    fprintf(stderr,"       export niikdir=/knakamura/kproj\n");
    fprintf(stderr,"       export niikdir=/Users/kunio/kproj\n");
    fprintf(stderr,"       export niikdir=/home/kunio/kproj\n");
    return NULL;
  }
  if(filename==NULL) {
    fprintf(stderr,"ERROR: filename is null\n");
    return NULL;
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return NULL;
  }
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg is null\n");
    return NULL;
  }
  sprintf(fname,"%s/data/CLADA/%s",NIIKDIR,filename);
  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask_get_CLADA_mask] reading %s\n",fname);
  if((outimg = nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: reading %s\n",fname);
    return NULL;
  }
  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask_get_CLADA_mask] warp image with %s interpolation\n",niik_interpolate_string(interp));
  if(!niik_image_apply_3d_warp_update(outimg,refimg,warpimg,NIIK_WARP_MAP_LOC,interp)) {
    fprintf(stderr,"ERROR: niik_image_3d_warp_update(outimg,img,warpimg,interp)\n");
    return NULL;
  }
  if(verbose>=1) niik_fc_display("niikcortex_gwi_mask_get_CLADA_mask",0);
  return outimg;
}

/************************************************
 *
 * niikcortex_nonctx
 *
 ************************************************/

int *niikcortex_non_cortex_label(int *ctx_label,nifti_image *nctx_img,kobj *ics, kobj *ocs,double dran,double delta,int iter)
/* nonctx is zero and
 * ctx is nonzero
 */
{
  int
    i,
    m,n,num;
  kvert
    **vlist,
    *vi,*vo;
  double
    cth,dval;
  int verbose=niik_verbose();
  const char *fcname="niikcortex_non_cortex_label";

  if(verbose>=1) {
    niik_fc_display("niikcortex_non_cortex_label",1);
    fprintf(stdout,"  non cortex mask         %s\n",nctx_img->fname);
    fprintf(stdout,"  white surface           %s\n",ics->fname);
    fprintf(stdout,"  pial surface            %s\n",ocs->fname);
    fprintf(stdout,"  range                   %-3.4f\n",dran);
    fprintf(stdout,"  delta                   %-3.4f\n",delta);
    fprintf(stdout,"  iter                    %i\n",iter);
  }
  NIIK_RET0(((ctx_label=(int *)realloc(ctx_label,ics->nvert*sizeof(int)))==NULL),fcname,"calloc for ctx_label");
  NIIK_RET0(((vlist = (kvert **)calloc(ics->nvert,sizeof(kvert *)))==NULL),fcname,"calloc for vlist");
  if(verbose>=2) fprintf(stdout,"[%s] interp\n",fcname);
  for(vi=ics->vert,vo=ocs->vert,num=n=0; vi!=NULL; vi=vi->next,vo=vo->next,n++) {
    ctx_label[n] = 1;
    if     (niik_image_interpolate_3d(nctx_img,vi->v,NIIK_INTERP_NN)>0) ctx_label[n] = 0;
    else if(niik_image_interpolate_3d(nctx_img,vo->v,NIIK_INTERP_NN)>0) ctx_label[n] = 0;
    else {
      cth = niikpt_distance(vi->v,vo->v);
      for(dval=-dran; dval<cth+dran; dval+=delta) {
        if(niik_image_interpolate_3d(nctx_img,niikpt_move_normal(vi->v,niikpt_unit(niikpt_sub(vo->v,vi->v)),dval),NIIK_INTERP_NN)>0)
          ctx_label[n] = 0;
      }
    }
    if(ctx_label[n]==0)
      vlist[num++] = vi;
  } /* each vertex */

  /* erosion */
  if(verbose>=2) fprintf(stdout,"[%s] erosion\n",fcname);
  for(i=0; i<iter; i++) {
    for(n=0; n<num; n++) {
      for(m=0; m<vlist[n]->nei; m++) {
        if(ctx_label[vlist[n]->neivert[m]->index-1]) break;
      }
      if(m<vlist[n]->nei) {
        vlist[n] = vlist[--num];
        n--;
      }
    }
  }
  /* dilation */
  if(verbose>=2) fprintf(stdout,"[%s] dilation\n",fcname);
  for(i=0; i<iter; i++) {
    num = off_get_local_kvert_check(vlist,ics->nvert,num);
  }
  /* set non-cortex */
  if(verbose>=2) fprintf(stdout,"[%s] set\n",fcname);
  for(i=0; i<ics->nvert; i++) {
    ctx_label[i]=1;
  }
  for(i=0; i<num; i++) {
    ctx_label[vlist[i]->index-1]=0;
  }
  /* free memory */
  free(vlist);
  if(verbose>=1) niik_fc_display(fcname,0);
  return ctx_label;
} /* niikcortex_non_cortex_label */

/* find non-cortex area */


/************************************************
 *
 * niikcortex_gwi_mask
 *
 ************************************************/

nifti_image *niikcortex_gwi_mask(nifti_image *img,
                                 nifti_image *brain_mask,
                                 nifti_image *ven_mask, double ven_dilate,
                                 nifti_image *wm_mask,
                                 double dgm_dilate,
                                 nifti_image *avoid_mask,
                                 nifti_image *lesion_mask,
                                 nifti_image *warpimg,double vessel_radius,
                                 double median_radius)
/* create gray matter - white matter interface mask
 *
 */

{
  nifti_image
  *msv_mask=NULL, /* midsagittal fill voxels around the lateral and third ventricles */
   *vsl_mask=NULL, /* bright vessels in sulci */
    *opw_mask=NULL, /* areas to open around the hippocampal folding */
     *dgm_mask=NULL, /* deep gray matter mask -dilated*/
      *dgm_fill=NULL, /* deep gray matter fill area -not dilated */
       *tmpimg=NULL,
        *outimg=NULL;
  unsigned char *bimg;
  const char *NIIKDIR=NULL;
  char fname[512];
  double openwm_radius=1.0;
  float *fimg;
  int oidx=1, i, verbose=niik_verbose();

  verbose = (g_niik_cortex_debug>0)?NIIK_IMAX(verbose,1):verbose;

  if(verbose>=1) niik_fc_display("niikcortex_gwi_mask",1);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(brain_mask==NULL) {
    fprintf(stderr,"ERROR: brain_mask is null\n");
    return NULL;
  }
  if(ven_mask==NULL) {
    fprintf(stderr,"ERROR: ven_mask is null\n");
    return NULL;
  }
  if(wm_mask==NULL) {
    fprintf(stderr,"ERROR: wm_mask is null\n");
    return NULL;
  }
  if(avoid_mask==NULL) {
    fprintf(stderr,"ERROR: avoid_mask is null\n");
    return NULL;
  }
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg is null\n");
    return NULL;
  }

  /* check image dimensions */
  if(niik_image_cmp_dim(img,brain_mask)) {
    fprintf(stderr,"ERROR: brain_mask dim\n");
    return NULL;
  }
  if(niik_image_cmp_dim(img,ven_mask)) {
    fprintf(stderr,"ERROR: ven_mask dim\n");
    return NULL;
  }
  if(niik_image_cmp_dim(img,wm_mask)) {
    fprintf(stderr,"ERROR: wm_mask dim\n");
    return NULL;
  }
  if(niik_image_cmp_dim(img,avoid_mask)) {
    fprintf(stderr,"ERROR: avoid_mask dim\n");
    return NULL;
  }
  if(lesion_mask!=NULL) {
    if(niik_image_cmp_dim(img,lesion_mask)) {
      fprintf(stderr,"ERROR: lesion_mask dim\n");
      return NULL;
    }
  }

  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    fprintf(stderr,"       for example,\n");
    fprintf(stderr,"       export niikdir=/knakamura/kproj\n");
    fprintf(stderr,"       export niikdir=/Users/kunio/kproj\n");
    fprintf(stderr,"       export niikdir=/home/kunio/kproj\n");
    return NULL;
  }

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] prepare output image\n");
  if((outimg = niik_image_copy(wm_mask))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy(wm_mask)\n");
    return NULL;
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  bimg = (unsigned char *)outimg->data;

  if((tmpimg = niik_image_copy(ven_mask))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy(wm_mask)\n");
    return NULL;
  }


  /****************************************
   *
   * 1. dilate/add ventricle mask
   *
   ****************************************/

  if(verbose  ) fprintf(stdout,"[niikcortex_gwi_mask] 1. ventricle mask\n");
  if(verbose>1) fprintf(stdout,"[niikcortex_gwi_mask] 1.a. dilate ventricle mask %6.2f\n",ven_dilate);
  if(ven_dilate>0) {
    if(!niik_image_morph_3d_radius_mask(tmpimg,NULL,NIIK_MORPH_DILATE,ven_dilate)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask(tmpimg,NULL,NIIK_MORPH_DILATE,ven_dilate)\n");
      return NULL;
    }
  }

  if(verbose>1) fprintf(stdout,"[niikcortex_gwi_mask] 1.b. add ventricle mask\n");
  for(i=0; i<img->nvox; i++) {
    bimg [i] = (bimg[i] || niik_image_get_voxel(tmpimg,i)>0);
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_addven.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }

  if(verbose>1) fprintf(stdout,"[niikcortex_gwi_mask] 1.c. close holes\n");
  if(!niik_image_close_holes(outimg)) {
    fprintf(stderr,"ERROR: niik_image_close_holes\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_venclose.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }


  /****************************************
   *
   * 2. dilate/add deep GM mask
   *
   * 2.a. dgm_mask   to be dilated
   * 2.b. dgm_fill   not to be dilated
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 2. deep gray matter mask\n");
  if((dgm_mask = niikcortex_gwi_mask_get_CLADA_mask("MNI152_T1_1mm_DGM_mask.nii.gz",img,warpimg,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niikcortex_gwi_mask_get_CLADA_mask(\"MNI152_T1_1mm_DGM_mask.nii.gz\",img,warpimg,NIIK_INTERP_NN)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_dgm_pre_dilate.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,dgm_mask);
  }

  if(verbose>1) fprintf(stdout,"[niikcortex_gwi_mask] 2.a. dilate deep GM mask %6.2f\n",dgm_dilate);
  if(dgm_dilate>0) {
    if(!niik_image_morph_3d_radius_mask(dgm_mask,NULL,NIIK_MORPH_DILATE,dgm_dilate)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask(dgm_mask,NULL,NIIK_MORPH_DILATE,dgm_dilate)\n");
      return NULL;
    }
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_dgm_post_dilate.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,dgm_mask);
  }

  if(verbose>1) fprintf(stdout,"[niikcortex_gwi_mask] 2.b. add deep GM mask\n");
  for(i=0; i<img->nvox; i++) {
    bimg [i] = (bimg[i] || niik_image_get_voxel(dgm_mask,i)>0);
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_adddgm.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }
  dgm_mask = niik_image_free(dgm_mask);

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 2.c. add deep GM mask no dilation\n");
  if((dgm_fill = niikcortex_gwi_mask_get_CLADA_mask("MNI152_T1_1mm_DGM_fill.nii.gz",img,warpimg,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niikcortex_gwi_mask_get_CLADA_mask(\"MNI152_T1_1mm_DGM_fill.nii.gz\",img,warpimg,NIIK_INTERP_NN)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_dgm_no_dilate.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,dgm_fill);
  }

  for(i=0; i<img->nvox; i++) {
    bimg [i] = (bimg[i] || niik_image_get_voxel(dgm_fill,i)>0);
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_adddgm.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }
  dgm_fill = niik_image_free(dgm_fill);


  /****************************************
   *
   * 3. remove cerebellum & brainstem AKA avoid mask
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 3. remove cerebellum & brainstem mask\n");
  for(i=0; i<img->nvox; i++) {
    bimg [i] = (bimg[i] && niik_image_get_voxel(avoid_mask,i)==0);
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_rmavoid.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }


  /****************************************
   *
   * 5. seed fill
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 5. seed fill\n");
  if(!niik_image_copy_data(ven_mask,tmpimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(ven_mask,tmpimg)\n");
    return NULL;
  }
  if(!niik_image_seed_fill(outimg,tmpimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill(outimg,tmpimg,0)\n");
    return NULL;
  }
  if(!niik_image_copy_data(tmpimg,outimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(tmpimg,outimg)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_seedfill.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }


  /****************************************
   *
   * 7. remove bright blood vessels
   *
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 7. remove bright blood vessels\n");
  if((vsl_mask = niikcortex_gwi_mask_get_CLADA_mask("MNI152_T1_1mm_bright_vessel_roi.nii.gz",img,warpimg,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niikcortex_gwi_mask_get_CLADA_mask(\"MNI152_T1_1mm_bright_vessel_roi.nii.gz\",img,warpimg,NIIK_INTERP_NN)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_bright_roi.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,vsl_mask);
  }

  niik_image_type_convert(vsl_mask,NIFTI_TYPE_FLOAT32);
  fimg=(float *)vsl_mask->data;
  for(i=0; i<vsl_mask->nvox; i++) {
    fimg[i] *= vessel_radius;
  }
  if(!niik_image_filter_gaussian_update(vsl_mask,9,vessel_radius)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }

  if(!niik_image_morph_3d_radius_map(outimg,vsl_mask,NIIK_MORPH_ERODE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(outimg,vsl_mask,NIIK_MORPH_ERODE)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_verode.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }

  /* OPENING AROUND THE HIPPOCAMPAL SULCI */
  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 7b. remove white matter from hippocampal folding\n");
  if((opw_mask = niikcortex_gwi_mask_get_CLADA_mask("MNI152_T1_1mm_DGM_open_ROI.nii.gz",img,warpimg,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niikcortex_gwi_mask_get_CLADA_mask(\"MNI152_T1_1mm_DGM_open_ROI.nii.gz\",img,warpimg,NIIK_INTERP_NN)\n");
    return NULL;
  }
  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_open_dgm_roi.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,vsl_mask);
  }
  niik_image_type_convert(opw_mask,NIFTI_TYPE_FLOAT32);
  fimg=(float *)opw_mask->data;
  for(i=0; i<opw_mask->nvox; i++) {
    fimg[i] *= openwm_radius;
  }
  if(!niik_image_filter_gaussian_update(opw_mask,9,openwm_radius)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  if(!niik_image_morph_3d_radius_map(outimg,opw_mask,NIIK_MORPH_ERODE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(outimg,opw_mask,NIIK_MORPH_ERODE)\n");
    return NULL;
  }
  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_remove_hippo_wm.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }
  opw_mask=niik_image_free(opw_mask);

  if(!niik_image_copy_data(ven_mask,tmpimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(ven_mask,tmpimg)\n");
    return NULL;
  }
  if(!niik_image_seed_fill(outimg,tmpimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill(outimg,tmpimg,0)\n");
    return NULL;
  }
  if(!niik_image_copy_data(tmpimg,outimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(tmpimg,outimg)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_vseedfill.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }

  if(!niik_image_morph_3d_radius_map(outimg,vsl_mask,NIIK_MORPH_DILATE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(outimg,vsl_mask,NIIK_MORPH_ERODE)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_vdilate.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }

  vsl_mask = niik_image_free(vsl_mask);


  /****************************************
   *
   * 8. add mid-sagittal voxels
   *
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 8. add mid-sagittal voxels\n");
  if((msv_mask = niikcortex_gwi_mask_get_CLADA_mask("MNI152_T1_2mm_mid_fill_mask.nii.gz",
                 img,warpimg,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niikcortex_gwi_mask_get_CLADA_mask(\"MNI152_T1_2mm_mid_fill_mask.nii.gz\",img,warpimg,NIIK_INTERP_NN)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_add_midsag.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,msv_mask);
  }


  if(!niik_image_type_convert(msv_mask,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(msv_mask,NIFTI_TYPE_UINT8)\n");
    return NULL;
  }

  bimg = outimg->data;
  for(i=0; i<img->nvox; i++) {
    bimg[i] = (bimg[i]>0 || niik_image_get_voxel(msv_mask,i)>0);
  }

  msv_mask = niik_image_free(msv_mask);


  /****************************************
   *
   * 9. median filter
   *
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 9. median filter %6.2f\n",median_radius);

  if(!niik_image_filter_median_radius(outimg,NULL,median_radius)) {
    fprintf(stderr,"ERROR: niik_image_filter_median_radius\n");
    return NULL;
  }

  if(!niik_image_copy_data(ven_mask,tmpimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(ven_mask,tmpimg)\n");
    return NULL;
  }
  if(!niik_image_seed_fill(outimg,tmpimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill(outimg,tmpimg,0)\n");
    return NULL;
  }
  if(!niik_image_copy_data(tmpimg,outimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data(tmpimg,outimg)\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_fin.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }


  /****************************************
   *
   * 10. close holes
   *
   *
   ****************************************/

  if(verbose) fprintf(stdout,"[niikcortex_gwi_mask] 10. close holes\n");
  if(!niik_image_close_holes(outimg)) {
    fprintf(stderr,"ERROR: niik_image_close_holes\n");
    return NULL;
  }

  if(verbose && g_niik_cortex_debug) {
    sprintf(fname,"tmp_niikcortex_gwi_mask%i_fin.nii.gz",oidx++);
    fprintf(stdout,"[niikcortex_gwi_mask] writing %s\n",fname);
    niik_image_write(fname,outimg);
  }

  tmpimg = niik_image_free(tmpimg);
  if(verbose>=1) niik_fc_display("niikcortex_gwi_mask",0);
  return outimg;
} /* nifti_image *niikcortex_gwi_mask */



nifti_image *niikcortex_modify_image(nifti_image *img,
                                     nifti_image *brain_mask,
                                     nifti_image *ven_mask, double ven_dilate, double ven_blur,
                                     nifti_image *avoid_mask,
                                     nifti_image *dgm_mask, double dgm_dilate, double dgm_blur,
                                     nifti_image *lesion_mask,
                                     double fillval)
/*
 * images include:
 *   t1g, brainmask, ventricle mask, avoid mask, deep GM mask
 */
{
  nifti_image
  *tmpimg,
  *outimg=NULL;
  niikpt
  p;
  int
  verbose=niik_verbose(),
  i,j,k,n;
  float *fimg;

  if(verbose>=1) niik_fc_display("niikcortex_modify_image",1);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  fimg = outimg->data;

  /* ventricle mask */
  if(verbose) fprintf(stdout,"[niikcortex_modify_image] ventricle mask %5.2f %5.2f\n",ven_dilate,ven_blur);
  if((tmpimg=niik_image_copy(ven_mask))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_morph_3d_radius(tmpimg,NIIK_MORPH_DILATE,ven_dilate)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
    return NULL;
  }
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(tmpimg,ven_blur*2,ven_blur)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  for(k=n=0; k<img->nz; k++) {
    p.z = k * img->dz;
    for(j=0; j<img->ny; j++) {
      p.y = j * img->dy;
      for(i=0; i<img->nx; n++,i++) {
        p.x = i * img->dx;
        fimg[n] = NIIK_FMAX(fimg[n],niik_image_interpolate_3d(tmpimg,p,NIIK_INTERP_LINEAR) * fillval);
      }
    }
  }
  tmpimg = niik_image_free(tmpimg);

  /* deep gray matter */
  if(verbose) fprintf(stdout,"[niikcortex_modify_image] deep gray matter mask %5.2f %5.2f\n",dgm_dilate,dgm_blur);
  if((tmpimg=niik_image_copy(dgm_mask))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_morph_3d_radius(tmpimg,NIIK_MORPH_DILATE,dgm_dilate)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
    return NULL;
  }
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(tmpimg,dgm_blur*2,dgm_blur)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  for(k=n=0; k<img->nz; k++) {
    p.z = k * img->dz;
    for(j=0; j<img->ny; j++) {
      p.y = j * img->dy;
      for(i=0; i<img->nx; n++,i++) {
        p.x = i * img->dx;
        fimg[n] = NIIK_FMAX(fimg[n],niik_image_interpolate_3d(tmpimg,p,NIIK_INTERP_LINEAR) * fillval);
      }
    }
  }
  tmpimg = niik_image_free(tmpimg);

  /* everything */
  for(k=n=0; k<img->nz; k++) {
    p.z = k * img->dz;
    for(j=0; j<img->ny; j++) {
      p.y = j * img->dy;
      for(i=0; i<img->nx; n++,i++) {
        p.x = i * img->dx;
        fimg[n] = NIIK_FMAX(fimg[n],niik_image_interpolate_3d(ven_mask,p,NIIK_INTERP_LINEAR) * fillval);
        fimg[n] = NIIK_FMAX(fimg[n],niik_image_interpolate_3d(dgm_mask,p,NIIK_INTERP_LINEAR) * fillval);
        fimg[n] = fimg[n] * (1.0 - niik_image_interpolate_3d(avoid_mask,p,NIIK_INTERP_LINEAR));
      }
    }
  }
  if(verbose>=1) niik_fc_display("niikcortex_modify_image",0);
  return outimg;
} /* nifti_image *niikcortex_modify_image */


nifti_image *niikcortex_wm_mask(nifti_image *img,nifti_image *maskimg,nifti_image *warpimg,double mag,double thresh)
/* white mater segmentation for niikcortex
 * 1. threshold image
 *    -lothresh area is pre-defined
 * 2. mask with brain
 */
{
  nifti_image
  *lo_roi=NULL,
   *wmmask=NULL;
  const char * NIIKDIR=NULL;
  char
  fname[4096],
        fcname[64]="niikcortex_wm_mask";
  int
  i,
  verbose=niik_verbose();
  double sFWHM=3.5;
  double t2;

  if(verbose>=1) niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg is null\n");
    return NULL;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return NULL;
  }
  if(niik_image_cmp_dim(img,maskimg)) {
    fprintf(stderr,"ERROR: niik_image_cmp_dim\n");
    return NULL;
  }

  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv NIIKDIR\n",fcname);
    return NULL;
  }

  sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_gwi_lothresh_area.nii.gz",NIIKDIR);
  if(verbose>=1) fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
  if((lo_roi=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: reading %s\n",fname);
    return NULL;
  }

  if(verbose>=1) fprintf(stdout,"[%s] warp low thresh image\n",fcname);
  if(!niik_image_apply_3d_warp_update(lo_roi,img,warpimg,NIIK_WARP_MAP_LOC,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_apply_3d_warp_update\n");
    return NULL;
  }

  if(verbose>=1) fprintf(stdout,"[%s] gaussian filter %5.2f\n",fcname,sFWHM);
  if(!niik_image_type_convert(lo_roi,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(lo_roi,9,sFWHM)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  if(0) niik_image_write("tmp_lotarea.nii.gz",lo_roi);

  if(verbose>=1) fprintf(stdout,"[%s] threshold image  factor %5.2f\n",fcname,mag);
  if((wmmask=niik_image_threshold_new(img,thresh))==NULL) {
    fprintf(stderr,"ERROR: niik_image_threshold_new\n");
    return NULL;
  }

  for(i=0; i<img->nvox; i++) {
    t2 = (1.0-mag*niik_image_get_voxel(lo_roi,i)) * thresh;
    if(niik_image_get_voxel(img,i)>t2)
      niik_image_set_voxel(wmmask,i,1);
    else
      niik_image_set_voxel(wmmask,i,0);
  }

  if(!niik_image_mask(wmmask,maskimg)) {
    fprintf(stderr,"ERROR: niik_image_mask\n");
    return NULL;
  }

  if(0) niik_image_write("tmp_wmmask.nii.gz",wmmask);

  if(verbose>=1) niik_fc_display(fcname,0);
  return wmmask;
} /* niikcortex_wm_mask */



/*******************************************************
 *
 * FOV-corrected averaging
 *
 *
 *******************************************************/

int niik_image_average_with_fov(nifti_image *refimg,nifti_image **imglist,niikmat **matlist,int num,int interp)
/*
 * average images with field-of-view weighting
 */
{
  nifti_image
  *fovimg,
  *tmpimg;
  niikmat
  *afmat;
  niikpt
  pfov,p,q;
  int
  m,n,i,j,k;
  double *dimg,*dfov;
  char fcname[64]="niik_image_average_with_fov";
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return 0;
  }
  if(matlist==NULL) {
    fprintf(stderr,"ERROR: matlist is null\n");
    return 0;
  }
  fprintf(stdout,"[niik_image_average_with_fov] %s \n",niik_interpolate_string(interp));
  NIIK_RET0(((dimg = (double *)calloc(refimg->nvox,sizeof(double)))==NULL),fcname,"calloc for dimg");
  NIIK_RET0(((dfov = (double *)calloc(refimg->nvox,sizeof(double)))==NULL),fcname,"calloc for dfov");
  niik_fc_display(fcname,1);
  for(n=0; n<num; n++) {
    if((afmat = niikmat_inverse(matlist[n]))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse, %i \n",n);
      return 0;
    }
    if((fovimg = niik_image_copy(imglist[n]))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(!niik_image_one(fovimg)) {
      fprintf(stderr,"ERROR: niik_image_one\n");
      return 0;
    }
    if(interp==NIIK_INTERP_BSPLINE) {
      if((tmpimg = niik_image_copy(imglist[n]))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy\n");
        return 0;
      }
      if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT64)) {
        fprintf(stderr,"ERROR: niik_image_type_convert, %i\n",n);
        return 0;
      }
      if(!niik_image_interpolate_convert_3d_bspline_coeff(tmpimg)) {
        fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff, %i\n",n);
        return 0;
      }
    } else {
      tmpimg = imglist[n];
    }
    pfov.x=(imglist[n]->nx-1)*imglist[n]->dx;
    pfov.y=(imglist[n]->ny-1)*imglist[n]->dy;
    pfov.z=(imglist[n]->nz-1)*imglist[n]->dz;

    for(k=0; k<refimg->nz; k++) {
      p.z = k * refimg->dz;
      m=k*refimg->nx*refimg->ny;
      for(j=0; j<refimg->ny; j++) {
        p.y = j * refimg->dy;
        for(i=0; i<refimg->nx; m++,i++) {
          p.x = i * refimg->dx;
          q = niikpt_affine_transform(afmat,p);
          if(interp==NIIK_INTERP_BSPLINE) {
            if     (q.x<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.y<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.z<1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.x>pfov.x-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.y>pfov.y-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else if(q.z>pfov.z-1) {
              dimg[m] += niik_image_interpolate_3d(imglist[n],q,NIIK_INTERP_LINEAR);
            } else {
              dimg[m] += niik_image_interpolate_3d(tmpimg,q,interp);
            }
          } else {
            dimg[m] += niik_image_interpolate_3d(tmpimg,q,interp);
          }
          dfov[m] += niik_image_interpolate_3d(fovimg,q,NIIK_INTERP_LINEAR);
          /*if(i==83 && j==201 && k==150) fprintf(stdout,"[%3i %3i %3i, %3i] %8.4f %5.2f \n",i,j,k,n,dimg[m],dfov[m]);*/
        }
      }
    }
    if(NIIK_INTERP_BSPLINE) tmpimg = niik_image_free(tmpimg);
    fovimg = niik_image_free(fovimg);
  }
  for(i=0; i<refimg->nvox; i++) {
    if(dfov[i]<0.25) dimg[i]=0;
    else dimg[i] = dimg[i] / dfov[i];
  }
  if(!niik_image_set_voxels_from_double_vector(refimg,dfov)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  /* niik_image_write("tmp_fov.nii.gz",refimg); */
  if(!niik_image_set_voxels_from_double_vector(refimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(dimg);
  free(dfov);
  niik_fc_display(fcname,0);
  return 1;
}


/*******************************************************
 *
 * cortical thickness calculation and visualization
 *
 *
 *******************************************************/

double niikcortex_calc_area_weighted_thickness(kobj *ics, kobj *ocs, double *thk,  int filter_type, int filter_size)
/* -calculates/returns the area-weighted cortical thickness
 * -calculates and updates thk
 * -uses niikcortex_calc_thickness
 */
{
  double
  *fai,*fao,
  sa=0,sc=0;
  kface *fi,*fo;
  kvert *vi,*vo;
  int n,fidx;
  if(ics==NULL) {
    fprintf(stderr,"ERROR: ics is null\n");
    return 0;
  }
  if(ocs==NULL) {
    fprintf(stderr,"ERROR: ocs is null\n");
    return 0;
  }
  if(thk==NULL) {
    fprintf(stderr,"ERROR: thk is null\n");
    return 0;
  }
  if(ics->nvert!=ocs->nvert) {
    fprintf(stderr,"ERROR: nvert is different %i %i\n",ics->nvert,ocs->nvert);
    return 0;
  }
  if(ics->nface!=ocs->nface) {
    fprintf(stderr,"ERROR: nface is different %i %i\n",ics->nface,ocs->nface);
    return 0;
  }
  fprintf(stdout,"[niikcortex_calc_area_weighted_thickness] vertex-wise measure thickness\n");
  if(!niikcortex_calc_thickness(ics, ocs, thk, NULL, NULL, filter_type, filter_size)) { /*TODO: add psi and the*/
    fprintf(stderr,"ERROR: niikcortex_calc_thickness(ics, ocs, thk, filter_type, filter_size)\n");
    return NIIKMAX;
  }
  fprintf(stdout,"[niikcortex_calc_area_weighted_thickness] calculate area\n");
  fai = niik_calloc_double_vector(ics->nface);
  fao = niik_calloc_double_vector(ocs->nface);
  for(fi=ics->face,fo=ocs->face,fidx=0; fi!=NULL; fi=fi->next,fo=fo->next,fidx++) {
    fai[fidx] = niikpt_area2(fi->vert[0]->v,fi->vert[1]->v,fi->vert[2]->v);
    fao[fidx] = niikpt_area2(fo->vert[0]->v,fo->vert[1]->v,fo->vert[2]->v);
  }
  fprintf(stdout,"[niikcortex_calc_area_weighted_thickness] calculate weighted average\n");
  for(vi=ics->vert,vo=ocs->vert; vi!=NULL; vi=vi->next,vo=vo->next) {
    for(n=0; n<vi->nei; n++) {
      sa +=  fai[vi->neiface[n]->index-1] + fao[vo->neiface[n]->index-1];
      sc += (fai[vi->neiface[n]->index-1] + fao[vo->neiface[n]->index-1]) * thk[vi->index-1];
    }
  }
  free(fai);
  free(fao);
  return sc/sa;
}

int niikcortex_calc_thickness(kobj *ics, kobj *ocs, double *thk, double *psi, double *the, int filter_type, int filter_size)
/* -calculates vertex-wise cortical thickness
 * -filter_type :  median or mean
 * -filter_size :  size of regional dilation
 * -thk         : ics->nvert-length vector; updated with a list of cortical thickness
 * -ics/ocs :  inner (white) and outer (pial) surfaces
 */
{
  kvert *vi,*vo,**vlist;
  int m,num,vidx;
  double *dlist,*dlist_phi,*dlist_the;
  int spherecoo=1;
  if(ics==NULL) {
    fprintf(stderr,"ERROR: ics is null\n");
    return 0;
  }
  if(ocs==NULL) {
    fprintf(stderr,"ERROR: ocs is null\n");
    return 0;
  }
  if(thk==NULL) {
    fprintf(stderr,"ERROR: thk is null\n");
    return 0;
  }
  if(ocs->spherecoo==0 || ics->spherecoo==0 || thk==NULL || psi==NULL)
    spherecoo=0;

  if(ics->nvert!=ocs->nvert) {
    fprintf(stderr,"ERROR: nvert is different %i %i\n",ics->nvert,ocs->nvert);
    return 0;
  }

  for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
    thk[vidx] = niikpt_distance(vi->v,vo->v);
  }
  if(spherecoo) {
    /*TODO: check for consistency between ics and ocs?*/
    for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
      psi[vidx] = vi->sph.psi;
      the[vidx] = vi->sph.the;
    }
  }

  vlist = (kvert **)calloc(ics->nvert,sizeof(kvert *));
  dlist = (double *)calloc(ics->nvert,sizeof(double));
//   dlist_phi = (double *)calloc(ics->nvert,sizeof(double));
//   dlist_the = (double *)calloc(ics->nvert,sizeof(double));

  switch(filter_type) {
  case 2: /* median filter */

    for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
      vlist[0]=vi;
      for(m=0,num=1; m<filter_size; m++) {
        num = off_get_local_kvert_check(vlist,ics->nvert,num);
      }

      for(m=0; m<num; m++) {
        dlist[m] = thk[vlist[m]->index-1];
      }
      vi->v.w = niik_median_quicksort_double(dlist,num);
    }

    for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
      thk[vidx] = vi->v.w;
    }

    if(spherecoo) { /*TODO: this is not optimal, I am redoing neigbour search two times extra*/
      for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
        vlist[0]=vi;
        for(m=0,num=1; m<filter_size; m++) {
          num = off_get_local_kvert_check(vlist,ics->nvert,num);
        }

        for(m=0; m<num; m++) {
          dlist[m] = psi[vlist[m]->index-1];
        }
        vi->v.w = niik_median_quicksort_double(dlist,num);
      }

      for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
        psi[vidx] = vi->v.w;
      }

      for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
        vlist[0]=vi;
        for(m=0,num=1; m<filter_size; m++) {
          num = off_get_local_kvert_check(vlist,ics->nvert,num);
        }

        for(m=0; m<num; m++) {
          dlist[m] = the[vlist[m]->index-1];
        }
        vi->v.w = niik_median_quicksort_double(dlist,num);
      }

      for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
        the[vidx] = vi->v.w;
      }
    }


    break;
  case 1: /* average filter */
    for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
      vlist[0]=vi;
      for(m=0,num=1; m<filter_size; m++) {
        num = off_get_local_kvert_check(vlist,ics->nvert,num);
      }
      for(m=0; m<num; m++) {
        dlist[m] = thk[vlist[m]->index-1];
      }
      vi->v.w = niik_get_mean_from_double_vector(dlist,num);
    }
    for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
      thk[vidx] = vi->v.w;
    }

    if(spherecoo) {      /*TODO: this is not optimal, I am redoing neigbour search two times extra*/
      for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
        vlist[0]=vi;
        for(m=0,num=1; m<filter_size; m++) {
          num = off_get_local_kvert_check(vlist,ics->nvert,num);
        }

        for(m=0; m<num; m++) {
          dlist[m] = psi[vlist[m]->index-1];
        }
        vi->v.w = niik_get_mean_from_double_vector(dlist,num);
      }

      for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
        psi[vidx] = vi->v.w;
      }

      for(vi=ics->vert,vo=ocs->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
        vlist[0]=vi;
        for(m=0,num=1; m<filter_size; m++) {
          num = off_get_local_kvert_check(vlist,ics->nvert,num);
        }

        for(m=0; m<num; m++) {
          dlist[m] = the[vlist[m]->index-1];
        }
        vi->v.w = niik_get_mean_from_double_vector(dlist,num);
      }
      for(vi=ics->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
        the[vidx] = vi->v.w;
      }
    }

    break;
  case 0: /* no filter at all*/
  default:
    break;
  }
  free(vlist);
  free(dlist);
//   free(dlist_phi);
//   free(dlist_the);
  return 1;
}


int niikcortex_add_thickness_color(kobj *ics, kobj *ocs, double *thk, double omin, double omax)
/* put color according to cortical thickness */
{
  kface *fi,*fo;
  double oran;
  int n,num=50;
  niikmat *cm;
  cm = niik_colormap_get(NIIK_COLORMAP_SPECTRAL,num);
  /*niikmat_display(cm);*/
  oran = (omax-omin)/(num-1.0);
  off_kobj_add_color(ics);
  off_kobj_add_color(ocs);
  for(fi=ics->face,fo=ocs->face; fi!=NULL; fi=fi->next,fo=fo->next) {
    n = ((thk[fo->vert[0]->index-1] + thk[fo->vert[1]->index-1] + thk[fo->vert[2]->index-1]) / 3.0 - omin) / oran;
    if(n<0) n=0;
    else if(n>=num) n=num-1;
    fo->color[0] = fi->color[0] = cm->m[n][0];
    fo->color[1] = fi->color[1] = cm->m[n][1];
    fo->color[2] = fi->color[2] = cm->m[n][2];
    fo->color[3] = fi->color[3] = 0;
  }
  cm = niikmat_free(cm);
  return 1;
}


int niikcortex_add_color(kobj *obj, double *var, double omin, double omax, int color_map_type,int color_levels)
/* put color according to cortical thickness */
{
  kface *f;
  double oran;
  int n;
  niikmat *cm;
  cm = niik_colormap_get(color_map_type, color_levels);
  oran = (omax-omin)/(color_levels-1.0);
  off_kobj_add_color(obj);
  for(f=obj->face; f!=NULL; f=f->next) {
    n = floor( ((var[f->vert[0]->index-1] + var[f->vert[1]->index-1] + var[f->vert[2]->index-1]) / 3.0 - omin) / oran + 0.5);
    if(n<0) n=0;
    else if(n>=color_levels) n=color_levels-1;
    f->color[0] = cm->m[n][0];
    f->color[1] = cm->m[n][1];
    f->color[2] = cm->m[n][2];
    f->color[3] = 0;
  }
  cm = niikmat_free(cm);
  return 1;
}

static int compare_dbl (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
    else return 0;
}

int niikcortex_add_color_discrete(kobj *obj, double *var, double omin, double omax, int color_map_type,int color_levels)
/* put color according to a value (without interpolating) */
{
  kface *f;
  double oran;
  int n;
  niikmat *cm;
  cm = niik_colormap_get(color_map_type, color_levels);
  oran = (omax-omin)/(color_levels-1.0);
  off_kobj_add_color(obj);
  for(f=obj->face; f!=NULL; f=f->next) {
    /*n = ((var[f->vert[0]->index-1] + var[f->vert[1]->index-1] + var[f->vert[2]->index-1]) / 3.0 - omin) / oran;*/
    double v[3]={var[f->vert[0]->index-1], var[f->vert[1]->index-1], var[f->vert[2]->index-1]};
    qsort(v,3,sizeof(double),compare_dbl);
   
    n=floor( ( (v[0]==v[1]?v[0]:v[1]==v[2]?v[2]:v[0])-omin) / oran + 0.5);

    if(n<0) n=0;
    else if(n>=color_levels) n=color_levels-1;
    f->color[0] = cm->m[n][0];
    f->color[1] = cm->m[n][1];
    f->color[2] = cm->m[n][2];
    f->color[3] = 0;
  }
  cm = niikmat_free(cm);
  return 1;
}



int niikcortex_make_fuzzy(nifti_image *input,nifti_image *output,
                          double value,double range) {
  int i;
  for(i=0; i<input->nvox; i++) {
    niik_image_set_voxel(output,i, NIIK_Heaviside11(niik_image_get_voxel(input,i)-value, range));
  }
  return 1;
}



/* end of 'cortical thickness calculation and visualization' */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/