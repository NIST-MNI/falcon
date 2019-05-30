/* FILENAME:     niik_nls_lesion.c
 * DESCRIPTION:  Kunio's nifti1 implementation
 *               using nonlocal segmentation (NLS) technique
 * AUTHOR:       Kunio Nakamura
 * DATE:         December 09, 2014
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

#define it1w  0 // t1w index and so on
#define iflr  1
#define it2w  2
#define ipdw  3
#define igad  4
#define iout  5
#define ncon    (5)       // # contrast 
#define nlibmax (4096)

#define NX 158

#define NIIK_NLS_FUSION_METHOD_MIN 0
#define NIIK_NLS_FUSION_METHOD_AVG 1
#define NIIK_NLS_FUSION_METHOD_MEDIAN 2

nifti_image *niik_nls_lesion2_patch_case(nifti_image **imgs,int nimg,int xpatch,int ypatch,int zpatch);

typedef struct {
  unsigned char ***tree;
  int *numlist;
  int ndim;
  int *dim;
  double *step;
  int num;
} kdtree;

kdtree *kdtree_free(kdtree *v) {
  int i;
  if(v==NULL) return NULL;
  if(v->tree!=NULL) {
    for(i=0; i<v->num; i++) {
      if(v->tree[i]!=NULL) free(v->tree[i]);
    }
    if(v->tree!=NULL) free(v->tree);
  }
  if(v->numlist!=NULL) free(v->numlist);
  if(v->dim!=NULL) free(v->dim);
  if(v->step!=NULL) free(v->step);
  free(v);
  return NULL;
}

typedef struct {
  nifti_image *flr;
  nifti_image *t1w;
  nifti_image *mask;
  nifti_image *patch;
  nifti_image **lib;
  nifti_image *pwm;
  nifti_image *pgm;
  int nlib;
  int ndim;   // = 150, 151, 154...
  char *lib_list_file;
  float *vtols; // list of tolerance values
  int *itols;   // index of tolerance values, bit confusing but correct (usually plus 1 to account for [0=lesion]
  int ntols;    // length of vtols and itols
  int xpatch;
  int ypatch;
  int zpatch;
  int topd;   // #top candidates to consider
  int verbose;
  int verbose_idx[4];
  int fusion_method;
  kdtree *tree;
} niik_nls_lesion2; // class

niik_nls_lesion2 *niik_nls_lesion2_init() {
  niik_nls_lesion2 *v=NULL;
  NIIK_RET0(((v=(niik_nls_lesion2 *)calloc(1,sizeof(niik_nls_lesion2)))==NULL),__func__,"calloc");
  v->fusion_method=NIIK_NLS_FUSION_METHOD_AVG;
  v->flr=v->t1w=v->mask=NULL;
  v->patch=NULL;
  v->lib=NULL;
  v->nlib=0;
  v->ndim=150;
  v->verbose=2;
  v->verbose_idx[0]=-1;
  v->lib_list_file=NULL;
  v->xpatch=v->ypatch=2;
  v->zpatch=1;
  v->topd=5;
  v->ntols=5;
  v->vtols=(float *)calloc(v->ntols,sizeof(float));
  v->itols=(int *)calloc(v->ntols,sizeof(int));
  // xyz coordinates
  v->vtols[0]=v->vtols[1]=v->vtols[2]=35;
  v->itols[0]=1;
  v->itols[1]=2;
  v->itols[2]=3;
  // mean values
  v->vtols[3]=v->vtols[4]=5;
  v->itols[3]=4;
  v->itols[4]=5;
  /* WM/GM probabilities
  v->vtols[5]=v->vtols[6]=25;
  v->itols[5]=6;
  v->itols[6]=7;*/
  if(0) v->tree=NULL;
  else {
    v->tree=(kdtree *)calloc(1,sizeof(kdtree));
    v->tree->ndim=v->ntols;
    v->tree->dim =(int    *)calloc(v->tree->ndim,sizeof(int));
    v->tree->step=(double *)calloc(v->tree->ndim,sizeof(double));
    v->tree->dim[0]=v->tree->dim[1]=51;
    v->tree->step[0]=v->tree->step[1]=5;
    v->tree->num=v->tree->dim[0]*v->tree->dim[1];
    v->tree->tree=NULL;
    v->tree->numlist=NULL;
  }
  return v;
}

niik_nls_lesion2 *niik_nls_lesion2_free(niik_nls_lesion2 *v) {
  int n;
  if(v==NULL) return NULL;
  if(v->patch!=NULL) v->patch=niik_image_free(v->patch);
  for(n=0; n<v->nlib; n++) {
    v->lib[n]=niik_image_free(v->lib[n]);
  }
  if(v->t1w!=NULL) v->t1w=niik_image_free(v->t1w);
  if(v->flr!=NULL) v->flr=niik_image_free(v->flr);
  if(v->mask!=NULL) v->mask=niik_image_free(v->mask);
  free(v);
  return NULL;
}

int niik_nls_lesion2_read_lib(niik_nls_lesion2 *v) {
  FILE *fp=NULL;
  int i=0,j=0,k=0;
  double dsum;
  char line[4096],*cptr;
  nifti_image *cimg=NULL;
  NIIK_RET0((v==NULL),__func__,"v is null");
  NIIK_RET0((v->lib_list_file==NULL),__func__,"lib_list_file is null");
  // get images
  v->nlib=0;
  NIIK_RET0(((fp=fopen(v->lib_list_file,"r"))==NULL),__func__,"fopen file");
  while( fgets(line,65536,fp) ) {
    v->nlib++;
  }
  fclose(fp);
  fprintf(stdout,"[%s] #lib = %i\n",__func__,v->nlib);
  NIIK_RET0(((v->lib=(nifti_image **)calloc(v->nlib,sizeof(nifti_image *)))==NULL),__func__,"calloc v->lib");
  NIIK_RET0(((fp=fopen(v->lib_list_file,"r"))==NULL),__func__,"fopen file");
  while( fgets(line,65536,fp) ) {
    if((cptr = strchr(line,'\n'))!=NULL) {
      *cptr = '\0';
    }
    if(v->verbose>2) fprintf(stdout,"[%s] reading %4i %s\n",__func__,i,line);
    NIIK_RET0(((v->lib[i++]=niik_image_read(line))==NULL),__func__,"reading patch lib");
    cimg=v->lib[i-1];
    NIIK_RET0((cimg->nx!=NX),__func__,"incorrect ny");
    if(v->verbose>2) fprintf(stdout,"[%s] %i %i %i\n",__func__,cimg->nx,cimg->ny,cimg->nz);
    if(cimg->nz>1) {
      cimg->ny=cimg->dim[2]=cimg->nvox/cimg->nx;
      cimg->nz=cimg->dim[3]=1;
      cimg->ndim=cimg->dim[0]=2;
      if(v->verbose>2) {
        fprintf(stdout,"[%s] fix dim %i,%i,%i\n",__func__,cimg->nx,cimg->ny,cimg->nz);
      }
      for(j=cimg->ny-1; j>0; j--) {
        for(k=0,dsum=0; k<cimg->nx; k++) {
          dsum+=niik_image_get_voxel(cimg,j*cimg->nx+k);
          if(dsum>0) break;
        }
        if(dsum>0.0) break;
      }
      cimg->ny=cimg->dim[2]=j+1;
      cimg->nvox=cimg->nx*cimg->ny;
      if(v->verbose>2) {
        fprintf(stdout,"[%s] fix dim %i,%i,%i *\n",__func__,cimg->nx,cimg->ny,cimg->nz);
      }
    } //fix dimensions
  }
  fclose(fp);
  return 1;
}


int niik_nls_lesion2_tree_idx(int *nt,kdtree *v) {
  int n,nn,k;
  for(n=0,nn=0,k=1; n<v->ndim; n++) {
    nn+=nt[n]*k;
    k*=v->dim[n];
  }
  return nn;
}

int niik_nls_lesion2_update_kdtree(niik_nls_lesion2 *v) {
  int i,nn,ni,nj,nx,nt[9];
  unsigned char *bimg;
  if(v->verbose>1) {
    niik_fc_display((char*)__func__,1);
    fprintf(stdout,"[%s] ndim %i\n",__func__,v->tree->ndim);
  }
  nx=v->lib[0]->nx;
  for(i=0; i<v->tree->ndim; i++) {
    v->tree->step[i]=v->vtols[i];
    v->tree->dim[i]=(255/v->tree->step[i]+1.5);
  }
  for(i=0,v->tree->num=1; i<v->tree->ndim; i++) {
    v->tree->num*=v->tree->dim[i];
  }
  if(v->verbose>1) {
    for(i=0; i<v->tree->ndim; i++) {
      fprintf(stdout,"\t%i %i %f\n",i,v->tree->dim[i],v->tree->step[i]);
    }
    fprintf(stdout,"total size %i\n",v->tree->num);
  }
  NIIK_RET0(((v->tree->tree=(unsigned char ***)calloc(v->tree->num,sizeof(char **)))==NULL),__func__,"calloc");
  NIIK_RET0(((v->tree->numlist=(int *)calloc(v->tree->num,sizeof(int)))==NULL),__func__,"calloc numlist");
  for(i=0; i<v->tree->num; i++) {
    v->tree->tree[i]=NULL;
    v->tree->numlist[i]=0;
  }
  for(nn=0; nn<v->nlib; nn++) {
    bimg=(unsigned char *)v->lib[nn]->data;
    for(nj=0; nj<v->lib[nn]->ny; nj++) {
      for(i=0; i<v->tree->ndim; i++) {
        nt[i]=floor(bimg[nj*nx+i+1]/v->tree->step[i]);
      }
      ni=niik_nls_lesion2_tree_idx(nt,v->tree);
      NIIK_RET0((ni<0),__func__,"ni < 0");
      NIIK_RET0((ni>=v->tree->num),__func__,"ni > num");
      v->tree->numlist[ni]++;
    }
  }
  for(i=0; i<v->tree->num; i++) {
    v->tree->tree[i]=(unsigned char **)calloc(v->tree->numlist[i],sizeof(char *));
    v->tree->numlist[i]=0;
  }
  for(nn=0; nn<v->nlib; nn++) {
    bimg=(unsigned char *)v->lib[nn]->data;
    for(nj=0; nj<v->lib[nn]->ny; nj++) {
      for(i=0; i<v->tree->ndim; i++) {
        nt[i]=floor(bimg[nj*nx+i+1]/v->tree->step[i]);
      }
      ni=niik_nls_lesion2_tree_idx(nt,v->tree);
      NIIK_RET0((ni<0),__func__,"ni < 0");
      NIIK_RET0((ni>=v->tree->num),__func__,"ni > num");
      v->tree->tree[ni][v->tree->numlist[ni]++]=bimg+nj*nx;
    }
  }
  if(v->verbose>1) niik_fc_display((char*)__func__,0);
  return 1;
}

int niik_nls_lesion2_patch_voxel(nifti_image **imgs,int nimg,nifti_image *mask,int i,int j,int k,int pi,int pj,int pk,unsigned char *bout) {
  int
  nn=0,
  ni,nj,nk,
  mi,mj,mk,mm,m,n;
  double sum[9];
  n=i+j*mask->nx+k*mask->nx*mask->ny;
  bout[nn++]=niik_image_get_voxel(imgs[0],n);
  bout[nn++]=(unsigned char)floor(NIIK_DMINMAX(i*mask->dx+0.5,0,255));
  bout[nn++]=(unsigned char)floor(NIIK_DMINMAX(j*mask->dy+0.5,0,255));
  bout[nn++]=(unsigned char)floor(NIIK_DMINMAX(k*mask->dz+0.5,0,255));
  nn+=2; // 4 and 5
  sum[1]=sum[2]=0;
  bout[nn++]=niik_image_get_voxel(imgs[4],n)*100;
  bout[nn++]=niik_image_get_voxel(imgs[5],n)*100;
  for(mk=-pk; mk<=pk; mk++) {
    nk=mk+k;
    if(nk<0) nk=0;
    else if(nk>=mask->nz) nk=mask->nz-1;
    for(mj=-pj; mj<=pj; mj++) {
      nj=mj+j;
      if(nj<0) nj=0;
      else if(nj>=mask->ny) nj=mask->ny-1;
      for(mi=-pi; mi<=pi; mi++) {
        ni=mi+i;
        if(ni<0) ni=0;
        else if(ni>=mask->nx) ni=mask->nx-1;
        m=ni+nj*mask->nx+nk*mask->nx*mask->ny;
        for(mm=1; mm<nimg; mm++) {
          bout[nn++]=niik_image_get_voxel(imgs[mm],m);
          sum[mm]+=(double)bout[nn-1];
        }
      }
    }
  }
  bout[4]=(unsigned char)floor(NIIK_DMINMAX(sum[1]/(pi*2+1)/(pj*2+1)/(pk*2+1)+0.5,0,255));
  bout[5]=(unsigned char)floor(NIIK_DMINMAX(sum[2]/(pi*2+1)/(pj*2+1)/(pk*2+1)+0.5,0,255));
  return 1;
}

nifti_image *niik_nls_lesion2_patch_case(nifti_image **imgs,int nimg,int xpatch,int ypatch,int zpatch)
/* assuming 1x1x3 input files
 * imgs[0] = lesion
 * imgs[1] = flr
 * imgs[2] = t1w
 * imgs[3] = mask
 * imgs[4] = pwm
 * imgs[5] = pgm
 * nimg must be 3
 * patch dim is 5x5x3
 * output image dim is nx=151, ny=sum(mask>0) or 32767, whichever is smaller
 * nx=151  -> [0]=lesion, [1-]=flr,t1w values 5x5x3
 * nx=154  -> [0]=lesion, [1-3]=x,y,z coordinates, [4-]=flr,t1w values 5x5x3
 * nx=156  -> [0]=lesion, [1-3]=x,y,z coordinates, [4-5]=flr,t1w mean values 5x5x3, [6-]=flr,t1w values 5x5x3
 * nx=158  -> [0]=lesion, [1-3]=x,y,z coordinates, [4-5]=flr,t1w mean values 5x5x3, [6,7]=pwm,pgm, [8-]=flr,t1w values 5x5x3
 * output's nx=0 is the lesion values
 */
{
  nifti_image *patch=NULL,*mask=NULL;
  int
  nx=NX,
  nn,nm,ny,nz,
  i,j,k,n;
  unsigned char *bimg=NULL,*bout=NULL;
  NIIK_RET0((imgs==NULL),__func__,"imgs is null");
  mask=imgs[3];
  NIIK_RET0((mask->datatype!=NIFTI_TYPE_UINT8),__func__,"mask is not uint8");
  nm=niik_image_count_mask(mask);
  fprintf(stdout,"[%s] mask size %i\n",__func__,nm);
  if(nm>32767) {
    nz=nm/32767+1;
    ny=32767;
    fprintf(stdout,"[%s] zdim resized %i\n",__func__,nz);
  } else {
    ny=nm;
    nz=1;
  }
  NIIK_RET0(((patch=niik_image_init(nx,ny,nz,1,1,1,1,1,1,1,1,1,1,1,NIFTI_TYPE_UINT8))==NULL),__func__,"niik_image_init");
  bimg=(unsigned char *)mask->data;
  bout=(unsigned char *)patch->data;
  for(k=n=nn=0; k<mask->nz; k++) {
    for(j=0; j<mask->ny; j++) {
      for(i=0; i<mask->nx; i++,n++) {
        if(bimg[n]==0) continue;
        if(nn>=patch->nvox) {
          i=j=k=32767;
          break;
        }
        niik_nls_lesion2_patch_voxel(imgs,nimg,mask,i,j,k,xpatch,ypatch,zpatch,bout+nn);
        nn+=nx;
      }
    }
  }
  return patch;
}

typedef struct {
  int imgid;  //img id
  int yid;    // y idx
  int treeid;
  double e;   // error
  float v;    // output value
} NIIKNLS_t;

double niiknls_tv_median(NIIKNLS_t *t,int num) {
  double *d=NULL,v;
  int i;
  d=(double *)calloc(num,sizeof(double));
  for(i=0; i<num; i++) d[i]=t[i].v;
  v=niik_median_quicksort_double(d,num);
  free(d);
  return v;
}

int niik_nls_lesion2_proc_voxel(niik_nls_lesion2 *v,nifti_image **imgs,int nimg,int n,int skip_tol) {
  int ni,nj,nn,nx;
  int i,j,k;
  int nt[9],nc[256],nnc;
  double e=0;
  unsigned char *bimg,*bout;
  unsigned char *mask;
  NIIKNLS_t *t,tt;
  int topd1;
  int tcount[9];

  mask=(unsigned char *)v->mask->data;
  if(mask[n]==0) return 1;
  //fprintf(stdout,"skiptol = %i\n",skip_tol);
  nx=v->lib[0]->nx;
  bout=(unsigned char *)calloc(nx,sizeof(double));
  niik_nls_lesion2_patch_voxel(imgs,nimg,v->mask,n%v->mask->nx,(n/v->mask->nx)%v->mask->ny,n/v->mask->nx/v->mask->ny,v->xpatch,v->ypatch,v->zpatch,bout);
  NIIK_RET0(((t=(NIIKNLS_t *)calloc(v->topd,sizeof(NIIKNLS_t)))==NULL),__func__,"calloc");
  for(i=0; i<v->topd; i++) {
    t[i].e=1e31;
    t[i].imgid=-1;
  }
  topd1=v->topd-1;

  for(i=0; i<9; i++)tcount[i]=0;
  if(v->tree==NULL || skip_tol) {
    for(nn=0; nn<v->nlib; nn++) {
      bimg=(unsigned char *)v->lib[nn]->data;
      for(nj=0; nj<v->lib[nn]->ny; nj++) {
        tcount[0]++;
        if(v->verbose>2) fprintf(stdout,"[%s:%i:%s] nn=%i nj=%i\n",__FILE__,__LINE__,__func__,nn,nj);
        if(!skip_tol) {
          for(ni=0; ni<v->ntols; ni++) {
            if(v->verbose>2)fprintf(stdout,"ni=%i, itol=%i, vtol=%.0f   %i, %i\n",ni,v->itols[ni],v->vtols[ni],bimg[v->itols[ni]+nj*nx],bout[v->itols[ni]]);
            if(fabs((double)bimg[v->itols[ni]+nj*nx]-(double)bout[v->itols[ni]])>v->vtols[ni]) {
              ni=-32;
              break;
            }
          }
          if(ni==-32) continue;
        }
        tcount[1]++;
        for(ni=4,e=0; ni<nx; ni++) {
          e+=fabs((double)bimg[ni+nj*nx] - (double)bout[ni]);
        }
        if(v->verbose>2) fprintf(stdout,"nn=%i nj=%i, e=%.0f %.0f\n",nn,nj,e/(nx-1.0),e);
        if(t[topd1].e>e) {
          tcount[2]++;
          if(v->verbose>2) fprintf(stdout,"nn=%i nj=%i, e=%.0f %.0f add\n",nn,nj,e/(nx-1.0),e);
          // update the top list
          t[topd1].e=e;
          t[topd1].imgid=nn;
          t[topd1].yid=nj;
          t[topd1].v=(float)bimg[nj*nx];
          // sort the top list
          for(i=0; i<v->topd; i++) {
            for(j=i+1; j<v->topd; j++) {
              if(t[i].e>t[j].e) {
                tt=t[i];
                t[i]=t[j];
                t[j]=tt; // swap
              }
            }
          }
        } // add to the top list
      } // each feature
    } // each lib
  } else {

    for(i=0; i<v->tree->ndim; i++) {
      nt[i]=floor(bout[i+1]/v->tree->step[i]-0.5);
      nt[i]=NIIK_IMINMAX(nt[i],0,v->tree->dim[i]-1);
    }
    nnc=0;
    nc[nnc++]=niik_nls_lesion2_tree_idx(nt,v->tree);
    for(i=0,ni=1; i<v->tree->ndim; i++) {
      k=nnc;
      for(j=0; j<k; j++)
        nc[nnc++]=nc[j]+ni;
      ni*=v->tree->dim[i];
    }
    /*for(i=0;i<v->tree->ndim;i++) {
      fprintf(stdout,"%3i %i\n",i,v->tree->dim[i]);}
      for(i=0;i<nnc;i++){
      fprintf(stdout,"%5i %12i %5i %5i %5i %5i %5i\n",i,nc[i],nc[i]%v->tree->dim[0],(nc[i]/v->tree->dim[0])%v->tree->dim[1],(nc[i]/v->tree->dim[0]/v->tree->dim[1])%v->tree->dim[2],
      (nc[i]/v->tree->dim[0]/v->tree->dim[1]/v->tree->dim[2])%v->tree->dim[3],
      (nc[i]/v->tree->dim[0]/v->tree->dim[1]/v->tree->dim[2]/v->tree->dim[3])%v->tree->dim[4]);} */

    for(ni=0; ni<nnc; ni++) {
      nn=nc[ni];
      for(k=0; k<v->tree->numlist[nn]; k++) {
        tcount[0]++;
        bimg=v->tree->tree[nn][k];
        //fprintf(stdout,"%5i [%i %i] [%i %i]\n",k,bout[4],bout[5],bimg[4],bimg[5]);
        for(i=0; i<v->ntols; i++) {
          if(fabs((double)bimg[v->itols[i]]-(double)bout[v->itols[i]])>v->vtols[i]) {
            i=-32;
            break;
          }
        }
        if(i==-32) continue;
        tcount[1]++;
        for(i=4,e=0; i<nx; i++) {
          e+=fabs((double)bimg[i] - (double)bout[i]);
        }
        if(t[topd1].e>e) {
          tcount[2]++;
          // update the top list
          t[topd1].e=e;
          t[topd1].imgid=nn;
          t[topd1].yid=k;
          t[topd1].v=(float)bimg[0];
          // sort the top list
          for(i=0; i<v->topd; i++) {
            for(j=i+1; j<v->topd; j++) {
              if(t[i].e>t[j].e) {
                tt=t[i];
                t[i]=t[j];
                t[j]=tt; // swap
              }
            }
          }
        } // add to the top list

      } // patch list
    } // tree
  } // using tree

  if(t[topd1].imgid<0) {  // did not find optimal patch
    free(t);
    free(bout);
    if(skip_tol) {
      fprintf(stderr,"[%s] ERROR: could not find optimal patch in second try\n",__func__);
      exit(0);
      return 0;
    }
    if(v->verbose>2)fprintf(stdout,"[%s] could not find optimal patch, try again\n",__func__);
    if(!niik_nls_lesion2_proc_voxel(v,imgs,nimg,n,1)) {
      fprintf(stderr,"[%s] ERROR: could not find optimal patch\n",__func__);
      return 0;
    }
    if(v->verbose>2)fprintf(stdout,"[%s]   found patch out = %i\n",__func__,mask[n]);
    return 1;
  }

  switch(v->fusion_method) {
  case NIIK_NLS_FUSION_METHOD_MIN:
    mask[n]=(unsigned char)floor(NIIK_DMINMAX(t[0].v+0.5,0,255));
    break;
  case NIIK_NLS_FUSION_METHOD_AVG:
    for(i=0,e=0; i<v->topd; i++) {
      e+=t[i].v;
    }
    mask[n]=(unsigned char)floor(NIIK_DMINMAX(e/v->topd+0.5,0,255));
    break;
  case NIIK_NLS_FUSION_METHOD_MEDIAN:
    e=niiknls_tv_median(t,v->topd);
    mask[n]=(unsigned char)floor(NIIK_DMINMAX(e+0.5,0,255));
    break;
  }

  if(v->verbose_idx[0]>=0) {
    /*fprintf(stdout,"test counts %i %i %i %i %i\n",
      testcount[0],testcount[1],testcount[2],testcount[3],testcount[4]);*/
    for(i=0; i<v->topd; i++) {
      if(v->tree==NULL || skip_tol) {
        bimg=v->lib[t[i].imgid]->data;
        bimg=bimg+t[i].yid*nx;
      } else {
        bimg=v->tree->tree[t[i].imgid][t[i].yid];
      }
      fprintf(stdout,"%3i %3i, %3i %3i | e=%-7.3f n=%-4i y=%-6i  o=%.0f ^\n",
              bout[4],bimg[4],bout[5],bimg[5],
              t[i].e/(nx-4.0),t[i].imgid,t[i].yid,t[i].v);
    }
  }

  if(v->verbose>1) {
    fprintf(stdout,"count = %i %i %i\n",tcount[0],tcount[1],tcount[2]);
    if(v->tree==NULL || skip_tol) {
      bimg=v->lib[t[0].imgid]->data;
      bimg=bimg+t[0].yid*nx;
    } else {
      bimg=v->tree->tree[t[0].imgid][t[0].yid];
    }
    fprintf(stdout,"[%3i %3i %3i, %9i] | %3i %3i, %3i %3i | e=%-7.3f n=%-3i y=%-6i  o=%i ^\n",
            n%v->mask->nx,(n/v->mask->nx)%v->mask->ny,n/v->mask->nx/v->mask->ny,n,
            bout[4],bimg[4],bout[5],bimg[5],
            t[0].e/(nx-4.0),t[0].imgid,t[0].yid,mask[n]);
  }

  free(t);
  free(bout);
  //  niik_fc_display((char *)__func__,0);
  return 1;
}

int niik_nls_lesion2_proc(niik_nls_lesion2 *v) {
  nifti_image *imgs[9];
  int const nimg=3;
  int n,nn;
  unsigned char *mask;

  NIIK_RET0((v==NULL),__func__,"v is null");
  NIIK_RET0((v->nlib<=0),__func__,"v->nlib is bad (<=0)");
  NIIK_RET0((v->t1w==NULL),__func__,"missing t1w");
  NIIK_RET0((v->flr==NULL),__func__,"missing flr");
  NIIK_RET0((v->mask==NULL),__func__,"missing mask");
  niik_fc_display((char *)__func__,1);

  for(n=nn=0; n<v->nlib; n++) {
    nn+=v->lib[n]->ny;
  }

  fprintf(stdout,"[%s] #lib        = %i\n",__func__,v->nlib);
  fprintf(stdout,"[%s] #patch cmp  = %i\n",__func__,nn);
  fprintf(stdout,"[%s] patch size  = %i\n",__func__,v->lib[0]->nx-1);
  fprintf(stdout,"[%s] top         = %i\n",__func__,v->topd);
  for(n=0; n<v->ntols; n++) {
    fprintf(stdout,"[%s] tolerance[%i] = %3i  [%i]\n",__func__,n,(int)v->vtols[n],v->itols[n]);
  }

  fprintf(stdout,"[%s] convert to uint8\n",__func__);
  niik_image_type_convert(v->mask,NIFTI_TYPE_UINT8);
  niik_image_type_convert(v->t1w,NIFTI_TYPE_UINT8);
  niik_image_type_convert(v->flr,NIFTI_TYPE_UINT8);
  for(n=0; n<v->nlib; n++) {
    niik_image_type_convert(v->lib[n],NIFTI_TYPE_UINT8);
  }

  imgs[0]=v->t1w;
  imgs[1]=v->flr;
  imgs[2]=v->t1w;
  imgs[3]=v->mask;
  imgs[4]=v->pwm;
  imgs[5]=v->pgm;
  NIIK_RET0((!niik_image_type_convert(v->mask,NIFTI_TYPE_UINT8)),__func__,"convert mask");

  if(v->tree!=NULL)
    NIIK_RET0((!niik_nls_lesion2_update_kdtree(v)),__func__,"niik_nls_lesion2_update_kdtree");

  if(v->verbose_idx[0]>=0) {
    niik_nls_lesion2_proc_voxel(v,imgs,nimg,
                                v->verbose_idx[1]+
                                v->verbose_idx[2]*v->mask->nx+
                                v->verbose_idx[3]*v->mask->nx*v->mask->ny,
                                0);
    exit(0);
  }

  mask=v->mask->data;

  fprintf(stderr,"[%s] loop start   ",__func__);
  #pragma omp parallel for ordered schedule(dynamic,512)
  for(n=0; n<v->mask->nvox; n++) {
    if(mask[n]==0) continue;
    if(v->verbose==1)
      fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b[%3i,%3i,%3i]",n%v->mask->nx,(n/v->mask->nx)%v->mask->ny,(n/v->mask->nx)/v->mask->ny);
    niik_nls_lesion2_proc_voxel(v,imgs,nimg,n,0);
  }
  fprintf(stderr,"\n");

  niik_fc_display((char *)__func__,0);
  return 1;
}

static char *prog_nls_lesion_version[] = {
  "  niik_nls_lesion history\n"
  "\n  0.0.0  December 09, 2014, Kunio Nakamura, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
  "\n  0.1.0  January 09, 2014, Kunio Nakamura, knakamura@mrs.bic.mcgill.ca\n"
  "  -second version\n"
  "  -complete remake\n"
};

static char *prog_nls_lesion_help[] = {
  "  niik_nls_lesion:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niik_nls_lesion\n");
  fprintf(stdout,"  usage: [options] <t1w.nii> <flr.nii> <mask.ni> <out.nii>\n\n");
  fprintf(stdout,"\n");
}


int main(int argc,char *argv[],char *envp[]) {
  char
  fcname[32]="niik_nls_lesion_main";
  char fname[1024];
  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  int
  i,*idx,
  nc=1,sc=1;
  niik_nls_lesion2 *nls=NULL;
  nifti_image *imgs[9];

  if(argc==1) {
    usage();
    exit(1);
  }

  niik_fc_display(fcname,1);
  srand(time(NULL));
  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  fprintf(stdout,"[%s] version %i.%i.%i   exec @ %s\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION,tmstr);

  NIIK_EXIT(((nls=niik_nls_lesion2_init())==NULL),__func__,"niik_nls_lesion2_init",1);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*prog_nls_lesion_version);
        exit(1);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*prog_nls_lesion_help);
        exit(1);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*prog_nls_lesion_help);
        exit(1);
      }

      else if(!strncmp(argv[nc],"-make-patch",11)) {
        NIIK_EXIT((argc<=nc+1),fcname,"missing value",1);
        nc++;
        sprintf(fname,"%s_ct2f_reg2t1p.mnc",argv[nc]);
        imgs[0]=niik_image_read(fname);
        sprintf(fname,"%s_flr_reg2t1p.mnc",argv[nc]);
        imgs[1]=niik_image_read(fname);
        sprintf(fname,"%s_t1p.mnc",argv[nc]);
        imgs[2]=niik_image_read(fname);
        sprintf(fname,"%s_t1p_brain_mask.mnc",argv[nc]);
        imgs[3]=niik_image_read(fname);
        sprintf(fname,"%s_pwm_reg2t1p.mnc",argv[nc]);
        imgs[4]=niik_image_read(fname);
        sprintf(fname,"%s_pgm_reg2t1p.mnc",argv[nc]);
        imgs[5]=niik_image_read(fname);
        for(i=0; i<6; i++) NIIK_EXIT((!niik_image_type_convert(imgs[i],NIFTI_TYPE_UINT8)),fcname,"type convert",1);
        imgs[0]=niik_nls_lesion2_patch_case(imgs,3,nls->xpatch,nls->ypatch,nls->zpatch);
        sprintf(fname,"%s_patch_flr-t1w.nii.gz",argv[nc]);
        fprintf(stdout,"[%s] writing output %s\n",__func__,fname);
        niik_image_write(fname,imgs[0]);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-tol-avg-t1w=",13)) {
        nls->vtols[4]=atof(argv[nc]+13);
        fprintf(stdout,"[%s] avg-t1w tolerance = %i\n",fcname,(int)nls->vtols[4]);
      } else if(!strncmp(argv[nc],"-tol-avg-flr=",13)) {
        nls->vtols[3]=atof(argv[nc]+13);
        fprintf(stdout,"[%s] avg-flr tolerance = %i\n",fcname,(int)nls->vtols[3]);
      }

      else if(!strncmp(argv[nc],"-lib=",5)) {
        nls->lib_list_file=argv[nc]+5;
        fprintf(stdout,"[%s] library list file = %s\n",fcname,nls->lib_list_file);
      }

      else if(!strncmp(argv[nc],"-notree",7)) {
        nls->tree=kdtree_free(nls->tree);
        fprintf(stdout,"[%s] turned off kd-tree\n",fcname);
      } else if(!strncmp(argv[nc],"-tree",5)) {
        fprintf(stdout,"[%s] turned on kd-tree\n",fcname);
      }

      else if(!strncmp(argv[nc],"-tol=",5)) {
        NIIK_EXIT(((idx=niik_csstring_to_int_list(argv[nc]+5,&i))==NULL),fcname,"niik_csstring_to_int_list",1);
        NIIK_EXIT((i!=2),fcname,"2 comma seprated values required",1);
        NIIK_EXIT((idx[0]>nls->ntols),fcname,"idx[0] too large",1);
        nls->vtols[idx[0]]=idx[1];
        fprintf(stdout,"[%s] tolerance = %i [%i]\n",fcname,(int)nls->vtols[idx[0]],idx[0]);
        free(idx);
      }


      else if(!strncmp(argv[nc],"-top=",5)) {
        nls->topd=atoi(argv[nc]+5);
        fprintf(stdout,"[%s] #top list = %i\n",fcname,nls->topd);
      }

      else if(!strncmp(argv[nc],"-idx=",5)) {
        nls->verbose_idx[0]=0;
        NIIK_EXIT(((idx=niik_csstring_to_int_list(argv[nc]+5,&i))==NULL),fcname,"niik_csstring_to_int_list",1);
        NIIK_EXIT((i!=3),fcname,"3 comma seprated values required",1);
        nls->verbose_idx[1]=idx[0];
        nls->verbose_idx[2]=idx[1];
        nls->verbose_idx[3]=idx[2];
        free(idx);
        fprintf(stdout,"[%s] idx = %i %i %i\n",fcname,nls->verbose_idx[1],nls->verbose_idx[2],nls->verbose_idx[3]);
      }


      else if(!strncmp(argv[nc],"-v=",3)) {
        nls->verbose=atoi(argv[nc]+3);
      } else if(!strncmp(argv[nc],"-u",2)) {
        usage();
        exit(1);

      } else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(1);
      }

      nc++;
    }

    else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  NIIK_EXIT((argc<7),fcname,"too few argments",1);
  NIIK_EXIT((argc>7),fcname,"too many argments",1);

  fprintf(stdout,"[%s] reading t1w   %s\n",fcname,argv[1]);
  NIIK_EXIT(((nls->t1w=niik_image_read(argv[1]))==NULL),__func__,"reading t1w",1);
  fprintf(stdout,"[%s] reading flr   %s\n",fcname,argv[2]);
  NIIK_EXIT(((nls->flr=niik_image_read(argv[2]))==NULL),__func__,"reading flr",1);
  fprintf(stdout,"[%s] reading mask  %s\n",fcname,argv[3]);
  NIIK_EXIT(((nls->mask=niik_image_read(argv[3]))==NULL),__func__,"reading mask",1);
  fprintf(stdout,"[%s] reading pwm   %s\n",fcname,argv[4]);
  NIIK_EXIT(((nls->pwm=niik_image_read(argv[4]))==NULL),__func__,"reading mask",1);
  fprintf(stdout,"[%s] reading pgm   %s\n",fcname,argv[5]);
  NIIK_EXIT(((nls->pgm=niik_image_read(argv[5]))==NULL),__func__,"reading mask",1);

  NIIK_EXIT((!niik_nls_lesion2_read_lib(nls)),__func__,"niik_nls_lesion2_read_lib",1);
  NIIK_EXIT((!niik_nls_lesion2_proc(nls)),__func__,"niik_nls_lesion2_proc",1);

  fprintf(stdout,"[%s] writing out   %s\n",fcname,argv[6]);
  NIIK_EXIT((!niik_image_write(argv[6],nls->mask)),__func__,"writing output",1);

  nls=niik_nls_lesion2_free(nls);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_nls_lesion */
