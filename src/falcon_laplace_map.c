/* Filename:     nifti1_kunio_laplace_map.c
 * Description:  function to calculate laplace map
 * Author:       Kunio Nakamura
 * Date:         March 4, 2012
 */

#ifndef _FALCON_LAPLACE_MAP_C_
#define _FALCON_LAPLACE_MAP_C_

#include "falcon.h"
#include "falcon_morph.h"

/********************************************************************************
 * calculates the laplacian map from Simg and Limg
 *
 * -Simg is the smaller obj
 * -Limg is the larger obj
 ********************************************************************************/

nifti_image *niik_image_laplace_map(nifti_image *Simg, nifti_image *Limg) {
  const char *fcname="niik_image_laplace_map";
  nifti_image *outimg=NULL;
  int
  xdim,dim12,
       i,j,k,m,n,p;
  unsigned char *simg,*limg;
  unsigned long *oimg;
  float *fimg;
  int verbose=niik_verbose();
  double LAPMAX=1000,LAPMIN=0;
  double *sa,*x,*b,**A=NULL;
  int *ija,N,nsa;

  if(verbose>0) niik_fc_display(fcname,1);
  if(Simg==NULL) {
    fprintf(stderr,"ERROR: Simg is a null pointer\n");
    return NULL;
  }
  if(Limg==NULL) {
    fprintf(stderr,"ERROR: Limg is a null pointer\n");
    return NULL;
  }

  if(verbose>2) fprintf(stdout,"[niik_image_laplace_map] niik_image_cmp_dim\n");
  if(niik_image_cmp_dim(Simg,Limg)!=0) {
    fprintf(stderr,"ERROR: niik_image_cmp_dim\n");
    return NULL;
  }

  if(verbose>2) fprintf(stdout,"[niik_image_laplace_map] check datatype \n");
  if(Simg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"ERROR: Simg is not uint8 \n");
    return NULL;
  } else if(Limg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"ERROR: Simg is not uint8 \n");
    return NULL;
  }
  simg=Simg->data;
  limg=Limg->data;
  xdim=Limg->dim[1];
  dim12=xdim*Limg->dim[2];

  if(verbose>2) fprintf(stdout,"[niik_image_laplace_map] niik_image_copy\n");
  if((outimg=niik_image_copy(Simg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return (0);
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT64)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(%s,NIFTI_TYPE_UINT8)\n",outimg->fname);
    return NULL;
  }
  oimg=outimg->data;
  for(i=0; i<Limg->nvox; i++) oimg[i]=0;


  /* set up sparse matrix (system of linear equations)
   *   at the outer border of Limg value = LAPMAX
   *   inside the Simg value = LAPMIN
   * we will estimate inbetween these masks
   *
   * for Laplacian values
   * -OUTER EDGE OF LARGE -> constant (only diagonal)
   * -INNER EDGE OF SMALL -> constant (only diagonal)
   * -IN LARGE & NOT IN SMALL -> variable (diagonal + 6 more off-diagonal terms)
   */

  if(verbose>2) fprintf(stdout,"[niik_image_laplace_map] determine the # voxels \n");
  N=nsa=0;
  for(k=n=0,m=1; k<Limg->dim[3]; k++) {
    for(j=0; j<Limg->dim[2]; j++) {
      for(i=0; i<Limg->dim[1]; n++,i++) {
        if(limg[n]==0) { /* check if it's an edge */
          if(i>0) if(limg[n-    1]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(j>0) if(limg[n- xdim]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(k>0) if(limg[n-dim12]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }

          if(i<Limg->dim[1]-1) if(limg[n+    1]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(j<Limg->dim[2]-1) if(limg[n+ xdim]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(k<Limg->dim[3]-1) if(limg[n+dim12]>0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
        }
        if(simg[n]>0) { /* check if it's an edge */
          if(i>0) if(simg[n-    1]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(j>0) if(simg[n- xdim]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(k>0) if(simg[n-dim12]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }

          if(i<Simg->dim[1]-1) if(simg[n+    1]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(j<Simg->dim[2]-1) if(simg[n+ xdim]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
          if(k<Simg->dim[3]-1) if(simg[n+dim12]==0) {
              N++;  /* constant */
              nsa++;
              oimg[n]=m++;
              continue;
            }
        }
        if(limg[n]>0 && simg[n]==0) {
          if(i==0) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }
          if(j==0) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }
          if(k==0) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }

          if(i==Limg->dim[1]-1) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }
          if(j==Limg->dim[2]-1) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }
          if(k==Limg->dim[3]-1) {
            N++;  /* constant */
            nsa++;
            oimg[n]=m++;
            continue;
          }
          N++;
          nsa+=7;
          oimg[n]=m++;
          continue;
        } /* not a constant */
      }
    }
  }

  if(verbose>2) {
    fprintf(stdout,"[niik_image_laplace_map] Matrix A is %i by %i\n",N,N);
    fprintf(stdout,"                         size of ija and sa is %i\n",nsa);
  }

  if(N==0) {
    fprintf(stderr,"ERROR: matrix A is %i by %i\n",N,N);
    return NULL;
  }

  if(nsa==0) {
    fprintf(stderr,"ERROR: size of ija or sa is %i\n",nsa);
    return NULL;
  }

  NIIK_RET0(((ija = (int    *)calloc(2*nsa,sizeof(int   )))==NULL),fcname,"calloc for ija");
  NIIK_RET0(((sa  = (double *)calloc(2*nsa,sizeof(double)))==NULL),fcname,"calloc for sa");
  NIIK_RET0(((x   = (double *)calloc(  2*N,sizeof(double)))==NULL),fcname,"calloc for x");
  NIIK_RET0(((b   = (double *)calloc(  2*N,sizeof(double)))==NULL),fcname,"calloc for b");

  /*A=(double **)calloc(N,sizeof(double *));
    for(n=0;n<N;n++){ A[n]=(double *)calloc(N,sizeof(double)); }*/

  if(verbose>2) fprintf(stdout,"[niik_image_laplace_map] prepare for conj grad\n");

  for(k=n=0,m=0,p=N+1; k<Limg->nz; k++) {
    for(j=0; j<Limg->ny; j++) {
      for(i=0; i<Limg->nx; n++,i++) {
        ija[m]=p;
        if(limg[n]>0 && simg[n]==0) {
          if(i==0) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          if(j==0) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          if(k==0) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          if(i==Limg->dim[1]-1) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          if(j==Limg->dim[2]-1) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          if(k==Limg->dim[3]-1) {
            x[m]=b[m]=LAPMAX;  /* constant */
            sa[m++]=1;
            continue;
          }
          b[m]=x[m]=0;
          sa[m]=-6;
          m++;
          sa[p]=1;
          ija[p]=oimg[n-1]-1;
          p++;
          sa[p]=1;
          ija[p]=oimg[n+1]-1;
          p++;
          sa[p]=1;
          ija[p]=oimg[n-xdim]-1;
          p++;
          sa[p]=1;
          ija[p]=oimg[n+xdim]-1;
          p++;
          sa[p]=1;
          ija[p]=oimg[n-dim12]-1;
          p++;
          sa[p]=1;
          ija[p]=oimg[n+dim12]-1;
          p++;
          continue;
        } /* not a constant */
        if(limg[n]==0) { /* check if it's an edge; on large image  */
          if(i>0) if(limg[n-    1]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(j>0) if(limg[n- xdim]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(k>0) if(limg[n-dim12]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(i<Limg->dim[1]-1) if(limg[n+    1]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(j<Limg->dim[2]-1) if(limg[n+ xdim]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(k<Limg->dim[3]-1) if(limg[n+dim12]>0) {
              x[m]=b[m]=LAPMAX;  /* constant */
              sa[m++]=1;
              continue;
            }
        }
        if(simg[n]>0) { /* check if it's an edge; on smaller image */
          if(i>0) if(simg[n-    1]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(j>0) if(simg[n- xdim]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(k>0) if(simg[n-dim12]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(i<Simg->dim[1]-1) if(simg[n+    1]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(j<Simg->dim[2]-1) if(simg[n+ xdim]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
          if(k<Simg->dim[3]-1) if(simg[n+dim12]==0) {
              x[m]=b[m]=LAPMIN;  /* constant */
              sa[m++]=1;
              continue;
            }
        }
      }
    }
  }

  if(A!=NULL) {
    niik_spmat_reverse(A,N,sa,ija);
    fprintf(stdout,"A   = ");
    niik_spmat_display_dig(A,N,N,0);
    fprintf(stdout,"ija = ");
    niik_display_int_vector(ija,nsa);
    fprintf(stdout,"sa  = ");
    niik_display_double_vector(sa,nsa);
    fprintf(stdout,"x   = ");
    niik_display_double_vector(x,N);
    fprintf(stdout,"b   = ");
    niik_display_double_vector(b,N);
  }

  niik_spmat_conj_grad(sa,ija,N,x,b);

  if(verbose>3) {
    fprintf(stdout,"\tx   = ");
    niik_display_double_vector(x,N);
  }

  if(verbose>2) fprintf(stdout,"\tfree vectors\n");
  free(sa);
  free(ija);
  free(b);

  NIIK_RET0((!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
  if(verbose>2) fprintf(stdout,"[%s] update image\n",fcname);

  fimg=outimg->data;
  for(k=n=0; k<Limg->dim[3]; k++) {
    for(j=0; j<Limg->dim[2]; j++) {
      for(i=0; i<Limg->dim[1]; n++,i++) {
        if((int)fimg[n]==0) continue;
        fimg[n]=x[(int)(fimg[n]-1)];
      }
    }
  }

  free(x);
  if(verbose>0) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_laplace_map */



int niik_image_correct_topology_for_cuberille(nifti_image *maskimg)
/* toplogy correction by connecting tissues
 * -useful in cuberille
 */
{
  nifti_image
  *dilimg=NULL,
   *lapmap=NULL;
  char fcname[64]="niik_image_correct_topology_for_cuberille";
  int i,j,k,n,xdim,area;
  double radius=15;
  float *fimg=NULL;
  unsigned char *bimg=NULL;

  niik_fc_display(fcname,1);

  if(!niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }

  fprintf(stdout,"[%s] dilation %.4f\n",fcname,radius);
  if((dilimg=niik_image_copy(maskimg))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
    return 0;
  }
  if(!niik_image_morph_3d_radius_mask(dilimg,NULL,NIIK_MORPH_DILATE,radius)) {
    fprintf(stderr,"[%s] ERROR: niik_image_morph_3d_radius_mask\n",fcname);
    return 0;
  }

  fprintf(stdout,"[%s] laplacian map\n",fcname);
  if((lapmap = niik_image_laplace_map(maskimg,dilimg))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_laplace_map\n",fcname);
    return 0;
  }
  fprintf(stdout,"[%s] laplacian map done\n",fcname);

  if(!niik_image_type_convert(lapmap,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }
  niik_image_write("tmp1.nii.gz",lapmap);
  exit(0);


  fimg = lapmap->data;
  bimg = maskimg->data;
  xdim = area = maskimg->nx;
  area *= maskimg->ny;

  for(k=1; k<maskimg->nz-1; k++) {
    for(j=1; j<maskimg->ny-1; j++) {
      for(i=1; i<maskimg->nx-1; i++) {
        n=i+j*xdim+k*area;
        if(fimg[n]>30) continue;
        if(fimg[n-1]>fimg[n]) {
          if(fimg[n+1]>fimg[n]) {
            if(bimg[n]==0) {
              fprintf(stdout,"[%s] adding [%3i,%3i,%3i] x\n",fcname,i,j,k);
              bimg[n]=1;
            }
          }
        }
        if(fimg[n-xdim]>fimg[n]) {
          if(fimg[n+xdim]>fimg[n]) {
            if(bimg[n]==0) {
              fprintf(stdout,"[%s] adding [%3i,%3i,%3i] y\n",fcname,i,j,k);
              bimg[n]=1;
            }
          }
        }
        if(fimg[n-area]>fimg[n]) {
          if(fimg[n+area]>fimg[n]) {
            if(bimg[n]==0) {
              fprintf(stdout,"[%s] adding [%3i,%3i,%3i] z\n",fcname,i,j,k);
              bimg[n]=1;
            }
          }
        }
      }
    }
  }

  dilimg=niik_image_free(dilimg);
  lapmap=niik_image_free(lapmap);
  niik_fc_display(fcname,0);
  return 1;
}


#endif /* _FALCON_LAPLACE_MAP_C_ */

/*
 kate: space-indent on; indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
 */