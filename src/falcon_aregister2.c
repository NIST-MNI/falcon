/* Filename:     nifti1_kunio_aregister2.c
 * Description:  affine registration function 2
 *               To replace nifti1_kunio_aregister.c
 * Author:       Kunio Nakamura
 * Date:         April 26, 2012
 *
 *
 * -new affine registration
 *
 * int niik_image_aregister2_multilevel(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,int nlevel,niikvec *FWHM, niikvec *delta, niikvec *tol,int *nseed,int *maxiter,niikmat *daffpar,niikmat *seed)
 *
 * Revision:    September 26, 2012, Kunio Nakamura
 *              -add function names for all error messages
 *              -niik_image_aregister2_test1 works with areg_cost = NIIK_REGISTER_NMI; I had not checked it yet...
 *
 */

#ifndef _FALCON_AREGISTER2_C_
#define _FALCON_AREGISTER2_C_

#include "falcon.h"

nifti_image *g_niik_aregister2_refimg;
nifti_image *g_niik_aregister2_refseg;
nifti_image *g_niik_aregister2_movimg;
niikmat     *g_niik_aregister2_imat=NULL;
double       g_niik_aregister2_param[20];
int          g_niik_aregister2_fixpar[20];
int          g_niik_aregister2_iter=0;
int          g_niik_aregister2_invmat=0;
int          g_niik_aregister2_method=0;
nmi_obj     *g_niik_aregister2_nmiobj=NULL;
int          g_niik_aregister2_verbose=0;



void niik_aregister2_set_verbose( int verbose ) {
  g_niik_aregister2_verbose=verbose;
}
int  niik_aregister2_get_verbose() {
  return g_niik_aregister2_verbose;
}

char *niik_aregister2_method_string( int reg_method ) {
  switch(reg_method) {
  case NIIK_REGISTER_SAD:
    return "Sum of absolute difference";
  case NIIK_REGISTER_SSQ:
    return "Sum of squared difference";
  case NIIK_REGISTER_CC:
    return "Correlation coefficient";
  case NIIK_REGISTER_NMI:
    return "Normalized mutual information";
  default:
    return "[niik_aregister2_method_string] ERROR: unknown registration method";
  }
  return "[niik_aregister2_method_string] ERROR: unknown registration method";
}


double niik_image_aregister2_obj_func(double *v)
/* cost function for registration */
{
  const char *fcname="niik_image_aregister2_obj_func";
  static double val_optim=1e5;
  niikmat *afmat = NULL;
  niikpt
  movfov,
  p,q;
  int
  verbose=0,
  m,n,num,
  ni[2],nj[2],
  xdim,ydim,area;

  static float *tmpdata1=NULL;
  static unsigned char *tmpdata2=NULL;
  float
  *fimg;
  unsigned char
  *bimg;
  double
  dint1,dint2,fi[2],dd[2],
        affpar[25],
        errval=0,
        wval,
        cc,nmi,
        xsum,ysum,xysum,xssq,yssq,wsum;
  verbose=g_niik_aregister2_verbose-1;
  /* initialize when inputs are null */
  if(v==NULL) {
    if(verbose>1) fprintf(stdout,"[niik_image_aregister_obj_func] initialization\n");
    val_optim = 1e6;
    g_niik_aregister2_iter = 0;
    if(g_niik_aregister2_refimg==NULL) {
      fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_refseg==NULL) {
      fprintf(stderr,"[%s] ERROR: refseg is null\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_movimg==NULL) {
      fprintf(stderr,"[%s] ERROR: movimg is null\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_refimg->datatype!=NIFTI_TYPE_FLOAT32) {
      fprintf(stderr,"[%s] ERROR: refimg is not float32\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_refseg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"[%s] ERROR: refseg is not uint8\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_movimg->datatype!=NIFTI_TYPE_FLOAT32) {
      fprintf(stderr,"[%s] ERROR: movimg is not float32\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_nmiobj==NULL && g_niik_aregister2_method==NIIK_REGISTER_NMI) {
      fprintf(stderr,"[%s] ERROR: nmiobj is null\n",fcname);
      return 0;
    }
    if(g_niik_aregister2_method==NIIK_REGISTER_NMI) { /* NMI registration */
      if(tmpdata1!=NULL) {
        free(tmpdata1);
        tmpdata1=NULL;
      }
      if(tmpdata2!=NULL) {
        free(tmpdata2);
        tmpdata2=NULL;
      }
      tmpdata1=(        float *)calloc(g_niik_aregister2_refimg->nvox,sizeof(float));
      tmpdata2=(unsigned char *)calloc(g_niik_aregister2_refimg->nvox,sizeof(char ));
    }
    if(verbose>1) fprintf(stdout,"[niik_image_aregister_obj_func] initialization success\n");
    return 1;
  } /* initialization */
  /*
   * create a matrix
   */
  if(verbose>2) fprintf(stdout,"[niik_image_aregister_obj_func] create matrix\n");
  afmat=niikmat_init(4,4);
  for(n=0; n<=17; n++) {
    affpar[n]=0;
    if(g_niik_aregister2_fixpar[n])
      affpar[n] = g_niik_aregister2_param[n];
    else
      affpar[n] = v[n];
  }
  if(verbose>2) {
    fprintf(stdout,"[niik_image_aregister_obj_func] affpar: ");
    niik_display_double_vector(affpar,17);
  }
  /* updating matrix in the objective function */
  if(verbose>2) fprintf(stdout,"[niik_image_aregister_obj_func]   create matrix from affpar\n");
  if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
    fprintf(stderr,"[%s] ERROR: niik_aregister_matrix_from_affpar_update\n",fcname);
    return NIIKMAX;
  }
  /* updating matrix with initial matrix */
  if(g_niik_aregister2_imat!=NULL) {
    if(verbose>2) fprintf(stdout,"[niik_image_aregister_obj_func]   init matrix\n");
    if(!niikmat_multiply_mat2(g_niik_aregister2_imat,afmat)) {
      fprintf(stderr,"[%s] ERROR: niikmat_multiply_mat2\n",fcname);
      return NIIKMAX;
    }
  }
  if(verbose>2) {
    /* fprintf(stdout,"  method = %i\n",g_reg_method);*/
    fprintf(stdout,"%4i | %5.1f %5.1f %5.1f | %3.0f %3.0f %3.0f | %4.2f %4.2f %4.2f | %4.2f %4.2f %4.2f | %3.0f %3.0f %3.0f \n",g_niik_aregister2_iter,
            affpar[1],affpar[2],affpar[3],
            affpar[4],affpar[5],affpar[6],
            affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
            affpar[11],affpar[12],affpar[13],
            affpar[14],affpar[15],affpar[16]);
    if(verbose>3) niikmat_display(afmat);
  }
  g_niik_aregister2_iter++;
  /* normally inverse b/c g_invmat is zero by default */
  if(!g_niik_aregister2_invmat) {
    if(!niikmat_inverse_update(afmat)) {
      fprintf(stderr,"[%s] ERROR: niikmat_inverse_update\n",fcname);
      return NIIKMAX;
    }
  }
  if(verbose>2) niikmat_display(afmat);
  /* initialize variables */
  xdim = g_niik_aregister2_refimg->nx;
  ydim = g_niik_aregister2_refimg->ny;
  area = xdim * ydim;
  bimg = (unsigned char *)g_niik_aregister2_refseg->data;
  fimg = (float         *)g_niik_aregister2_refimg->data;
  num = 0;
  movfov.x = (g_niik_aregister2_refimg->dx) * (g_niik_aregister2_refimg->nx-1);
  movfov.y = (g_niik_aregister2_refimg->dy) * (g_niik_aregister2_refimg->ny-1);
  movfov.z = (g_niik_aregister2_refimg->dz) * (g_niik_aregister2_refimg->nz-1);
  xsum = ysum = xysum = xssq = yssq = wsum = wval = 0;

  /*
   * for different registration method
   */
  switch(g_niik_aregister2_method) {

  case NIIK_REGISTER_SAD:      /* Sum of absolute difference  */
    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      if(bimg[m]==0) continue;
      p.x = (m%xdim)        * g_niik_aregister2_refimg->dx;
      p.y = ((m/xdim)%ydim) * g_niik_aregister2_refimg->dy;
      p.z = floor(m/area)   * g_niik_aregister2_refimg->dz;
      p.w = q.w = 0;
      num++;
      /* transform and check with fov, each direction separately */
      q.x = afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0) continue;
      if(q.x >= movfov.x) continue;
      q.y = afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0) continue;
      if(q.y >= movfov.y) continue;
      q.z = afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0) continue;
      if(q.z >= movfov.z) continue;
      dint1 = fimg[m];
      dint2 = niik_image_interpolate_float_image_3d_linear(g_niik_aregister2_movimg,q);
      xsum += fabs(dint2-dint1);
    }
    afmat = niikmat_free(afmat);
    return xsum;

  case NIIK_REGISTER_SSQ:      /* Sum of squared difference */
    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      if(bimg[m]==0) continue;
      p.x = (m%xdim)        * g_niik_aregister2_refimg->dx;
      p.y = ((m/xdim)%ydim) * g_niik_aregister2_refimg->dy;
      p.z = floor(m/area)   * g_niik_aregister2_refimg->dz;
      p.w = q.w = 0;
      num++;
      /* transform and check with fov, each direction separately */
      q.x = afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0) continue;
      if(q.x >= movfov.x) continue;
      q.y = afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0) continue;
      if(q.y >= movfov.y) continue;
      q.z = afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0) continue;
      if(q.z >= movfov.z) continue;
      dint1 = fimg[m];
      dint2 = niik_image_interpolate_float_image_3d_linear(g_niik_aregister2_movimg,q);
      xssq += NIIK_SQ(dint2-dint1);
    }
    afmat = niikmat_free(afmat);
    return xssq;

  case NIIK_REGISTER_NMI:    /* NMI registration */
    /* initialize histograms */
    if(verbose>1) fprintf(stdout,"[niik_image_aregister_obj_func] start NMI\n");
    for(n=0; n<2; n++)
      for(m=0; m<g_niik_aregister2_nmiobj->hnum[1]; m++)
        g_niik_aregister2_nmiobj->h1[n][m]=0;
    for(n=0; n<g_niik_aregister2_nmiobj->hnum[0]; n++)
      for(m=0; m<g_niik_aregister2_nmiobj->hnum[1]; m++)
        g_niik_aregister2_nmiobj->h2[n][m]=0;
    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      tmpdata2[m]=0;
    }
    num=0;
    #pragma omp parallel for private(p,q) reduction(+:num)
    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      if(bimg[m]==0) continue;
      num++;
      p.x = (m%xdim)        * g_niik_aregister2_refimg->dx;
      p.y = ((m/xdim)%ydim) * g_niik_aregister2_refimg->dy;
      p.z = floor(m/area)   * g_niik_aregister2_refimg->dz;
      p.w = q.w = 0;
      /* transform and check with fov, each direction separately */
      q.x = afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0) continue;
      if(q.x >= movfov.x) continue;
      q.y = afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0) continue;
      if(q.y >= movfov.y) continue;
      q.z = afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0) continue;
      if(q.z >= movfov.z) continue;
      tmpdata1[m] = niik_image_interpolate_float_image_3d_linear(g_niik_aregister2_movimg,q);
      tmpdata2[m] = 1;
    }

    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      if(tmpdata2[m]==0) continue;
      /* histogram positions */
      fi[0] = fimg[m]     / g_niik_aregister2_nmiobj->hran[0];
      fi[1] = tmpdata1[m] / g_niik_aregister2_nmiobj->hran[1];
      for(n=0; n<2; n++) {
        ni[n] = floor(fi[n]);
        nj[n] = ni[n]+1;
        ni[n] = (ni[n]>=g_niik_aregister2_nmiobj->hnum[n])?-1:ni[n];
        nj[n] = (nj[n]>=g_niik_aregister2_nmiobj->hnum[n])?-1:nj[n];
        dd[n] = fi[n] - ni[n];
        if(ni[n]>=0) {
          g_niik_aregister2_nmiobj->h1[n][ni[n]] += (1.0 - dd[n]);
        }
        if(nj[n]>=0) {
          g_niik_aregister2_nmiobj->h1[n][nj[n]] += dd[n];
        }
      }
      if(ni[0]>=0) {
        if(ni[1]>=0) {
          g_niik_aregister2_nmiobj->h2[ni[0]][ni[1]] += (1.0-dd[0]) * (1.0-dd[1]);
        }
        if(nj[1]>=0) {
          g_niik_aregister2_nmiobj->h2[ni[0]][nj[1]] += (1.0-dd[0]) * dd[1];
        }
      }
      if(nj[0]>=0) {
        if(ni[1]>=0) {
          g_niik_aregister2_nmiobj->h2[nj[0]][ni[1]] += dd[0] * (1.0-dd[1]);
        }
        if(nj[1]>=0) {
          g_niik_aregister2_nmiobj->h2[nj[0]][nj[1]] += dd[0] * dd[1];
        }
      }
      wsum+=1.0;
    } /* each sample point */
    /* entropy calculation */
    if(verbose>=3) fprintf(stdout,"[niik_image_aregister_obj_func] entropy\n");
    fi[0]=fi[1]=0;
    for(n=0; n<2; n++)
      for(m=0; m<g_niik_aregister2_nmiobj->hnum[n]; m++) {
        g_niik_aregister2_nmiobj->h1[n][m]/=num;
        if(g_niik_aregister2_nmiobj->h1[n][m]<1e-5) continue;
        fi[0]-=g_niik_aregister2_nmiobj->h1[n][m]*log10(g_niik_aregister2_nmiobj->h1[n][m]);
      }
    for(n=0; n<g_niik_aregister2_nmiobj->hnum[0]; n++)
      for(m=0; m<g_niik_aregister2_nmiobj->hnum[1]; m++) {
        g_niik_aregister2_nmiobj->h2[n][m]/=num;
        if(g_niik_aregister2_nmiobj->h2[n][m]<1e-5) continue;
        fi[1]-=g_niik_aregister2_nmiobj->h2[n][m]*log10(g_niik_aregister2_nmiobj->h2[n][m]);
      }
    /* nmi calculation */
    nmi=fi[0]/(fi[1]+1e-6)/2.0;
    if(fi[1]<0.2) nmi=0;
    if(verbose>=3) fprintf(stdout,"[niik_image_aregister_obj_func] entropy %.3f <- %.3f %.3f\n",nmi,fi[0],fi[1]);
    errval =  3.0 * NIIK_Heaviside(0.3 - wsum / num,0.1);
    errval += (1.0 - nmi);
    if(verbose>=3) fprintf(stdout,"[niik_image_aregister_obj_func] errval %.3f | wsum %.3f\n",errval,wsum);
    if(verbose>=2) fprintf(stdout,"e %6.3f | nmi %6.3f | %4i | %8i | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.3f %5.3f %5.3f | %5.2f %5.2f %5.2f | %5.1f %5.1f %5.1f\n",
                             2.0 * NIIK_Heaviside(0.3 - wsum / wval,0.1),nmi,g_niik_aregister2_iter,num,
                             affpar[1],affpar[2],affpar[3],
                             affpar[4],affpar[5],affpar[6],
                             affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                             affpar[11],affpar[12],affpar[13],
                             affpar[14],affpar[15],affpar[16]);
    if(val_optim > errval) {
      val_optim = errval;
      if(verbose>=1) {
        fprintf(stdout,"e %6.3f | nmi %6.3f | %4i | %8i | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %5.3f %5.3f %5.3f | %5.2f %5.2f %5.2f | %5.1f %5.1f %5.1f *\n",
                2.0 * NIIK_Heaviside(0.3 - wsum / wval,0.1),nmi,g_niik_aregister2_iter,num,
                affpar[1],affpar[2],affpar[3],
                affpar[4],affpar[5],affpar[6],
                affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                affpar[11],affpar[12],affpar[13],
                affpar[14],affpar[15],affpar[16]);
      }
    }
    afmat = niikmat_free(afmat);
    return errval;

  case NIIK_REGISTER_CC:    /* CC registration, aka normalized correlation (Jenkinson 2002 NeuroImage) */
    #pragma omp parallel for private(p,q,dint1,dint2,wval) reduction(+:xsum) reduction(+:ysum) reduction(+:xssq) reduction(+:yssq) reduction(+:xysum) reduction(+:wsum)
    for(m=0; m<g_niik_aregister2_refimg->nvox; m++) {
      if(bimg[m]==0) continue;
      p.x = (m%xdim)        * g_niik_aregister2_refimg->dx;
      p.y = ((m/xdim)%ydim) * g_niik_aregister2_refimg->dy;
      p.z = floor(m/area)   * g_niik_aregister2_refimg->dz;
      p.w = q.w = 0;
      num++;
      /* transform and check with fov, each direction separately */
      q.x = afmat->m[0][0]*p.x + afmat->m[0][1]*p.y + afmat->m[0][2]*p.z + afmat->m[0][3];
      if(q.x < 0) continue;
      if(q.x >= movfov.x) continue;
      q.y = afmat->m[1][0]*p.x + afmat->m[1][1]*p.y + afmat->m[1][2]*p.z + afmat->m[1][3];
      if(q.y < 0) continue;
      if(q.y >= movfov.y) continue;
      q.z = afmat->m[2][0]*p.x + afmat->m[2][1]*p.y + afmat->m[2][2]*p.z + afmat->m[2][3];
      if(q.z < 0) continue;
      if(q.z >= movfov.z) continue;
      dint1 = fimg[m];
      dint2 = niik_image_interpolate_float_image_3d_linear(g_niik_aregister2_movimg,q);
      /* apodization
       * as in Jenkinson 2002 NeuroImage */
      wval = 1; /*niik_image_aregister_obj_func_apodization(g_niik_aregister2_movimg,q,apodization_thresh);  */
      if(fabs(wval)<1e-5) continue;
      xsum  += dint1 * wval;
      ysum  += dint2 * wval;
      xssq  += dint1*dint1 * wval;
      yssq  += dint2*dint2 * wval;
      xysum += dint1*dint2 * wval;
      wsum  += wval;
    } /* each voxel */
    if(verbose>=2) {
      fprintf(stdout,"  xsum  = %12.6f\n",xsum);
      fprintf(stdout,"  ysum  = %12.6f\n",ysum);
      fprintf(stdout,"  xssq  = %12.6f\n",xssq);
      fprintf(stdout,"  yssq  = %12.6f\n",yssq);
      fprintf(stdout,"  xysum = %12.6f\n",xysum);
      fprintf(stdout,"  wsum  = %12.6f\n",wsum);
      fprintf(stdout,"  nsum  = %12i\n",num);
    }
    if(wsum>5) {
      if(wsum*xssq < xsum*xsum+1e-5) cc=-1e3;
      else if(wsum*yssq < ysum*ysum+1e-5) cc=-1e3;
      else
        cc = (wsum * xysum - xsum * ysum) / sqrt(wsum*xssq-xsum*xsum) / sqrt(wsum*yssq-ysum*ysum);
    } else {
      cc = -1e3;
    }
    errval =  3.0 * NIIK_Heaviside(0.3 - wsum / num,0.2);
    errval += (1.0 - cc);
    if(verbose>=1) {
      fprintf(stdout,"   err %-6.4f cc %-6.4f | %4i | %6i | %6.1f %6.1f %6.1f | %6.2f %6.2f %6.2f | %5.3f %5.3f %5.3f | %6.3f %6.3f %6.3f | %5.1f %5.1f %5.1f\n",
              3.0*NIIK_Heaviside(0.3 - wsum / num,0.2),cc,g_niik_aregister2_iter,(int)wsum,
              affpar[1],affpar[2],affpar[3],
              affpar[4],affpar[5],affpar[6],
              affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
              affpar[11],affpar[12],affpar[13],
              affpar[14],affpar[15],affpar[16]);
    }
    if(val_optim > errval) {
      val_optim = errval;
      if(verbose>=1) {
        fprintf(stdout,"   err %-6.4f cc %-6.4f | %4i | %6i | %6.1f %6.1f %6.1f | %6.2f %6.2f %6.2f | %5.3f %5.3f %5.3f | %6.3f %6.3f %6.3f | %5.1f %5.1f %5.1f *\n",
                3.0*NIIK_Heaviside(0.3 - wsum / num,0.2),cc,g_niik_aregister2_iter,(int)wsum,
                affpar[1],affpar[2],affpar[3],
                affpar[4],affpar[5],affpar[6],
                affpar[7]*affpar[8],affpar[7]*affpar[9],affpar[7]*affpar[10],
                affpar[11],affpar[12],affpar[13],
                affpar[14],affpar[15],affpar[16]);
      }
      /*niik_image_write("tmp_movimg.nii.gz",g_niik_aregister2_movimg);
        niik_image_write("tmp_refimg.nii.gz",g_niik_aregister2_refimg);
        exit(0);*/
    }
    afmat = niikmat_free(afmat);
    return errval;
  default:
    fprintf(stderr,"[%s] ERROR: unkown registration type, %s\n",niik_aregister_method_string(g_niik_aregister2_method),fcname);
    return NIIKMAX;
  }
  return NIIKMAX;
} /* niik_image_aregister2_obj_func */



int niik_image_aregister2(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,niikmat *imat,
                          double *affpar,double *daffpar,int maxiter, double tol,int register_method)
/*
 * niik_image_aregister2
 *   refimg     -reference target image
 *   refseg     -reference target mask image
 *   movimg     -moving image
 *   imat       -initial matrix
 *   affpar     -initial affine registration parameters
 *   daffpar    -amount of perturbation for each affine parameter
 *              -if zero, then fixed
 *   register_method
 *              -registration method
 *
 * -sampling distance and size is defined by refimg
 *
 */
{
  char fcname[50]="niik_image_aregister2";
  int
  ndim=20,
  m,n,nvar,
  cmaxiter,
  verbose=1;
  niikmat
  *p;
  double
  atol,
  plim=1e-7,
  (* pfn)();
  verbose=g_niik_aregister2_verbose;
  if(verbose>1) fprintf(stdout,"[%s] start\n",fcname);
  /* make sure images are present */
  if(refimg==NULL)  {
    fprintf(stderr,"[%s] ERROR: refimg is a null pointer\n",fcname);
    return 0;
  }
  if(movimg==NULL)  {
    fprintf(stderr,"[%s] ERROR: movimg is a null pointer\n",fcname);
    return 0;
  }
  if(affpar==NULL)  {
    fprintf(stderr,"[%s] ERROR: affpar is a null pointer\n",fcname);
    return 0;
  }
  if(daffpar==NULL) {
    fprintf(stderr,"[%s] ERROR: d'affpar is a null pointer\n",fcname);
    return 0;
  }
  if(refseg!=NULL) {
    if(refseg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"[%s] ERROR: refseg is not NIFTI_TYPE_UINT8\n",fcname);
      return 0;
    }
  }
  if(verbose>1) fprintf(stdout,"  niik_affine_register nvar\n");
  for(n=1,nvar=0; n<=16; n++)
    if(fabs(daffpar[n])>plim) nvar++;
  /* check number of variables */
  if(nvar==0) {
    fprintf(stderr,"[%s] ERROR: no variable\n",fcname);
    return 0;
  }
  /* check nmi object */
  if(register_method==NIIK_REGISTER_NMI)  {
    if(verbose>1) fprintf(stdout,"  niik_affine_register check g_niik_aregister2_nmiobj\n");
    if(g_niik_aregister2_nmiobj==NULL) {
      fprintf(stderr,"[%s] ERROR: g_niik_aregister2_nmiobj is undefined\n",fcname);
      return 0;
    }
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] niik_affine_register \n",fcname);
    fprintf(stdout,"    DOF             %2i: ",nvar);
    if(fabs(daffpar[ 1])>plim) fprintf(stdout,"Rx ");
    if(fabs(daffpar[ 2])>plim) fprintf(stdout,"Ry ");
    if(fabs(daffpar[ 3])>plim) fprintf(stdout,"Rz ");
    if(fabs(daffpar[ 4])>plim) fprintf(stdout,"Tx ");
    if(fabs(daffpar[ 5])>plim) fprintf(stdout,"Ty ");
    if(fabs(daffpar[ 6])>plim) fprintf(stdout,"Tz ");
    if(fabs(daffpar[ 7])>plim) fprintf(stdout,"Sg ");
    if(fabs(daffpar[ 8])>plim) fprintf(stdout,"Sx ");
    if(fabs(daffpar[ 9])>plim) fprintf(stdout,"Sy ");
    if(fabs(daffpar[10])>plim) fprintf(stdout,"Sz ");
    if(fabs(daffpar[11])>plim) fprintf(stdout,"Kx ");
    if(fabs(daffpar[12])>plim) fprintf(stdout,"Ky ");
    if(fabs(daffpar[13])>plim) fprintf(stdout,"Kz ");
    if(fabs(daffpar[14])>plim) fprintf(stdout,"Cx ");
    if(fabs(daffpar[15])>plim) fprintf(stdout,"Cy ");
    if(fabs(daffpar[16])>plim) fprintf(stdout,"Cz ");
    fprintf(stdout,"\n");
    if(imat==NULL) {
      fprintf(stdout,"    no initial matrix\n");
    } else {
      fprintf(stdout,"    using initial matrix\n");
    }
    fprintf(stdout,"    cost function   %s\n",niik_aregister_method_string(register_method));
    if(refseg==NULL) {
      fprintf(stdout,"    no mask\n");
    } else {
      fprintf(stdout,"    mask            %s\n",refseg->fname);
    }
    if(register_method==NIIK_REGISTER_NMI)  {
      fprintf(stdout,"    NMI hist        %5.1f %5.1f %6.1f %-5i (ref)\n",
              g_niik_aregister2_nmiobj->hmin[0],g_niik_aregister2_nmiobj->hran[0],g_niik_aregister2_nmiobj->hmax[0],g_niik_aregister2_nmiobj->hnum[0]);
      fprintf(stdout,"                    %5.1f %5.1f %6.1f %-5i (mov)\n",
              g_niik_aregister2_nmiobj->hmin[1],g_niik_aregister2_nmiobj->hran[1],g_niik_aregister2_nmiobj->hmax[1],g_niik_aregister2_nmiobj->hnum[1]);
    }
    niik_aregister_display_affine(affpar);
  } /* show info */

  /* set other global variables */
  g_niik_aregister2_refimg = refimg;
  g_niik_aregister2_refseg = refseg;
  g_niik_aregister2_movimg = movimg;

  /* update the global variables */
  g_niik_aregister2_method = register_method;
  for(n=0; n<18; n++) {
    g_niik_aregister2_param[n] = affpar[n];
  }

  for(n=1; n<=16; n++) {
    if(fabs(daffpar[n])<=1e-7)  /* no daffparation -> don't change */
      g_niik_aregister2_fixpar[n] = 1;
    else
      g_niik_aregister2_fixpar[n] = 0;
  }

  /* prepare for simplex optimization variables */
  if(verbose>1) fprintf(stdout,"    initialize\n");
  pfn=niik_image_aregister2_obj_func;
  if((pfn(NULL))==0) {
    fprintf(stderr,"[%s] ERROR: initialization for niik_image_aregister2_obj_func\n",fcname);
    return 0;
  }
  g_niik_aregister2_imat = imat;

  /* update with initial values */
  if(verbose>1) fprintf(stdout,"    initialize NM matrix\n");
  p=niikmat_init(ndim+1,ndim);
  for(m=0; m<=ndim; m++) {
    for(n=1; n<17; n++) {
      p->m[m][n] = affpar[n];
    }
  }
  for(n=1,m=0; n<17; n++) {
    p->m[m][n] += (fabs(daffpar[n])>plim) * NIIK_DMAX(0.01,p->m[m][n] * niik_get_rand() * 0.01);
  }
  for(m=1; m<17; m++) {
    if(fabs(daffpar[m])>plim) {
      p->m[m][m] += daffpar[m];
    }
  }
  for(m=17; m<=ndim; m++) {
    for(n=0; n<=16; n++) {
      p->m[m][n] += (niik_get_rand()-0.5)*2.0 * daffpar[n];
    }
  }

  /* for debug
  fprintf(stdout,"\naffpar = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7.2f ",affpar[m]); }
  fprintf(stdout,"\npeturb = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7.2f ",daffpar[m]); }
  fprintf(stdout,"\nfixpar = ");
  for(m=0;m<=16;m++) { fprintf(stdout,"%7i ",g_niik_aregister2_fixpar[m]); }
  fprintf(stdout,"\n");
  niikmat_display(p);  */

  atol = tol;
  cmaxiter = maxiter;
  niik_image_aregister2_obj_func(NULL);

  if(verbose>1) fprintf(stdout,"    start nelder mead\n");
  if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&cmaxiter)) {
    fprintf(stderr,"[%s] ERROR: nifti_k_nelder_mead\n",fcname);
    return 0;
  }

  for(n=0; n<=16; n++) {
    p->m[ndim][n] = affpar[n];
  }

  atol = tol;
  cmaxiter = maxiter;
  if(verbose>1) fprintf(stdout,"    start nelder mead (2)\n");
  if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&cmaxiter)) {
    fprintf(stderr,"[%s] ERROR: nifti_k_nelder_mead\n",fcname);
    return 0;
  }

  /* put the registration paramters into the output vector */
  if(verbose>1) fprintf(stdout,"    output vector\n");
  for(n=0; n<=16; n++) {
    affpar[n] = p->m[0][n];
  }
  p = niikmat_free(p);

  /* update error (saved at index=0) */
  affpar[0] = atol;
  if(verbose>1) {
    fprintf(stdout,"    optimized %9.5f  %i\n",atol,g_niik_aregister2_iter);
    niik_aregister_display_affine(affpar);
  }

  return 1;
} /* niik_image_aregister2 */



/*******************************************************
 *
 * nmi object
 *
 *
 ********************************************************/

nmi_obj *niik_aregister2_nmi_obj_alloc() {
  nmi_obj *obj;
  obj = (nmi_obj *)calloc(1,sizeof(nmi_obj));
  obj->nmi = 0;
  obj->hmin[0]=obj->hmin[1]=0;
  obj->hmax[0]=obj->hmax[1]=0;
  obj->hran[0]=obj->hran[1]=0;
  obj->h1=obj->h2=NULL;
  obj->hnum[0]=obj->hnum[1]=0;
  return obj;
}

void niik_aregister2_nmi_obj_free(nmi_obj *obj) {
  int n;
  if(obj==NULL) return;
  if(obj->h1!=NULL) {
    free(obj->h1[0]);
    free(obj->h1[1]);
    free(obj->h1);
  }
  if(obj->h2!=NULL) {
    for(n=0; n<obj->hnum[0]; n++) {
      free(obj->h2[n]);
    }
    free(obj->h2);
  }
  free(obj);
}

nmi_obj *niik_aregister2_nmi_obj_init(double min1,double max1,double min2,double max2,int num1,int num2) {
  nmi_obj *obj;
  int n;
  obj = niik_aregister_nmi_obj_alloc();
  obj->hnum[0]=num1;
  obj->hnum[1]=num2;
  obj->hmin[0]=min1;
  obj->hmin[1]=min2;
  obj->hmax[0]=max1;
  obj->hmax[1]=max2;
  obj->hran[0]=(max1-min1)/(num1-1);
  obj->hran[1]=(max2-min2)/(num2-1);
  obj->h1=(double **)calloc(   2,sizeof(double *));
  obj->h2=(double **)calloc(num1,sizeof(double *));
  for(n=0; n<   2; n++) {
    obj->h1[n]=(double *)calloc(obj->hnum[n],sizeof(double));  /* not sure */
  }
  for(n=0; n<num1; n++) {
    obj->h2[n]=(double *)calloc(num2,sizeof(double));
  }
  return obj;
}

void niik_image_aregister_set_nmi_obj(nmi_obj *obj) {
  g_niik_aregister2_nmiobj=obj;
}
nmi_obj *niik_image_aregister_get_nmi_obj() {
  return g_niik_aregister2_nmiobj;
}


/******************** end of nmi object *****************/


/************************************************************
 *
 * generic cost functions
 *
 ************************************************************/

int niik_image_aregister2_cost_cc(nifti_image *refimg,nifti_image *refmask,nifti_image *img,double *cc) {
  int
  i,j,k,n,
  verbose=0;
  niikpt p;
  double
  ii,ri,
  xsum,ysum,xssq,yssq,xysum,wsum;
  char fcname[64]="niik_image_aregister2_cost_cc";
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  xsum=ysum=xssq=yssq=xysum=wsum=0;
  for(k=0; k<refimg->nz; k++) {
    p.z=k*refimg->dz;
    n=k*refimg->nx*refimg->ny;
    for(j=0; j<refimg->ny; j++) {
      p.y=j*refimg->dy;
      for(i=0; i<refimg->nx; n++,i++) {
        p.x=i*refimg->dx;
        if(refmask!=NULL) {
          if((int)niik_image_get_voxel(refmask,n)==0)
            continue;
        }
        ri = niik_image_get_voxel(refimg,n);
        ii = niik_image_interpolate_3d(img,p,NIIK_INTERP_LINEAR);
        xsum += ri;
        ysum += ii;
        xssq += ri*ri;
        yssq += ii*ii;
        xysum += ii*ri;
        wsum += 1.0;
      }
    }
  }
  *cc = 0;
  if(wsum<5) {
    fprintf(stderr,"[%s] ERROR: wsum %f\n",fcname,wsum);
    return 0;
  }
  if(wsum*yssq-ysum*ysum<1e-3) {
    fprintf(stderr,"[%s] ERROR: yssq %f\n",fcname,wsum*yssq-ysum*ysum);
    return 0;
  }
  if(wsum*xssq-xsum*xsum<1e-3) {
    fprintf(stderr,"[%s] ERROR: xssq %f\n",fcname,wsum*xssq-xsum*xsum);
    return 0;
  }
  * cc = (wsum * xysum - xsum * ysum) / sqrt(wsum*xssq-xsum*xsum) / sqrt(wsum*yssq-ysum*ysum);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
} /* niik_image_aregister2_cost_cc */


int niik_image_aregister2_cost_nmi(nifti_image *refimg,nifti_image *refmask,nifti_image *img,nmi_obj *nmi) {
  int
  num=0,
  hi,hj,hk,
  ni[2],nj[2],
  i,j,k,n,m,
  verbose=0;
  niikpt p;
  double fi[2],dd[4];
  char fcname[64]="niik_image_aregister2_cost_nmi";
  if(verbose>=1) niik_fc_display(fcname,1);
  NIIK_RET0((refimg==NULL),fcname,"refimg is null");
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((nmi==NULL),fcname,"nmi is null");
  for(hk=0; hk<2; hk++) {
    NIIK_RET0((nmi->hmax[hk]<=nmi->hmin[hk]),fcname,"bad min max");
    nmi->hran[hk]=(nmi->hmax[hk]-nmi->hmin[hk])/(nmi->hnum[hk]+1.0);
    if(verbose>0) fprintf(stdout,"[%s] img %i   %9.3f %9.3f\n",fcname,hk,nmi->hmin[hk],nmi->hmax[hk]);
    for(hi=0; hi<nmi->hnum[hk]; hi++) {
      nmi->h1[hk][hi]=0;
    }
  }
  for(hi=0; hi<nmi->hnum[hk]; hi++) {
    for(hj=0; hj<nmi->hnum[1]; hj++) {
      nmi->h2[hi][hj]=0;
    }
  }
  if(verbose>1) fprintf(stdout,"[%s:%i:%s] nmi calc\n",__FILE__,__LINE__,__func__);
  p.w=0;
  for(k=n=0; k<refimg->nz; k++) {
    p.z=k*refimg->dz;
    for(j=0; j<refimg->ny; j++) {
      p.y=j*refimg->dy;
      for(i=0; i<refimg->nx; n++,i++) {
        p.x=i*refimg->dx;
        if(refmask!=NULL) {
          if((int)niik_image_get_voxel(refmask,n)==0)
            continue;
        }
        num++;
        //fi[0] = niik_image_get_voxel(img,n);
        fi[1] = niik_image_get_voxel(refimg,n);
        fi[0] = niik_image_interpolate_3d(img,p,NIIK_INTERP_LINEAR);
        /* histogram positions */
        for(hk=0; hk<2; hk++) {
          ni[hk] = floor(fi[hk]);
          nj[hk] = ni[hk]+1;
          ni[hk] = (ni[hk]>=nmi->hnum[hk])?-1:ni[hk];
          nj[hk] = (nj[hk]>=nmi->hnum[hk])?-1:nj[hk];
          dd[hk] = fi[hk] - ni[hk];
          if(ni[hk]>=0) {
            nmi->h1[hk][ni[hk]] += (1.0 - dd[hk]);
          }
          if(nj[hk]>=0) {
            nmi->h1[hk][nj[hk]] += dd[hk];
          }
        }
        if(ni[0]>=0) {
          if(ni[1]>=0) {
            nmi->h2[ni[0]][ni[1]] += (1.0-dd[0]) * (1.0-dd[1]);
          }
          if(nj[1]>=0) {
            nmi->h2[ni[0]][nj[1]] += (1.0-dd[0]) * dd[1];
          }
        }
        if(nj[0]>=0) {
          if(ni[1]>=0) {
            nmi->h2[nj[0]][ni[1]] += dd[0] * (1.0-dd[1]);
          }
          if(nj[1]>=0) {
            nmi->h2[nj[0]][nj[1]] += dd[0] * dd[1];
          }
        }
      }
    }
  } /* each sample point */
  /* entropy calculation */
  if(verbose>=3) fprintf(stdout,"[%s:%i:%s] entropy\n",__FILE__,__LINE__,__func__);
  fi[0]=fi[1]=0;
  for(hk=0; hk<2; hk++)
    for(m=0; m<nmi->hnum[hk]; m++) {
      nmi->h1[hk][m]/=num;
      if(nmi->h1[hk][m]<1e-5) continue;
      fi[0]-=nmi->h1[hk][m]*log10(nmi->h1[hk][m]);
    }
  for(hk=0; hk<nmi->hnum[0]; hk++)
    for(m=0; m<nmi->hnum[1]; m++) {
      nmi->h2[hk][m]/=num;
      if(nmi->h2[hk][m]<1e-5) continue;
      fi[1]-=nmi->h2[hk][m]*log10(nmi->h2[hk][m]);
    }
  /* nmi calculation */
  nmi->nmi=fi[0]/(fi[1]+1e-6)/2.0;
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_aregister2_cost_nmi */



/************************************************************
 *
 * registration functions
 *
 *
 * niik_image_aregister2_multilevel
 * -generic function for multilevel registration
 * -need to pre-define various parameters
 *
 * niik_image_aregister2_multiseed
 * -not well validated...
 *
 ************************************************************/

int niik_image_aregister2_test1(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int dof,int areg_cost,nmi_obj *nmiobj) {
  char fcname[48]="niik_image_aregister2_test1";
  niikvec
  *FWHM=NULL,*delta=NULL,*tol=NULL;
  niikmat *daffpar=NULL,*seed=NULL;
  int
  i,j,
  nlevel,ns,
  *nseed=NULL,*maxiter=NULL;

  if(refimg==NULL) {
    fprintf(stderr,"[%s] refimg is null\n",fcname);
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"[%s] img is null\n",fcname);
    return 0;
  }
  if(areg_cost == NIIK_REGISTER_NMI && nmiobj==NULL ) {
    fprintf(stderr,"[%s] nmiobj is null\n",fcname);
    return 0;
  }

  /* parameters */
  nlevel=5;
  nseed=(int *)calloc(nlevel,sizeof(int));
  maxiter=(int *)calloc(nlevel,sizeof(int));
  FWHM=niikvec_init(nlevel);
  delta=niikvec_init(nlevel);
  tol=niikvec_init(nlevel);
  daffpar=niikmat_init(nlevel,20);

  switch(dof) {
  case 6:
    nseed[0]=100;
    nseed[1]=25;
    nseed[2]=12;
    nseed[3] = 6;
    nseed[4] = 2;
    maxiter[0]=40;
    maxiter[1]=80;
    maxiter[2]=120;
    maxiter[3]=150;
    maxiter[4]=200;
    FWHM->v[0]=16;
    FWHM->v[1]=12;
    FWHM->v[2]=8;
    FWHM->v[3]=4;
    FWHM->v[4]=2;
    delta->v[0]=12.1;
    delta->v[1]=8.1;
    delta->v[2]=6.4;
    delta->v[3]=3.9;
    delta->v[4]=1.9;
    tol->v[0]=1e-4;
    tol->v[1]=51e-5;
    tol->v[2]=1e-5;
    tol->v[3]=5e-6;
    tol->v[4]=1e-6;
    i=0;
    daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=30;
    daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=40;
    for(i=1; i<nlevel; i++) {
      daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=daffpar->m[i-1][1]/2.0;
      daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=daffpar->m[i-1][4]/2.0;
    }
    ns = niik_get_max_from_int_vector(nseed,nlevel);
    seed = niikmat_init(ns,17);
    for(i=0; i<ns; i++) {
      for(j=0; j<17; j++)
        seed->m[i][j]=affpar[j];
    }
    for(i=1; i<ns/2; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
    }
    for(; i<ns; i++) { /*VF: ?*/
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
    }
    /*fprintf(stdout,"\tseed display\n");
      niikmat_display(seed); */
    break;

  case 7:
    nseed  [0]=151;
    nseed[1]=35;
    nseed[2]=12;
    nseed[3] = 6;
    nseed[4] = 2;
    maxiter[0]=40;
    maxiter[1]=80;
    maxiter[2]=120;
    maxiter[3]=150;
    maxiter[4]=200;
    FWHM->v[0]=16;
    FWHM->v[1]=12;
    FWHM->v[2]=8;
    FWHM->v[3]=4;
    FWHM->v[4]=2;
    delta->v[0]=12.1;
    delta->v[1]=8.1;
    delta->v[2]=6.4;
    delta->v[3]=3.9;
    delta->v[4]=1.9;
    tol->v[0]=1e-4;
    tol->v[1]=51e-5;
    tol->v[2]=1e-5;
    tol->v[3]=5e-6;
    tol->v[4]=1e-6;
    i=0;
    daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=30;
    daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=40;
    daffpar->m[i][7]=0.2;
    for(i=1; i<nlevel; i++) {
      daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=daffpar->m[i-1][1]/2.0;
      daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=daffpar->m[i-1][4]/2.0;
      daffpar->m[i][7]=daffpar->m[i-1][7]/1.2;
    }
    ns = niik_get_max_from_int_vector(nseed,nlevel);
    seed = niikmat_init(ns,17);
    for(i=0; i<ns; i++) {
      for(j=0; j<17; j++)
        seed->m[i][j]=affpar[j];
    }
    for(i=1; i<ns*0.5; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][7] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for(; i<ns*0.6; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
    }
    for(; i<ns*0.7; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
    }
    for(; i<ns*0.8; i++) {
      seed->m[i][7] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for(; i<ns*0.9; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][7] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for(; i<ns; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][7] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    break;

  case 9:
    nseed  [0]=163;
    nseed[1]=35;
    nseed[2]=12;
    nseed[3] = 6;
    nseed[4] = 2;
    maxiter[0]=40;
    maxiter[1]=80;
    maxiter[2]=120;
    maxiter[3]=150;
    maxiter[4]=200;
    FWHM->v[0]=16;
    FWHM->v[1]=12;
    FWHM->v[2]=8;
    FWHM->v[3]=4;
    FWHM->v[4]=2;
    delta->v[0]=12.1;
    delta->v[1]=8.1;
    delta->v[2]=6.4;
    delta->v[3]=3.9;
    delta->v[4]=1.9;
    tol->v[0]=1e-4;
    tol->v[1]=51e-5;
    tol->v[2]=1e-5;
    tol->v[3]=5e-6;
    tol->v[4]=1e-6;
    i=0;
    daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=30;
    daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=40;
    daffpar->m[i][8]=daffpar->m[i][9]=daffpar->m[i][10]=0.4;
    for(i=1; i<nlevel; i++) {
      daffpar->m[i][1]=daffpar->m[i][2]=daffpar->m[i][3]=daffpar->m[i-1][1]/2.0;
      daffpar->m[i][4]=daffpar->m[i][5]=daffpar->m[i][6]=daffpar->m[i-1][4]/2.0;
      daffpar->m[i][8]=daffpar->m[i][9]=daffpar->m[i][10]=daffpar->m[i-1][8]/1.2;
    }
    ns = niik_get_max_from_int_vector(nseed,nlevel);
    seed = niikmat_init(ns,17);
    for(i=0; i<ns; i++) {
      for(j=0; j<17; j++)
        seed->m[i][j]=affpar[j];
    }
    for(i=1; i<ns*0.5; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.6; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
    }
    for( ; i<ns*0.7; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
    }
    for( ; i<ns*0.8; i++) {
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.9; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    break;

  case 12:
    nseed  [0]=215;
    nseed[1]=35;
    nseed[2]=12;
    nseed[3] = 6;
    nseed[4] = 2;
    maxiter[0]=40;
    maxiter[1]=80;
    maxiter[2]=120;
    maxiter[3]=150;
    maxiter[4]=200;
    FWHM->v[0]=16;
    FWHM->v[1]=12;
    FWHM->v[2]=8;
    FWHM->v[3]=4;
    FWHM->v[4]=2;
    delta->v[0]=12.1;
    delta->v[1]=8.1;
    delta->v[2]=6.4;
    delta->v[3]=3.9;
    delta->v[4]=1.9;
    tol->v[0]=1e-4;
    tol->v[1]=51e-5;
    tol->v[2]=1e-5;
    tol->v[3]=5e-6;
    tol->v[4]=1e-6;
    i=0;
    daffpar->m[i][ 1]=daffpar->m[i][ 2]=daffpar->m[i][ 3]=30;
    daffpar->m[i][ 4]=daffpar->m[i][ 5]=daffpar->m[i][ 6]=40;
    daffpar->m[i][ 8]=daffpar->m[i][ 9]=daffpar->m[i][10]=0.4;
    daffpar->m[i][11]=daffpar->m[i][12]=daffpar->m[i][13]=0.4;
    for(i=1; i<nlevel; i++) {
      daffpar->m[i][ 1]=daffpar->m[i][ 2]=daffpar->m[i][ 3]=daffpar->m[i-1][ 1]/2.0;
      daffpar->m[i][ 4]=daffpar->m[i][ 5]=daffpar->m[i][ 6]=daffpar->m[i-1][ 4]/2.0;
      daffpar->m[i][ 8]=daffpar->m[i][ 9]=daffpar->m[i][10]=daffpar->m[i-1][ 8]/1.2;
      daffpar->m[i][11]=daffpar->m[i][12]=daffpar->m[i][13]=daffpar->m[i-1][11]/1.2;
    }
    ns = niik_get_max_from_int_vector(nseed,nlevel);
    seed = niikmat_init(ns,17);
    for(i=0; i<ns; i++) {
      for(j=0; j<17; j++)
        seed->m[i][j]=affpar[j];
    }
    for(i=1; i<ns*0.5; i++) {
      seed->m[i][ 1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
      niik_get_rand();
    }
    for(niik_get_rand(); i<ns*0.55; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
    }
    for( ; i<ns*0.60; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
    }
    for( ; i<ns*0.65; i++) {
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.70; i++) {
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns*0.75; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( niik_get_rand(); i<ns-0.80; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.85; i++) {
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for(; i<ns*0.85; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns*0.90; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( ; i<ns*0.95; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns; i++) {
      seed->m[i][ 1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    break;

  case 15:
    nseed  [0]=450;
    nseed[1]=35;
    nseed[2]=12;
    nseed[3] = 6;
    nseed[4] = 2;
    maxiter[0]=40;
    maxiter[1]=80;
    maxiter[2]=120;
    maxiter[3]=150;
    maxiter[4]=200;
    FWHM->v[0]=16;
    FWHM->v[1]=12;
    FWHM->v[2]=8;
    FWHM->v[3]=4;
    FWHM->v[4]=2;
    delta->v[0]=12.1;
    delta->v[1]=8.1;
    delta->v[2]=6.4;
    delta->v[3]=3.9;
    delta->v[4]=1.9;
    tol->v[0]=1e-4;
    tol->v[1]=51e-5;
    tol->v[2]=1e-5;
    tol->v[3]=5e-6;
    tol->v[4]=1e-6;
    i=0;
    daffpar->m[i][ 1]=daffpar->m[i][ 2]=daffpar->m[i][ 3]=30;
    daffpar->m[i][ 4]=daffpar->m[i][ 5]=daffpar->m[i][ 6]=40;
    daffpar->m[i][ 8]=daffpar->m[i][ 9]=daffpar->m[i][10]=0.4;
    daffpar->m[i][11]=daffpar->m[i][12]=daffpar->m[i][13]=0.4;
    for(i=1; i<nlevel; i++) {
      daffpar->m[i][ 1]=daffpar->m[i][ 2]=daffpar->m[i][ 3]=daffpar->m[i-1][ 1]/2.0;
      daffpar->m[i][ 4]=daffpar->m[i][ 5]=daffpar->m[i][ 6]=daffpar->m[i-1][ 4]/2.0;
      daffpar->m[i][ 8]=daffpar->m[i][ 9]=daffpar->m[i][10]=daffpar->m[i-1][ 8]/1.2;
      daffpar->m[i][11]=daffpar->m[i][12]=daffpar->m[i][13]=daffpar->m[i-1][11]/1.2;
    }
    ns = niik_get_max_from_int_vector(nseed,nlevel);
    seed = niikmat_init(ns,17);
    for(i=0; i<ns; i++) {
      for(j=0; j<17; j++)
        seed->m[i][j]=affpar[j];
    }
    for(i=1; i<ns*0.5; i++) {
      seed->m[i][ 1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
      niik_get_rand();
    }
    for( niik_get_rand(); i<ns*0.55; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
    }
    for( ; i<ns*0.60; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
    }
    for( ; i<ns*0.65; i++) {
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.70; i++) {
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns*0.75; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( niik_get_rand(); i<ns-0.80; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
    }
    for( ; i<ns*0.85; i++) {
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( ; i<ns*0.85; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns*0.90; i++) {
      seed->m[i][4] = (niik_get_rand()-0.5) * 150;
      seed->m[i][5] = (niik_get_rand()-0.5) * 150;
      seed->m[i][6] = (niik_get_rand()-0.5) * 150;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( ; i<ns*0.95; i++) {
      seed->m[i][1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for( niik_get_rand(); i<ns; i++) {
      seed->m[i][ 1] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 2] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 3] = (niik_get_rand()-0.5) * 360;
      seed->m[i][ 8] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][ 9] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][10] = (niik_get_rand()-0.5) * 0.6 + 1.0;
      seed->m[i][11] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][12] = (niik_get_rand()-0.5) * 0.6;
      seed->m[i][13] = (niik_get_rand()-0.5) * 0.6;
    }
    for(i=1; i<ns; i++) {
      switch(i%3) {
      case 0:
        seed->m[i][ 8] *= -1;
        break;
      case 1:
        seed->m[i][ 9] *= -1;
        break;
      case 2:
        seed->m[i][10] *= -1;
        break;
      }
    }
    break;

  default:
    fprintf(stdout,"[niik_image_aregister2_test1] ERROR: unkown dof %i\n",dof);
    return 0;
  }

  if(!niik_image_aregister2_multilevel(refimg,refseg,img,afmat,affpar,areg_cost,nmiobj,nlevel,FWHM,delta,tol,nseed,maxiter,daffpar,seed)) {
    fprintf(stderr,"[niik_image_aregister2_test1] ERROR: niik_image_aregister2_multilevel\n");
    return 0;
  }

  free(nseed);
  free(maxiter);
  FWHM=niikvec_free(FWHM);
  delta=niikvec_free(delta);
  tol=niikvec_free(tol);
  daffpar=niikmat_free(daffpar);
  return 1;
} /* niik_image_aregister2_test1 */


int niik_image_aregister2_test2(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,
                                int nlevel,double fFWHM,double fmaxiter,int fnseed,double fdelta,niikvec *idaffpar,nmi_obj *nmiobj) {
  niikvec *range=NULL;
  char fcname[32]="niik_image_aregister2_test2";
  int verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);
  range=niikvec_init(20);
  range->v[1]=range->v[2]=range->v[3]=30.0;
  range->v[4]=range->v[5]=range->v[6]=30.0;
  range->v[7]=range->v[8]=range->v[9]=range->v[10]=0.2;
  NIIK_RET0((!niik_image_aregister2_test_range(refimg,refseg,img,afmat,affpar,areg_cost,nlevel,fFWHM,fmaxiter,fnseed,fdelta,range,idaffpar,nmiobj)),
            fcname,
            "niik_image_aregister2_test_range");
  range=niikvec_free(range);
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
}

int niik_image_aregister2_test_range(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,
                                     int nlevel,double fFWHM,double fmaxiter,int fnseed,double fdelta,niikvec *range,niikvec *idaffpar,nmi_obj *nmiobj) {
  char fcname[64]="niik_image_aregister2_test2";
  niikvec
  *FWHM,*delta,*tol;
  niikmat *daffpar,*seed;
  int
  i,j,
  ns,
  *nseed,*maxiter;

  if(refimg==NULL) {
    fprintf(stderr,"[%s] refimg is null\n",fcname);
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"[%s] img is null\n",fcname);
    return 0;
  }
  if(areg_cost==NIIK_REGISTER_NMI && nmiobj==NULL ) {
    fprintf(stderr,"[%s] nmiobj is null\n",fcname);
    return 0;
  }

  /* parameters */
  nseed=(int *)calloc(nlevel,sizeof(int));
  maxiter=(int *)calloc(nlevel,sizeof(int));
  FWHM=niikvec_init(nlevel);
  delta=niikvec_init(nlevel);
  tol=niikvec_init(nlevel);
  daffpar=niikmat_init(nlevel,20);

  nseed[nlevel-1]=fnseed;
  maxiter[nlevel-1]=fmaxiter;
  FWHM->v[nlevel-1]=fFWHM;
  delta->v[nlevel-1]=fdelta;

  for(i=nlevel-2,j=1; i>=0; i--,j++) {
    nseed[i]=nseed[i+1]*pow(5,j);
    maxiter[i]=NIIK_IMAX(20,maxiter[nlevel-1]/2);
    FWHM->v[i]=FWHM->v[i+1]*2.0;
    delta->v[i]=delta->v[i+1]*2.0;
  }
  for(i=0; i<nlevel; i++)
    tol->v[i]=pow(10.0,-i-4.0);

  i=0;
  for(j=0; j<17; j++) daffpar->m[i][j]=idaffpar->v[j];
  for(i=1; i<nlevel; i++) {
    daffpar->m[i][ 1]=daffpar->m[i][ 2]=daffpar->m[i][ 3]=daffpar->m[i-1][ 1]/2.0;
    daffpar->m[i][ 4]=daffpar->m[i][ 5]=daffpar->m[i][ 6]=daffpar->m[i-1][ 4]/2.0;
    daffpar->m[i][ 8]=daffpar->m[i][ 9]=daffpar->m[i][10]=daffpar->m[i-1][ 8]/2.0;
    daffpar->m[i][11]=daffpar->m[i][12]=daffpar->m[i][13]=daffpar->m[i-1][11]/2.0;
  }
  ns = niik_get_max_from_int_vector(nseed,nlevel);
  seed = niikmat_init(ns,17);
  for(i=0; i<ns; i++) {
    for(j=0; j<17; j++)
      seed->m[i][j]=affpar[j];
  }
  for(i=1; i<ns*0.5; i++) {
    seed->m[i][ 1] = (niik_get_rand()-0.5) * range->v[1];
    seed->m[i][ 2] = (niik_get_rand()-0.5) * range->v[2];
    seed->m[i][ 3] = (niik_get_rand()-0.5) * range->v[3];
    seed->m[i][ 4] = (niik_get_rand()-0.5) * range->v[4];
    seed->m[i][ 5] = (niik_get_rand()-0.5) * range->v[5];
    seed->m[i][ 6] = (niik_get_rand()-0.5) * range->v[6];
    seed->m[i][ 8] = (niik_get_rand()-0.5) * range->v[8]  + 1.0;
    seed->m[i][ 9] = (niik_get_rand()-0.5) * range->v[9]  + 1.0;
    seed->m[i][10] = (niik_get_rand()-0.5) * range->v[10] + 1.0;
  }
  for( ; i<ns*0.6; i++) {
    seed->m[i][1] = (niik_get_rand()-0.5) * range->v[1];
    seed->m[i][2] = (niik_get_rand()-0.5) * range->v[2];
    seed->m[i][3] = (niik_get_rand()-0.5) * range->v[3];
  }
  for( ; i<ns*0.7; i++) {
    seed->m[i][4] = (niik_get_rand()-0.5) * range->v[4];
    seed->m[i][5] = (niik_get_rand()-0.5) * range->v[5];
    seed->m[i][6] = (niik_get_rand()-0.5) * range->v[6];
  }
  for( ; i<ns*0.8; i++) {
    seed->m[i][ 8] = (niik_get_rand()-0.5) * range->v[8]  + 1.0;
    seed->m[i][ 9] = (niik_get_rand()-0.5) * range->v[9]  + 1.0;
    seed->m[i][10] = (niik_get_rand()-0.5) * range->v[10] + 1.0;
  }
  for( ; i<ns*0.9; i++) {
    seed->m[i][1] = (niik_get_rand()-0.5) * range->v[1];
    seed->m[i][2] = (niik_get_rand()-0.5) * range->v[2];
    seed->m[i][3] = (niik_get_rand()-0.5) * range->v[3];
    seed->m[i][ 8] = (niik_get_rand()-0.5) * range->v[ 8] + 1.0;
    seed->m[i][ 9] = (niik_get_rand()-0.5) * range->v[ 9] + 1.0;
    seed->m[i][10] = (niik_get_rand()-0.5) * range->v[10] + 1.0;
  }
  for( ; i<ns; i++) {
    seed->m[i][4] = (niik_get_rand()-0.5) * range->v[4];
    seed->m[i][5] = (niik_get_rand()-0.5) * range->v[5];
    seed->m[i][6] = (niik_get_rand()-0.5) * range->v[6];
    seed->m[i][ 8] = (niik_get_rand()-0.5) * range->v[ 8] + 1.0;
    seed->m[i][ 9] = (niik_get_rand()-0.5) * range->v[ 9] + 1.0;
    seed->m[i][10] = (niik_get_rand()-0.5) * range->v[10] + 1.0;
  }
  /*niikmat_display(seed);*/

  if(!niik_image_aregister2_multilevel(refimg,refseg,img,afmat,affpar,areg_cost,nmiobj,nlevel,FWHM,delta,tol,nseed,maxiter,daffpar,seed)) {
    fprintf(stderr,"[%s] ERROR: niik_image_aregister2_multilevel\n",fcname);
    return 0;
  }

  free(nseed);
  free(maxiter);
  FWHM=niikvec_free(FWHM);
  delta=niikvec_free(delta);
  tol=niikvec_free(tol);
  daffpar=niikmat_free(daffpar);

  return 1;
}


int niik_image_aregister2_multilevel(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,
                                     int areg_cost,nmi_obj *nmiobj,
                                     int nlevel,niikvec *FWHM, niikvec *delta, niikvec *tol,int *nseed,int *maxiter,niikmat *daffpar,niikmat *seed) {
  char fcname[128]="niik_image_aregister2_multilevel";
  nifti_image **refseg_list=NULL;
  int n;
  if((refseg_list=(nifti_image **)calloc(nlevel,sizeof(nifti_image *)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc\n",fcname);
    return 0;
  }
  for(n=0; n<nlevel; n++) {
    if(refseg==NULL) {
      refseg_list[n]=NULL;
      continue;
    }
    if((refseg_list[n] = niik_image_copy_as_type(refseg,NIFTI_TYPE_UINT8))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
      return 0;
    }
  }
  if(!niik_image_aregister2_multilevel_multimask(refimg,refseg_list,img,afmat,affpar,areg_cost,nmiobj,nlevel,FWHM,delta,tol,nseed,maxiter,daffpar,seed)) {
    fprintf(stderr,"[%s] ERROR: niik_image_aregister2_multilevel_multimask\n",fcname);
    return 0;
  }
  for(n=0; n<nlevel; n++) {
    refseg_list[n] = niik_image_free(refseg_list[n]);
  }
  free(refseg_list);
  return 1;
} /* niik_image_aregister2_multilevel */


int niik_image_aregister2_multilevel_multimask(nifti_image *refimg,nifti_image **refseg,nifti_image *img,niikmat *afmat,double *affpar,
    int areg_cost,nmi_obj *nmiobj,
    int nlevel,niikvec *FWHM, niikvec *delta, niikvec *tol,int *nseed,int *maxiter,
    niikmat *daffpar,niikmat *seed)
/* main registration function for multilevel with multiple masks
 *
 * refimg      : reference grayscale image
 * refseg      : reference's mask for each level
 * img         : moving image
 * afmat       : initial registration matrix
 * affpar      : registration parameters (results replace here)
 * areg_cost   : affine registration cost method
 * nlevel      : number of levels
 * FWHM        : nlevel-long vector of blurring kernel (FWHM)
 * delta       : nlevel-long vector of sampling distance (in mm)
 * tol         : nlevel-long vector of registration tolerance
 * nseed       : nlevel-long vector of number of seed points per level
 * nseed       : nlevel-long vector of max number of iteration per level
 * daffpar     : nlevel-by-16 matrix containing the magnitude of
 *               perturbations per iteration for each parameter
 * seed        : max-of-nseed[*]-by-16 matrix containing the initial parameters
 *             : these parameters are updated at each level
 */
{
  char fcname[64]="niik_image_aregister2_multilevel_multimask";
  nifti_image
  *cimg=NULL,
   *crefseg=NULL,
    *crefimg=NULL;
  int
  verbose=1,
  cmaxiter,sumiter=0,
           i,j,k,
           nl;
  niikmat
  *tmpmat;
  double
  atol;
  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  srand(time(NULL));
  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  if(verbose) fprintf(stdout,"[niik_image_aregister2_multilevel] start at %s\n",tmstr);

  /* check inputs */
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(nlevel!=FWHM->num) {
    fprintf(stderr,"[%s] ERROR: FWHM size (%i) is not %i\n",fcname,FWHM->num,nlevel);
    return 0;
  }
  if(nlevel!=delta->num) {
    fprintf(stderr,"[%s] ERROR: delta size (%i) is not %i\n",fcname,delta->num,nlevel);
    return 0;
  }
  if(nlevel!=tol->num) {
    fprintf(stderr,"[%s] ERROR: tol size (%i) is not %i\n",fcname,tol->num,nlevel);
    return 0;
  }
  if(nlevel!=daffpar->row) {
    fprintf(stderr,"[%s] ERROR: daffpar row (%i) is not %i\n",fcname,daffpar->row,nlevel);
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"\tseed display\n");
    niikmat_display(seed);
  }

  fprintf(stdout,"[niik_image_aregister2_multilevel] parameters\n");
  fprintf(stdout,"  #level      %i\n",nlevel);
  fprintf(stdout,"  FWHM        ");
  for(nl=0; nl<nlevel; nl++) fprintf(stdout,"%8.5f ",FWHM->v[nl]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  delta       ");
  for(nl=0; nl<nlevel; nl++) fprintf(stdout,"%8.5f ",delta->v[nl]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  tol         ");
  for(nl=0; nl<nlevel; nl++) fprintf(stdout,"%8.1e ",tol->v[nl]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  maxiter     ");
  for(nl=0; nl<nlevel; nl++) fprintf(stdout,"%8i ",maxiter[nl]);
  fprintf(stdout,"\n");
  fprintf(stdout,"  #seed       ");
  for(nl=0; nl<nlevel; nl++) fprintf(stdout,"%8i ",nseed[nl]);
  fprintf(stdout,"\n");
  for(nl=0; nl<nlevel; nl++) {
    fprintf(stdout,"  daffpar %3i  ",nl);
    for(i=1; i<17; i++) fprintf(stdout,"%8.5f ",daffpar->m[nl][i]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"  cost        %s, %i\n",niik_aregister_method_string(areg_cost),areg_cost);

  if(areg_cost == NIIK_REGISTER_NMI) {
    if(nmiobj==NULL) {
      fprintf(stderr,"[%s] ERROR: nmi obj is missing\n",fcname);
      return 0;
    }
    g_niik_aregister2_nmiobj=nmiobj;
  } else {
    g_niik_aregister2_nmiobj=NULL;
  }

  if(afmat!=NULL) {
    fprintf(stdout,"  initial matrix\n");
    niikmat_display(afmat);
  }

  for(nl=0; nl<nlevel; nl++) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    sumiter=0;

    fprintf(stdout,"[niik_image_aregister2_multilevel] level %-2i %s  FWHM=%-7.3f delta=%-7.3f #seed=%-4i tol=%5.1e maxiter=%-4i\n",nl+1,tmstr,FWHM->v[nl],delta->v[nl],nseed[nl],tol->v[nl],maxiter[nl]);

    /* prepare input images
     * -blur moving and reference images
     * -ref image is resampled (linearly)
     * -if ref mask is not used, create it
     * -else ref mask is also resampled (nearest neighbor)
     */
    if(cimg==NULL) {
      if((cimg = niik_image_filter_gaussian(img,FWHM->v[nl]*2,FWHM->v[nl]))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
        return 0;
      }
      if(!niik_image_type_convert(cimg,NIFTI_TYPE_FLOAT32)) {
        fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
        return 0;
      }
    }
    if(crefimg==NULL) {
      if(verbose>=2) {
        fprintf(stdout,"[%s] running niik_image_filter_gaussian\n",fcname);
      }
      if((crefimg = niik_image_filter_gaussian(refimg,FWHM->v[nl]*2,FWHM->v[nl]))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
        return 0;
      }
      if(verbose>=2) {
        fprintf(stdout,"[%s] running niik_image_resample_3d_update\n",fcname);
      }
      if(!niik_image_resample_3d_update(crefimg,delta->v[nl],delta->v[nl],delta->v[nl],-1,-1,-1,NIIK_INTERP_LINEAR)) {
        fprintf(stderr,"[%s] ERROR: niik_image_resample_3d_update\n",fcname);
        return 0;
      }
      if(!niik_image_type_convert(crefimg,NIFTI_TYPE_FLOAT32)) {
        fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
        return 0;
      }
      if(refseg[nl]==NULL) { /* no refseg */
        if(verbose>=2) {
          fprintf(stdout,"[%s] making refseg\n",fcname);
        }
        if((crefseg=niik_image_copy(refimg))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
          return 0;
        }
        if(!niik_image_type_convert(crefseg,NIFTI_TYPE_UINT8)) {
          fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
          return 0;
        }
        if(!niik_image_one(crefseg)) {
          fprintf(stderr,"[%s] ERROR: niik_image_one\n",fcname);
          return 0;
        }
        free(crefseg->fname);
        crefseg->fname=(char *)calloc(10,sizeof(char));
        sprintf(crefseg->fname,"ones.nii");
      } /* no refseg */
      else if((crefseg=niik_image_affine_transform_3d(refseg[nl],crefimg,NULL,NIIK_INTERP_NN))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d\n",fcname);
        return 0;
      } /* with refseg */
    } /* crefimg is NULL */

    if(1) {
      for(i=0; i<nseed[nl]; i++) {
        atol = tol->v[nl];
        cmaxiter = maxiter[nl];
        if(verbose>=2) {
          fprintf(stdout,"[%s] running niik_image_aregister2, seed ->\n",fcname);
          niik_display_double_vector(seed->m[i],17);
          niik_aregister_display_affine(seed->m[i]);
        }
        if(!niik_image_aregister2(crefimg,crefseg,cimg,afmat,seed->m[i],daffpar->m[nl],cmaxiter,atol,areg_cost)) {
          fprintf(stderr,"[%s] ERROR: niik_image_aregister2\n",fcname);
          return 0;
        }
        sumiter+=cmaxiter;
        if(verbose>=1) {
          fprintf(stdout,"[niik_image_aregister2_multilevel] level %i seed %-4i / %-4i %12.8f\n",nl+1,i,nseed[nl],1.0-seed->m[i][0]);
          if(verbose>=0) {
            fprintf(stdout,"[niik_image_aregister2_multilevel] level %i seed %-4i        %12.8f | ",nl+1,i,1.0-seed->m[i][0]);
            niik_display_double_vector(seed->m[i]+1,16);
          }
        }
        for(j=1; j<=3; j++) {
          seed->m[i][j]=(seed->m[i][j]> 180)?seed->m[i][j]-360:seed->m[i][j];
          seed->m[i][j]=(seed->m[i][j]<-180)?seed->m[i][j]+360:seed->m[i][j];
        }
      } /* each seed */
      /* niikmat_display(seed); */
    } else {
      if(!niik_image_aregister2_multiseed(crefimg,crefseg,cimg,afmat,affpar,areg_cost,nmiobj,0,0,tol->v[nl],maxiter[nl],daffpar->m[nl],seed)) {
        fprintf(stderr,"[%s] ERROR: niik_image_aregister2_multiseed\n",fcname);
        return 0;
      }
    }

    /* re-order the matrix */
    for(i=0; i<nseed[nl]; i++) {
      for(j=i+1; j<nseed[nl]; j++) {
        if(seed->m[i][0]>seed->m[j][0]) {
          for(k=0; k<17; k++)
            NIIK_DSWAP(&seed->m[i][k],&seed->m[j][k]);
        }
      }
    }

    /*fprintf(stdout,"free %i %0.4e\n",nseed[nl],tol->v[nl]);*/
    if(nl+1==nlevel) {
      cimg    = niik_image_free(cimg);
      crefimg = niik_image_free(crefimg);
      crefseg = niik_image_free(crefseg);
    } else if(FWHM->v[nl]==FWHM->v[nl+1]) {
      fprintf(stdout,"[%s]   not freeing image for next level\n",fcname);
    } else {
      cimg    = niik_image_free(cimg);
      crefimg = niik_image_free(crefimg);
      crefseg = niik_image_free(crefseg);
    }

    fprintf(stdout,"[%s] optimized cost %.5f  %6i\n",fcname,seed->m[0][0],sumiter);
    niik_aregister_display_affine(seed->m[0]);
  }  /* each level */

  for(k=0; k<17; k++) {
    affpar[k]=seed->m[0][k];
  }
  niik_aregister_display_affine(affpar);

  tmpmat = niik_aregister_matrix_from_affpar(affpar);
  niikmat_display(tmpmat);
  tmpmat=niikmat_free(tmpmat);

  return 1;
} /* niik_image_aregister2_multilevel */



int niik_image_aregister2_multiseed(nifti_image *refimg,nifti_image *refseg,
                                    nifti_image *img,niikmat *afmat,
                                    double *affpar,int areg_cost,nmi_obj *nmiobj,
                                    double FWHM, double delta, double tol,
                                    int maxiter, double *daffpar, niikmat *seed) {
  nifti_image
  *cimg,
  *crefseg,
  *crefimg;
  double atol,plim=1e-7;
  int
  nseed,
  nvar,
  i,j,k,cmaxiter,
  verbose=1;
  char fcname[64]="niik_image_aregister2_multiseed";

  if(verbose>=1) niik_fc_display(fcname,1);
  if(FWHM>0) {
    if(verbose>=2) fprintf(stdout,"[%s] gaussian blurring FWHM = %-9.3f\n",fcname,FWHM);
    if((crefimg = niik_image_filter_gaussian(refimg,FWHM*2,FWHM))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
      return 0;
    }
    if((cimg = niik_image_filter_gaussian(img,FWHM*2,FWHM))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian\n",fcname);
      return 0;
    }
  } else {
    if(verbose>=2) fprintf(stdout,"[%s] no gaussian blurring\n",fcname);
    if((crefimg = niik_image_copy(refimg))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
      return 0;
    }
    if((cimg = niik_image_copy(img))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
      return 0;
    }
  }
  if(delta>0) {
    if(!niik_image_resample_3d_update(crefimg,delta,delta,delta,-1,-1,-1,NIIK_INTERP_LINEAR)) {
      fprintf(stderr,"[%s] ERROR: niik_image_resample_3d_update\n",fcname);
      return 0;
    }
  }
  if(!niik_image_type_convert(cimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }
  if(!niik_image_type_convert(crefimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }

  if(refseg==NULL) { /* no refseg */
    if((crefseg=niik_image_copy_as_type(refimg,NIFTI_TYPE_UINT8))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
      return 0;
    }
    if(!niik_image_one(crefseg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_one\n",fcname);
      return 0;
    }
    free(crefseg->fname);
    crefseg->fname=(char *)calloc(10,sizeof(char));
    sprintf(crefseg->fname,"ones.nii");
  } /* no refseg */
  else if((crefseg=niik_image_affine_transform_3d(refseg,crefimg,NULL,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d\n",fcname);
    return 0;
  } /* with refseg */

  if(verbose>=1) fprintf(stdout,"[%s]   prepare %i %i\n",fcname,seed->row,seed->col);
  if(areg_cost == NIIK_REGISTER_NMI) {
    g_niik_aregister2_nmiobj=nmiobj;
  } else {
    g_niik_aregister2_nmiobj=NULL;
  }

  for(i=1,nvar=0; i<=16; i++)
    if(fabs(daffpar[i])>plim) nvar++;
  fprintf(stdout,"    DOF             %2i: ",nvar);
  if(fabs(daffpar[ 1])>plim) fprintf(stdout,"Rx ");
  if(fabs(daffpar[ 2])>plim) fprintf(stdout,"Ry ");
  if(fabs(daffpar[ 3])>plim) fprintf(stdout,"Rz ");
  if(fabs(daffpar[ 4])>plim) fprintf(stdout,"Tx ");
  if(fabs(daffpar[ 5])>plim) fprintf(stdout,"Ty ");
  if(fabs(daffpar[ 6])>plim) fprintf(stdout,"Tz ");
  if(fabs(daffpar[ 7])>plim) fprintf(stdout,"Sg ");
  if(fabs(daffpar[ 8])>plim) fprintf(stdout,"Sx ");
  if(fabs(daffpar[ 9])>plim) fprintf(stdout,"Sy ");
  if(fabs(daffpar[10])>plim) fprintf(stdout,"Sz ");
  if(fabs(daffpar[11])>plim) fprintf(stdout,"Kx ");
  if(fabs(daffpar[12])>plim) fprintf(stdout,"Ky ");
  if(fabs(daffpar[13])>plim) fprintf(stdout,"Kz ");
  fprintf(stdout,"\n");

  nseed = seed->row;
  for(i=0; i<nseed; i++) {
    atol = tol;
    cmaxiter = maxiter;
    if(!niik_image_aregister2(crefimg,crefseg,cimg,afmat,seed->m[i],daffpar,cmaxiter,atol,areg_cost)) {
      fprintf(stderr,"[%s] ERROR: niik_image_aregister2\n",fcname);
      return 0;
    }
    for(j=0; j<i; j++) {
      if(seed->m[i][0]<seed->m[j][0]) break;
    }
    if(j<i) {
      for(k=0; k<17; k++) {
        NIIK_DSWAP(&seed->m[i][k],&seed->m[j][k]);
      }
    }
    if(j==0) {
      fprintf(stdout,"[%s]   %3i %12.9f\n",fcname,i,seed->m[i][0]);
      niik_display_double_vector(seed->m[0]+1,16);
      niik_aregister_display_affine(seed->m[0]);
    }
    for(j=1; j<=3; j++) { /* correct rotations */
      seed->m[i][j]=(seed->m[i][j]> 180)?seed->m[i][j]-360:seed->m[i][j];
      seed->m[i][j]=(seed->m[i][j]<-180)?seed->m[i][j]+360:seed->m[i][j];
    }
  } /* each seed */

  /*fprintf(stdout,"free %i %0.4e\n",nseed[nl],tol->v[nl]);*/
  cimg    = niik_image_free(cimg);
  crefimg = niik_image_free(crefimg);
  crefseg = niik_image_free(crefseg);

  for(k=0; k<17; k++) {
    affpar[k] = seed->m[0][k];
  }

  g_niik_aregister2_nmiobj=NULL;

  if(verbose>=1) fprintf(stdout,"[%s] success\n",fcname);
  return 1;
} /* niik_image_aregister2_multiseed */


/******************** end of  registration functions *******************************/


int niik_image_update_sform(nifti_image *img,niikmat *afmni)
/* 2012-11-07, Kunio Nakamura
 * -updates sform in img with afmni
 */
{
  nifti_image *mni_img=NULL;
  char
  *FSLDIR=NULL,
   fname[4096],
   fcname[32]="niik_image_update_sform";
  niikmat *afmat=NULL;
  int m,n;
  fprintf(stdout,"[%s] reading mni img   %s\n",fcname,fname);
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[%s] FSLDIR is not defined\n",fcname);
    return 0;
  }
  sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
  if((mni_img=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
    exit(0);
  }
  afmat=niikmat_copy(afmni);
  niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(img->dx,img->dy,img->dz));
  niikmat_multiply_mat2_free1(niikmat_mat44_matrix(mni_img->sto_xyz),afmat);
  mni_img=niik_image_free(mni_img);
  img->sform_code=NIFTI_XFORM_MNI_152;
  for(m=0; m<4; m++) /* forward matrix */
    for(n=0; n<4; n++)
      img->sto_xyz.m[m][n]=afmat->m[m][n];
  if(!niikmat_inverse_update(afmat)) {
    fprintf(stderr,"[%s] ERROR: niikmat_inverse_update\n",fcname);
    return 0;
  }
  for(m=0; m<4; m++) /* inverse matrix */
    for(n=0; n<4; n++)
      img->sto_ijk.m[m][n]=afmat->m[m][n];
  afmat=niikmat_free(afmat);
  return 1;
} /* niik_image_update_sform */




/*************************************************************************
 *
 * IMAGE LEFT-RIGHT SYMMETRY OPTIMIZATION
 *
 *************************************************************************/

niikpt g_niik_image_aregister2_hemispheric_symmetry_optimization_center;
nifti_image *g_niik_image_aregister2_hemispheric_symmetry_optimization_image=NULL;

double niik_image_aregister2_hemispheric_symmetry_optimization_obj_func(double *v) {
  static int iter=0;
  niikmat *afmat=NULL;
  double
  err=0,
  phi,
  dphi = 0.05,
  rad,
  drad = 0.4,
  dist,
  dd = 1.0;
  niikpt
  center,
  pp[3];
  char fcname[128]="niik_image_aregister2_hemispheric_symmetry_optimization_obj_func";
  int
  num=0,
  verbose=0;
  /* initialization */
  if(v==NULL) {
    iter = 0;
    return 0;
  }
  iter++;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if((afmat = niikmat_init(4,4))==NULL) {
    fprintf(stderr,"[%s] ERROR: niikmat_init\n",fcname);
    fprintf(stderr,"[%s] exiting\n",fcname);
    exit(0);
  }
  if(verbose>2) {
    fprintf(stdout,"[%s] identity\n",fcname);
    niikmat_display(afmat);
  }
  niikmat_rotate_matrix_update(afmat,0,-v[0],-v[1]);
  if(verbose>2) {
    fprintf(stdout,"[%s] rotate\n",fcname);
    niikmat_display(afmat);
  }
  center = g_niik_image_aregister2_hemispheric_symmetry_optimization_center;
  niikmat_multiply_mat2_free1(niikmat_translate_matrix(center.x,center.y,center.z),
                              afmat);
  if(verbose>2) {
    fprintf(stdout,"[%s] translate\n",fcname);
    niikmat_display(afmat);
  }
  niikmat_multiply_mat2_free1(niikmat_scale_matrix   (1.0/g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dx,
                              1.0/g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dy,
                              1.0/g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dz),
                              afmat);
  if(verbose>2) {
    fprintf(stdout,"[%s] scale (pixel size)\n",fcname);
    niikmat_display(afmat);
  }
  if(verbose>1) {
    fprintf(stdout,"[%s] R = %7.3f %7.3f\n",fcname,v[0],v[1]);
    niikmat_display(afmat);
  }
  phi = 0;
  rad = 0;
  for(;;) {
    phi += dphi;
    rad += drad;
    pp[0].x =0;
    pp[0].y = cos(phi) * rad;
    pp[0].z = sin(phi) * rad;
    pp[1] = niikpt_affine_transform(afmat,pp[0]);
    if(verbose>2) {
      fprintf(stdout,"\npp0  %6.1f %6.1f %6.1f\n",pp[0].x,pp[0].y,pp[0].z);
      fprintf(stdout,"pp1  %6.1f %6.1f %6.1f\n",pp[1].x,pp[1].y,pp[1].z);
    }
    if(pp[1].x<0) break;
    if(pp[1].y<0) break;
    if(pp[1].z<0) break;
    if(pp[1].x>g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dim[1]) break;
    if(pp[1].y>g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dim[2]) break;
    if(pp[1].z>g_niik_image_aregister2_hemispheric_symmetry_optimization_image->dim[3]) break;
    for(dist=0; dist<100; dist+=dd) {
      /* Right direction */
      pp[0].x = dist;
      pp[1] = niikpt_affine_transform(afmat,pp[0]);
      if(verbose>3) fprintf(stdout,"  %6.1f %6.1f %6.1f\n",pp[1].x,pp[1].y,pp[1].z);
      if(pp[1].x<0) continue;
      if(pp[1].y<0) continue;
      if(pp[1].z<0) continue;
      if(pp[1].x >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->nx) continue;
      if(pp[1].y >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->ny) continue;
      if(pp[1].z >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->nz) continue;
      /* Left direction */
      pp[0].x = -dist;
      pp[2] = niikpt_affine_transform(afmat,pp[0]);
      if(verbose>2) fprintf(stdout,"  %6.1f %6.1f %6.1f\n",pp[2].x,pp[2].y,pp[2].z);
      if(pp[2].x<0) continue;
      if(pp[2].y<0) continue;
      if(pp[2].z<0) continue;
      if(pp[2].x >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->nx) continue;
      if(pp[2].y >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->ny) continue;
      if(pp[2].z >= g_niik_image_aregister2_hemispheric_symmetry_optimization_image->nz) continue;
      /* interpolate */
      pp[1].w = niik_image_interpolate_float_image_3d_linear_ijk(g_niik_image_aregister2_hemispheric_symmetry_optimization_image,pp[1]);
      pp[2].w = niik_image_interpolate_float_image_3d_linear_ijk(g_niik_image_aregister2_hemispheric_symmetry_optimization_image,pp[2]);
      if(verbose>3) {
        fprintf(stdout,"  %6.1f %6.1f %6.1f | %5.1f\n",pp[1].x,pp[1].y,pp[1].z,pp[1].w);
        fprintf(stdout,"  %6.1f %6.1f %6.1f | %5.1f\n",pp[2].x,pp[2].y,pp[2].z,pp[2].w);
      }
      /* absolute difference is the cost function */
      err+=fabs(pp[1].w - pp[2].w);
      num++;
    } /* symmetric interpolation about mid-plane */
    if(rad>250) break;
  } /* spiral loop on the mid-plane */
  if(verbose) fprintf(stdout,"  error = %9.3f   num = %8i   r[1] = %6.1f r[2] = %6.1f \n",
                        err/num,num,
                        v[0],v[1]);
  if(num==0) exit(0);
  niikmat_free(afmat);
  return err / num;
} /* niik_image_aregister2_hemispheric_symmetry_optimization_obj_func */

int niik_image_aregister2_hemispheric_symmetry_optiization(nifti_image *img,double FWHM,niikpt imgctr,double *ry, double *rz) {
  char fcname[128]="niik_image_aregister2_hemispheric_symmetry_optimization";
  int
  verbose=0,
  m,n,
  nm_maxiter,
  ndim;
  niikmat *p=NULL;
  double
  tol=1-4,
  (* pfn)(); /* symmetry-maximization function */
  if(verbose>=1) niik_fc_display(fcname,1);
  g_niik_image_aregister2_hemispheric_symmetry_optimization_center = imgctr;
  niik_image_aregister2_hemispheric_symmetry_optimization_obj_func (NULL);
  ndim=12;
  p=niikmat_init(ndim+1,ndim);
  for(m=1; m<ndim+1; m++) {
    for(n=0; n<ndim; n++) {
      p->m[m][n] = niik_get_rand() * 60.0 - 30.0;
    }
  }
  pfn=niik_image_aregister2_hemispheric_symmetry_optimization_obj_func;
  nm_maxiter = 700;
  if(FWHM<=0) {
    if(verbose>=1) fprintf(stdout,"[%s] no gaussian filter\n",fcname);
    if((g_niik_image_aregister2_hemispheric_symmetry_optimization_image=niik_image_copy(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
  } else {
    if(verbose>=1) fprintf(stdout,"[%s] apply gaussian filter %5.2f\n",fcname,FWHM);
    if((g_niik_image_aregister2_hemispheric_symmetry_optimization_image=niik_image_filter_gaussian(img,2.5*FWHM,FWHM))==NULL) {
      fprintf(stderr,"ERROR: niik_image_filter_gaussian\n");
      return 0;
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] optimization function\n",fcname);
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&nm_maxiter)) {
    fprintf(stderr,"ERROR: niik_nelder_mead\n");
    return 0;
  }
  if(verbose>=1)
    fprintf(stdout,"[%s] symmetry-maximization: r = %7.3f %7.3f %7.3f | iter %i | error %9.5f\n",fcname,0.0f,p->m[0][0],p->m[0][1],nm_maxiter,tol);
  *ry=p->m[0][0];
  *rz=p->m[0][1];
  niikmat_free(p);
  return 1;
} /* niik_image_aregister2_hemispheric_symmetry_optiization */


#endif

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/