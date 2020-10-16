/* Filename:     nifti1_kunio_off_cuberille.c
 * Description:  object function for cuberille
 * Author:       Kunio Nakamura
 * Date:         March 1, 2012
 */


#ifndef _FALCON_OFF_CUBERILLES_C_
#define _FALCON_OFF_CUBERILLES_C_

#include "falcon.h"
#include "falcon_surfaces.h"


int off_cuberille_find_largest_obj(kobj *obj);



/**********************************************
 * cuberille function
 *
 * -creates mesh from segmented image
 * -brute force approach
 *  where we get vertices and triangles first
 *  then clean up to remove redundant vertices
 * -program does not check for multiple objects or shared vertices
 *
 *
 **********************************************/



int niik_image_correct_for_cuberille(nifti_image *img,int flag)
/* prepare for cuberille functoin
 * -removes shared edges and vertices
 */
{
  const char *fcname="niik_image_correct_for_cuberille";
  int
  xdim,area,
       n,m,i,j,k,
       sum=1,sumsum=0,
       iter=0;
  unsigned char *bimg=NULL,*btmp=NULL;

  niik_fc_display(fcname,1);

  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }

  /* remove shared edges */
  fprintf(stdout,"[%s] preprocessing -removing shared edges / and vertices\n",fcname);
  if((btmp=niik_image_get_voxels_as_uint8_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector (btmp)\n");
    return 0;
  }
  if((bimg=niik_image_get_voxels_as_uint8_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector (btmp)\n");
    return 0;
  }

  xdim = img->dim[1];
  area = img->dim[1] * img->dim[2];

  if(flag==0) {
    sum=1;
    iter=0;
    while(sum) {
      iter++;
      sum=0;

      /* correcting foreground */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            m =
              btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
              btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
            if(0 && abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i\n",i,j,k,m);
            if(m!=2) continue;
            /*if(abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i %i %i %i   %i %i %i %i\n",i,j,k,
              btmp[n],btmp[n+1],btmp[n+xdim],btmp[n+xdim+1],
              btmp[n+area],btmp[n+1+area],btmp[n+xdim+area],btmp[n+xdim+1+area]);*/
            if(btmp[n]>0 && btmp[n+area+xdim+1]>0) {
              bimg[n     ]=0;
              sum++;
            }
            if(btmp[n+1]>0 && btmp[n+area+xdim]>0) {
              bimg[n+1   ]=0;
              sum++;
            }
            if(btmp[n+xdim]>0 && btmp[n+area+1]>0) {
              bimg[n+xdim]=0;
              sum++;
            }
            if(btmp[n+area]>0 && btmp[n+xdim+1]>0) {
              bimg[n+area]=0;
              sum++;
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %7i voxels (edges+fg vertices)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }

      /* correcting background */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            m =
              btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
              btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
            if(m!=6) continue;
            if(btmp[n]==0 && btmp[n+area+xdim+1]==0) {
              bimg[n+1   ]=0;
              sum++;
            }
            if(btmp[n+1]==0 && btmp[n+area+xdim]==0) {
              bimg[n     ]=0;
              sum++;
            }
            if(btmp[n+xdim]==0 && btmp[n+area+1]==0) {
              bimg[n+area]=0;
              sum++;
            }
            if(btmp[n+area]==0 && btmp[n+xdim+1]==0) {
              bimg[n+xdim]=0;
              sum++;
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %7i voxels (edges+fg/bg vertices)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }

      if(1) {
        for(k=0; k<img->nz; k++) {
          for(j=0; j<img->ny-1; j++) {
            for(i=0; i<img->nx-1; n++,i++) {
              n=i+j*xdim+k*area;
              /* x-dir */
              m = btmp[n] + btmp[n+xdim] + btmp[n+area] + btmp[n+area+xdim];
              if(m==2) {
                if(btmp[n]>0 && btmp[n+area+xdim]>0) {
                  bimg[n     ]=0;
                  sum++;
                }
                if(btmp[n+xdim]>0 && btmp[n+area]>0) {
                  bimg[n+xdim]=0;
                  sum++;
                }
              }
              /* y-dir */
              m = btmp[n] + btmp[n+1] + btmp[n+area] + btmp[n+area+1];
              if(m==2) {
                if(btmp[n]>0 && btmp[n+area+1]>0) {
                  bimg[n  ]=0;
                  sum++;
                }
                if(btmp[n+1]>0 && btmp[n+area]>0) {
                  bimg[n+1]=0;
                  sum++;
                }
              }
              /* z-dir */
              m = btmp[n] + btmp[n+1] + btmp[n+xdim] + btmp[n+xdim+1];
              if(m==2) {
                if(btmp[n]>0 && btmp[n+xdim+1]>0) {
                  bimg[n  ]=0;
                  sum++;
                }
                if(btmp[n+1]>0 && btmp[n+xdim]>0) {
                  bimg[n+1]=0;
                  sum++;
                }
              }
            }
          }
        }
        fprintf(stdout,"[%s]   modified %7i voxels (edges+fg/bg vertices,2d)\n",fcname,sum);
        for(i=0; i<img->nvox; i++) {
          btmp[i]=bimg[i];
        }
      }

      sumsum+=sum;
    } /* iteration */
    fprintf(stdout,"[%s] modified %7i voxels overall\n",fcname,sumsum);
  }

  else { /* ADD */
    sum=1;
    iter=0;
    while(sum) {
      iter++;
      sum=0;

      /* correcting foreground */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            m =
              btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
              btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
            /*if(abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i\n",i,j,k,m);*/
            if(m!=2) continue;
            /*if(abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i %i %i %i   %i %i %i %i\n",i,j,k,
              btmp[n],btmp[n+1],btmp[n+xdim],btmp[n+xdim+1],
              btmp[n+area],btmp[n+1+area],btmp[n+xdim+area],btmp[n+xdim+1+area]);*/
            if(btmp[n]>0 && btmp[n+area+xdim+1]>0) {
              bimg[n+1   ]=1;
              sum++;
            }
            if(btmp[n+1]>0 && btmp[n+area+xdim]>0) {
              bimg[n+xdim]=1;
              sum++;
            }
            if(btmp[n+xdim]>0 && btmp[n+area+1]>0) {
              bimg[n+area]=1;
              sum++;
            }
            if(btmp[n+area]>0 && btmp[n+xdim+1]>0) {
              bimg[n+xdim]=1;
              sum++;
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %i voxels (edges+fg vertices)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }

      /* correcting background */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            m =
              btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
              btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
            if(m!=6) continue;
            if(btmp[n]==0 && btmp[n+area+xdim+1]==0) {
              bimg[n     ]=1;
              sum++;
            }
            if(btmp[n+1]==0 && btmp[n+area+xdim]==0) {
              bimg[n+1   ]=1;
              sum++;
            }
            if(btmp[n+xdim]==0 && btmp[n+area+1]==0) {
              bimg[n+xdim]=1;
              sum++;
            }
            if(btmp[n+area]==0 && btmp[n+xdim+1]==0) {
              bimg[n+area]=1;
              sum++;
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %i voxels (edges+fg/bg vertices)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }

      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            /* x-dir */
            m = btmp[n] + btmp[n+xdim] + btmp[n+area] + btmp[n+area+xdim];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+area+xdim]>0) {
                bimg[n+xdim]=1;
                sum++;
              }
              if(btmp[n+xdim]>0 && btmp[n+area]>0) {
                bimg[n     ]=1;
                sum++;
              }
            }
            /* y-dir */
            m = btmp[n] + btmp[n+1] + btmp[n+area] + btmp[n+area+1];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+area+1]>0) {
                bimg[n+1]=1;
                sum++;
              }
              if(btmp[n+1]>0 && btmp[n+area]>0) {
                bimg[n  ]=1;
                sum++;
              }
            }
            /* z-dir */
            m = btmp[n] + btmp[n+1] + btmp[n+xdim] + btmp[n+xdim+1];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+xdim+1]>0) {
                bimg[n+1]=1;
                sum++;
              }
              if(btmp[n+1]>0 && btmp[n+xdim]>0) {
                bimg[n  ]=1;
                sum++;
              }
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %i voxels (edges+fg/bg vertices,2d)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }

      sumsum+=sum;
    } /* iteration */
    fprintf(stdout,"[%s] modified %i voxels overall\n",fcname,sumsum);
  } /* adding */

  if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_set_voxels_from_uint8_vector\n",fcname);
    return 0;
  }

  free(btmp);
  free(bimg);
  niik_fc_display(fcname,0);
  return 1;
}


kobj *off_cuberille_kobj(nifti_image *img,int flag_all_obj)
/* -img is a binary mask
 * -flag_all_obj is a flag to include all objects (nonzero) or
 *  only return the largest object
 */

{
  kobj *obj;
  kvert *v,*nv[4];
  kface *f;
  kedge *e,*et=NULL;
  char fname[512];
  const char *fcname="off_cuberille_kobj";
  unsigned char *bimg,*btmp;
  int
  iter=0,
  area,xdim,
  sum,sumsum=0,
      nn,m,n,i,j,k;
  niikvec *vec;

  niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if((bimg=niik_image_get_voxels_as_uint8_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector (bimg)\n");
    return NULL;
  }
  if((btmp=niik_image_get_voxels_as_uint8_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector (btmp)\n");
    return NULL;
  }

  obj=off_obj_init();
  obj->vert = v = off_vert_init();
  v->v=niikpt_problem();

  obj->face = f = off_face_init();

  xdim = img->dim[1];
  area = img->dim[1] * img->dim[2];


  /* remove shared edges */

  fprintf(stdout,"[%s] preprocessing -removing shared edges / and vertices\n",fcname);

  sum=1;
  iter=0;
  while(sum) {
    iter++;
    sum=0;
    if(0) {
      for(k=n=0; k<img->dim[3]; k++) {
        for(j=0; j<img->dim[2]; j++) {
          for(i=0; i<img->dim[1]; n++,i++) {
            if(btmp[n]==0) continue;
            if(i==0) continue;
            if(j==0) continue;
            if(k==0) continue;
            if(i==img->dim[1]-1) continue;
            if(j==img->dim[2]-1) continue;
            if(k==img->dim[3]-1) continue;

            if(btmp[n-area-xdim]>0)
              if(btmp[n-area]==0 &&
                  btmp[n-xdim]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k-1); */
                bimg[n-area]=1;
              }

            if(btmp[n-area-1]>0)
              if(btmp[n-area]==0 &&
                  btmp[n-   1]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k-1); */
                bimg[n-area]=1;
              }

            if(btmp[n-area+1]>0)
              if(btmp[n-area]==0 &&
                  btmp[n   +1]==0 ) {
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k-1); */
                bimg[n-area]=1;
              }

            if(btmp[n-area+xdim]>0)
              if(btmp[n-area]==0 &&
                  btmp[n+xdim]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k-1); */
                bimg[n-area]=1;
              }

            if(btmp[n-xdim-1]>0)
              if(btmp[n-xdim]==0 &&
                  btmp[n   -1]==0 ) {
                sum++;
                /* 	    fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j-1,k); */
                bimg[n-xdim]=1;
              }

            if(btmp[n-xdim+1]>0)
              if(btmp[n-xdim]==0 &&
                  btmp[n   +1]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j-1,k); */
                bimg[n-xdim]=1;
              }

            if(btmp[n+xdim+1]>0)
              if(btmp[n+xdim]==0 &&
                  btmp[n   +1]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j+1,k); */
                bimg[n+xdim]=1;
              }

            if(btmp[n+xdim-1]>0)
              if(btmp[n+xdim]==0 &&
                  btmp[n   -1]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j+1,k); */
                bimg[n+xdim]=1;
              }

            if(btmp[n+area-xdim]>0)
              if(btmp[n+area]==0 &&
                  btmp[n-xdim]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k+1); */
                bimg[n+area]=1;
              }

            if(btmp[n+area+xdim]>0)
              if(btmp[n+area]==0 &&
                  btmp[n+xdim]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k+1); */
                bimg[n+area]=1;
              }

            if(btmp[n+area-1]>0)
              if(btmp[n+area]==0 &&
                  btmp[n-1]==0 ) {
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k+1); */
                bimg[n+area]=1;
              }

            if(btmp[n+area+1]>0)
              if(btmp[n+area]==0 &&
                  btmp[n+1]==0 ) {
                sum++;
                /* fprintf(stdout,"      adding [%3i %3i %3i]\n",i,j,k+1); */
                bimg[n+area]=1;
              }

          }
        }
      }
      fprintf(stdout,"[%s]   modified %i voxels (edges)\n",fcname,sum);
    }

    /* correcting foreground */
    for(k=0; k<img->nz; k++) {
      for(j=0; j<img->ny-1; j++) {
        for(i=0; i<img->nx-1; n++,i++) {
          n=i+j*xdim+k*area;
          m =
            btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
            btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
          if(0 && abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i\n",i,j,k,m);
          if(m!=2) continue;
          if(0 && abs(i-70)<=2 && abs(j-49)<=2 && abs(k-61)<=2) fprintf(stdout,"[%3i %3i %3i] %i %i %i %i   %i %i %i %i\n",i,j,k,
                btmp[n],btmp[n+1],btmp[n+xdim],btmp[n+xdim+1],
                btmp[n+area],btmp[n+1+area],btmp[n+xdim+area],btmp[n+xdim+1+area]);
          if(btmp[n]>0 && btmp[n+area+xdim+1]>0) {
            bimg[n     ]=0;
            sum++;
          }
          if(btmp[n+1]>0 && btmp[n+area+xdim]>0) {
            bimg[n+1   ]=0;
            sum++;
          }
          if(btmp[n+xdim]>0 && btmp[n+area+1]>0) {
            bimg[n+xdim]=0;
            sum++;
          }
          if(btmp[n+area]>0 && btmp[n+xdim+1]>0) {
            bimg[n+area]=0;
            sum++;
          }
        }
      }
    }
    fprintf(stdout,"[%s]   modified %i voxels (edges+fg vertices)\n",fcname,sum);
    for(i=0; i<img->nvox; i++) {
      btmp[i]=bimg[i];
    }

    /* correcting background */
    for(k=0; k<img->nz; k++) {
      for(j=0; j<img->ny-1; j++) {
        for(i=0; i<img->nx-1; n++,i++) {
          n=i+j*xdim+k*area;
          m =
            btmp[n     ] + btmp[n+1     ] + btmp[n+xdim     ] + btmp[n+xdim+1     ] +
            btmp[n+area] + btmp[n+1+area] + btmp[n+xdim+area] + btmp[n+xdim+1+area] ;
          if(m!=6) continue;
          if(btmp[n]==0 && btmp[n+area+xdim+1]==0) {
            bimg[n+1   ]=0;
            sum++;
          }
          if(btmp[n+1]==0 && btmp[n+area+xdim]==0) {
            bimg[n     ]=0;
            sum++;
          }
          if(btmp[n+xdim]==0 && btmp[n+area+1]==0) {
            bimg[n+area]=0;
            sum++;
          }
          if(btmp[n+area]==0 && btmp[n+xdim+1]==0) {
            bimg[n+xdim]=0;
            sum++;
          }
        }
      }
    }
    fprintf(stdout,"[%s]   modified %i voxels (edges+fg/bg vertices)\n",fcname,sum);
    for(i=0; i<img->nvox; i++) {
      btmp[i]=bimg[i];
    }

    if(1) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny-1; j++) {
          for(i=0; i<img->nx-1; n++,i++) {
            n=i+j*xdim+k*area;
            /* x-dir */
            m = btmp[n] + btmp[n+xdim] + btmp[n+area] + btmp[n+area+xdim];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+area+xdim]>0) {
                bimg[n     ]=0;
                sum++;
              }
              if(btmp[n+xdim]>0 && btmp[n+area]>0) {
                bimg[n+xdim]=0;
                sum++;
              }
            }
            /* y-dir */
            m = btmp[n] + btmp[n+1] + btmp[n+area] + btmp[n+area+1];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+area+1]>0) {
                bimg[n  ]=0;
                sum++;
              }
              if(btmp[n+1]>0 && btmp[n+area]>0) {
                bimg[n+1]=0;
                sum++;
              }
            }
            /* z-dir */
            m = btmp[n] + btmp[n+1] + btmp[n+xdim] + btmp[n+xdim+1];
            if(m==2) {
              if(btmp[n]>0 && btmp[n+xdim+1]>0) {
                bimg[n  ]=0;
                sum++;
              }
              if(btmp[n+1]>0 && btmp[n+xdim]>0) {
                bimg[n+1]=0;
                sum++;
              }
            }
          }
        }
      }
      fprintf(stdout,"[%s]   modified %i voxels (edges+fg/bg vertices,2d)\n",fcname,sum);
      for(i=0; i<img->nvox; i++) {
        btmp[i]=bimg[i];
      }
      if(0) {
        niik_image_set_voxels_from_uint8_vector(img,bimg);
        sprintf(fname,"tmp_img-iter%i.nii.gz",iter);
        fprintf(stdout,"\twriting image:    %s\n",fname);
        niik_image_write(fname,img);
      }
    }

    sumsum+=sum;
  } /* iteration */

  fprintf(stdout,"[%s] modified %i voxels overall\n",fcname,sumsum);

  niik_image_set_voxels_from_uint8_vector(img,bimg);
  if(0) {
    sprintf(fname,"tmp_img2.nii.gz");
    fprintf(stdout,"\twriting image:    %s\n",fname);
    niik_image_write(fname,img);
    /*exit(0);*/
  } /* writing image */

  free(btmp);





  fprintf(stdout,"[%s] create vertices and faces\n",fcname);

  for(k=n=0; k<img->dim[3]; k++) {
    /*    for(tv=obj->vert,i=0;tv!=NULL;tv=tv->next){ i++; }
    fprintf(stdout,"  k = %i  nvert = %i\n",k,i);*/

    for(j=0; j<img->dim[2]; j++) {
      for(i=0; i<img->dim[1]; n++,i++) {
        if(bimg[n]==0) continue;

        /* LEFT */
        nn=0;
        if(i==0) nn=1;
        else if(bimg[n-1]==0) nn=1;
        if(nn) {
          /* add 4 new vertices */
          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i,j,k))==NULL) {
            nv[0] = v->next = off_vert_init_with_v       (i,j,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k))==NULL) {
            nv[1] = v->next = off_vert_init_with_v       (i,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k+1))==NULL) {
            nv[2] = v->next = off_vert_init_with_v       (i,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i,j,k+1))==NULL) {
            nv[3] = v->next = off_vert_init_with_v       (i,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on left side */
          /* fprintf(stdout,"  LEFT [%3i %3i %3i] FACE \n",i,j,k);*/
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[1];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[3];
          f -> vert[2] = nv[2];

        } /* LEFT */

        /* RIGHT */
        nn=0;
        if(i==img->dim[1]-1) nn=1;
        else if(bimg[n+1]==0) nn=1;
        if(nn) {
          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k))==NULL) {
            nv[0] = v->next = off_vert_init_with_v      ( i+1,j,k );
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k))==NULL) {
            nv[1] = v->next = off_vert_init_with_v      ( i+1,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k+1))==NULL) {
            nv[2] = v->next = off_vert_init_with_v  (     i+1,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k+1))==NULL) {
            nv[3] = v->next = off_vert_init_with_v     (  i+1,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on right side */
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[1];
          f -> vert[2] = nv[2];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[3];

        } /* RIGHT */


        /* ANTERIOR */
        nn=0;
        if(j==0) nn=1;
        else if(bimg[n-img->dim[1]]==0) nn=1;
        if(nn) {
          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i,j,k))==NULL) {
            nv[0] = v->next = off_vert_init_with_v       (i,j,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i,j,k+1))==NULL) {
            nv[1] = v->next = off_vert_init_with_v       (i,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k+1))==NULL) {
            nv[2] = v->next = off_vert_init_with_v       (i+1,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k))==NULL) {
            nv[3] = v->next = off_vert_init_with_v       (i+1,j,k);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on anterior side */
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[1];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[3];
          f -> vert[2] = nv[2];

        } /* ANTERIOR */


        /* POSTERIOR */
        nn=0;
        if(j==img->dim[2]-1) nn=1;
        else if(bimg[n+img->dim[1]]==0) nn=1;
        if(nn) {
          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k))==NULL) {
            nv[0] = v->next = off_vert_init_with_v       (i,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k+1))==NULL) {
            nv[1] = v->next = off_vert_init_with_v       (i,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k+1))==NULL) {
            nv[2] = v->next = off_vert_init_with_v       (i+1,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k))==NULL) {
            nv[3] = v->next = off_vert_init_with_v       (i+1,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on posterior side */
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[1];
          f -> vert[2] = nv[2];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[3];

        } /* POSTERIOR */


        /* UP */
        nn=0;
        if(k==0) nn=1;
        else if(bimg[n-area]==0) nn=1;
        if(nn) {

          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i,j,k))==NULL) {
            nv[0] = v->next = off_vert_init_with_v       (i,j,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k))==NULL) {
            nv[1] = v->next = off_vert_init_with_v       (i+1,j,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k))==NULL) {
            nv[2] = v->next = off_vert_init_with_v       (i+1,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k))==NULL) {
            nv[3] = v->next = off_vert_init_with_v       (i,j+1,k);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on up side */
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[1];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[3];
          f -> vert[2] = nv[2];

        } /* UP */

        /* DOWN */
        nn=0;
        if(k==img->dim[3]) nn=1;
        else if(bimg[n+area]==0) nn=1;
        if(nn) {
          if((nv[0]=off_cuberille_kobj_search_kvert(obj,v,i,j,k+1))==NULL) {
            nv[0] = v->next = off_vert_init_with_v       (i,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[1]=off_cuberille_kobj_search_kvert(obj,v,i+1,j,k+1))==NULL) {
            nv[1] = v->next = off_vert_init_with_v       (i+1,j,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[2]=off_cuberille_kobj_search_kvert(obj,v,i+1,j+1,k+1))==NULL) {
            nv[2] = v->next = off_vert_init_with_v       (i+1,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          if((nv[3]=off_cuberille_kobj_search_kvert(obj,v,i,j+1,k+1))==NULL) {
            nv[3] = v->next = off_vert_init_with_v       (i,j+1,k+1);
            v->next->prev = v;
            v=v->next;
          }

          /* add 2 new faces on down side */
          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[1];
          f -> vert[2] = nv[2];

          f -> next = off_face_init();
          f -> next -> prev = f;
          f = f -> next;
          f -> vert[0] = nv[0];
          f -> vert[1] = nv[2];
          f -> vert[2] = nv[3];

        } /* DOWN */

      }
    }
  }  /* ijk to get the vertex and surface */

  free(bimg);


  /**********************************************************
   *
   * clean up
   *
   * -remove the extra vert and face at the top of linked list
   *
   *
   ***********************************************************/

  fprintf(stdout,"[%s] clean up the empty top of the double linked list\n",fcname);
  v = obj->vert;
  obj->vert = obj->vert -> next;
  off_kvert_free(v);
  obj->vert->prev = NULL;

  f = obj->face;
  obj->face = obj->face -> next;
  off_kface_free(f);
  obj->face->prev = NULL;

  for(v=obj->vert,n=1; v!=NULL; v=v->next) {
    v->index = n++;
  }
  obj->nvert = n-1;

  for(f=obj->face,n=1; f!=NULL; f=f->next) {
    f->index = n++;
  }
  obj->nface = n-1;
  fprintf(stdout,"[%s] vfe %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  /* convert to world coordinates */
  fprintf(stdout,"[%s] voxel to world conversion\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) {
    niikpt p = v->v;
    niik_ijk_to_world(img, &p, &v->v);
  }

  /* make edges */
  fprintf(stdout,"[%s] create edges\n",fcname);
  if(!off_kobj_read_off_make_edge(obj)) {
    fprintf(stderr,"ERROR: off_kobj_read_off_make_edge\n");
    exit(0);
  }
  vec=niikvec_init(obj->nedge);
  for(e=obj->edge; e!=NULL; e=e->next) vec->v[e->index-1]=niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
  fprintf(stdout,"[%s]   elen %9.5f  %9.5f    %9.5f %9.5f\n",fcname,
          niik_get_mean_from_double_vector(vec->v,vec->num),
          niik_get_stdv_from_double_vector(vec->v,vec->num),
          niik_get_min_from_double_vector(vec->v,vec->num),
          niik_get_max_from_double_vector(vec->v,vec->num));
  fprintf(stdout,"[%s] vfe %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  /* very exhausive edge checking, do not run it unless object is small */
  if(0) {
    fprintf(stdout,"     edge checking\n");
    for(e=obj->edge; e!=NULL; e=e->next) {
      for(et=e->next; et!=NULL; et=et->next) {
        if(e->endpts[0] == et->endpts[0])
          if(e->endpts[1] == et->endpts[1]) {
            fprintf(stderr,"ERROR: same endpts on different edge \n");
            off_display_edge_info(e,10);
            off_display_edge_info(et,10);
            exit(0);
          }
        if(e->endpts[1] == et->endpts[0])
          if(e->endpts[0] == et->endpts[1]) {
            fprintf(stderr,"ERROR: same endpts on different edge \n");
            off_display_edge_info(e,10);
            off_display_edge_info(et,10);
            exit(0);
          }
      }
    }
  }


  /* update neighbors */
  fprintf(stdout,"[%s] update neighbors\n",fcname);
  off_kobj_update_num(obj);
  fprintf(stdout,"     make neighbors\n");
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->nei=0;
  }
  for(e=obj->edge; e!=NULL; e=e->next) {
    if(!e->endpts[0]->nei)
      off_kobj_update_vertex_nei(e->endpts[0],e);
    if(!e->endpts[1]->nei)
      off_kobj_update_vertex_nei(e->endpts[1],e);
  }


  /* 2012-08-23 Kunio
   * check if a vertex is used on the other side of surface */
  fprintf(stdout,"[%s] check for shared vertices, create duplicates\n",fcname);
  for(e=obj->edge; e!=NULL && 0; e=e->next) {
    if((n=off_kobj_test_common_neighbor(e))!=2) {
      fprintf(stdout,"\n\n[off_kobj_cuberille] need to correct here #hit=%i\n",n);
      off_display_edge_info(e,20);
      off_display_vert_info(e->endpts[0],20);
      off_display_vert_info(e->endpts[1],20);
      off_display_face_info(e->adjface[0],20);
      off_display_face_info(e->adjface[1],20);
      exit(0);
    }
  }

  off_kobj_update_num(obj);

  /* TESTING */
  if(!off_kobj_test(obj)) {
    fprintf(stderr,"ERROR: off_kobj_test\n");
    return NULL;
  }
  fprintf(stdout,"[%s] passed off_kobj_test\n",fcname);

  if(!flag_all_obj) {
    /* FIND THE LARGEST OBJECT */
    if(!off_cuberille_find_largest_obj(obj)) {
      fprintf(stderr,"ERROR: off_cuberille_find_largest_obj\n");
      return NULL;
    }
  } /* find the largest object */


  niik_fc_display(fcname,0);
  return obj;
}

kvert *off_cuberille_kobj_search_kvert(kobj *obj,kvert *vtail,int x,int y, int z)
/* fast search of a vertex to reduce redundancy
 * -it becames faster when we searched from the last vertex in the double linked list */
{
  kvert *v;
  for(v=vtail; v->prev!=NULL; v=v->prev) {
    if(v==NULL) break;
    if((int)v->v.z < z-5) return NULL;
    if((int)v->v.x == x)
      if((int)v->v.y == y)
        if((int)v->v.z == z)
          return v;
  }
  return NULL;
  for(v=obj->vert->next; v!=NULL; v=v->next) {
    if(v==NULL) break;
    if((int)v->v.x == x)
      if((int)v->v.y == y)
        if((int)v->v.z == z)
          return v;
  }
  /* fprintf(stdout,"  no vertex \n"); */
  return NULL;
} /* off_cuberlile_kobj_search_kvert */


/**********************************************************************************
 *
 * largest connected object
 *
 **********************************************************************************/

int off_cuberille_find_largest_obj_mark_neighbors(kvert *vert,double w) 
{
  
  kvert ** working_set = malloc(sizeof(kvert *));
  kvert ** next_working_set = NULL;
  int working_set_size=1;
  working_set[0] = vert;

  for(;;)
  {
    kvert **tmp;
    int next_working_set_size=0;
    int n;

    // find neighbors connected and not labelled
    for(n=0; n<working_set_size; n++) {
      int i;
      // maximum new entries
      next_working_set=realloc(next_working_set, (next_working_set_size + working_set[n]->nei)*sizeof(kvert *));
      for(i=0;i<working_set[n]->nei;++i)
      {
        if( (int)working_set[n]->neivert[i]->v.w == 0) //unassigned
        {
          working_set[n]->neivert[i]->v.w = w;
          next_working_set[next_working_set_size] = working_set[n]->neivert[i];
          next_working_set_size++; 
        }
      }
    }
    if(!next_working_set_size) break; // no more vertices to process
    //swap buffers
    tmp = working_set;
    working_set = next_working_set;
    working_set_size = next_working_set_size;
    next_working_set = tmp;
  }

  free(working_set);
  free(next_working_set);
  return 1;
}


int off_cuberille_find_largest_obj(kobj *obj)
/* -finds the largest object and
 *  removes points and edges of small objects
 * -uses v.w as the marker so this will be altered.
 */
{
  kvert *v,**vlist;
  kedge *e,**elist;
  kface *f,**flist;
  int
  verbose=1,
  m,n,w=1,fw=0,wsum=0,wmax=0;
  const char *fcname="off_cuberille_find_largest_obj";

  if(verbose>=1) niik_fc_display(fcname,1);
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }

  off_kobj_update_num(obj);
  fprintf(stdout,"[%s] vfe before %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  /* initialize */
  if(verbose>=1) fprintf(stdout,"[%s] initialization v.w\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) v->v.w=0;

  /* repeat... */
  if(verbose>=1) fprintf(stdout,"[%s] loop until all vertices are labeld\n",fcname);
  for(;;) {
    for(v=obj->vert; v!=NULL; v=v->next) {
      if((int)v->v.w==0) break;
    }
    if(v==NULL) break;
    v->v.w=w;
    off_cuberille_find_largest_obj_mark_neighbors(v,w);
    for(v=obj->vert,wsum=0; v!=NULL; v=v->next) {
      if((int)v->v.w==w)
        wsum++;
    }
    fprintf(stdout,"[%s]   %i    # %i\n",fcname,(int)w,(int)wsum);
    if(!wsum) break;
    if(wsum>wmax) {
      wmax=wsum;
      fw=w;
    }
    if(wmax > obj->nvert*0.51) break; //VF?
    w++;
  }

  if(verbose>=1) fprintf(stdout,"[%s] clean up edges\n",fcname);
  for(e=obj->edge,n=0; e!=NULL; e=e->next) {
    if((int)e->endpts[0]->v.w!=fw) {
      n++;
    }
  }
  elist=(kedge **)calloc(n,sizeof(kedge *));
  for(e=obj->edge,n=0; e!=NULL; e=e->next) {
    if((int)e->endpts[0]->v.w != fw) {
      elist[n++]=e;
    }
  }
  fprintf(stdout,"[%s] removing %d edges\n",fcname,n);
  for(m=0; m<n; m++) {
    /*fprintf(stdout,"[%s]   removing edge %i [%i]\n",fcname, elist[m]->index,m);*/
    off_kedge_remove(elist[m],obj);
  }
  free(elist);

  if(verbose>=1) fprintf(stdout,"[%s] clean up faces\n",fcname);
  for(f=obj->face,n=0; f!=NULL; f=f->next) {
    if((int)f->vert[0]->v.w != fw) {
      n++;
    }
  }
  flist=(kface **)calloc(n,sizeof(kface *));
  for(f=obj->face,n=0; f!=NULL; f=f->next) {
    if((int)f->vert[0]->v.w!=fw) {
      flist[n++]=f;
    }
  }
  fprintf(stdout,"[%s] removing %d faces\n",fcname,n);
  for(m=0; m<n; m++) {
    /*fprintf(stdout,"[%s]   removing edge %i [%i]\n",fcname,elist[m]->index,m);*/
    off_kface_remove(flist[m],obj);
  }
  free(flist);

  if(verbose>=1) fprintf(stdout,"[%s] clean up verts\n",fcname);
  for(v=obj->vert,n=0; v!=NULL; v=v->next) {
    if((int)v->v.w != fw) {
      n++;
    }
  }
  vlist=(kvert **)calloc(n,sizeof(kvert *));
  for(v=obj->vert,n=0; v!=NULL; v=v->next) {
    if((int)v->v.w != fw) {
      vlist[n++]=v;
    }
  }
  fprintf(stdout,"[%s] removing %d verts\n", fcname, n);

  for(m=0; m<n; m++) {
    /*fprintf(stdout,"[%s]   removing edge %i [%i]\n",fcname,elist[m]->index,m);*/
    off_kvert_remove(vlist[m],obj);
  }
  free(vlist);

  off_kobj_update_all_index(obj);
  fprintf(stdout,"[%s] vfe after %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}

#endif /* _FALCON_OFF_CUBERILLES_C_ */
