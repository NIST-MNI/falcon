/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)


int niik_image_patch_test1(nifti_image *img,nifti_image *lib,int i,int j,int k,int len,double *out,int method);


static char *prog_version[] = {
  "  program history\n"
  "  0.0.0, October 4, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
};

static char *prog_describe[] = {
  "  Description: None\n"
};

static char *prog_help[] = {
  "  program help:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"test program\n");
  fprintf(stdout,"  usage: [options] <img.mnc> <mask.mnc> <out.mnc> whatever\n\n");
}


int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *lib=NULL,
   *libg=NULL,
    *libw=NULL,
     *outimg=NULL,
      *mask=NULL,
       *img=NULL;

  int
  verbose=1,
  num_req_arg=4,
  ijk[4],test_ijk[4],
  n,nc,sc;
  char
  fcname[32]="niik_test";
  double
  vmedian,
  v;
  niikmat
  *m_classify=NULL,
   *m_feature=NULL;
  char
  fname[4096],
        libdir[4096],
        **csstr;
  int ncsstr=-1;

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;
  test_ijk[0]=0;

  sprintf(libdir,"/data/kproj/data/classify/bw/");

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--describe",10)) {
        fprintf(stdout,"%s",*prog_describe);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-libdir=",8)) {
        sprintf(libdir,"%s",argv[nc]+8);
      }

      else if(!strncmp(argv[nc],"-ijk=",5)) {
        csstr=niik_csstring_to_list(argv[nc]+5,&ncsstr);
        for(n=0; n<ncsstr; n++) {
          test_ijk[1+n]=atoi(csstr[n]);
          free(csstr[n]);
        }
        test_ijk[0]=1;
        fprintf(stdout,"[%s] testing voxel %i %i %i\n",fcname,test_ijk[1],test_ijk[2],test_ijk[3]);
        verbose=2;
        free(csstr);
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

  if(argc>num_req_arg) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(9);
  } else if(argc<num_req_arg) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(9);
  }

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"reading img, niik_image_read",9);

  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[2]);
  NIIK_EXIT(((mask=niik_image_read(argv[2]))==NULL),fcname,"reading mask, niik_image_read",9);

  sprintf(fname,"%s/lib.nii.gz",libdir);
  fprintf(stdout,"[%s] reading lib:          %s\n",fcname,fname);
  NIIK_EXIT(((lib=niik_image_read(fname))==NULL),fcname,"reading library, niik_image_read",9);

  sprintf(fname,"%s/libg.nii.gz",libdir);
  fprintf(stdout,"[%s] reading lib g:        %s\n",fcname,fname);
  NIIK_EXIT(((libg=niik_image_read(fname))==NULL),fcname,"reading library g, niik_image_read",9);

  sprintf(fname,"%s/libw.nii.gz",libdir);
  fprintf(stdout,"[%s] reading lib w:        %s\n",fcname,fname);
  NIIK_EXIT(((libw=niik_image_read(fname))==NULL),fcname,"reading library w, niik_image_read",9);

  NIIK_RET0(((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_UINT8))==NULL),fcname,"niik_image_copy_as_type");
  niik_image_clear(outimg);

  /* feature classification matrix */
  m_classify=niikmat_init(1,4);
  m_classify->m[0][0]=-0.001;
  m_classify->m[0][1]=-0.005/ 27.0;
  m_classify->m[0][2]=-0.010/125.0;
  m_classify->m[0][3]=1.000;

  if(test_ijk[0]>0)
    test_ijk[0]=test_ijk[1]+test_ijk[2]*img->nx+test_ijk[3]*img->nx*img->ny;

  fprintf(stdout,"[%s] #templates = %i\n",fcname,lib->nt);

  for(ijk[3]=0; ijk[3]<img->nz; ijk[3]++) {
    printf("[%d%%   %i of %i]",(int)(100.0*ijk[3]/img->nz),ijk[3]+1,img->nz);
    printf("\r");
    fflush(stdout);

    for(ijk[2]=0; ijk[2]<img->ny; ijk[2]++) {
      for(ijk[1]=0; ijk[1]<img->nx; ijk[0]++,ijk[1]++) {
        ijk[0]=ijk[1]+ijk[2]*img->nx+ijk[3]*img->nx*img->ny;

        if(niik_image_get_voxel(mask,ijk[0])<0.5) continue;
        if(test_ijk[0]>0) {
          if(test_ijk[0]!=ijk[0]) continue;
          else fprintf(stdout,"[%s] testing...\n",fcname);
        }

        m_feature=niikmat_init(4,lib->nt);

        niik_image_patch_test1(img,lib,ijk[1],ijk[2],ijk[3],0,m_feature->m[0],0);
        niik_image_patch_test1(img,lib,ijk[1],ijk[2],ijk[3],1,m_feature->m[1],1);
        niik_image_patch_test1(img,lib,ijk[1],ijk[2],ijk[3],2,m_feature->m[2],1);
        niik_image_patch_test1(img,lib,ijk[1],ijk[2],ijk[3],3,m_feature->m[3],2);
        if(verbose>1) {
          fprintf(stdout,"[%s] features:\n",fcname);
          niikmat_display(m_feature);
          fprintf(stdout,"[%s] classifier:\n",fcname);
          niikmat_display(m_classify);
        }

        niikmat_multiply_mat2(m_classify,m_feature);
        if(verbose>1) {
          fprintf(stdout,"[%s] multiplied:\n",fcname);
          niikmat_display(m_feature);
        }

        /* filtering (only use the good half) */
        vmedian = niik_get_median_from_sorted_double_vector(m_feature->m[0],m_feature->col);
        for(n=0; n<lib->nt; n++) {
          if(m_feature->m[0][n]<vmedian) m_feature->m[0][n]=0;
          else m_feature->m[0][n]=pow(m_feature->m[0][n],5);
        }
        if(verbose>1) {
          fprintf(stdout,"[%s] modulated:\n",fcname);
          niikmat_display(m_feature);
        }

        /* normalize */
        for(v=0,n=0; n<lib->nt; n++) {
          v+=m_feature->m[0][n];
        }
        for(n=0; n<lib->nt; n++) {
          m_feature->m[0][n]/=v;
        }

        if(verbose>1) {
          fprintf(stdout,"[%s] normalized probability:\n",fcname);
          niikmat_display(m_feature);
        }

        for(v=0,n=0; n<lib->nt; n++) {
          /*fprintf(stdout,"[%3i,%3i,%3i] %12.6f %12.6f\n",ijk[1],ijk[2],ijk[3],m_feature->m[0][n],niik_image_get_voxel(libg,ijk[0]+n*img->ny*img->nx*img->nz));*/
          v += m_feature->m[0][n] * niik_image_get_voxel(libg,ijk[0]+n*img->ny*img->nx*img->nz);
        }

        if(verbose>1) {
          fprintf(stdout,"[%3i,%3i,%3i]",ijk[1],ijk[2],ijk[3]);
          for(n=0; n<lib->nt; n++) {
            fprintf(stdout,"%15.9f ",niik_image_get_voxel(libg,ijk[0]+n*img->ny*img->nx*img->nz));
          }
          fprintf(stdout,"\n[%3i,%3i,%3i] %12.6f\n",ijk[1],ijk[2],ijk[3],v);
        } /* verbose */

        niik_image_set_voxel(outimg,ijk[0],v);

        m_feature=niikmat_free(m_feature);

      }
    }
  } /* each voxel */

  fprintf(stdout,"[%s] writing output img:   %s\n",fcname,argv[3]);
  NIIK_EXIT((!niik_image_write(argv[3],outimg)),fcname,"niik_image_write",9);

  m_feature=niikmat_free(m_feature);
  m_classify=niikmat_free(m_classify);
  img=niik_image_free(img);
  lib=niik_image_free(lib);
  libg=niik_image_free(libg);
  libw=niik_image_free(libw);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_test1 */



int niik_image_patch_test1(nifti_image *img,nifti_image *lib,int i,int j,int k,int len,double *out,int method)
/* -lib has data in t-direction
 * -len is the 3d distance
 * -method: 1=absdiff
 */
{
  int
  verbose=0,
  area,size=0;
  int
  ilo,jlo,klo,ihi,jhi,khi,
      m,n,
      ii,jj,kk,nn;
  char fcname[32]="niik_image_patch_test1";
  double x,y,sxy=0,sx=0,sy=0,sxx=0,syy=0;
  area=img->nx*img->ny;
  size=area*img->nz;
  for(n=0; n<lib->nt; n++) {
    switch(method) {
    default:
      fprintf(stderr,"[%s] ERROR: unknown method, %i\n",fcname,method);
      return 0;
    case 0: /* intensity difference */
      out[n]+=fabs(niik_image_get_voxel(img,i+j*img->nx+k*area) -
                   niik_image_get_voxel(lib,i+j*img->nx+k*area+n*size) );
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]);
      break;
    case 1: /* absolute difference */
      out[n]=0;
      for(kk=k-len; kk<=k+len; kk++) {
        for(jj=j-len; jj<=j+len; jj++) {
          for(ii=i-len,m=kk*area+jj*img->nx+ii; ii<=i+len; m++,ii++) {
            /* fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]); */
            out[n]+=fabs(niik_image_get_voxel(img,m) -
                         niik_image_get_voxel(lib,m+n*size) );
          }
        }
      }
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f\n",i,j,k,out[n]);
      break;
    case 2: /* correlation */
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
      for(kk=klo; kk<=khi; kk++) {
        for(jj=jlo; jj<=jhi; jj++) {
          for(ii=ilo,m=kk*area+jj*img->nx+ii; ii<=ihi; m++,ii++) {
            x=niik_image_get_voxel(img,m);
            y=niik_image_get_voxel(lib,m+n*size);
            sx+=x;
            sy+=y;
            sxx+=x*x;
            syy+=y*y;
            sxy+=x*y;
            nn++;
          }
        }
      }
      out[n] = (nn * sxy - sx * sy) / sqrt(nn*sxx-sx*sx) / sqrt(nn*syy-sy*sy);
      if(verbose>0)fprintf(stdout,"[%3i,%3i,%3i] %12.7f   | %6.0f %6.0f %6.0f %6.0f\n",i,j,k,out[n],sx,sy,sxy,(double)nn);
      break;
    } /* method */
  } /* lib img */
  return 1;
}
