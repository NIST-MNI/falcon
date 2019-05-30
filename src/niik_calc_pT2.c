/* FILENAME:     niik_calc_pT2.c
 * DESCRIPTION:  Kunio's nifti1 pseudo-T2 calculation
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  program history\n"
  "  0.0.0 May 5, 2013; knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_describe[] = {
  "  [niik_calc_pT2.c] description\n"
};

static char *prog_help[] = {
  "  niik_calc_pT2:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
  "  -thresh  <thresh>           : threshold for PDw\n"
  "  -uthresh  <uthresh>         : upper threshold for T2w\n"
  "  -out-PD  <p-pd.mnc>         : pseudo-PD image for output\n"
  "  -T1SE <T1.mnc> <pT1.mnc> <TE> <TR>\n"
  "                              : input T1 image, output T1 filename, and its TE and TR\n"
  "  -FLASH <T1.mnc> <pT1.mnc> <TE> <TR> <FA>\n"
  "                              : input T1 image, output T1 filename, and its TE, TR, FA\n"
};

void usage() {
  fprintf(stdout,"niik_calc_pT2\n");
  fprintf(stdout,"  usage: [options] <pdw.mnc> <pd_te> <t2w.mnc> <t2_te> <t2_tr> <out_pT2.mnc>\n\n");
  fprintf(stdout,"    te and tr in msec\n\n");
}


int main(int argc,char *argv[],char *envp[]) {
  nifti_image
    *pt1img=NULL,
    *t1wimg=NULL,
    *ppdimg=NULL,
    *pt2img=NULL,
    *t2wimg=NULL,
    *pdwimg=NULL;
  int
    j,i,nc,sc;
  char
    *PD_fname=NULL,
    *T1_fname=NULL,
    fcname[32]="niik_calc_pT2";

  double
    d,d1,d2,dmin,*pt1=NULL,dd,
                uthresh_T2=5000,
                athresh=-1,
                EPS=1e-4,
                sFA,cFA,
                uni,Lt1,Ht1,dt,t1,pd,
                *pt2,*ppd,*t1w,*t2w,*pdw;
                
  double TET1=-1,FAT1=-1,FAT1rad=0,TRT1=-1,TET2=-1,TRPD=-1,TEPD=-1;

  niik_fc_display(fcname,1);

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

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

      else if(!strncmp(argv[nc],"-uthresh",8)) {
        NIIK_RET0((++nc>=argc),fcname,"missing uthresh value");
        uthresh_T2=atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-thresh",7)) {
        NIIK_RET0((++nc>=argc),fcname,"missing thresh value");
        athresh=atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-out-PD",7)) {
        PD_fname=argv[++nc];
      }

      else if(!strncmp(argv[nc],"-FLASH",6)) {
        NIIK_RET0((nc+5>=argc),fcname,"missing args for '-FLASH'");
        fprintf(stdout,"[%s] reading FLASH image:  %s\n",fcname,argv[nc+1]);
        NIIK_RET0(((t1wimg=niik_image_read(argv[++nc]))==NULL),fcname,"reading t1w, niik_image_read");
        NIIK_RET0((!niik_image_type_convert_scl(t1wimg,NIFTI_TYPE_FLOAT64,1)),fcname,"niik_image_type_convert");
        T1_fname=argv[++nc];
        TET1=atof(argv[++nc]);
        fprintf(stdout,"[%s] TE [T1]:              %-9.4f msec\n",fcname,TET1);
        TRT1=atof(argv[++nc]);
        fprintf(stdout,"[%s] TR [T1]:              %-9.4f msec\n",fcname,TRT1);
        FAT1=atof(argv[++nc]);
        FAT1rad = FAT1 / 180.0 * NIIK_PI;
        fprintf(stdout,"[%s] FA [T1]:              %-9.4f degrees\n",fcname,FAT1);
      }

      else if(!strncmp(argv[nc],"-T1SE",5)) {
        NIIK_RET0((nc+4>=argc),fcname,"missing args for '-T1'");
        fprintf(stdout,"[%s] reading SE image:     %s\n",fcname,argv[nc+1]);
        NIIK_RET0(((t1wimg=niik_image_read(argv[++nc]))==NULL),fcname,"reading t1w, niik_image_read");
        NIIK_RET0((!niik_image_type_convert_scl(t1wimg,NIFTI_TYPE_FLOAT64,1)),fcname,"niik_image_type_convert");
        T1_fname=argv[++nc];
        TET1=atof(argv[++nc]);
        TRT1=atof(argv[++nc]);
        fprintf(stdout,"[%s] TE [T1]:              %-9.4f msec\n",fcname,TET1);
        fprintf(stdout,"[%s] TR [T1]:              %-9.4f msec\n",fcname,TRT1);
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

  if(argc>7) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(0);
  } else if(argc<7) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(0);
  }


  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[1]);
  NIIK_RET0(((pdwimg=niik_image_read(argv[1]))==NULL),fcname,"reading pdw, niik_image_read");
  NIIK_RET0((!niik_image_type_convert_scl(pdwimg,NIFTI_TYPE_FLOAT64,1)),fcname,"niik_image_type_convert");
  TEPD = atof(argv[2]);
  fprintf(stdout,"[%s] TE [PD]:              %-9.4f msec\n",fcname,TEPD);
  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[3]);
  NIIK_RET0(((t2wimg=niik_image_read(argv[3]))==NULL),fcname,"reading pdw, niik_image_read");
  NIIK_RET0((!niik_image_type_convert_scl(t2wimg,NIFTI_TYPE_FLOAT64,1)),fcname,"niik_image_type_convert");
  TET2 = atof(argv[4]);
  fprintf(stdout,"[%s] TE [T2]:              %-9.4f msec\n",fcname,TET2);
  TRPD = atof(argv[5]);
  fprintf(stdout,"[%s] TR [T2]:              %-9.4f msec\n",fcname,TRPD);

  NIIK_RET0((pdwimg->nvox!=t2wimg->nvox),fcname,"wrong image size");
  pdw=pdwimg->data;
  t2w=t2wimg->data;

  NIIK_RET0(((pt2img=niik_image_copy(t2wimg))==NULL),fcname,"niik_image_copy(t2wimg)");
  pt2=pt2img->data;

  if(athresh<0) {
    athresh=-athresh;
    NIIK_RET0((!niik_image_thresh_ridler(pdwimg,NULL,&athresh)),fcname,"niik_image_thresh_ridler");
    athresh/=2;
  }
  fprintf(stdout,"[%s] threshold = %8.3f\n",fcname,athresh);


  /* t2 calculation */
  fprintf(stdout,"[%s] calculate pt2\n",fcname);
  for(i=0; i<pdwimg->nvox; i++) {
    pt2[i]=0;
    if(pdw[i]<athresh) continue;
    if(t2w[i]<EPS*3)   continue;
    pt2[i] = (TET2 - TEPD) / log( pdw[i]/t2w[i] );
    if(t2w[i]>pdw[i]) {
      pt2[i]=uthresh_T2;
      continue;
    }
    if(pt2[i]>uthresh_T2) {
      pt2[i]=uthresh_T2;
    } else if(pt2[i]<0) {
      pt2[i]=0;
      fprintf(stdout,"[%3i %3i %3i] %.5f %.5f\n",i%pdwimg->nx,(i/pdwimg->nx)%pdwimg->ny,i/pdwimg->nx/pdwimg->ny,pdw[i],t2w[i]);
    }
  }

  NIIK_RET0(((ppdimg=niik_image_copy(t2wimg))==NULL),fcname,"niik_image_copy(ppdimg)");
  NIIK_RET0((!niik_image_clear(ppdimg)),fcname,"niik_image_clear(ppdimg)");
  ppd=ppdimg->data;

  if(t1wimg!=NULL) {
    fprintf(stdout,"[%s] calculate pT1\n",fcname);
    NIIK_RET0(((t1wimg->nvox!=pdwimg->nvox)),fcname,"wrong #vox (pdw vs t1w)");

    NIIK_RET0(((pt1img=niik_image_copy(t2wimg))==NULL),fcname,"niik_image_copy(t2wimg)");
    NIIK_RET0((!niik_image_clear(pt1img)),fcname,"niik_image_clear");
    NIIK_RET0((!niik_image_clear(ppdimg)),fcname,"niik_image_clear");
    t1w=(double *)t1wimg->data;
    t2w=(double *)t2wimg->data;
    pdw=(double *)pdwimg->data;
    pt2=(double *)pt2img->data;
    pt1=(double *)pt1img->data;
    ppd=(double *)ppdimg->data;

    if(FAT1<0) {
      fprintf(stdout,"[%s] calculate pT1 from spin echo\n",fcname);
      for(i=0; i<pdwimg->nvox; i++) {
        if(pt2[i]<1) continue;
        if(t1w[i]<1) continue;
        for(j=0,uni=10000,Lt1=0,Ht1=100000; j<10; j++) { /* iterative search */
          for(d=Lt1,dd=1e10,dmin=0; d<=Ht1; d+=uni) {
            dt = pdw[i] * (1.0 - exp(-TRT1/d)) * exp(-TET1/pt2[i]) / (1.0 - exp(-TRPD/d)) / exp(-TEPD/pt2[i]);
            dt += t2w[i] * (1.0 - exp(-TRT1/d)) * exp(-TET1/pt2[i]) / (1.0 - exp(-TRPD/d)) / exp(-TET2/pt2[i]);
            dt /= 2.0;
            if(fabs(dt-t1w[i])<dd) {
              dmin=d;
              dd=fabs(dt-t1w[i]);
            }
          }
          Lt1=dmin-uni; /* for next resolution */
          Ht1=dmin+uni;
          uni/=10;
          Lt1=(Lt1<0)?0:Lt1;
        } /* iterative search */
        t1=dmin;
        pd=(pdw[i]/exp(-TEPD/pt2[i])/(1.0-exp(-TRPD/t1)) +
            t2w[i]/exp(-TET2/pt2[i])/(1.0-exp(-TRPD/t1)) +
            t1w[i]/exp(-TET1/pt2[i])/(1.0-exp(-TRT1/t1)) ) / 3.0;
        pt1[i]=t1;
        ppd[i]=pd;
      } /* each voxel */

    } else { /* FLASH */

      fprintf(stdout,"[%s] calculate pT1 from FLASH (gradient echo)\n",fcname);

      sFA=sin(FAT1rad);
      cFA=cos(FAT1rad);
      for(i=0; i<pdwimg->nvox; i++) {
        if(pt2[i]<1) continue;
        if(t1w[i]<1) continue;
        for(j=0,uni=10000,Lt1=0,Ht1=100000; j<10; j++) { /* iterative search */
          for(d=Lt1,dd=1e10,dmin=0; d<=Ht1; d+=uni) {
            dt = pdw[i] * sFA * (1.0 - exp(-TRT1/d)) * exp(-TET1/pt2[i]) / (1.0 - exp(-TRPD/d)) / (1.0 - exp(-TRT1/d)*cFA) / exp(-TEPD/pt2[i]);
            if(fabs(dt-t1w[i])<dd) {
              dmin=d;
              dd=fabs(dt-t1w[i]);
              if(i==1938814000) {  /* debug */
                fprintf(stdout,"[%s] [%3i,%3i,%3i] %6.0f   e=%-.12f | %9.2f %9.2f\n",fcname,
                        i%pdwimg->nx,(i/pdwimg->nx)%pdwimg->ny,i/pdwimg->nx/pdwimg->ny,dmin,dd,dt,t1w[i]);
              } /* debug */
            }
          } /* each t1 */
          Lt1=dmin-uni; /* for next resolution */
          Ht1=dmin+uni;
          uni/=10;
          Lt1=(Lt1<0)?0:Lt1;
        } /* iterative search */
        t1=dmin;
        pd=(pdw[i]/exp(-TEPD/pt2[i])/(1.0-exp(-TRPD/t1)) +
            t2w[i]/exp(-TET2/pt2[i])/(1.0-exp(-TRPD/t1)) +
            t1w[i]/exp(-TET1/pt2[i])/(1.0-exp(-TRT1/t1))/sFA*(1.0-exp(-TRT1/t1)*cFA) ) / 3.0;
        pt1[i]=t1;
        ppd[i]=pd;

        if(i==1938814) {  /* debug */
          fprintf(stdout,"[%s] [%3i,%3i,%3i] %6.0f %6.0f %6.0f  ->  %4.0f   %6.0f %6.0f\n",fcname,
                  i%pdwimg->nx,(i/pdwimg->nx)%pdwimg->ny,i/pdwimg->nx/pdwimg->ny,
                  pdw[i],t2w[i],t1w[i],pt2[i],
                  pd,t1);
          fprintf(stdout,"[%s] [%3i,%3i,%3i] E2 = %12.7f\n",fcname,i%pdwimg->nx,(i/pdwimg->nx)%pdwimg->ny,i/pdwimg->nx/pdwimg->ny,exp(-TEPD/pt2[i]));
        } /* debug */

      } /*each voxel*/
    } /* FLASH */
  } else { /* do PD processing without t1 */
    fprintf(stdout,"[%s] calculate pPD (with t1-weighting)\n",fcname);
    /* S = k PD (1 - exp(-TR/T1)) exp(-TE/T2*)
     */
    for(i=0; i<pdwimg->nvox; i++) {
      if(pt2[i]<EPS) {
        continue;
      }
      d1 = exp(-TEPD/pt2[i]);
      if(d1 < EPS) {
        continue;
      }
      d2 = exp(-TET2/pt2[i]);
      if(d2 < EPS) {
        continue;
      }
      ppd[i] =
        ( pdw[i] / ( EPS + d1 ) +
          t2w[i] / ( EPS + d2 ) ) / 2;
    }
  }



  if(PD_fname!=NULL) {
    fprintf(stdout,"[%s] write pseudo-PD\n",fcname);
    NIIK_RET0((!niik_image_type_convert(ppdimg,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
    fprintf(stdout,"[%s] writing image:        %s\n",fcname,PD_fname);
    NIIK_RET0((!niik_image_write(PD_fname,ppdimg)),fcname,"niik_write_image");
  }

  if(T1_fname!=NULL) {
    fprintf(stdout,"[%s] write pseudo-T1\n",fcname);
    NIIK_RET0((!niik_image_type_convert(pt1img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
    fprintf(stdout,"[%s] writing image:        %s\n",fcname,T1_fname);
    NIIK_RET0((!niik_image_write(T1_fname,pt1img)),fcname,"niik_write_image");
  }

  fprintf(stdout,"[%s] write pseudo-T2\n",fcname);
  NIIK_RET0((!niik_image_type_convert(pt2img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
  fprintf(stdout,"[%s] writing image:        %s\n",fcname,argv[6]);
  NIIK_RET0((!niik_image_write(argv[6],pt2img)),fcname,"niik_write_image");

  pdwimg=niik_image_free(pdwimg);
  t2wimg=niik_image_free(t2wimg);
  t2wimg=niik_image_free(pt2img);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_calc_T2 */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/