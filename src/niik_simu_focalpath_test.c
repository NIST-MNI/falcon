/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_morph.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (1)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"\n");
  fprintf(stdout,"  creates simulated lesion for simlation MRI\n");
  fprintf(stdout,"  usage: <in.mnc> <les.mnc> <model.mnc> <pwm.mnc> <in_WM_val> <out_WM_val> <out.mnc> \n\n");
  fprintf(stdout,"  <in.mnc>          : real T1-weighted MRI where lesion was segmented\n");
  fprintf(stdout,"  <les.mnc>         : lesion mask\n");
  fprintf(stdout,"  <model.mnc>       : simulated MRI on which lesion will be created\n");
  fprintf(stdout,"  <pwd.mnc>         : simulation MRI's WM probability map\n");
  fprintf(stdout,"                      i.e. model.mnc's wm probability map\n");
  fprintf(stdout,"  <in_WM_val>       : WM intensity for <in.mnc>\n");
  fprintf(stdout,"  <out_WM_val>      : WM intensity for <model.mnc>\n");
  fprintf(stdout,"  <out.mnc>         : output image\n");
  fprintf(stdout,"\n  all images need to be in the same space with same dimensions\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optonal usage:\n");
  fprintf(stdout,"  -fwhm=<FWHM>       : gaussian filter for lesion mask [default=1.3]\n");
  fprintf(stdout,"  -rice=<s>          : rice distribution parameters [default=10]\n");
  fprintf(stdout,"                     : e.g., -rice=15,15\n");
  fprintf(stdout,"  -fill=<fill.mnc>   : filled image, the intensity is used to estimate NAWM\n");
  fprintf(stdout,"                     : ignores <in_WM_val>\n");
}

void version() {
  fprintf(stdout,"version \n");
  fprintf(stdout,"0.0.0, June 13, 2013, Kunio Nakamura, knakamura@mrs.mni.mcgill.ca\n");
  fprintf(stdout,"  -initial version\n\n");
  fprintf(stdout,"0.1.0, August 6, 2013, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n");
  fprintf(stdout,"  -added option to include filled image to estimate NAWM intensity\n\n");
}


int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *fill=NULL,
   *pwm=NULL,
    *out=NULL,
     *mask=NULL, *maskblur=NULL,
      *img=NULL;
  unsigned char *bimg=NULL;
  int
  i,n,
  nc=1,
  sc=1;
  char
  **strlist,
  fcname[32]="test";
  double
  srice=10,
  blurFWHM=1.3,
  dilateR=2.3,
  iGM=0,oGM=0,
  iWM=0,oWM=0,v,*dpwm=NULL;

  if(argc==1) {
    usage();
    exit(0);
  }

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--help",6)) {
        usage();
        exit(9);
      }

      else if(!strncmp(argv[nc],"-dilate=",8)) {
        dilateR=atof(argv[nc]+8);
      }

      else if(!strncmp(argv[nc],"-fill=",6)) {
        fprintf(stdout,"[%s] reading fill image:   %s\n",fcname,argv[nc]+6);
        NIIK_EXIT(((fill=niik_image_read(argv[nc]+6))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-fwhm=",6)) {
        blurFWHM=atof(argv[nc]+6);
      }

      else if(!strncmp(argv[nc],"-rice=",6)) {
        srice=atof(argv[nc]+6);
      }

      /* more options here */

      else if(!strncmp(argv[nc],"-gm=",4)) {
        strlist=niik_csstring_to_list(argv[nc]+4,&n);
        NIIK_EXIT((n!=2),fcname,"need 2 args for -gm=<iGM>,<oGM>",1);
        iGM=atof(strlist[0]);
        oGM=atof(strlist[1]);
        free(strlist[0]);
        free(strlist[1]);
        free(strlist);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        usage();
        exit(9);
      } else if(!strncmp(argv[nc],"-u",2)) {
        usage();
        exit(9);
      } else if(!strncmp(argv[nc],"-v",2)) {
        version();
        exit(9);
      } else {
        fprintf(stderr,"[%s] ERROR: unknown option, %s\n",fcname,argv[nc]);
        exit(9);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  }
  argc=sc;
  NIIK_EXIT((argc!=8),fcname,"wrong usage",9);

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading image:        %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"reading img, niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading lesion:       %s\n",fcname,argv[2]);
  NIIK_EXIT(((mask=niik_image_read(argv[2]))==NULL),fcname,"reading mask, niik_image_read",9);
  fprintf(stdout,"[%s]   gaussian filter, %.5f\n",fcname,blurFWHM);
  NIIK_EXIT(((maskblur=niik_image_copy_as_type(mask,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type",9);
  NIIK_EXIT((!niik_image_filter_gaussian_update(maskblur,9,blurFWHM)),fcname,"niik_image_filter_gaussian_update",9);
  fprintf(stdout,"[%s]   dilation, %.5f\n",fcname,dilateR);
  NIIK_EXIT((!niik_image_morph_gray_dilate(maskblur,dilateR,0.001)),fcname,"niik_image_morph_gray_dilate",9);
  niik_image_write("tmp_mask_blur.mnc",maskblur);
  NIIK_EXIT(((bimg=niik_image_get_voxels_as_uint8_vector(mask))==NULL),fcname,"niik_image_get_voxels_as_uint8_vector",9);
  for(i=0; i<img->nvox; i++) {
    if(bimg[i]>0) niik_image_set_voxel(maskblur,i,1.0);
  }

  fprintf(stdout,"[%s] reading model:        %s\n",fcname,argv[3]);
  NIIK_EXIT(((out=niik_image_read(argv[3]))==NULL),fcname,"reading model, niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert(out,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading wm map:       %s\n",fcname,argv[4]);
  NIIK_EXIT(((pwm=niik_image_read(argv[4]))==NULL),fcname,"reading pwm, niik_image_read",9);
  NIIK_EXIT(((dpwm=niik_image_get_voxels_as_double_vector(pwm))==NULL),fcname,"niik_image_get_voxels_as_double_vector",9);

  iWM = atof(argv[5]);
  oWM = atof(argv[6]);

  fprintf(stdout,"[%s] noise  %9.3f \n",fcname,srice);
  fprintf(stdout,"[%s] WM   %9.3f %9.3f \n",fcname,iWM,oWM);
  fprintf(stdout,"[%s] GM   %9.3f %9.3f \n",fcname,iGM,oGM);

  for(i=0; i<img->nvox; i++) {
    dpwm[i] *= niik_image_get_voxel(maskblur,i);
    if(dpwm[i]<0.0001) continue;
    /* this is good, but oWM is not constant over the entire volume ...
       v = niik_image_get_voxel(img,i) / iWM * oWM + NIIK_RiceRnd(vrice,srice); */
    if(i==136+170*img->nx+74*img->nx*img->ny) {
      if(fill==NULL)
        fprintf(stdout,"[%s] %3i %3i %3i  %7.2f  %7.2f %4.2f\n",fcname,
                i%img->nx,(i/img->nx)%img->ny,i/img->nx/img->ny,
                niik_image_get_voxel(img,i),
                niik_image_get_voxel(out,i),
                niik_image_get_voxel(pwm,i));
      else
        fprintf(stdout,"[%s] %3i %3i %3i  %7.2f  %7.2f  %7.2f  %7.2f  %4.2f\n",fcname,
                i%img->nx,(i/img->nx)%img->ny,i/img->nx/img->ny,
                niik_image_get_voxel(img,i),
                niik_image_get_voxel(fill,i),
                niik_image_get_voxel(out,i),
                dpwm[i],
                niik_image_get_voxel(pwm,i));
    }
    if(fill!=NULL) iWM=niik_image_get_voxel(fill,i);
    oWM = niik_image_get_voxel(out,i);
    v = (niik_image_get_voxel(img,i)-iGM) / (iWM-iGM) * (oWM-oGM) + oGM;
    v = NIIK_RiceRnd(v,srice);
    v = v * dpwm[i] + niik_image_get_voxel(out,i) * (1.0 - dpwm[i]);
    if(i==99+80*img->nx+90*img->nx*img->ny)
      fprintf(stdout,"[%s] %3i %3i %3i  %7.5f -> %7.5f\n",fcname,
              i%img->nx,(i/img->nx)%img->ny,i/img->nx/img->ny,
              niik_image_get_voxel(out,i),v);
    niik_image_set_voxel(out,i,v);
  } /* each voxel */

  fprintf(stdout,"[%s] writing image:        %s\n",fcname,argv[7]);
  NIIK_EXIT((!niik_image_write(argv[7],out)),fcname,"niik_image_write output",9);

  img=niik_image_free(img);
  out=niik_image_free(out);

  niik_fc_display(fcname,0);
  exit(0);
} /* main */
