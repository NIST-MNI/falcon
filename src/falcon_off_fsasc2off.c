/* Filename:     nifti1_kunio_off_fsasc2off.c
 * Description:  converts FreeSurfer's surface in asci format to OFF format
 * Author:       Kunio Nakamura
 * Date:         March 8, 2012
 */

#ifndef _FALCON_OFF_FSASC2OFF_C_
#define _FALCON_OFF_FSASC2OFF_C_

#include "falcon.h"

/*****************************************************
 *
 * off_fsasc2off
 *
 * -need to convert the freesurfer file to ascii format
 *  for example,
 * -this is a quick solution (not permanent)
 * -there is a 4th number for each vertex and
 *  for the triangles.
 * -according to FreeSurfer webpage
 *  http://surfer.nmr.mgh.harvard.edu/fswiki/FileFormats
 *  this is a ripflag but I don't know what that is...
 *
 * $ mris_convert lh.pial lh.test.asc
 * $ niikmath -fsasc=lh.test.asc -out=lh.test.off fsasc2off
 * $ geomview lh.test.off
 *
 *****************************************************/

int niik_off_fsasc2off(char *iname,char *oname) {
  FILE *ifp,*ofp;
  char line1[1024];
  int
  n,
  verbose=1,
  tri[4],
  nvert,nface;
  niikpt p;
  char fcname[32]="niik_off_fsasc2off";

  niik_fc_display(fcname,1);
  if(iname==NULL) {
    fprintf(stderr,"ERROR: iname is null\n");
    return 0;
  }
  if(oname==NULL) {
    fprintf(stderr,"ERROR: oname is null\n");
    return 0;
  }

  if((ifp=fopen(iname,"r"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",iname);
    return 0;
  }
  if((ofp=fopen(oname,"w"))==NULL) {
    fprintf(stderr,"ERROR: fopen %s\n",oname);
    return 0;
  }

  NIIK_RET0(((fgets(line1,1024,ifp))==NULL),fcname,"fgets problem");
  fprintf(stdout,"  first line:  %s#END#\n",line1);

  NIIK_RET0((!fscanf(ifp,"%i %i\n",&nvert,&nface)),fcname,"fscan problem 1");
  fprintf(ofp,"OFF\n%i %i 0\n",nvert,nface);
  if(verbose) fprintf(stdout,"[%s] v,f,e = %i %i 0\n",fcname,nvert,nface);

  for(n=0; n<nvert; n++) {
    NIIK_RET0((!fscanf(ifp,"%lf %lf %lf %lf\n",&p.x,&p.y,&p.z,&p.w)),fcname,"fscanf problem 2");
    fprintf(ofp,"%15.9f %15.9f %15.9f\n",p.x,p.y,p.z);
  }
  for(n=0; n<nface; n++) {
    NIIK_RET0((!fscanf(ifp,"%i %i %i %i\n",&tri[0],&tri[1],&tri[2],&tri[3])),fcname,"fscan problem 3");
    fprintf(ofp,"3 %7i %7i %7i\n",tri[0],tri[1],tri[2]);
  }

  fclose(ifp);
  fclose(ofp);

  niik_fc_display(fcname,0);
  return 1;
} /* int niik_off_fsasc2off(char *iname,char *oname) */


#endif /* _FALCON_OFF_FSASC2OFF_C_ */
