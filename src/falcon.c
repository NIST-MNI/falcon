/* Filename:     falcon.c
 * Description:  general nifti1 functions by Kunio
 * Author:       Kunio Nakamura
 * Date:         February 23, 2012
 */

#include <stdlib.h>
#include <pwd.h>
#include <unistd.h>
#include "falcon.h"

/***************************************************
 *
 * BASIC FUNCTIONS
 *
 ***************************************************/



double NIIK_RiceRnd(double v,double s) {
  double x,y;
  x = s * niik_get_rand_normal() + v;
  y = s * niik_get_rand_normal();
  return sqrt(x*x+y*y);
}

double NIIK_RiceRnd2(double v,double s) {
  double x,y;
  x = s * NIIK_BoxMuller((rand()%10000)/10000.0,(rand()%10000)/10000.0) + v;
  y = s * NIIK_BoxMuller((rand()%10000)/10000.0,(rand()%10000)/10000.0);
  /* fprintf(stdout,"[NIIK_RiceRnd2] %9.4f %9.4f\n",x,y); */
  return sqrt(x*x+y*y);
}

int niik_access_file(char *f) {
  if(access(f,F_OK)!=-1) return 1;
  return 0;
}


double NIIK_Heaviside(double val,double eps) {
  if     (val>= eps) return 1;
  else if(val<=-eps) return 0;
  val/=eps;
  val = 0.5*(1.0+val+sin(NIIK_PI*val)/NIIK_PI);
  return val;
}

double NIIK_Heaviside11(double val,double eps) {
  return NIIK_Heaviside(val,eps)*2.0-1.0;
}

double NIIK_DiracDelta(double val,double eps) {
  if(val<-eps) return 0;
  else if(val>eps) return 0;
  return 0.5/eps+0.5*cos(NIIK_PI*val/eps)/eps;
}

double NIIK_BoxMuller(double v1,double v2)
/* -v1 and v2 are [0,1)
 * -output is normally distribued (0,1) */
{
  if(v1==0) v1+=1e-5;
  return sqrt(-2.0*log(v1))*cos(NIIK_PI*2.0*v2);
}

double niik_linear_step(double val,double eps) {
  if     (val>= eps) return 1;
  else if(val<=-eps) return 0;
  return (val+eps)/eps/2.0;
}

double niik_pvc(double val,double m1, double m2) {
  if(val<=m1) return 1;
  if(val>=m2) return 0;
  return (val-m2)/(m1-m2);
}

double niik_pv(double val,double m1, double m2)
/* val is [0, 1]  (inclusive)
 * m1 is the output when val is 0
 * m2 is the output when val is 1
 * for the rest, output is linearly interpolated according to val */
{
  return val*m2 + (1.0-val)*m1;
}

int niik_triangular_number(int num)
/* output = 1 + 2 + 3 + ... up to num */
{
  if(num<=0) return 0;
  return num*(num+1)/2;
}

double niik_legendre_func(double x,int deg) {
  switch(deg) {
  case 0:
    return 1;
  case 1:
    return x;
  case 2:
    return (3.0 * x * x - 1.0)/2.0;
  case 3:
    return (5.0 * x * x * x - 3.0 * x)/2.0;
  case 4:
    x*=x;
    return (35.0 * x * x - 30.0 * x + 3.0)/8.0;
  case 5:
    return (63.0 * x*x*x*x*x - 70.0 * x*x*x + 15.0 * x)/8.0;
  case 6:
    x*=x;
    return (231.0 * x*x*x - 315.0 * x*x + 105.0 * x - 5.0)/16.0;
  case 7:
    return (429.0 * x*x*x*x*x*x*x - 693.0 * x*x*x*x*x + 315.0 * x*x*x - 35.0 * x)/16.0;
  case 8:
    x*=x;
    return (6435.0 * x*x*x*x - 12012.0 * x*x*x + 6930.0 * x*x - 1260.0 * x + 35.0)/128.0;
  case 9:
    return (12155.0 * x*x*x*x*x*x*x*x*x - 25740.0 * x*x*x*x*x*x*x + 18018.0 * x*x*x*x*x - 4620.0 * x*x*x + 315.0 * x)/128.0;
  case 10:
    x*=x;
    return (46189.0 * x*x*x*x*x - 109395.0 * x*x*x*x + 90090.0 * x*x*x - 30030.0 * x*x + 3465.0 * x - 63.0)/256.0;
    fprintf(stderr,"[niik_legendre_func] ERROR: unknown degree, %i\n",deg);
    exit(0);
  }
  return 0;
} /* niik_legendre_func */

int niik_legendre_func_with_param(double x,double *y,double *param,int deg) {
  int n;
  double d=0;
  if(param==NULL) {
    fprintf(stderr,"[niik_legendre_func_with_param]: ERROR: param is null\n");
    return 0;
  }
  if(deg>10) {
    fprintf(stderr,"[niik_legendre_func_with_param]: ERROR: deg is too large\n");
    return 0;
  }
  d=param[0];
  for(n=1; n<=deg; n++) {
    d += param[n] * niik_legendre_func(x,n);
  }
  *y=d;
  return 1;
} /* niik_legendre_func_with_param */

int niik_version_display(const char *fcname,int major,int minor,int micro) {
  fprintf(stdout,"  %s: version %i.%i.%i\n",fcname,major,minor,micro);
  return 1;
}

int niik_fc_display(const char *fcname,int flag_start) {
  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);
  switch(flag_start) {
  case 0:
    fprintf(stdout,"[%s] finished at %s\n",fcname,tmstr);
    break;
  case 1:
    fprintf(stdout,"[%s] start    at %s\n",fcname,tmstr);
    break;
  default:
    fprintf(stdout,"[%s] func\n",fcname);
    break;
  }
  return 1;
} /* niik_fc_display */

int niik_check_fsldir_exists() {
  char *FSLDIR=NULL;
  char fcname[32]="niik_check_fsdir_exists";
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv FSLDIR\n",fcname);
    fprintf(stderr,"  export FSLDIR=/usr/local/fsl\n");
    fprintf(stderr,"  export FSLDIR=/lab1/fsl/4.1.7\n");
    return 1;
  }
  return 0;
} /* niik_check_fsdir_exists */


unsigned short niik_float16(float h) {
  float *pf;
  unsigned long l,*pl;
  unsigned short o=0;
  int s,e,m;
  pf=&h;
  pl=(unsigned long *)pf;
  l=(*pl);
  if(h==0x0000000) return 0x0000;
  if(h==0x8000000) return 0x8000;
  if(h> 65504) return 0x7BFF;
  if(h<-65504) return 0xFBFF;
  if(fabs(h)<(1.0/(1<<24)))  {
    if(l&0x8000000) return 0x8000;
    return 0x0000;
  }
  if(fabs(h)<(1.0/(1<<14))) {
    s=(l & 0x80000000)>>31;
    m=(l & 0x007FFFFF);
    h*=(1<<14); /* account for 2^-14 expo term */
    h*=(1<<10); /* account for 2^10 mantissa term */
    m=(int)floor(fabs(h)+0.5); /* calculate mantissa */
    o=0;
    o = ((l&0x80000000)>>16) + (m&0x000003FF);  /* s & e term */
    return o;
  }
  s=(l & 0x80000000)>>31;
  e=(l & 0x7F800000)>>23;
  m=(l & 0x007FFFFF);
  if(e==255) {
    if(m==0) {
      if(s) return 0xFC00; /* negative infinity */
      else  return 0x7C00; /* positive infinity */
    } else {
      if(s) return 0xFC01; /* negative NaN */
      else  return 0x7C01; /* positive NaN */
    }
  }
  e=e-127+15;
  o = (e<<10); /* e term */
  o = o +
      (((l&0x80000000)>>31)<<15); /* s term */
  o = o +
      ((l&0x007FE000)>>13);       /* m term */
  if(l&0x00001000) {
    o+=1;
  }
  return o;
}

float niik_unfloat16(unsigned short h) {
  int s,e,m;
  float f;
  s=h>>15;
  e = (h & 0x7C00u) >> 10;
  m=(h & 0x03FF);
  if(h == 0x7c00) return  HUGE_VAL; /* positive infinity */
  if(h == 0xfc00) return -HUGE_VAL; /* negative infinity */
  if(h == 0x0000) return  0;  /* positive zero */
  if(h == 0x8000) return -0;  /* negative zero */
  if(e==0) {
    if(m==0) {
      if(s==0) return 0x00000000;
      else     return 0x80000000;
    }
    f=1.0/(1<<14);
    f=f*((float)m/(1<<10));
    if(s) return -f;
    return f;
  } else if(e==31) {
    if(m>0) return 0x7FFFFFFFul;
  }
  if(e==15) {
    f=1;
  } else if(e>15) {
    e=e-15;
    f=(1<<e);
  } else  {
    e=15-e;
    f=(1<<e);
    f=1.0/f;
  }
  f=f*(1.0+((float)m/(1<<10)));
  if(s) return -f;
  return f;
}



int niik_check_double_problem(double val) {
  if(val>1e300) return 1;
  return 0;
}

void niik_disp_exec_info(char *progname,int major,int minor,int micro) {
  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  char hostname[512];
  char *p=NULL;
  size_t len=512;
  register struct passwd *pw;
  register uid_t uid;
  uid = geteuid ();
  pw = getpwuid (uid);
  if (pw) {
    p=(char *)pw->pw_name;
  }
  gethostname(hostname,len);
  srand(time(NULL));
  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);
  fprintf(stdout,"[%s;%s@%s] version %i.%i.%i   exec @ %s\n",progname,p,hostname,major,minor,micro,tmstr);
  return;
}

int niik_csstring_count(char *CSstr) {
  char
  *prevstr,*currstr,
  *cptr;
  int
  slen[4096],
       n,N;
  if(CSstr==NULL) {
    fprintf(stderr,"ERROR: CSstr is null\n");
    return 0;
  }
  for(n=0; n<4096; n++) slen[n]=0;
  N=0;
  cptr=currstr=prevstr=CSstr;
  while((cptr=strchr(cptr,','))!=NULL) {
    /* cptr[0] = '\0'; */
    cptr++;
    currstr=cptr;
    while(prevstr!=currstr-1) {
      slen[N]++;
      prevstr++;
    }
    prevstr=currstr;
    N++;
  }
  while(prevstr[0]!=0 || prevstr==NULL || prevstr[0]==32) {
    slen[N]++;
    prevstr++;
  }
  N++;
  return N;
} /* niik_csstring_count */


void    niik_free_list(char **list,int num) {
  int i;
  for(i=0; i<num; i++)
    free(list[i]);
  free(list);
}

int *niik_csstring_to_int_list(char *CSstr,int *num) {
  char fcname[64]="niik_csstring_to_int_list";
  int *v=NULL,N,i;
  char **CSlist;
  if((CSlist=niik_csstring_to_list(CSstr,&N))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_csstring_to_list\n",__FILE__,__LINE__,fcname);
    return NULL;
  }
  *num=N;
  NIIK_RET0(((v=(int *)calloc(N,sizeof(int)))==NULL),fcname,"calloc, N");
  for(i=0; i<N; i++) {
    v[i]=atoi(CSlist[i]);
    free(CSlist[i]);
  }
  free(CSlist);
  return v;
} /* niik_csstring_to_int_list */

double *niik_csstring_to_double_list(char *CSstr,int *num) {
  double *v=NULL;
  char **CSlist;
  int N,i;
  NIIK_RET0(((CSlist=niik_csstring_to_list(CSstr,&N))==NULL),__func__,"niik_csstring_to_list");
  *num=N;
  v=(double *)calloc(N,sizeof(double));
  for(i=0; i<N; i++) {
    v[i]=atof(CSlist[i]);
    free(CSlist[i]);
  }
  free(CSlist);
  return v;
} /* niik_csstring_to_double_list */


char **niik_txtfile_to_list(char *filename,int *num) {
  FILE *fp=NULL;
  char **str;
  char line[65536];
  char *cptr;
  int nline=0;

  // count # lines
  NIIK_RETURN(((fp=fopen(filename,"r"))==NULL),"could not open file",NULL);
  while( fgets(line,65536,fp) ) {
    nline++;
  }
  str=(char **)calloc(nline,sizeof(char *));
  fclose(fp);
  *num=nline;
  nline=0;

  // copy the contents
  NIIK_RETURN(((fp=fopen(filename,"r"))==NULL),"could not open file",NULL);
  while( fgets(line,65536,fp) ) {
    if((cptr = strchr(line,'\n'))!=NULL) {
      *cptr = '\0';
    }
    str[nline]=(char *)calloc(strlen(line)+1,sizeof(char));
    sprintf(str[nline],"%s",line);
    nline++;
  }
  fclose(fp);
  return str;
}

char **niik_csstring_to_list(char *CSstr,int *num)
/* -replaces ',' with '\0'
 * -returns the list of char pointers to the character after \0 and the first one
 * -num is replaced with number of output (char **) elements
 */
{
  char
    *prevstr,*currstr,
    *cptr,
    **outstr;
  int
      slen[4096], m,n,N;
  if(CSstr==NULL) {
    fprintf(stderr,"ERROR: CSstr is null\n");
    return NULL;
  }
  for(n=0; n<4096; n++) slen[n]=0;
  N=0;
  cptr=currstr=prevstr=CSstr;
  while((cptr=strchr(cptr,','))!=NULL && N<4095) {
    /* cptr[0] = '\0'; */
    cptr++;
    currstr=cptr;
    while(prevstr!=currstr-1) {
      slen[N]++;
      prevstr++;
    }
    prevstr=currstr;
    N++;
  }
  while(prevstr[0]!=0 || prevstr==NULL || prevstr[0]==32) { /*why NULL*/
    slen[N]++;
    prevstr++;
  }
  N++;
  outstr=(char **)calloc(N,sizeof(char *));
  for(n=0,cptr=CSstr; n<N; n++,cptr++) {
    outstr[n]=(char *)calloc(slen[n]+1,sizeof(char));
    for(m=0; m<slen[n]; m++,cptr++)
      outstr[n][m]=*cptr;
    outstr[n][m]=0;
  }
  *num=N;
  return outstr;
} /* niik_csstring_to_list */

int niik_underscore_to_space(char *str) {
  char *c;
  /*fprintf(stdout,"     %s\n",str);*/
  if((c=strchr(str,'_'))==NULL) {
    return 1;
  }
  c[0]=' ';
  niik_underscore_to_space(c+1);
  return 1;
}




/***************************************************
 *
 * COLOR FUNCTIONS
 *
 ***************************************************/

niikmat *niik_colormap_get(int ctype, int num) {
  niikmat *mat;
  switch(ctype) {
  case NIIK_COLORMAP_ATROPHY:
    if((mat=niik_colormap_get_atrophy(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_atrophy\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_NEG_ATROPHY:
    if((mat=niik_colormap_get_neg_atrophy(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_neg_atrophy\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_POS_ATROPHY:
    if((mat=niik_colormap_get_pos_atrophy(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_pos_atrophy\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_JACOBIAN:
    if((mat=niik_colormap_get_jacobian(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_jacobian\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_SUMMER:
    if((mat=niik_colormap_get_summer(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_summer\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_SPECTRAL:
    if((mat=niik_colormap_get_spectral(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_spectral\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_GRAYSCALE:
    if((mat=niik_colormap_get_grayscale(num))==NULL) {
      fprintf(stderr,"[niik_colormap_get] ERROR: niik_colormap_get_grayscale\n");
      return NULL;
    }
    return mat;
  case NIIK_COLORMAP_UNKNOWN:
  default:
    fprintf(stderr,"[niik_colormap_get] ERROR: unknown type: %i\n",ctype);
    return NULL;
  }
  return NULL;
}

niikmat *niik_colormap_get_atrophy(int num) {
  niikmat *mat;
  int i;
  double d;
  double r[3]= {0,1,1};
  double g[3]= {0,1,0};
  double b[3]= {1,1,0};
  mat = niikmat_init(num,3);
  for(i=0; i<num; i++) {
    d=(double)i/(num-1.0)*2.0;
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,3,d);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,3,d);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,3,d);
  }
  return mat;
}

niikmat *niik_colormap_get_pos_atrophy(int num) {
  niikmat *mat;
  int i;
  double d;
  double r[3]= {0,1,1};
  double g[3]= {0,1,0};
  double b[3]= {1,1,0};
  mat = niikmat_init(num,3);

  for(i=0; i<num; i++) {
    d=(double)i/(num-1.0)+1.0;
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,3,d);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,3,d);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,3,d);
  }
  return mat;
}

niikmat *niik_colormap_get_neg_atrophy(int num) {
  niikmat *mat;
  int i;
  double d;
  double r[3]= {0,1,1};
  double g[3]= {0,1,0};
  double b[3]= {1,1,0};
  mat = niikmat_init(num,3);

  for(i=0; i<num; i++) {
    d=(double)i/(num-1.0);
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,3,d);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,3,d);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,3,d);
  }
  return mat;
}


niikmat *niik_colormap_get_jacobian(int num) {
  niikmat *mat;
  int i;
  double d;
  double r[5]= {0.8,0,0,1,1.0};
  double g[5]= {0.8,0,0,0,0.8};
  double b[5]= {1.0,1,0,0,0.8};
  mat = niikmat_init(num,3);
  for(i=0; i<num; i++) {
    d=(double)i/(num-1.0)*4.0;
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,5,d);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,5,d);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,5,d);
  }
  return mat;
}

niikmat *niik_colormap_get_summer(int num) {
  niikmat *mat;
  int i;
  mat = niikmat_init(num,3);
  for(i=0; i<num; i++) {
    mat->m[i][0]=1.0-i/(num-1.0);
    mat->m[i][1]=1.0-i/(num-1.0)*0.5;
    mat->m[i][2]=0.40;
  }
  return mat;
}

niikmat *niik_colormap_get_spectral(int num) {
  niikmat *mat;
  int i;
  double x;
  double r[21]= {0.0000,0.4667,0.5333,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.7333,0.9333,1.0000,1.0000,1.0000,0.8667,0.8000,0.8000};
  double g[21]= {0.0000,0.0000,0.0000,0.0000,0.0000,0.4667,0.6000,0.6667,0.6667,0.6000,0.7333,0.8667,1.0000,1.0000,0.9333,0.8000,0.6000,0.0000,0.0000,0.0000,0.8000};
  double b[21]= {0.0000,0.5333,0.6000,0.6667,0.8667,0.8667,0.8667,0.6667,0.5333,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.8000};
  if(num<21) {
    fprintf(stderr,"ERROR: num is too small for spectral\n");
    return NULL;
  }
  mat = niikmat_init(num,3);
  for(i=0; i<num; i++) {
    x = (double)i/(num-1.0)*20.0;
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,21,x);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,21,x);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,21,x);
  }
  return mat;
}

niikmat *niik_colormap_get_grayscale(int num) {
  niikmat *mat;
  int i;
  double x;
  double r[2]= {0.0000,1.0000};
  double g[2]= {0.0000,1.0000};
  double b[2]= {0.0000,1.0000};
  if(num<2) {
    fprintf(stderr,"ERROR: num is too small for grayscale\n");
    return NULL;
  }
  mat = niikmat_init(num,3);
  for(i=0; i<num; i++) {
    x = (double)i/(num-1.0);
    mat->m[i][0]=niik_interp1d_linear_in_double_vector(r,2,x);
    mat->m[i][1]=niik_interp1d_linear_in_double_vector(g,2,x);
    mat->m[i][2]=niik_interp1d_linear_in_double_vector(b,2,x);
  }
  return mat;
}


nifti_image *niik_image_map_color_overlay_mask(nifti_image *img,double imin,double imax,
    nifti_image *map,double map_min,double map_max,
    double bg_factor,
    nifti_image *maskimg,niikmat *cmap)
/* creates a 'nice' image for QC
 * -img is the grayscale image to be overlapped
 * -imin imax are the min/max for the img
 * -map is the image with colormap
 * -map_min map_max are the min/max for map image
 * -maskimg is the mask image
 * -cmap color map, should be jacobian atrophy type
 */
{
  nifti_image *outimg=NULL;
  char fcname[64]="niik_image_map_color_overlay_mask";
  double step,d,rgb[3];
  int i,j,m,n;
  int const verbose=1;

  if(verbose>0) niik_fc_display(fcname,1);
  step=(map_max-map_min)/(cmap->row-1.0);
  bg_factor = (bg_factor<0.0)?0.2:bg_factor;

  if(verbose>0) {
    fprintf(stdout,"[%s] image min/max   %12.6f %12.6f\n",fcname,imin,imax);
    fprintf(stdout,"[%s] cmap min/max    %12.6f %12.6f\n",fcname,map_min,map_max);
    fprintf(stdout,"[%s] bg factor       %12.6f\n",fcname,bg_factor);
  }

  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((map==NULL),fcname,"map is null");
  NIIK_RET0((cmap==NULL),fcname,"cmap is null");
  /* NIIK_RET0((img->nvox!=map->nvox),fcname,"nvox doesn't match, img,map"); */
  if(img->nx!=map->nx) {
    fprintf(stderr,"[%s] ERROR: nx didn't match, %i %i\n",fcname,img->nx,map->nx);
    return 0;
  }
  if(img->ny!=map->ny) {
    fprintf(stderr,"[%s] ERROR: ny didn't match, %i %i\n",fcname,img->ny,map->ny);
    return 0;
  }
  if(img->nz!=map->nz) {
    fprintf(stderr,"[%s] ERROR: nz didn't match, %i %i\n",fcname,img->nz,map->nz);
    return 0;
  }
  NIIK_RET0((niik_image_calc_nvox3d(img)!=niik_image_calc_nvox3d(map)),fcname,"nvox doesn't match, img,map");

  NIIK_RET0(((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET0((!niik_image_convert_to_color_image(outimg)),fcname,"niik_image_convert_to_color_image");
  for(i=0; i<img->nvox; i++) {
    d=160.0*(niik_image_get_voxel(img,i)-imin)/(imax-imin);
    n=floor((niik_image_get_voxel(map,i)-map_min)/step+0.5);
    if(n<0) n=0;
    else if(n>=cmap->row) n=cmap->row-1;
    j=i;
    for(m=0; m<3; m++) {
      rgb[m] = cmap->m[n][m] * d + d;
    }
    if(maskimg!=NULL) {
      if(niik_image_get_voxel(maskimg,i)==0) {
        for(m=0; m<3; m++) rgb[m]*=0.2;
      }
    }
    niik_image_set_voxel(outimg,j,rgb[0]);
    j+=img->nvox;
    niik_image_set_voxel(outimg,j,rgb[1]);
    j+=img->nvox;
    niik_image_set_voxel(outimg,j,rgb[2]);
  } /* each voxel */

  if(verbose>0) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_map_color_overlay_map */

nifti_image *niik_image_map_color(nifti_image *img,double imin,double imax,niikmat *cmap) {
  char fcname[48]="niik_image_map_color";
  int const verbose=1;
  nifti_image *outimg=NULL;
  int i,j,n;
  double d,step;
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((cmap==NULL),fcname,"cmap is null");
  step=(imax-imin)/(cmap->row-1.0);
  if(verbose>0) {
    niik_fc_display(fcname,1);
    fprintf(stdout,"[%s] min  = %f\n",fcname,imin);
    fprintf(stdout,"[%s] max  = %f\n",fcname,imax);
    fprintf(stdout,"[%s] step = %f\n",fcname,step);
  }
  NIIK_RET0(((outimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  NIIK_RET0((!niik_image_convert_to_color_image(outimg)),fcname,"niik_image_convert_to_color_image");
  for(i=0; i<img->nvox; i++) {
    d=niik_image_get_voxel(img,i);
    n=floor((d-imin)/step+0.5);
    //fprintf(stdout,"%i --> %i %f\n",i,n,d);
    if(n<0) n=0;
    else if(n>=cmap->row) n=cmap->row-1;
    j=i;
    niik_image_set_voxel(outimg,j,cmap->m[n][0]);
    j+=img->nvox;
    niik_image_set_voxel(outimg,j,cmap->m[n][1]);
    j+=img->nvox;
    niik_image_set_voxel(outimg,j,cmap->m[n][2]);
  } /* each voxel */
  if(verbose>0) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_map_color */


/***************************************************
 *
 * BASIC NIFTI_IMAGE OPERATIONS
 *
 * image type converters
 *   int niik_image_type_convert(nifti_image * img,int datatype);
 *   int niik_image_type_convert_scl(nifti_image * img,int datatype,int scl_flag);
 *
 * image copy
 *   nifti_image *niik_image_copy(nifti_image * src);
 *
 ***************************************************/

int niik_image_type_convert(nifti_image * img,int datatype) {
  const char* fcname="niik_image_type_convert";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(!niik_image_type_convert_scl(img,datatype,0)) {
    fprintf(stderr,"ERROR: niik_image_type_convert_scl\n");
    return 0;
  }
  return 1;
}

int niik_image_type_convert_scl(nifti_image * img,int datatype,int scl_flag)
/* converts the image type
 * if scl_flag is on, then applies the intensity scaling from scl_slope and scl_interp */
{
  double *dimg;
  const char *fcname="niik_image_type_convert_scl";
  int i;
  int const verbose=0;
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) start\n");
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) {
    fprintf(stderr,"ERROR: complex or rgb are not supported \n");
    return 0;
  }
  /* skip if the new type is the same as old one */
  if(scl_flag==0) if(datatype==img->datatype) return 1;
  if(verbose) {
    fprintf(stdout,"-v (niik_image_type_convert_scl) checked new/old types\n");
    fprintf(stdout,"-v (niik_image_type_convert_scl) old = %i %s\n",img->datatype,nifti_datatype_string(img->datatype));
    fprintf(stdout,"-v (niik_image_type_convert_scl) new = %i %s\n",datatype,nifti_datatype_string(datatype));
  }
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) copy image\n");
  if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  /* update the data type */
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) update datatype\n");
  img->datatype=datatype;
  nifti_datatype_sizes(img->datatype,&(img->nbyper),&(img->swapsize));
  /* don't scale if scl_slope is zero */
  if(img->scl_slope==0.0) {
    scl_flag = 0;
  }
  /* create a new memory space */
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) allocate memory\n");
  free(img->data);
  if((img->data=(void *)malloc(img->nvox*img->nbyper))==NULL) {
    fprintf(stderr,"ERROR: malloc (img->data)\n");
    return 0;
  }
  /* update each voxel intensity */
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) update voxel intensity\n");
  if(scl_flag) { /* scaling is used */
    if(fabs(img->scl_slope)>=1e-5) {
      for(i=0; i<img->nvox; i++) {
        dimg[i] = dimg[i] * img->scl_slope + img->scl_inter;
      }
      fprintf(stdout,"[%s] reset scale slope = %.5f   intercept =  %.5f  %s\n",fcname,img->scl_slope,img->scl_inter,img->fname);
      img->scl_slope = img->scl_inter = 0;
    }
  }
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) set voxels\n");
  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(dimg);
  if(verbose) fprintf(stdout,"-v (niik_image_type_convert_scl) done\n");
  return 1;
}

int niik_image_auto_type_convert(nifti_image * img)
/* automatically convert image to smaller size */
{
  double *dimg,dmin,dmax;
  int i,n,datatype=-1;
  char fcname[32]="niik_image_auto_type_convert";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if((dimg = niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  dmax = dimg[niik_get_max_index_double_vector(dimg,img->nvox)];
  dmin = dimg[niik_get_min_index_double_vector(dimg,img->nvox)];
  for(i=0,n=1; i<img->nvox; i++) {
    if(dimg[i]!=floor(dimg[i])) n=0;
  }
  if(n) {
    if(dmin<0) { /* signed type */
      dmax = NIIK_DMAX(fabs(dmax),fabs(dmin));
      if(dmax>NIIK_VAL_LMAX) {}
      else if(dmax>NIIK_VAL_IMAX) { /* long int type, this is 64-bit */
      } else if(dmax>NIIK_VAL_SMAX) { /* int type, this is 32-bit */
        datatype=NIFTI_TYPE_INT32;
      } else if(dmax>NIIK_VAL_CMAX) { /* short type, this is 16-bit */
        datatype=NIFTI_TYPE_INT16;
      } else { /* char type, this is 16-bit */
        datatype=NIFTI_TYPE_INT8;
      }
    } else { /* unsigned type */
      if(dmax>NIIK_VAL_ULMAX) {}
      else if(dmax>NIIK_VAL_UIMAX) { /* unsigned long int type, this is 64-bit */
      } else if(dmax>NIIK_VAL_USMAX) { /* unsigned int type, this is 32-bit */
        datatype=NIFTI_TYPE_UINT32;
      } else if(dmax>NIIK_VAL_UCMAX) { /* unsigned short type, this is 16-bit */
        datatype=NIFTI_TYPE_UINT16;
      } else { /* unsigned char type, this is 16-bit */
        datatype=NIFTI_TYPE_UINT8;
      }
    }
  }
  if(datatype<0) return 1;
  if(!niik_image_type_convert_scl(img,datatype,0)) {
    fprintf(stderr,"ERROR: niik_image_type_convert_scl\n");
    return 0;
  }
  return 1;
}


nifti_image *niik_image_copy(nifti_image * src)
/* copies the image */
{
  nifti_image *dest = NULL;
  int const verbose=0;
  if(src==NULL)       {
    fprintf(stderr,"ERROR: src is a null pointer\n");
    return NULL;
  }
  if(src->data==NULL) {
    fprintf(stderr,"ERROR: src->data is a null pointer\n");
    return NULL;
  }

  if(verbose) fprintf(stdout,"-v (niik_image_copy) nifti_copy_nim_info: %s\n",src->fname);

  if((dest=nifti_copy_nim_info(src))==NULL) {
    fprintf(stderr,"ERROR: nifti_copy_nim_info \n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_copy) check src/dest\n");
  if(src==dest) {
    fprintf(stderr,"ERROR: src == dest \n");
    return NULL;
  }

  if(verbose) fprintf(stdout,"-v (niik_image_copy) check nvox\n");
  if(dest->nvox != src->nvox) {
    fprintf(stderr,"ERROR: nvox did not match \n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_copy) calloc data\n");
  if((dest->data = (void *)calloc(dest->nvox,dest->nbyper))==NULL) {
    fprintf(stderr,"ERROR: malloc for dest (%s) \n",dest->iname);
    return NULL;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_copy) %i %i  -->   %i %i \n",src->nvox,src->nbyper,dest->nvox,dest->nbyper);
  if(verbose) fprintf(stdout,"-v (niik_image_copy) copy data\n");
  memcpy(dest->data, src->data, (size_t)(dest->nvox * dest->nbyper));
  if(verbose) fprintf(stdout,"-v (niik_image_copy) return data\n");
  return dest;
}

nifti_image *niik_image_copy_as_type(nifti_image *img,int datatype) {
  nifti_image *outimg=NULL;
  const char *fcname="niik_image_copy_as_type";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }

  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_type_convert(outimg,datatype)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  return outimg;
}

nifti_image **niik_image_array(int num) {
  nifti_image **img=NULL;
  img=(nifti_image **)calloc(num,sizeof(nifti_image *));
  return img;
}

nifti_image *niik_image_init(int x,int y,int z,int t,int u,int v,int w,
                             double dx,double dy,double dz,double dt,double du,double dv,double dw,int datatype) {
  nifti_image *img=NULL;
  char fcname[32]="niik_image_init";
  int
  m,n;
  int const verbose=0;
  if(verbose>0) {
    fprintf(stdout,"[%s] %i %i %i %i %i %i %i\n",fcname,x,y,z,t,u,v,w);
    fprintf(stdout,"[%s] %f %f %f %f %f %f %f\n",fcname,dx,dy,dz,dt,du,dv,dw);
  }
  img = (nifti_image *)calloc(1,sizeof(nifti_image));
  img->qform_code = img->sform_code = NIFTI_XFORM_UNKNOWN;
  for(m=0; m<4; m++) for(n=0; n<4; n++) {
      img->qto_xyz.m[m][n]=(m==n);
      img->sto_xyz.m[m][n]=(m==n);
    }
  img->nx = img->dim[1] = x;
  img->ny = img->dim[2] = y;
  img->nz = img->dim[3] = z;
  img->nt = img->dim[4] = t;
  img->nu = img->dim[5] = u;
  img->nv = img->dim[6] = v;
  img->nw = img->dim[7] = w;
  img->dx = img->pixdim[1] = dx;
  img->dy = img->pixdim[2] = dy;
  img->dz = img->pixdim[3] = dz;
  img->dt = img->pixdim[4] = dt;
  img->du = img->pixdim[5] = du;
  img->dv = img->pixdim[6] = dv;
  img->dw = img->pixdim[7] = dw;
  for(img->ndim=7; img->ndim>0; img->ndim--) {
    if(img->dim[img->ndim]>1)
      break;
  }
  img->dim[0] = img->ndim;
  if(verbose>0) fprintf(stdout,"[%s] dim = %i\n",fcname,img->ndim);
  for(n=img->nvox=1; n<=img->ndim; n++) {
    img->nvox*=img->dim[n];
  }
  if(verbose>0) fprintf(stdout,"[%s] nvox = %i\n",fcname,img->nvox);
  img->datatype=datatype;
  nifti_datatype_sizes(datatype,&img->nbyper,&img->swapsize);
  if(verbose>0) fprintf(stdout,"[%s] datatype = %i %i %s\n",fcname,img->nbyper,img->swapsize,nifti_datatype_string(img->datatype));
  if((img->data=(void *)calloc(img->nvox,img->nbyper))==NULL) {
    fprintf(stderr,"ERROR: calloc for img->data\n");
    return NULL;
  }
  if(verbose>0) fprintf(stdout,"[%s] finish\n",fcname);
  return img;
} /* niik_image_init */

int niik_image_iscale(nifti_image *img,double imin,double imax,double omin,double omax) {
  double d,df;
  int i;
  char fcname[32]="niik_image_iscale";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,NULL);
  if(niik_check_double_problem(imax)) imax=niik_image_get_max(img,NULL);
  if(niik_check_double_problem(omin)) omin=0;
  if(niik_check_double_problem(omax)) omin=140;
  fprintf(stdout,"\tiscale %8.3f %8.3f -> %8.3f %8.3f\n",imin,imax,omin,omax);
  df=(omax-omin)/(imax-imin);
  for(i=0; i<img->nvox; i++) {
    d=niik_image_get_voxel(img,i);
    d=(d-imin)*df+omin;
    niik_image_set_voxel(img,i,d);
  }
  /*img->minc_history = NULL;*/
  return 1;
}

/***************************************************
 *
 * BASIC NIFTI_IMAGE OPERATIONS
 *
 ***************************************************/

char *niik_nifti1_xform_string( int code ) {
  switch( code ) {
  case NIFTI_XFORM_UNKNOWN:
    return "nifti_xform_unknown";
  case NIFTI_XFORM_SCANNER_ANAT:
    return "nifti_xform_scanner_anat";
  case NIFTI_XFORM_ALIGNED_ANAT:
    return "nifti_xform_aligned_anat";
  case NIFTI_XFORM_TALAIRACH:
    return "nifti_xform_talairach";
  case NIFTI_XFORM_MNI_152:
    return "nifti_xform_mni_152";
  }
  return "unknown";
}


char *niik_image_dim_string( int code ) {
  switch( code ) {
  case 0:
    return "unknown";
  case 1:
    return "x-dim";
  case 2:
    return "y-dim";
  case 3:
    return "z-dim";
  case 4:
    return "t-dim";
  case 5:
    return "u-dim";
  case 6:
    return "v-dim";
  case 7:
    return "w-dim";
  }
  return "unknown";
}


int niik_image_cmp_dim(nifti_image *img1,nifti_image *img2)
/* -returns zero for equal dimensions
 * -returns negative number for errors
 * -returns the dimension index where the dim is different
 */
{
  int n;
  char fcname[32]="niik_image_cmp_dim";
  if(img1==NULL) {
    fprintf(stderr,"[%s] ERROR: img1 is null\n",fcname);
    return -1;
  }
  if(img2==NULL) {
    fprintf(stderr,"[%s] ERROR: img2 is null\n",fcname);
    return -1;
  }
  if(img1->ndim!=img2->ndim) {
    fprintf(stderr,"[%s] ERROR: ndim is diff %i %i\n",fcname,img1->ndim,img2->ndim);
    return -2;
  }
  for(n=1; n<img1->ndim; n++)
    if(img1->dim[n]!=img2->dim[n]) {
      fprintf(stderr,"[%s] ERROR: dim[%i] is diff %i %i\n",fcname,n,img1->dim[n],img2->dim[n]);
      fprintf(stderr,"       img %-5i %s\n",img1->dim[n],img1->fname);
      fprintf(stderr,"       img %-5i %s\n",img2->dim[n],img2->fname);
      return n;
    }
  return 0;
}

double niik_image_get_voxel_size(nifti_image *img) {
  char fcname[32]="niik_image_get_voxel_size";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NIIKMAX;
  }
  switch(img->ndim) {
  case 2:
    return img->dx*img->dy;
  case 3:
    return img->dx*img->dy*img->dz;
  }
  return NIIKMAX;
}

double niik_image_get_mask_vol(nifti_image *maskimg) {
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return NIIKMAX;
  }
  return niik_image_get_voxel_size(maskimg) * niik_image_count_mask(maskimg) * 0.001;
}

int niik_image_get_index(nifti_image *img,int x,int y,int z)
/* if error, returns -1 */
{
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return -1;
  }
  if(x<0) return -1;
  if(y<0) return -1;
  if(z<0) return -1;
  if(x>=img->dim[1]) return -1;
  if(y>=img->dim[2]) return -1;
  if(z>=img->dim[3]) return -1;
  return x + y*img->dim[1] + z*img->dim[1]*img->dim[2];
}

niikpt niik_image_get_pt_from_index(nifti_image *img,int idx,int *x,int *y,int *z) {
  niikpt p,pw;
  *x = idx%img->nx;
  *y = (idx/img->nx) % img->ny;
  *z = idx/img->nx/img->ny;

  p.x=(*x);
  p.y=(*y);
  p.z=(*z);

  p.w=1;
  niik_ijk_to_world(img,&p,&pw);
  return pw;
} /* niik_image_get_pt_from_index */


int niik_image_get_index_niikpt(nifti_image *img,niikpt pw) {
  int n;
  niikpt p;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  niik_world_to_ijk(img,&pw,&p);

  n=niik_image_get_index(img,
                         (int)floor(p.x+0.5),
                         (int)floor(p.y+0.5),
                         (int)floor(p.z+0.5));
  if(n<0) {
    fprintf(stderr,"ERROR: niik_image_get_index_niikpt %7.3f %7.3f %7.3f\n",pw.x,pw.y,pw.z);
    return -1;
  }
  return n;
}

int niik_image_check_index_niikpt(nifti_image *img,niikpt pw)
/* -same as niik_image_get_index_niikpt
 *  but does not display the error message
 * -returns the index position if valid
 * -returns -1 if outside the field of view
 * -returns -2 if error
 */
{
  int n;
  niikpt p;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return -2;
  }
  niik_world_to_ijk(img,&pw,&p);

  n=niik_image_get_index(img,
                         (int)floor(p.x+0.5),
                         (int)floor(p.y+0.5),
                         (int)floor(p.z+0.5));
  if(n<0) {
    /*fprintf(stderr,"ERROR: niik_image_get_index_niikpt %7.3f %7.3f %7.3f\n",p.x,p.y,p.z); */
    return -1;
  }
  return n;
}

int niik_image_count_mask(nifti_image *maskimg)
/* count the number of voxels and returns it */
{
  unsigned char *bimg;
  int i,n;
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return -1;
  }
  if((bimg = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
    return 0;
  }
  for(i=n=0; i<maskimg->nvox; i++) {
    n+=(bimg[i]>0);
  }
  free(bimg);
  return n;
}

int niik_image_write(const char *fname,nifti_image *img) {
  const char *fcname="niik_image_write";
  int const verbose=niik_verbose();
  if(verbose>1) fprintf(stdout,"[%s] writing output  %s\n",fcname,fname);
  NIIK_RET0((img  ==NULL),fcname,"img is a null pointer");
  NIIK_RET0((fname==NULL),fcname,"fname is a null pointer");
  if(!strncmp(fname+(strlen(fname)-7),".nii.gz",7)) {
    if(verbose>2) {
      fprintf(stdout,"[%s] writing nii.gz  %s\n",fcname,fname);
    }
    if(nifti_set_filenames(img,fname,0,img->byteorder)) {
      fprintf(stderr,"[%s] ERROR: nifti_set_filenames %s\n",fcname,fname);
      return 0;
    }
    if(verbose>2) {
      fprintf(stdout,"[%s] write\n",fcname);
    }
    nifti_image_write(img);
    if(nifti_io_get_problem_code_kunio()>0) {
      fprintf(stderr,"[%s] ERROR: nifti_image_write\n",fcname);
      return 0;
    }
  }
#if NIIK_MINC_SUPPORT==TRUE
  else if(  !strncmp(fname+(strlen(fname)-4),".mnc",4) ) {
    if(verbose>2) {
      fprintf(stdout,"[%s] writing minc    %s\n",fcname,fname);
    }
    if(!niik_image_write_minc2((char *)fname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write_minc %s\n",fcname,fname);
      return 0;
    }
  }
#endif
  else if(!strncmp(fname+(strlen(fname)-4),".nii",4)) {
    if(nifti_set_filenames(img,fname,0,img->byteorder)) {
      fprintf(stderr,"[%s] ERROR: nifti_set_filenames %s\n",fcname,fname);
      return 0;
    }
    nifti_image_write(img);
    if(nifti_io_get_problem_code_kunio()>0) {
      fprintf(stderr,"[%s] ERROR: nifti_image_write\n",fcname);
      return 0;
    }
  } else {
    fprintf(stderr,"[%s] ERROR: unknown output image format\n",fcname);
    return 0;
  }
  return 1;
} /* niik_image_write */

nifti_image **niik_image_read_multiple_from_file(char *fname,int *num) {
  char fcname[64]="niik_image_read_multiple_from_file";
  FILE *fp=NULL;
  nifti_image **imglist=NULL;
  char *cptr,filelist[65536];
  int i=0,ii;
  int const verbose=niik_verbose();
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  for(i=0; i<65536; i++) filelist[i]=0;
  NIIK_RET0(((fp=fopen(fname,"r"))==NULL),fcname,"fopen file");
  i=0;
  while(!feof(fp)) {
    filelist[i++]=fgetc(fp);
    /* fprintf(stdout,"%i %c\n",i,filelist[i-1]); */
  }
  fclose(fp);
  ii=i;
  /* fprintf(stdout,"i = %i\n",i); */
  for(i=0; i<ii; i++) {
    // fprintf(stdout,"---%c\n",filelist[i]);
    if(filelist[i]=='\n') filelist[i]=',';
  }
  cptr=strrchr(filelist,',');
  *cptr=0;
  // fprintf(stdout,"+++%s\n",filelist);
  NIIK_RET0(((imglist=niik_image_read_multiple(filelist,num))==NULL),fcname,"niik_image_read_multiple");
  return imglist;
}

nifti_image **niik_image_read_multiple(char *fname,int *num) {
  char fcname[64]="niik_image_read_multiple";
  int n,N=-1;
  nifti_image **imglist=NULL;
  char **CSlist=NULL;
  if(fname==NULL) {
    fprintf(stderr,"[%s] ERROR: fname is null\n",fcname);
    return NULL;
  }
  if((CSlist=niik_csstring_to_list(fname,&N))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_csstring_to_list\n",fcname);
    return NULL;
  }
  fprintf(stdout,"[%s] # imglist = %i\n",fcname,N);
  if((imglist = (nifti_image **)calloc(N,sizeof(nifti_image *)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc\n",fcname);
    return NULL;
  }
  for(n=0; n<N; n++) {
    fprintf(stdout,"[%s] reading image:    %s\n",fcname,CSlist[n]);
    if((imglist[n]=niik_image_read(CSlist[n]))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,CSlist[n]);
      exit(0);
    }
    free(CSlist[n]);
  }
  free(CSlist);
  *num = N;
  return imglist;
} /* niik_image_read_multiple */


nifti_image *niik_image_read(const char *fname) {
  nifti_image *img=NULL;
  char fcname[64]="niik_image_read";
  /*fprintf(stdout,"[%s] reading image  %s\n",fcname,fname);*/
  if(!strncmp(fname+(strlen(fname)-7),".nii.gz",7)) {
    if((img=nifti_image_read(fname,1))==NULL) {
      fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
      return NULL;
    }
  }
#if NIIK_MINC_SUPPORT==TRUE
  else if(!strncmp(fname+(strlen(fname)-7),".mnc.gz",7)) {
    if((img=niik_image_read_minc2(fname))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_read_minc %s\n",fcname,fname);
      return NULL;
    }
  } else if(!strncmp(fname+(strlen(fname)-4),".mnc",4)) {
    if((img=niik_image_read_minc2(fname))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_read_minc %s\n",fcname,fname);
      return NULL;
    }
  }
#endif
  else if(!strncmp(fname+(strlen(fname)-4),".nii",4)) {
    if((img=nifti_image_read(fname,1))==NULL) {
      fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
      return NULL;
    }
  } else {
    fprintf(stderr,"[%s:%i:%s] ERROR: unknown format\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  switch(img->ndim) {
  case 1:
    img->ny=img->dim[2]=1;
  case 2:
    img->nz=img->dim[3]=1;
  case 3:
    img->nt=img->dim[4]=1;
  case 4:
    img->nu=img->dim[5]=1;
  case 5:
    img->nv=img->dim[6]=1;
  case 6:
    img->nw=img->dim[7]=1;
    break;
  default:
    fprintf(stdout,"[%s] ERROR: unknown dimensions, %i\n",fcname,img->ndim);
    return NULL;
  }
  return img;
}

nifti_image *niik_image_free(nifti_image *img) {
  if(img==NULL) return NULL;
  nifti_image_free(img);
  return NULL;
}

int niik_image_clear(nifti_image *img)
/* clear the image = set all to zero */
{
  int i;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_clear] ERROR: img is null\n");
    return 0;
  }
  for(i=0; i<img->nvox; i++) {
    if(!niik_image_set_voxel(img,i,0)) {
      fprintf(stderr,"ERROR: niik_image_set_voxel\n");
      return 0;
    }
  }
  return 1;
}

int niik_image_one(nifti_image *img)
/* set all voxels to one */
{
  int i;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_one] ERROR: img is null\n");
    return 0;
  }
  for(i=0; i<img->nvox; i++) {
    if(!niik_image_set_voxel(img,i,1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxel\n");
      return 0;
    }
  }
  return 1;
}

int niik_image_calc_nvox3d(nifti_image *img) {
  int nvox=0;
  char fcname[32]="niik_image_update_nvox";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return -1;
  }
  nvox = img->nx;
  if(img->ndim>1) nvox*=img->ny;
  if(img->ndim>2) nvox*=img->nz;
  return nvox;
}

int niik_image_calc_nvox(nifti_image *img) {
  int nvox=0;
  char fcname[32]="niik_image_update_nvox";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return -1;
  }
  nvox = img->nx;
  if(img->ndim>1) nvox*=img->ny;
  if(img->ndim>2) nvox*=img->nz;
  if(img->ndim>3) nvox*=img->nt;
  if(img->ndim>4) nvox*=img->nu;
  if(img->ndim>5) nvox*=img->nv;
  if(img->ndim>6) nvox*=img->nw;
  return nvox;
}

int niik_image_val(nifti_image *img,double val)
/* set all voxels to val */
{
  int i;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_val] ERROR: img is null\n");
    return 0;
  }
  for(i=0; i<img->nvox; i++) {
    if(!niik_image_set_voxel(img,i,val)) {
      fprintf(stderr,"ERROR: niik_image_set_voxel\n");
      return 0;
    }
  }
  return 1;
}

int niik_image_flip_bgfg(nifti_image *maskimg) {
  int i;
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: input maskimg is null\n");
    return 0;
  }
  for(i=0; i<maskimg->nvox; i++) {
    if(niik_image_voxel_get(maskimg,i)>0)
      niik_image_set_voxel(maskimg,i,0);
    else
      niik_image_set_voxel(maskimg,i,1);
  }
  return 1;
}

int niik_image_gaussian_noise(nifti_image *img,double FWHM, double mean) {
  int i;
  fprintf(stdout,"[niik_image_gaussian_noise] FWHM %8.3f   mean %8.3f\n",FWHM,mean);
  if(img==NULL) {
    fprintf(stderr,"[niik_image_gaussian_noise] ERROR: img is null\n");
    return 0;
  }
  if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"error: niik_image_gaussian_noise\n");
    return 0;
  }
  for(i=0; i<img->nvox; i++) {
    if(!niik_image_set_voxel(img,i,mean+FWHM*NIIK_BoxMuller((rand()%NIIK_VAL_USMAX)/(double)NIIK_VAL_USMAX,
                             (rand()%NIIK_VAL_USMAX)/(double)NIIK_VAL_USMAX))) {
      fprintf(stderr,"ERROR: niik_image_set_voxel\n");
      return 0;
    }
  }
  return 1;
}

int niik_image_mask(nifti_image *img,nifti_image *maskimg)
/* mask img with maskimg */
{
  int i;
  NIIK_RET0((img==NULL),__func__,"img is null");
  NIIK_RET0((maskimg==NULL),__func__,"maskimg is a null pointer");
  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(maskimg,i)>0) {
      continue;
    }
    if(!niik_image_set_voxel(img,i,0)) {
      fprintf(stderr,"ERROR: niik_image_set_voxel\n");
      return 0;
    }
  }
  return 1;
}

int niik_image_maskout(nifti_image *img,nifti_image *maskimg)
/* mask-out img with maskimg */
{
  int i;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_maskout] ERROR: img is null\n");
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is a null pointer \n");
    return 0;
  }
  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(maskimg,i)>0) {
      if(!niik_image_set_voxel(img,i,0)) {
        fprintf(stderr,"ERROR: niik_image_set_voxel\n");
        return 0;
      }
    }
  }
  return 1;
}

int niik_image_add_masks(nifti_image *modified,nifti_image *unchange)
/* add 2 masks:
 * -modified is the modified mask
 * -unchange is the unchanged mask */
{
  int i;
  char fcname[32]="niik_image_add_masks";
  if(modified==NULL) {
    fprintf(stderr,"[%s] ERROR: modified mask is null\n",fcname);
    return 0;
  }
  if(unchange==NULL) {
    fprintf(stderr,"[%s] ERROR: unchanged mask is null\n",fcname);
    return 0;
  }
  if(modified->nvox!=unchange->nvox) {
    fprintf(stderr,"[%s] ERROR: nvox is different %i %i\n",fcname,modified->nvox,unchange->nvox);
    return 0;
  }
  for(i=0; i<modified->nvox; i++) {
    niik_image_add_voxel(modified,i,niik_image_get_voxel(unchange,i));
  }
  return 1;
}

int niik_image_display_stats_mask(nifti_image *img,nifti_image *maskimg) {
  const char *fcname="niik_iamge_display_stats_mask";
  int i,minidx,maxidx;
  double
  wavg,ww,g,
       dmean=0,dstdv=0,derr=0,
       dmin=0,dmax=0,dran=0;
  niikvec *v;
  niikpt ctr,pmin,pmax;
  if( img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: mask is a null pointer\n",fcname);
    return 0;
  }
  if(img->nvox!=maskimg->nvox) {
    fprintf(stderr,"[%s] ERROR: nvox is different %i, %i\n",fcname,img->nvox,maskimg->nvox);
    return 0;
  }
  minidx = niik_image_get_min_index(img,maskimg);
  if(minidx<0) {
    fprintf(stderr,"[%s] ERROR: niik_image_get_min_index\n",fcname);
    return 0;
  }
  maxidx = niik_image_get_max_index(img,maskimg);
  dmin = niik_image_get_voxel(img,minidx);
  dmax = niik_image_get_voxel(img,maxidx);
  dran = dmax - dmin;
  fprintf(stdout,"          sum         %15.9f\n",niik_image_get_sum(img,maskimg));
  fprintf(stdout,"         mean         %15.9f\n",niik_image_get_mean(img,maskimg));
  /* 2013-03-26, Kunio Nakamura, knakamura@mrs.mni.mcgill.ca
   * -added weighted average if the input mask is not binary
   * -it still does everything else the same way... it should really check if the mask is binary...
   */
  if( (maskimg->datatype==NIFTI_TYPE_FLOAT32) ||
      (maskimg->datatype==NIFTI_TYPE_FLOAT32) ) {
    for(i=0,wavg=ww=0; i<maskimg->nvox; i++) {
      g = niik_image_get_voxel(maskimg,i);
      ww += g;
      wavg += niik_image_get_voxel(img,i) * g;
    }
    fprintf(stdout," weighted mean        %15.9f\n",wavg/ww);
  }
  fprintf(stdout,"         stdv         %15.9f\n",sqrt(niik_image_get_var(img,maskimg)));
  fprintf(stdout,"         skew         %15.9f\n",niik_image_get_skew(img,maskimg));
  fprintf(stdout,"         min          %15.9f  at %5i %5i %5i\n",dmin,minidx%img->nx,(minidx/img->nx)%img->ny,minidx/(img->nx*img->ny));
  fprintf(stdout,"         max          %15.9f  at %5i %5i %5i\n",dmax,maxidx%img->nx,(maxidx/img->nx)%img->ny,maxidx/(img->nx*img->ny));
  if((v=niik_image_get_voxels_as_double_vector_within_mask(img,maskimg))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_double_vector_within_mask\n",fcname);
    return 0;
  }
  if(!niikvec_sort(v)) {
    fprintf(stderr,"[%s] ERROR: niikvec_sort\n",fcname);
    return 0;
  }
  fprintf(stdout,"         median       %15.9f\n",niik_get_median_from_sorted_double_vector(v->v,v->num));
  fprintf(stdout,"         25%%          %15.9f\n",v->v[(int)floor(v->num*0.25+0.5)]);
  fprintf(stdout,"         75%%          %15.9f\n",v->v[(int)floor(v->num*0.75+0.5)]);
  v=niikvec_free(v);
  /*num=(dmax-dmin)/niik_image_histogram_optim_bin_size(img,maskimg,2);
    fprintf(stdout,"         mode         %15.9f\n",niik_image_get_mode(img,maskimg,dmin,dmax,num,3));*/
  if((int)dran<1e-2) { /* 2012-06-15, Kunio: check for num (or dran) */
  } else
    fprintf(stdout,"         mode         %15.9f\n",niik_image_get_mode(img,maskimg,dmin,dmax,(int)dran,dran/20));
  fprintf(stdout,"         count        %15i\n",niik_image_count_mask(maskimg));
  if(!niik_image_fit_gaussian(img,maskimg,-1,&dmean,&dstdv,&derr)) {
    fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian\n",fcname);
    return 0;
  }
  fprintf(stdout,"  gaussian fit        %15.9f +/- %15.9f    e = %15.9f\n",dmean,dstdv,derr);
  ctr=niikpt_image_get_centroid(img,maskimg);
  if(niik_check_double_problem(ctr.x)) {
    fprintf(stderr,"[%s] ERROR: niikpt_image_get_centroid\n",fcname);
    return 0;
  }
  fprintf(stdout,"        center        %15.9f %15.9f %15.9f\n",ctr.x,ctr.y,ctr.z);
  if(!niikpt_image_get_max_ijk_position(maskimg,&pmax)) {
    fprintf(stderr,"[%s] ERROR: niikpt_image_get_max_ijk_position\n",fcname);
    return 0;
  }
  if(!niikpt_image_get_min_ijk_position(maskimg,&pmin)) {
    fprintf(stderr,"[%s] ERROR: niikpt_image_get_min_ijk_position\n",fcname);
    return 0;
  }
  fprintf(stdout,"        min position  %5i %5i %5i    %12.6f %12.6f %12.6f\n",(int)pmin.x,(int)pmin.y,(int)pmin.z,
          pmin.x*maskimg->dx,pmin.y*maskimg->dy,pmin.z*maskimg->dz);
  fprintf(stdout,"        max position  %5i %5i %5i    %12.6f %12.6f %12.6f\n",(int)pmax.x,(int)pmax.y,(int)pmax.z,
          pmax.x*maskimg->dx,pmax.y*maskimg->dy,pmax.z*maskimg->dz);
  return 1;
} /* with maskimg */


int niik_image_display_stats(nifti_image *img,nifti_image *maskimg) {
  char fcname[64]="niik_iamge_display_stats";
  niikpt ctr;
  niikvec *v=NULL;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  /* global stats */
  fprintf(stdout,"  global mean         %15.9f\n",niik_image_get_mean(img,NULL));
  fprintf(stdout,"  global stdv         %15.9f\n",sqrt(niik_image_get_var(img,NULL)));
  fprintf(stdout,"  global skew         %15.9f\n",niik_image_get_skew(img,NULL));
  fprintf(stdout,"  global min          %15.9f\n",niik_image_get_min(img,NULL));
  fprintf(stdout,"  global max          %15.9f\n",niik_image_get_max(img,NULL));
  v=(niikvec *) calloc(1,sizeof(niikvec));
  if((v->v=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_double_vector\n",fcname);
    return 0;
  }
  v->num=img->nvox;
  if(!niikvec_sort(v)) {
    fprintf(stderr,"[%s] ERROR: niikvec_sort\n",fcname);
    return 0;
  }
  fprintf(stdout,"  global median       %15.9f\n",niik_get_median_from_sorted_double_vector(v->v,v->num));
  fprintf(stdout,"  global 25%%          %15.9f\n",v->v[(int)floor(v->num*0.25+0.5)]);
  fprintf(stdout,"  global 75%%          %15.9f\n",v->v[(int)floor(v->num*0.75+0.5)]);
  fprintf(stdout,"  global count        %15i\n",img->nvox);
  ctr=niikpt_image_get_centroid(img,NULL);
  if(niik_check_double_problem(ctr.x)) {
    fprintf(stderr,"[%s] ERROR: niikpt_image_get_centroid\n",fcname);
    return 0;
  }
  fprintf(stdout,"  global center       %15.9f %15.9f %15.9f\n",ctr.x,ctr.y,ctr.z);
  v=niikvec_free(v);
  if(maskimg!=NULL) {
    if(!niik_image_display_stats_mask(img,maskimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_display_stats_mask\n",fcname);
      return 0;
    }
  }
  return 1;
} /* niik_image_display_stats */

int niik_image_label_display_stats(nifti_image *img,nifti_image *labelimg,int label) {
  char fcname[64]="niik_iamge_label_display_stats";
  nifti_image *maskimg=NULL;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(labelimg!=NULL) {
    if((maskimg = niik_image_label_to_mask(labelimg,label))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_label_to_mask\n",fcname);
      return 0;
    }
  }
  if(!niik_image_display_stats_mask(img,maskimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_display_stats\n",fcname);
    return 0;
  }
  maskimg=niik_image_free(maskimg);
  return 1;
} /* niik_image_label_display_stats */


int niik_image_dim_display(nifti_image *img,char *s) {
  int n;
  char str[4096];
  if(img==NULL) {
    fprintf(stderr,"[niik_image_dim_display] ERROR: img is null\n");
    return 0;
  }
  if(s==NULL) {
    fprintf(stderr,"[niik_image_dim_display] ERROR: string is null\n");
    return 0;
  }
  sprintf(str,"image dimension %s: ",img->fname);
  for(n=1; n<=img->ndim; n++) {
    sprintf(s,"%s %3i ",str,img->dim[n]);
    sprintf(str,"%s",s);
  }
  return 1;
}

int niik_image_pad_3d(nifti_image *img,int xdim,int ydim,int zdim)
/* image padding with zeros
 * -sform and qform are updated (January 14, 2014)
 */
{
  nifti_image *tmpimg;
  int i,j,k,n;
  char fcname[32]="niik_image_pad_3d";
  niikmat *m=NULL;
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim>3),fcname,"not a 3d image");
  NIIK_RET0(((tmpimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  img->dim[1]=img->nx=img->nx+xdim*2;
  img->dim[2]=img->ny=img->ny+ydim*2;
  img->dim[3]=img->nz=img->nz+zdim*2;
  img->nvox=(img->nx)*(img->ny)*(img->nz);
  free(img->data);
  img->data=(void *)calloc(img->nvox,img->nbyper);
  for(k=0,n=0; k<tmpimg->nz; k++) {
    for(j=0; j<tmpimg->ny; j++) {
      for(i=0; i<tmpimg->nx; n++,i++) {
        niik_image_set_voxel(img,(i+xdim)+(j+ydim)*img->nx+(k+zdim)*img->nx*img->ny,niik_image_get_voxel(tmpimg,n));
      }
    }
  }
  tmpimg=niik_image_free(tmpimg);
  if(img->sform_code>0) {
    NIIK_RET0(((m=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix (img->sto_xyz)");
    fprintf(stdout,"[%s] init s matrix\n",fcname);
    niikmat_display(m);
    niikmat_multiply_mat2_free1(niikmat_translate_matrix(-xdim*img->dx,-ydim*img->dy,-zdim*img->dz),
                                m);
    fprintf(stdout,"[%s] new  s matrix\n",fcname);
    niikmat_display(m);
    img->sto_xyz=niikmat_make_mat44(m);
    niikmat_inverse_update(m);
    img->sto_ijk=niikmat_make_mat44(m);
    m=niikmat_free(m);
  }
  if(img->qform_code>0) {
    NIIK_RET0(((m=niikmat_mat44_matrix(img->qto_xyz))==NULL),fcname,"niikmat_mat44_matrix (img->qto_xyz)");
    niikmat_display(m);
    niikmat_multiply_mat2_free1(niikmat_translate_matrix(-xdim*img->dx,-ydim*img->dy,-zdim*img->dz),
                                m);
    niikmat_display(m);
    img->qto_xyz=niikmat_make_mat44(m);
    niikmat_inverse_update(m);
    img->qto_ijk=niikmat_make_mat44(m);
    m=niikmat_free(m);
  }

  return 1;
} /* niik_image_pad */

int niik_image_crop(nifti_image *img,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax) {
  char fcname[32]="niik_image_crop";
  double *dimg;
  int
  i,j,k,n;
  niikmat
  *afmat,*cropmat;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(img->ndim!=3) {
    fprintf(stderr,"[%s] ERROR: not a 3d image\n",fcname);
    return 0;
  }
  if(xmin>xmax) NIIK_ISWAP(&xmin,&xmax);
  if(ymin>ymax) NIIK_ISWAP(&ymin,&ymax);
  if(zmin>zmax) NIIK_ISWAP(&zmin,&zmax);
  if((dimg = niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  for(k=zmin,n=0; k<=zmax; k++) {
    for(j=ymin; j<=ymax; j++) {
      for(i=xmin; i<=xmax; n++,i++) {
        niik_image_set_voxel(img,n,dimg[i+j*img->nx+k*img->nx*img->ny]);
      }
    }
  }
  free(dimg);
  if(img->sform_code) {
    afmat = niikmat_mat44_matrix(img->sto_xyz);
    fprintf(stdout,"[%s] orig sto_xyz\n",fcname);
    niikmat_display(afmat);
    cropmat = niikmat_translate_matrix(xmin*img->dx,ymin*img->dy,zmin*img->dz);
    niikmat_multiply_mat2_free1(cropmat,afmat);
    fprintf(stdout,"[%s] new sto_xyz\n",fcname);
    niikmat_display(afmat);
    img->sto_xyz=niikmat_make_mat44(afmat);
    niikmat_inverse_update(afmat);
    img->sto_ijk=niikmat_make_mat44(afmat);
    afmat=niikmat_free(afmat);
  }
  if(img->qform_code) {
    afmat = niikmat_mat44_matrix(img->qto_xyz);
    cropmat = niikmat_translate_matrix(xmin*img->dx,ymin*img->dy,zmin*img->dz);
    niikmat_multiply_mat1_free2(afmat,cropmat);
    img->qto_xyz=niikmat_make_mat44(afmat);
    niikmat_inverse_update(afmat);
    img->qto_ijk=niikmat_make_mat44(afmat);
    afmat=niikmat_free(afmat);
  }
  img->nx=img->dim[1]=xmax-xmin+1;
  img->ny=img->dim[2]=ymax-ymin+1;
  img->nz=img->dim[3]=zmax-zmin+1;
  img->nvox=img->nx*img->ny*img->nz;
  return 1;
}

int niik_image_get_mask_crop_roi(nifti_image *maskimg,int *xmin,int *xmax,int *ymin,int *ymax,int *zmin,int *zmax) {
  unsigned char *bimg;
  int i,j,k,n;
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return 0;
  }
  if((bimg=niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
    return 0;
  }
  *xmin=*ymin=*zmin=1e9;
  *xmax=*ymax=*zmax=-1e9;
  for(k=n=0; k<maskimg->nz; k++) {
    for(j=0; j<maskimg->ny; j++) {
      for(i=0; i<maskimg->nx; n++,i++) {
        if(!bimg[n]) continue;
        if(*xmin>i) *xmin=i;
        if(*xmax<i) *xmax=i;
        if(*ymin>j) *ymin=j;
        if(*ymax<j) *ymax=j;
        if(*zmin>k) *zmin=k;
        if(*zmax<k) *zmax=k;
      }
    }
  }
  free(bimg);
  return 1;
}


/***************************************************
 *
 * VOXEL INFO
 *
 *
 *
 ***************************************************/

int niik_image_copy_data(nifti_image *src,nifti_image *dest)
/* copies image data from src to dest */
{
  double *dimg;
  if(src==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  if(dest==NULL) {
    fprintf(stderr,"ERROR: dest is null\n");
    return 0;
  }
  if(src->nvox!=dest->nvox) {
    fprintf(stderr,"ERROR: #voxels did not match %i %i\n",src->nvox,dest->nvox);
    return 0;
  }
  if((dimg = niik_image_get_voxels_as_double_vector(src))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  if(!niik_image_set_voxels_from_double_vector(dest,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(dimg);
  return 1;
}

int niik_image_copy_ref_info(nifti_image *src,nifti_image *dest)
/* copy some header info
 * -useful for resampling functions
 */
{
  if(dest==NULL) {
    fprintf(stderr,"ERROR: dest is null\n");
    return 0;
  }
  if(src==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  dest->ndim = dest->dim[0] = src->ndim;
  dest->nx = dest->dim[1] = src->nx;
  dest->ny = dest->dim[2] = src->ny;
  dest->nz = dest->dim[3] = src->nz;
  dest->dx = dest->pixdim[1] = src->dx;
  dest->dy = dest->pixdim[2] = src->dy;
  dest->dz = dest->pixdim[3] = src->dz;
  dest->nvox = dest->nx*dest->ny*dest->nz*dest->nt*dest->nu*dest->nv*dest->nw;
  dest->sform_code = src->sform_code;
  dest->qform_code = src->qform_code;
  /*dest-> freq_dim = src-> freq_dim;
    dest->phase_dim = src->phase_dim;
    dest->slice_dim = src->slice_dim;
    dest->slice_start = src->slice_start;
    dest->slice_end = src->slice_end;
    dest->slice_duration = src->slice_duration;*/
  if(dest->qform_code) {
    dest->quatern_b = src->quatern_b;
    dest->quatern_c = src->quatern_c;
    dest->quatern_d = src->quatern_d;
    dest->qoffset_x = src->qoffset_x;
    dest->qoffset_x = src->qoffset_y;
    dest->qoffset_z = src->qoffset_z;
    dest->qfac      = src->qfac;
    dest->qto_xyz   = src->qto_xyz;
    dest->qto_ijk   = src->qto_ijk;
  }
  if(dest->sform_code) {
    dest->sto_xyz   = src->sto_xyz;
    dest->sto_ijk   = src->sto_ijk;
  }
  dest->toffset    = src->toffset;
  dest->xyz_units  = src->xyz_units;
  dest->time_units = src->time_units;
  return 1;
}

double *niik_image_get_voxels_as_double_vector(nifti_image *img) {
  double *v;
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector] ERROR: img is null\n");
    return NULL;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector] ERROR: img->data is null\n");
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NULL;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) calloc %i\n",img->nvox);
  if((v=(double *)calloc(img->nvox,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)cptr[i];
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)sptr[i];
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)iptr[i];
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)lptr[i];
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ucptr[i];
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)usptr[i];
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)uiptr[i];
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ulptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)fptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)dptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ldptr[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX128:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX256:
    free(v);
    return NULL;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"ERROR: unknown datatype\n");
    free(v);
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) return\n");
  return v;
}

int niik_image_get_voxels_as_double_vector_update(nifti_image *img,double *v) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector_update] ERROR: img is null\n");
    return 0;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector_update] ERROR: img->data is null\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)cptr[i];
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)sptr[i];
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)iptr[i];
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)lptr[i];
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ucptr[i];
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)usptr[i];
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)uiptr[i];
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ulptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)fptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)dptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(double)ldptr[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    free(v);
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    free(v);
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    free(v);
    return 0;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"ERROR: unknown datatype\n");
    free(v);
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) return\n");
  return 1;
}

niikvec *niik_image_get_voxels_as_double_vector_within_mask(nifti_image *img,nifti_image *maskimg) {
  niikvec *v;
  unsigned char *ucptr,*bimg;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,j,verbose=0;
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector_within_mask] ERROR: img is null\n");
    return NULL;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_double_vector_within_mask] ERROR: img->data is null\n");
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector_within_mask) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NULL;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector_within_mask) create niikvec\n");
  if((v=(niikvec *)calloc(1,sizeof(niikvec)))==NULL) {
    fprintf(stderr,"ERROR: niikvec calloc\n");
    return NULL;
  }
  v->num=niik_image_count_mask(maskimg);
  if((v->v=(double *)calloc(v->num,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc v->v\n");
    return NULL;
  }
  if((bimg=niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
    v=niikvec_free(v);
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector_within_mask) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)cptr[i];
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)sptr[i];
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)iptr[i];
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)lptr[i];
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)ucptr[i];
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)usptr[i];
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)uiptr[i];
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)ulptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)fptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)dptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=j=0; i<img->nvox; i++) {
      if(bimg[i]) v->v[j++]=(double)ldptr[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    niikvec_free(v);
    free(bimg);
    return NULL;
  case NIFTI_TYPE_COMPLEX128:
    niikvec_free(v);
    free(bimg);
    return NULL;
  case NIFTI_TYPE_COMPLEX256:
    niikvec_free(v);
    free(bimg);
    return NULL;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"ERROR: unknown datatype\n");
    niikvec_free(v);
    free(bimg);
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector_within_mask) return\n");
  return v;
}

float *niik_image_get_voxels_as_float_vector(nifti_image *img) {
  float *v;
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_float_vector] ERROR: img is null\n");
    return NULL;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_get_voxels_as_float_vector] ERROR: img->data is null\n");
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NULL;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) calloc %i\n",img->nvox);
  if((v=(float *)calloc(img->nvox,sizeof(float)))==NULL) {
    fprintf(stderr,"ERROR: calloc for data\n");
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)cptr[i];
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)sptr[i];
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)iptr[i];
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)lptr[i];
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)ucptr[i];
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)usptr[i];
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)uiptr[i];
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)ulptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)fptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)dptr[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(float)ldptr[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX128:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX256:
    free(v);
    return NULL;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"ERROR: unknown datatype\n");
    free(v);
    return NULL;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) return\n");
  return v;
}

unsigned char *niik_image_get_voxels_as_uint8_vector(nifti_image *img) {
  unsigned char *v;
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  char fcname[64]="niik_image_get_voxels_as_uint8_vector";
  NIIK_RET0((img      ==NULL),fcname,"img is nlul");
  NIIK_RET0((img->data==NULL),fcname,"img->data is null");
  if(verbose) niik_fc_display(fcname,1);
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NULL;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) calloc %i\n",img->nvox);
  NIIK_RET0(((v=(unsigned char *)calloc(img->nvox,sizeof(char)))==NULL),fcname,"calloc");
  if(verbose) fprintf(stderr,"[%s] switch\n",fcname);
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(cptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(sptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(iptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(lptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(ucptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(usptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(uiptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(ulptr[i],0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float       *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(fptr[i]+0.5,0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double      *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(dptr[i]+0.5,0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(unsigned char)NIIK_DMINMAX(ldptr[i]+0.5,0,NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX128:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX256:
    free(v);
    return NULL;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"[%s] ERROR: unknown datatype\n",fcname);
    free(v);
    return NULL;
  }
  if(verbose) niik_fc_display(fcname,0);
  return v;
}

long *niik_image_get_voxels_as_long_vector(nifti_image *img) {
  long *v;
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  char fcname[64]="niik_image_get_voxels_as_long_vector";
  NIIK_RET0((img      ==NULL),fcname,"img is nlul");
  NIIK_RET0((img->data==NULL),fcname,"img->data is null");
  if(verbose) niik_fc_display(fcname,1);
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NULL;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxels_as_double_vector) calloc %i\n",img->nvox);
  NIIK_RET0(((v=(long *)calloc(img->nvox,sizeof(long)))==NULL),fcname,"calloc");
  if(verbose) fprintf(stderr,"[%s] switch\n",fcname);
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(cptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(sptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(iptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(lptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(ucptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(usptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(uiptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(ulptr[i],-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float       *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(fptr[i]+0.5,-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double      *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(dptr[i]+0.5,-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      v[i]=(long)NIIK_DMINMAX(ldptr[i]+0.5,-NIIK_VAL_LMAX,NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX128:
    free(v);
    return NULL;
  case NIFTI_TYPE_COMPLEX256:
    free(v);
    return NULL;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
    fprintf(stderr,"[%s] ERROR: unknown datatype\n",fcname);
    free(v);
    return NULL;
  }
  if(verbose) niik_fc_display(fcname,0);
  return v;
}


int niik_image_set_voxels_from_double_vector(nifti_image *img,double *v) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(  v==NULL) {
    fprintf(stderr,"ERROR: v is a null pointer\n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_set_voxels_from_double_vector) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=img->data;
    for(i=0; i<img->nvox; i++) {
      cptr[i]=(char)  NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_CMAX);
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=img->data;
    for(i=0; i<img->nvox; i++) {
      sptr[i]=(short) NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_SMAX);
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=img->data;
    for(i=0; i<img->nvox; i++) {
      iptr[i]=(int)   NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_IMAX);
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=img->data;
    for(i=0; i<img->nvox; i++) {
      lptr[i]=(long)  NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_LMAX);
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=img->data;
    for(i=0; i<img->nvox; i++) {
      ucptr[i]=(unsigned char) NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_UCMAX);
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=img->data;
    for(i=0; i<img->nvox; i++) {
      usptr[i]=(unsigned short)NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_USMAX);
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=img->data;
    for(i=0; i<img->nvox; i++) {
      uiptr[i]=(unsigned int)  NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_UIMAX);
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=img->data;
    for(i=0; i<img->nvox; i++) {
      ulptr[i]=(unsigned long) NIIK_DMINMAX(floor(v[i]+0.5),0.0,(double)NIIK_VAL_ULMAX);
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=img->data;
    for(i=0; i<img->nvox; i++) {
      fptr[i]=(float) v[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=img->data;
    for(i=0; i<img->nvox; i++) {
      dptr[i]=(double)v[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=img->data;
    for(i=0; i<img->nvox; i++) {
      ldptr[i]=(long double) v[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    return 0;
  case NIFTI_TYPE_RGB24:
    return 0;
  case NIFTI_TYPE_RGBA32:
    return 0;
  default:
    fprintf(stderr,"ERROR: unknown datatype\n");
    return 0;
  }
  return 1;
}

int niik_image_set_voxels_from_uint8_vector(nifti_image *img,unsigned char *v) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  int i,verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(  v==NULL) {
    fprintf(stderr,"ERROR: v is a null pointer\n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_set_voxels_from_double_vector) switch\n");
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char*)img->data;
    for(i=0; i<img->nvox; i++) {
      cptr[i]=(char)  v[i];
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short*)img->data;
    for(i=0; i<img->nvox; i++) {
      sptr[i]=(short) v[i];
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int*)img->data;
    for(i=0; i<img->nvox; i++) {
      iptr[i]=(int)   v[i];
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long*)img->data;
    for(i=0; i<img->nvox; i++) {
      lptr[i]=(long)  v[i];
    }
    break;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char *)img->data;
    for(i=0; i<img->nvox; i++) {
      ucptr[i]=(unsigned char) v[i];
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(i=0; i<img->nvox; i++) {
      usptr[i]=(unsigned short)v[i];
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int *)img->data;
    for(i=0; i<img->nvox; i++) {
      uiptr[i]=(unsigned int)  v[i];
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long *)img->data;
    for(i=0; i<img->nvox; i++) {
      ulptr[i]=(unsigned long) v[i];
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float *)img->data;
    for(i=0; i<img->nvox; i++) {
      fptr[i]=(float) v[i];
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(i=0; i<img->nvox; i++) {
      dptr[i]=(double)v[i];
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(i=0; i<img->nvox; i++) {
      ldptr[i]=(long double) v[i];
    }
    break;
  case NIFTI_TYPE_COMPLEX64:
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    return 0;
  case NIFTI_TYPE_RGB24:
    return 0;
  case NIFTI_TYPE_RGBA32:
    return 0;
  default:
    fprintf(stderr,"[%s:%i:%s] ERROR: unknown datatype\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  return 1;
}

int niik_image_set_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val)
/* -set voxel values to be val within ROI [xmin,xmax] ... inclusive
 * -that is, IMG[xmin,ymin,zmin] will be val
 */
{
  int i,j,k,n;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_set_voxels_ROI] ERROR: img is null\n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) {
    fprintf(stderr,"ERROR: complex or rgb are invalid\n");
    return 0;
  }
  /* check for min/max */
  if(xmax<xmin) NIIK_ISWAP(&xmin,&xmax);
  if(ymax<ymin) NIIK_ISWAP(&ymin,&ymax);
  if(zmax<zmin) NIIK_ISWAP(&zmin,&zmax);
  /* check for extremes */
  xmin=NIIK_IMAX(0,xmin);
  ymin=NIIK_IMAX(0,ymin);
  zmin=NIIK_IMAX(0,zmin);
  xmax=NIIK_IMIN(img->dim[1]-1,xmax);
  ymax=NIIK_IMIN(img->dim[2]-1,ymax);
  zmax=NIIK_IMIN(img->dim[3]-1,zmax);
  /* put the values */
  for(k=zmin; k<=zmax; k++) {
    for(j=ymin; j<=ymax; j++) {
      for(i=xmin; i<=xmax; i++) {
        n=i+j*img->nx+k*img->nx*img->ny;
        if(!niik_image_set_voxel(img,n,val)) {
          fprintf(stderr,"ERROR: niik_image_set_voxel(img,%i,%f)\n",n,val);
          return 0;
        }
      }
    }
  }
  return 1;
} /* int niik_image_set_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val) */


int niik_image_add_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val)
/* add value (val) to image voxels defined by min/max */
{
  int i,j,k,n;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_add_voxels_ROI] ERROR: img is null\n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) {
    fprintf(stderr,"ERROR: complex or rgb are invalid\n");
    return 0;
  }
  /* check for min/max */
  if(xmax<xmin) NIIK_ISWAP(&xmin,&xmax);
  if(ymax<ymin) NIIK_ISWAP(&ymin,&ymax);
  if(zmax<zmin) NIIK_ISWAP(&zmin,&zmax);
  /* check for extremes */
  xmin=NIIK_IMAX(0,xmin);
  ymin=NIIK_IMAX(0,ymin);
  zmin=NIIK_IMAX(0,zmin);
  xmax=NIIK_IMIN(img->dim[1]-1,xmax);
  ymax=NIIK_IMIN(img->dim[2]-1,ymax);
  zmax=NIIK_IMIN(img->dim[3]-1,zmax);
  /* add the values */
  for(k=zmin; k<=zmax; k++) {
    for(j=ymin; j<=ymax; j++) {
      for(i=xmin; i<=xmax; i++) {
        n=i+j*img->nx+k*img->nx*img->ny;
        if(!niik_image_add_voxel(img,n,val)) {
          fprintf(stderr,"ERROR: niik_image_set_voxel(img,%i,%f)\n",n,val);
          return 0;
        }
      }
    }
  }
  return 1;
} /* int niik_image_add_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val) */

int niik_image_mul_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val)
/* multiply value (val) to image voxels defined by min/max */
{
  int i,j,k,n;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_mul_voxels_ROI] ERROR: img is null\n");
    return 0;
  }
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGB24 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) {
    fprintf(stderr,"ERROR: complex or rgb are invalid\n");
    return 0;
  }
  /* check for min/max */
  if(xmax<xmin) NIIK_ISWAP(&xmin,&xmax);
  if(ymax<ymin) NIIK_ISWAP(&ymin,&ymax);
  if(zmax<zmin) NIIK_ISWAP(&zmin,&zmax);
  /* check for extremes */
  xmin=NIIK_IMAX(0,xmin);
  ymin=NIIK_IMAX(0,ymin);
  zmin=NIIK_IMAX(0,zmin);
  xmax=NIIK_IMIN(img->dim[1]-1,xmax);
  ymax=NIIK_IMIN(img->dim[2]-1,ymax);
  zmax=NIIK_IMIN(img->dim[3]-1,zmax);
  /* multiply the values */
  for(k=zmin; k<=zmax; k++) {
    for(j=ymin; j<=ymax; j++) {
      for(i=xmin; i<=xmax; i++) {
        n=i+j*img->nx+k*img->nx*img->ny;
        if(!niik_image_mul_voxel(img,n,val)) {
          fprintf(stderr,"ERROR: niik_image_set_voxel(img,%i,%f)\n",n,val);
          return 0;
        }
      }
    }
  }
  return 1;
} /* int niik_image_add_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val) */


nifti_image *niik_image_3d_average_filter(nifti_image *img,int size) {
  char fcname[64]="niik_image_average_filter";
  nifti_image *tmp;
  int
  size21,
  xy,
  i,j,k,n,m,mm,
  verbose=0;
  double *tmp1,*tmp2;

  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((size<1),fcname,"size is too small");
  //NIIK_RET0(((tmp1=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector");
  NIIK_RET0(((tmp1=(double *)calloc(img->nvox,sizeof(double)))==NULL),fcname,"calloc tmp1");
  NIIK_RET0(((tmp2=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector tmp2");
  size21=size*2+1;
  xy=img->nx*img->ny;

  // x-dir
  if(verbose>0) fprintf(stdout,"[%s] x-dir\n",fcname);
  /*#pragma omp parallel for private(i,j,n,m,mm)*/
  for(k=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        n=i+j*img->nx+k*xy;
        tmp1[n]=0;
        for(m=-size; m<=size; m++) {
          mm=i+m;
          if(mm<1)             tmp1[n]+=tmp2[n];
          else if(mm>=img->nx) tmp1[n]+=tmp2[n];
          else                 tmp1[n]+=tmp2[n+m];
        }
        tmp1[n]/=size*2+1;
      }
    }
  }

  // y-dir
  if(verbose>0) fprintf(stdout,"[%s] y-dir\n",fcname);
  /*#pragma omp parallel for private(i,j,n,m,mm)*/
  for(k=0; k<img->nz; k++) {
    for(i=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        n=i+j*img->nx+k*xy;
        tmp2[n]=0;
        for(m=-size; m<=size; m++) {
          mm=j+m;
          if(mm<1)             tmp2[n]+=tmp1[n];
          else if(mm>=img->ny) tmp2[n]+=tmp1[n];
          else                 tmp2[n]+=tmp1[n+m*img->nx];
        }
        tmp2[n]/=size21;
      }
    }
  }

  // z-dir
  if(verbose>0) fprintf(stdout,"[%s] z-dir\n",fcname);
  /*#pragma omp parallel for private(i,j,n,m,mm)*/
  for(k=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        n=i+j*img->nx+k*xy;
        tmp1[n]=0;
        for(m=-size; m<=size; m++) {
          mm=k+m;
          if(mm<1)             tmp1[n]+=tmp2[n];
          else if(mm>=img->nz) tmp1[n]+=tmp2[n];
          else                 tmp1[n]+=tmp2[n+m*xy];
        }
        tmp1[n]/=size21;
      }
    }
  }

  free(tmp2);
  if(verbose>0) fprintf(stdout,"[%s] copy\n",fcname);
  NIIK_RET0(((tmp=niik_image_copy(img))==NULL),fcname,"niik_image_copy, tmp");
  for(i=0; i<img->nvox; i++)
    niik_image_set_voxel(tmp,i,tmp1[i]);
  free(tmp1);

  if(verbose>0) niik_fc_display(fcname,0);
  return tmp;
}

/***************************************************
 *
 * GAUSSIAN FILTER
 *
 * int niik_image_filter_gaussian_update(nifti_image *img,int kdim,double FWHM);
 *   -main function for 3d gaussian filter
 *   -works up to 5 dimensions by the filter is in 3d
 *
 * nifti_image *niik_image_filter_gaussian(nifti_image *img,int kdim,double FWHM);
 *   -wrapper for gaussian filter
 *   -probably too much memory allocation (can be reduced with coding)
 *
 ***************************************************/

nifti_image *niik_image_filter_gaussian(nifti_image *img,int kdim,double FWHM) {
  nifti_image *outimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(outimg,kdim,FWHM)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  return outimg;
}

int niik_image_filter_gaussian_update(nifti_image *img,int kdim,double FWHM)
/* 3d gaussian filter
 * -img is filtered
 * -kdim is the filter kernel size radius is like (kdim-1)/2
 * -FWHM of filter
 */
{
  char fcname[64]="niik_image_filter_gaussian_update";
  double
  *dimg1,
  *dimg,
  *dv,
  dsum,
  sigma;
  niikmat *tmat;
  int
  verbose=0,
  i,j,k,m,n,u,nu,
  ii,jj,kk,
  kmid,
  maxdim,
  area,size;
  /* error checking */
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer \n",fcname);
    return 0;
  }
  if(kdim<=0)   {
    fprintf(stderr,"[%s] ERROR: kdim is invalid %i\n",fcname,kdim);
    return 0;
  }
  if(kdim<0) kdim=-kdim; /* i could use this as a flag */
  if(!(kdim%2)) {
    kdim++;
  }
  if(kdim<=0) {
    fprintf(stderr,"[%s] ERROR: kdim is invalid\n",fcname);
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) start %s\n",img->fname);
  /* prepare variables */
  if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) prepare variables\n");
  if((dv = (double *)calloc(kdim,sizeof(double)))==NULL) {
    fprintf(stderr,"[niik_image_filter_gaussian_update] ERROR: calloc for dv\n");
    return 0;
  }
  kmid = (kdim-1)/2;
  area  = img->nx*img->ny;
  size  = img->nz*area;
  sigma = NIIK_FWHM2SIGMA(FWHM);
  maxdim = NIIK_IMAX(img->nx,NIIK_IMAX(img->ny,img->nz));
  dimg1 = dimg = niik_image_get_voxels_as_double_vector(img);
  tmat=niikmat_init(maxdim,maxdim);
  if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) more than 3d?\n");
  if(img->ndim>5) {
    fprintf(stderr,"ERROR: can't handle more than 5 dimensions\n");
    return 0;
  } else if(img->ndim<=3) nu=1;
  else if(img->ndim==4) nu=img->nt;
  else if(img->ndim==5) nu=img->nu * img->nt;
  else {
    fprintf(stderr,"ERROR: please check the image dimensions %i\n",img->ndim);
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) main outer loop %i\n",nu);
  for(u=0; u<nu; u++,dimg+=size) { /* each 3d image */
    /* x-dir */
    if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) x-dir\n");
    for(n=0,dsum=0; n<kdim; n++) {
      dv[n]=NIIK_GaussPDF((n-kmid)*img->pixdim[1],sigma);
      dsum+=dv[n];
    }
    for(n=0; n<kdim; n++) {
      dv[n]/=dsum;
    }
    /*#pragma omp parallel for private(i,j,n,m,ii)*/
    for(k=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        n = j*img->nx + k*area;
        for(i=0; i<img->nx; n++,i++) {
          tmat->m[k][i] = 0;
          for(m=0; m<kdim; m++) { /* get the appropriate value (edge or non-edge) */
            ii=i+m-kmid;
            if     (ii<0)        tmat->m[k][i] += dv[m] * dimg[n];
            else if(ii>=img->nx) tmat->m[k][i] += dv[m] * dimg[n];
            else                 tmat->m[k][i] += dv[m] * dimg[n+m-kmid];
          }
        } /* x-dir */
        n = j*img->nx + k*area;
        for(i=0; i<img->nx; i++,n++) {
          dimg[n] = tmat->m[k][i];
        }
      }
    } /* along x for y+z plane */
    /* y-dir */
    if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) y-dir\n");
    for(n=0,dsum=0; n<kdim; n++) {
      dv[n]=NIIK_GaussPDF((n-kmid)*img->pixdim[2],sigma);
      dsum+=dv[n];
    }
    for(n=0; n<kdim; n++) {
      dv[n]/=dsum;
    }
    /*#pragma omp parallel for private(i,j,n,m,jj)*/
    for(k=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        n=i+k*area;
        for(j=0; j<img->ny; n+=img->nx,j++) {
          tmat->m[k][j]=0;
          for(m=0; m<kdim; m++) { /* get the appropriate value (edge or non-edge) */
            jj=j+m-kmid;
            if     (jj<0)        tmat->m[k][j] += dv[m] * dimg[n];
            else if(jj>=img->ny) tmat->m[k][j] += dv[m] * dimg[n];
            else                 tmat->m[k][j] += dv[m] * dimg[n+(m-kmid)*img->nx];
          }
        }
        n = i + k*area;
        for(j=0; j<img->ny; j++,n+=img->nx) {
          dimg[n] = tmat->m[k][j];
        }
      }
    }
    /* z-dir */
    if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) z-dir\n");
    for(n=0,dsum=0; n<kdim; n++) {
      dv[n]=NIIK_GaussPDF((n-kmid)*img->pixdim[3],sigma);
      dsum+=dv[n];
    }
    for(n=0; n<kdim; n++) {
      dv[n]/=dsum;
    }

    /*#pragma omp parallel for private(j,k,n,m,kk)*/
    for(i=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        n=i+j*img->nx;
        for(k=0; k<img->nz; n+=area,k++) {
          tmat->m[i][k]=0;
          for(m=0; m<kdim; m++) { /* get the appropriate value (edge or non-edge) */
            kk=k+m-kmid;
            if     (kk<0)        tmat->m[i][k] += dv[m] * dimg[n];
            else if(kk>=img->nz) tmat->m[i][k] += dv[m] * dimg[n];
            else                 tmat->m[i][k] += dv[m] * dimg[n+(m-kmid)*area];
          }
        }
        n=i+j*img->nx;
        for(k=0; k<img->nz; n+=area,k++) {
          dimg[n] = tmat->m[i][k];
        }
      }
    }
  } /* larger than 3 dimensions */
  if(verbose) fprintf(stdout,"-d (niik_image_filter_gaussian_update) set the image values\n");
  if(!niik_image_set_voxels_from_double_vector(img,dimg1)) {
    fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector \n");
    return 0;
  }
  free(dv);
  free(dimg1);
  tmat=niikmat_free(tmat);
  return 1;
}




/***************************************************
 *
 * STATS FUNCTIONS
 *
 *
 * double niik_image_min(nifti_image *img,nifti_image *maskimg);
 * double niik_image_max(nifti_image *img,nifti_image *maskimg);
 * double niik_image_mean(nifti_image *img,nifti_image *maskimg);
 * double niik_image_var(nifti_image *img,nifti_image *maskimg);
 *
 *
 ***************************************************/

double niik_image_get_min(nifti_image *img,nifti_image *maskimg) {
  double out;
  int m,n;
  char fcname[32]="niik_image_get_min";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    out=niik_image_get_voxel(img,0);
    for(n=1; n<img->nvox; n++) {
      out=NIIK_DMIN(out,niik_image_get_voxel(img,n));
    }
  } else {
    if(img->nvox!=maskimg->nvox) {
      fprintf(stderr,"[%s] ERROR: img and maskimg size did not match, %i %i\n",fcname,img->nvox,maskimg->nvox);
      return NIIKMAX;
    }
    for(n=m=0,out=NIIKMAX; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        out=NIIK_DMIN(out,niik_image_get_voxel(img,n));
      }
    }
  }
  return out;
}

double niik_image_get_max(nifti_image *img,nifti_image *maskimg) {
  double out;
  int m,n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    out=niik_image_get_voxel(img,0);
    for(n=1; n<img->nvox; n++) {
      out=NIIK_DMAX(out,niik_image_get_voxel(img,n));
    }
  } else {
    for(n=m=0,out=-NIIKMAX; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        out=NIIK_DMAX(out,niik_image_get_voxel(img,n));
      }
    }
  }
  return out;
}

double niik_image_get_mean_label(nifti_image *img,nifti_image *labelimg,int label) {
  char fcname[32]="niik_image_get_mean_label";
  double out;
  int m,n;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NIIKMAX;
  }
  if(labelimg==NULL) {
    out = niik_image_get_mean(img,NULL);
    if(niik_check_double_problem(out)) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_mean\n",fcname);
      return NIIKMAX;
    }
  } else {
    for(n=m=0,out=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(labelimg,n)==label)  {
        out+=niik_image_get_voxel(img,n);
        m++;
      }
    }
    out/=m;
  }
  return out;
} /* niik_image_get_mean_label */

double niik_image_get_stdv_label(nifti_image *img,nifti_image *labelimg,int label) {
  const char *fcname="niik_image_get_stdv_label";
  double out,mu;
  int m,n;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NIIKMAX;
  }
  if(labelimg==NULL) {
    out = niik_image_get_stdv(img,NULL);
    if(niik_check_double_problem(out)) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_mean\n",fcname);
      return NIIKMAX;
    }
  } else {
    for(n=m=0,mu=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(labelimg,n)==label)  {
        mu+=niik_image_get_voxel(img,n);
        m++;
      }
    }
    mu/=m;
    for(n=0,out=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(labelimg,n)==label)  {
        out+=NIIK_SQ(niik_image_get_voxel(img,n)-mu);
      }
    }
    out=sqrt(out/m);
  }
  return out;
} /* niik_image_get_stdv_label */


double niik_image_get_mean(nifti_image *img,nifti_image *maskimg) {
  double out=0.0;
  int m,n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
#if _OPENMP>=201307
    #pragma omp simd
#endif
    for(n=0; n<img->nvox; n++) {
      out+=niik_image_get_voxel(img,n);
    }
    out/=img->nvox;
  } else {
#if _OPENMP>=201307
    #pragma omp simd
#endif
    for(n=m=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        out+=niik_image_get_voxel(img,n);
        m++;
      }
    }
    out/=m;
  }
  return out;
}

double niik_image_get_stdv(nifti_image *img,nifti_image *maskimg) {
  return sqrt(niik_image_get_var(img,maskimg));
}

double niik_image_get_var(nifti_image *img,nifti_image *maskimg) {
  double dsum,dssq,d;
  int m,n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    for(n=0,dsum=dssq=0; n<img->nvox; n++) {
      d=niik_image_get_voxel(img,n);
      dsum+=d;
      dssq+=d*d;
    }
    return (dssq/img->nvox)-NIIK_SQ(dsum/img->nvox);
  } else {
    for(n=m=0,dsum=dssq=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        m++;
        d=niik_image_get_voxel(img,n);
        dsum+=d;
        dssq+=d*d;
      }
    }
  }
  return dssq/m-NIIK_SQ(dsum/m);
}

double niik_image_get_skew(nifti_image *img, nifti_image *maskimg) {
  const char *fcname="niik_image_get_skew";
  unsigned char *bimg=NULL;
  double skew=0,dmean=0,dstdv=0,d;
  int i,n;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NIIKMAX;
  }
  dmean = niik_image_get_mean(img,maskimg);
  dstdv = sqrt(niik_image_get_var(img,maskimg));
  if(maskimg!=NULL) {
    if((bimg=niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
      fprintf(stderr,"[%s] niik_image_get_voxels_as_uint8_vector\n",fcname);
      return NIIKMAX;
    }
  } else {
    bimg=(unsigned char *)calloc(img->nvox, sizeof(unsigned char));
    for(i=0; i<img->nvox; i++) {
      bimg[i]=255;
    }
  }
  for(i=n=0; i<img->nvox; i++) {
    if(bimg[i]>0) {
      d=(niik_image_get_voxel(img,i) - dmean)/dstdv;
      skew += d*d*d;
      n++;
    }
  }
  free(bimg);
  if(n==0) {
    fprintf(stderr,"[%s] ERROR: missing mask voxels\n",fcname);
    return NIIKMAX;
  }
  skew /= n;
  return skew;
} /* niik_image_get_skew */


double niik_image_get_sum(nifti_image *img,nifti_image *maskimg) {
  double dsum;
  int n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    for(n=0,dsum=0; n<img->nvox; n++) {
      dsum+=niik_image_get_voxel(img,n);
    }
    return dsum;
  } else {
    if(img->nvox!=maskimg->nvox) {
      fprintf(stderr,"[niik_image_get_sum] ERROR: img->nvox and maskimg->nvox are different: %i %i\n",img->nvox,maskimg->nvox);
      return NIIKMAX;
    }
    for(n=0,dsum=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        dsum+=niik_image_get_voxel(img,n);
      }
    }
  }
  return dsum;
}


double niik_image_get_percentile(nifti_image *img,nifti_image *maskimg, double percent) {
  double *dimg,d;
  int m,n,num;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    num=img->nvox;
    if(!niik_sort_double_vector(dimg,num)) {
      fprintf(stderr,"ERROR: niik_sort_double_vector\n");
      return NIIKMAX;
    }
    n=(int)NIIK_DMINMAX(num*percent,0,img->nvox-1);
  } else {
    for(n=m=0; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        dimg[m++]=dimg[n];
      }
    }
    num=m;
    if(!niik_sort_double_vector(dimg,num)) {
      fprintf(stderr,"ERROR: niik_sort_double_vector\n");
      return NIIKMAX;
    }
    n=(int)NIIK_DMINMAX(num*percent,0,num-1);
  }
  d=dimg[n];
  free(dimg);
  return d;
}

nifti_image *niik_image_label_to_mask(nifti_image *labelimg,double label) {
  nifti_image *maskimg=NULL;
  char fcname[53]="niik_image_label_to_mask";
  int i;
  if(labelimg==NULL) {
    fprintf(stderr,"[%s] ERROR: labelimg is a null pointer\n",fcname);
    return NULL;
  }
  if((maskimg=niik_image_copy_as_type(labelimg,NIFTI_TYPE_UINT8))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  niik_image_clear(maskimg);
  for(i=0; i<labelimg->nvox; i++) {
    if(niik_image_get_voxel(labelimg,i)==label)
      niik_image_set_voxel(maskimg,i,1);
  }
  /*niik_fc_display(fcname,0);*/
  return maskimg;
}

double niik_image_get_mode_label(nifti_image *img,nifti_image *labelimg,int label,double dmin,double dmax,int num,int avgnum) {
  nifti_image *maskimg=NULL;
  const char *fcname="niik_image_get_mode_label";
  double out;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  maskimg = niik_image_label_to_mask(labelimg,(double)label);
  out = niik_image_get_mode(img,maskimg,dmin,dmax,num,avgnum);
  maskimg=niik_image_free(maskimg);
  /*niik_fc_display(fcname,0);*/
  return out;
}

double niik_image_get_mode(nifti_image *img,nifti_image *maskimg,double dmin,double dmax,int num,int avgnum) {
  double
  dx,d,
  *hx,*hy;
  int n;
  const int verbose=0;

  const char *fcname="niik_image_get_mode";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(num<=0) {
    fprintf(stderr,"[%s] ERROR: num is invalid %i\n",fcname,num);
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_get_mode) memory alloc + init %i\n",num);
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
  if(verbose) fprintf(stdout,"-d (niik_image_get_mode) histogram %8.3f %8.3f %8.3f %i\n",dmin,dx,dmax,num);
  if(!niik_image_histogram_limits(img,maskimg,hx,hy,num)) {
    fprintf(stderr,"[%s] ERROR: niik_image_histogram\n",fcname);
    return NIIKMAX;
  }
  /* 2012-03-14 Kunio
   * -check if histogram exists */
  for(n=0,d=0.0; n<num; n++) {
    d+=hy[n];  /*VF: this is dubious sometimes*/
  }
  if(d<1e-3) {
    fprintf(stderr,"[%s] ERROR: no histogram\n",fcname);
    free(hx);
    free(hy);
    return 0.0;
  }
  if(verbose>1) {
    fprintf(stdout,"-d (niik_image_get_mode) writing tmp_x.txt\n");
    niik_write_double_vector("tmp_x.txt",hx,num);
    niik_write_double_vector("tmp_y0.txt",hy,num);
  }
  if(verbose) fprintf(stdout,"-d (niik_image_get_mode) average: vec %i   avg %i\n",num,avgnum);
  if(!niik_runavg_double_vector(hy,num,avgnum)) {
    fprintf(stderr,"[%s] ERROR: niik_runavg_double_vector\n",fcname);
    free(hx);
    free(hy);
    return NIIKMAX;
  }
  if(verbose>1) {
    fprintf(stdout,"-d (niik_image_get_mode) writing tmp_y1.txt\n");
    niik_write_double_vector("tmp_y1.txt",hy,num);
  }
  if((n=niik_get_max_index_double_vector(hy,num))<0) {
    fprintf(stderr,"[%s] ERROR: niik_get_min_index_double_vector\n",fcname);
    free(hx);
    free(hy);
    return NIIKMAX;
  }
  d=hx[n];
  free(hx);
  free(hy);
  return d;
} /* niik_image_get_mode */

double niik_image_get_median(nifti_image *img,nifti_image *maskimg) {
  niikvec *v;
  double
  out;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(maskimg==NULL) {
    v=(niikvec *)calloc(1,sizeof(niikvec));
    if((v->v=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
      return 0;
    }
    v->num=img->nvox;
  } else {
    if((v=niik_image_get_voxels_as_double_vector_within_mask(img,maskimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector_within_mask\n");
      return 0;
    }
  }
  if(!niikvec_sort(v)) {
    fprintf(stderr,"ERROR: niikvec_sort\n");
    return 0;
  }
  out = niik_get_median_from_sorted_double_vector(v->v,v->num);
  v = niikvec_free(v);
  return out;
} /* niik_image_get_median */

int niik_image_get_min_index(nifti_image *img,nifti_image *maskimg) {
  double out,dval;
  int m,n,k=-1;
  char fcname[32]="niik_image_get_min_index";
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return -1;
  }
  if(maskimg==NULL) {
    for(n=0,out=1e99; n<img->nvox; n++) {
      dval=niik_image_get_voxel(img,n);
      if(out>dval) {
        out=dval;
        k=n;
      }
    }
  } else {
    if(img->nvox!=maskimg->nvox) {
      fprintf(stderr,"[%s] ERROR: nvox is different for img %i and maskimg %i\n",fcname,img->nvox,maskimg->nvox);
      return -1;
    }
    for(n=m=0,out=1e99; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        dval=niik_image_get_voxel(img,n);
        if(out>dval) {
          /* fprintf(stdout,"%i %f\n",n,dval); */
          out=dval;
          k=n;
        }
      }
    }
  }
  return k;
}

int niik_image_get_max_index(nifti_image *img,nifti_image *maskimg) {
  double out,dval;
  int m,n,k=-1;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return -1;
  }
  if(maskimg==NULL) {
    for(n=0,out=-1e99; n<img->nvox; n++) {
      dval=niik_image_get_voxel(img,n);
      if(out<dval) {
        out=dval;
        k=n;
      }
    }
  } else {
    for(n=m=0,out=-1e99; n<img->nvox; n++) {
      if(niik_image_get_voxel(maskimg,n)>0)  {
        dval=niik_image_get_voxel(img,n);
        if(out<dval) {
          out=dval;
          k=n;
        }
      }
    }
  }
  return k;
}



/**** end of stat functions ****/







/***************************************************
 *
 * TIFF NIFTI_IMAGE OPERATIONS
 *
 * int niik_tiff_write_xslice(const char *fname,nifti_image *img,double dmin,double dmax,int xslice);
 * int niik_tiff_write_yslice(const char *fname,nifti_image *img,double dmin,double dmax,int yslice);
 * int niik_tiff_write_zslice(const char *fname,nifti_image *img,double dmin,double dmax,int zslice);
 *
 ***************************************************/

int niik_tiff_write_xmontage(const char *fname,nifti_image *img,double dmin,double dmax,int *xslice,int num_xslice) {

  return 1;
}

#if 0
int niik_tiff_write_xslice(const char *fname,nifti_image *img,double dmin,double dmax,int xslice) {
  TIFF *tifimg;
  int
  i,j,k,m,
  isize,
  area,size;
  double
  thismin,thismax,
          d,dran;
  unsigned char *bout;
  int verbose=0;
  char fcname[64]="niik_tiff_write_xslice";
  NIIK_RET0((img==NULL),fcname,"img is null");
  if(verbose>0) {
    fprintf(stdout,"  sagittal slice for tiff image %i %s\n",xslice,fname);
  }
  size=area=img->ny*img->nz;
  isize=img->nx * img->ny * img->nz * img->nt * img->nu;
  dran=dmax - dmin;
  if(img->nv==3) {
    bout = (unsigned char *)calloc(3*size,sizeof(char));
    if(verbose>0) {
      fprintf(stdout,"[%s] color image\n",fcname);
    }
  } else
    bout = (unsigned char *)calloc(size,sizeof(char));
  area=img->nx*img->ny;
  fprintf(stdout,"  min/max  %8.3f %8.3f\n",dmin,dmax);
  for(k=0; k<img->nx; k++) {
    if(k!=xslice) continue;
    if((tifimg = TIFFOpen(fname, "w")) == NULL) {
      fprintf(stderr,"ERROR: could not open for writing, %i\n",k);
      return 0;
    }
    TIFFSetField( tifimg, TIFFTAG_IMAGEWIDTH,  img->ny );
    TIFFSetField( tifimg, TIFFTAG_IMAGELENGTH, img->nz );
    TIFFSetField( tifimg, TIFFTAG_BITSPERSAMPLE,   8);
    if(img->nv==3) {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 3);
    } else {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 1);
    }
    TIFFSetField( tifimg, TIFFTAG_ROWSPERSTRIP, img->nz);
    TIFFSetField( tifimg, TIFFTAG_COMPRESSION,  1);
    if(img->nv==3) {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_RGB);
    } else {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_MINISBLACK);
    }
    TIFFSetField( tifimg, TIFFTAG_FILLORDER,    FILLORDER_LSB2MSB);
    TIFFSetField( tifimg, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField( tifimg, TIFFTAG_XRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_YRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
    TIFFSetField( tifimg, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);

    if(img->nv==3) {
      for(j=m=0; j<img->nz; j++) {
        for(i=0; i<img->ny; i++) {
          d = niik_image_get_voxel(img,k+i*img->nx+j*area);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,k+i*img->nx+j*area+isize);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,k+i*img->nx+j*area+isize*2);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
        }
      } /* each voxel on a slice */
    } else {
      for(j=m=0,thismin=dmax,thismax=dmin; j<img->nz; j++) {
        for(i=0; i<img->ny; m++,i++) {
          d = niik_image_get_voxel(img,k+i*img->nx+j*area);
          bout[m] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          if(d>thismax) thismax=d;
          if(d<thismin) thismin=d;
        }
      } /* each voxel on a slice */
      /* fprintf(stdout,"\t  %i  %7.2f %7.2f\n",k,thismin,thismax); */
    }
    TIFFSetField(tifimg, TIFFTAG_BITSPERSAMPLE,  8);
    if(img->nv==3) {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->ny * img->nz * 3);
    } else {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->ny * img->nz );
    }
    TIFFClose(tifimg);
    fprintf(stdout,"    writing %s\n",fname);
  } /* shift in x */
  if(bout!=NULL) {
    free(bout);
    bout=NULL;
  }
  return 1;
} /* doing sagittal pictures */




int niik_tiff_write_yslice(const char *fname,nifti_image *img,double dmin,double dmax,int yslice) {
  TIFF *tifimg;
  int
  i,j,k,m,
  isize,
  area,size;
  double
  thismin,thismax,
          d,dran;
  unsigned char *bout;
  if(img==NULL) {
    fprintf(stderr,"ERROR: niik_tiff_write_yslice\n");
    return 0;
  }
  fprintf(stdout,"  coronal slice for tiff image %i %s\n",yslice,fname);
  size=area=img->nx*img->nz;
  isize=img->nx * img->ny * img->nz * img->nt * img->nu;
  dran=dmax - dmin;
  if(img->nv==3)
    bout = (unsigned char *)calloc(3*size,sizeof(char));
  else
    bout = (unsigned char *)calloc(size,sizeof(char));
  area=img->nx*img->ny;
  fprintf(stdout,"  min/max  %8.3f %8.3f\n",dmin,dmax);
  for(k=0; k<img->ny; k++) {
    if(k!=yslice) continue;
    if((tifimg = TIFFOpen(fname, "w")) == NULL) {
      fprintf(stderr,"ERROR: could not open for writing, %i\n",k);
      return 0;
    }
    TIFFSetField( tifimg, TIFFTAG_IMAGEWIDTH,  img->nx );
    TIFFSetField( tifimg, TIFFTAG_IMAGELENGTH, img->nz );
    TIFFSetField( tifimg, TIFFTAG_BITSPERSAMPLE,   8);
    if(img->nv==3) {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 3);
    } else {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 1);
    }
    TIFFSetField( tifimg, TIFFTAG_ROWSPERSTRIP, img->nz);
    TIFFSetField( tifimg, TIFFTAG_COMPRESSION,  1);
    if(img->nv==3) {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_RGB);
    } else {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_MINISBLACK);
    }
    TIFFSetField( tifimg, TIFFTAG_FILLORDER,    FILLORDER_LSB2MSB);
    TIFFSetField( tifimg, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField( tifimg, TIFFTAG_XRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_YRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
    TIFFSetField( tifimg, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    if(img->nv==3) {
      for(j=m=0; j<img->nz; j++) {
        for(i=0; i<img->nx; i++) {
          d = niik_image_get_voxel(img,i+k*img->nx+j*area);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,i+k*img->nx+j*area+isize);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,i+k*img->nx+j*area+isize*2);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
        }
      } /* each voxel on a slice */
    } else {
      for(j=m=0,thismin=dmax,thismax=dmin; j<img->nz; j++) {
        for(i=0; i<img->nx; m++,i++) {
          d = niik_image_get_voxel(img,i+k*img->nx+j*area);
          bout[m] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          if(d>thismax) thismax=d;
          if(d<thismin) thismin=d;
        }
      } /* each voxel on a slice */
      /* fprintf(stdout,"\t  %i  %7.2f %7.2f\n",k,thismin,thismax); */
    }
    TIFFSetField(tifimg, TIFFTAG_BITSPERSAMPLE,  8);
    if(img->nv==3) {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->nx * img->nz * 3);
    } else {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->nx * img->nz );
    }
    TIFFClose(tifimg);
    fprintf(stdout,"    writing %s\n",fname);
  } /* shift in x */
  if(bout!=NULL) {
    free(bout);
    bout=NULL;
  }
  return 1;
} /* doing coronal pictures */

int niik_tiff_write_zslice(const char *fname,nifti_image *img,double dmin,double dmax,int zslice) {
  TIFF *tifimg;
  int
  i,j,k,m,
  isize,
  area,size;
  double
  thismin,thismax,
          d,dran;
  unsigned char *bout;
  if(img==NULL) {
    fprintf(stderr,"ERROR: niik_tiff_write_zslice\n");
    return 0;
  }
  fprintf(stdout,"  axial slice for tiff image %i %s\n",zslice,fname);
  size=area=img->nx*img->ny;
  isize=img->nx * img->ny * img->nz * img->nt * img->nu;
  dran=dmax-dmin;
  if(img->nv==3)
    bout = (unsigned char *)calloc(3*size,sizeof(char));
  else
    bout = (unsigned char *)calloc(size,sizeof(char));
  area=img->nx*img->ny;
  fprintf(stdout,"  min/max  %8.3f %8.3f\n",dmin,dmax);
  for(k=0; k<img->nz; k++) {
    if(k!=zslice) continue;
    if((tifimg = TIFFOpen(fname, "w")) == NULL) {
      fprintf(stderr,"ERROR: could not open for writing, %i\n",k);
      return 0;
    }
    TIFFSetField( tifimg, TIFFTAG_IMAGEWIDTH,  img->nx );
    TIFFSetField( tifimg, TIFFTAG_IMAGELENGTH, img->ny );
    TIFFSetField( tifimg, TIFFTAG_BITSPERSAMPLE,   8);
    if(img->nv==3) {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 3);
    } else {
      TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 1);
    }
    TIFFSetField( tifimg, TIFFTAG_ROWSPERSTRIP, img->ny);
    TIFFSetField( tifimg, TIFFTAG_COMPRESSION,  1);
    if(img->nv==3) {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_RGB);
    } else {
      TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_MINISBLACK);
    }
    TIFFSetField( tifimg, TIFFTAG_FILLORDER,    FILLORDER_LSB2MSB);
    TIFFSetField( tifimg, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField( tifimg, TIFFTAG_XRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_YRESOLUTION, 150.0);
    TIFFSetField( tifimg, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
    TIFFSetField( tifimg, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    if(img->nv==3) {
      for(j=m=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          d = niik_image_get_voxel(img,i+j*img->nx+k*area);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,i+j*img->nx+k*area+isize);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          d = niik_image_get_voxel(img,i+j*img->nx+k*area+isize*2);
          bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
        }
      } /* each voxel on a slice */
    } else {
      for(j=m=0,thismin=dmax,thismax=dmin; j<img->ny; j++) {
        for(i=0; i<img->nx; m++,i++) {
          d = niik_image_get_voxel(img,i+j*img->nx+k*area);
          bout[m] = (unsigned char)NIIK_DMINMAX(255.0*(d-dmin)/dran,0,255);
          if(d>thismax) thismax=d;
          if(d<thismin) thismin=d;
        }
      } /* each voxel on a slice */
      /* fprintf(stdout,"\t  %i  %7.2f %7.2f\n",k,thismin,thismax); */
    }
    TIFFSetField(tifimg, TIFFTAG_BITSPERSAMPLE,  8);
    if(img->nv==3) {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->nx * img->ny * 3);
    } else {
      TIFFWriteEncodedStrip(tifimg, 0, bout, img->nx * img->ny );
    }
    TIFFClose(tifimg);
    fprintf(stdout,"    writing %s\n",fname);
  } /* shift in x */
  if(bout!=NULL) {
    free(bout);
    bout=NULL;
  }
  return 1;
} /* doing axial pictures */

#endif
/**** end of tiff processing ****/




/***************************************************
 *
 * IMAGE OPERATION IN MORE THAN 3D
 *
 ***************************************************/

nifti_image *niik_image_combine(nifti_image **imglist,int nimg,int dimnum,double scale_max)
/* image combination
 * -imglist is a list of pointers to differnet (or can be the same) image
 * -individual image in imglist should have the same number of voxels with same dimensions
 * -dimnum can be 4 for time series or 5 for some other series
 * -scale_max is the output target maximum intensity (not actual)
 *  -image's 90 percentile will be scaled to this value
 *  -if scale_max is zero, then no intensity scaling (normalization)
 * -not sure what would happen if individual images had 4th or 5th dimensions on input
 * -updated January 21, 2014, Kunio Nakamura
 * --does not convert from double
 */
{
  nifti_image *outimg=NULL;
  double
  imax,fac;
  int
  nvox,
  n,m,i;
  char fcname[64]="niik_image_combine";
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is a null pointer\n");
    return NULL;
  }
  /* check for x,y,z dimensions; we may add others */
  for(n=1; n<nimg; n++) {
    if(imglist[0]->nx!=imglist[n]->nx) {
      fprintf(stderr,"ERROR: different nx %i %i for img %i\n",imglist[0]->nx,imglist[n]->nx,n);
      return NULL;
    }
    if(imglist[0]->ny!=imglist[n]->ny) {
      fprintf(stderr,"ERROR: different ny %i %i for img %i\n",imglist[0]->ny,imglist[n]->ny,n);
      return NULL;
    }
    if(imglist[0]->nz!=imglist[n]->nz) {
      fprintf(stderr,"ERROR: different nz %i %i for img %i\n",imglist[0]->nz,imglist[n]->nz,n);
      return NULL;
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  switch(dimnum) {
  case 6:  /* v series */
    for(n=1; n<nimg; n++) {
      if(imglist[0]->nu!=imglist[n]->nu) {
        fprintf(stderr,"ERROR: different nu %i %i for img %i\n",imglist[0]->nu,imglist[n]->nu,n);
        return NULL;
      }
    }
  case 5:  /* u series */
    for(n=1; n<nimg; n++) {
      if(imglist[0]->nt!=imglist[n]->nt) {
        fprintf(stderr,"ERROR: different nt %i %i for img %i\n",imglist[0]->nt,imglist[n]->nt,n);
        return NULL;
      }
    }
  case 4:  /* time series */
    break;
  default:
    fprintf(stderr,"[%s] ERROR: unknown dimension\n",fcname);
    return NULL;
  }
  for(n=nvox=0; n<nimg; n++) {
    nvox += imglist[n]->nvox;
  }
  free(outimg->data);
  outimg->data=NULL;
  NIIK_RET0(((outimg->data=(void *)calloc(nvox,outimg->nbyper))==NULL),fcname,"calloc for outimg->data");
  for(n=m=0; n<nimg; n++) {
    if(fabs(scale_max)<1e-5)
      fac = 1;
    else {
      imax = niik_image_get_percentile(imglist[n],NULL,0.9);
      if(niik_check_double_problem(imax)) {
        fprintf(stderr,"ERROR: niik_check_double_problem\n");
        outimg = niik_image_free(outimg);
        return NULL;
      }
      fac = scale_max / imax;
    }
    /* fprintf(stdout,"      %i factor = %6.3f | %6.3f\n",n,fac,scale_max); */
    for(i=0; i<imglist[n]->nvox; i++) {
      niik_image_set_voxel(outimg,m++,niik_image_get_voxel(imglist[n],i) * fac);
    }
  }
  outimg->ndim=0;
  outimg->nx=outimg->dim[1]=imglist[0]->nx;
  outimg->ny=outimg->dim[2]=imglist[0]->ny;
  outimg->nz=outimg->dim[3]=imglist[0]->nz;
  outimg->nt=outimg->dim[4]=1;
  outimg->nu=outimg->dim[5]=1;
  outimg->nv=outimg->dim[6]=1;
  outimg->nw=outimg->dim[7]=1;
  if(dimnum==4) { /* time series */
    for(n=m=0; n<nimg; n++) m+=(imglist[n]->nt==0)?1:imglist[n]->nt;
    outimg -> ndim = outimg->dim[0] = NIIK_IMAX(4,outimg->ndim);
    outimg -> nt   = outimg->dim[4] = m;
    outimg -> dt   = outimg->pixdim[4] = 1;
  } else if(dimnum==5) { /* vector image */
    for(n=m=0; n<nimg; n++) m+=imglist[n]->nu;
    outimg -> ndim = outimg->dim[0] = 5;
    outimg -> nt   = outimg->dim[4] = imglist[0]->nt;
    outimg -> nu   = outimg->dim[5] = m;
    outimg -> dt   = outimg->pixdim[4] = imglist[0]->dt;
    outimg -> du   = outimg->pixdim[5] = imglist[0]->du;
  } else { /* assume v for now */
    for(n=m=0; n<nimg; n++) m+=imglist[n]->nv;
    outimg -> ndim = outimg->dim[0] = 6;
    outimg -> nt   = outimg->dim[4] = imglist[0]->nt;
    outimg -> nu   = outimg->dim[5] = imglist[0]->nu;
    outimg -> nv   = outimg->dim[6] = m;
    outimg -> dt   = outimg->pixdim[4] = imglist[0]->dt;
    outimg -> du   = outimg->pixdim[5] = imglist[0]->du;
    outimg -> dv   = outimg->pixdim[6] = 1; /*imglist[0]->dv;*/
  }
  outimg -> nvox = 1;
  for(n=1; n<=outimg->ndim; n++) {
    outimg->nvox *= outimg -> dim[n];
  }
  outimg -> xyz_units  = NIFTI_UNITS_MM;
  outimg -> time_units = NIFTI_UNITS_SEC;
  /*outimg -> data = dimg;
    outimg -> datatype = NIFTI_TYPE_FLOAT64;
    outimg -> nbyper = sizeof(double);
    if(!niik_image_type_convert(outimg,imglist[0]->datatype)){
    fprintf(stderr,"[%s] ERROR: nifti_k_convert_image\n",fcname);
    return NULL; } */
  return outimg;
} /* niik_image_combine */

nifti_image **niik_image_unmerge3d(nifti_image *img) {
  nifti_image **imglist=NULL;
  int num,n,i,j;
  char fcname[64]="niik_image_unmerge3d";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  num = img->nt * img->nu * img->nv * img->nw;
  imglist=(nifti_image **)calloc(num,sizeof(nifti_image *));
  for(n=i=0; n<num; n++) {
    if(n==0) {
      imglist[n]=niik_image_copy(img);
      free(imglist[n]->data);
      imglist[n]->dim[0]=imglist[n]->ndim=3;
      imglist[n]->nvox=img->nx*img->ny*img->nz;
      imglist[n]->data=(void *)malloc(imglist[n]->nvox*imglist[n]->nbyper);
    } else {
      if((imglist[n]=niik_image_copy(imglist[0]))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
        return NULL;
      }
    }
    for(j=0; j<imglist[n]->nvox; i++,j++) {
      niik_image_set_voxel(imglist[n],j,niik_image_get_voxel(img,i));
    }
  }
  return imglist;
} /* niik_image_unmerge3d */

nifti_image *niik_image_unmerge3d_extract_one(nifti_image *img,int n)
/* n'th 'frame'
 * the order of higher dimensions (>3) is t, u, v, and w
 */
{
  nifti_image *outimg=NULL;
  int i,j;
  NIIK_RET0((img==NULL),__func__,"img is null");
  NIIK_RET0(((outimg=nifti_copy_nim_info(img))==NULL),__func__,"nifti_copy_nim_info");
  outimg->ndim=outimg->dim[0]=3;
  outimg->nt=outimg->dim[4]=1;
  outimg->nu=outimg->dim[5]=1;
  outimg->nv=outimg->dim[6]=1;
  outimg->nw=outimg->dim[7]=1;
  outimg->nvox=outimg->nx*outimg->ny*outimg->nz;
  NIIK_RET0(((outimg->data=(void *)calloc(outimg->nvox,outimg->nbyper))==NULL),__func__,"calloc for data");
  for(i=0,j=n*outimg->nvox; i<outimg->nvox; i++,j++) {
    niik_image_set_voxel(outimg,i,niik_image_get_voxel(img,j));
  }
  return outimg;
} /* niik_image_unmerge3d */


int niik_image_convert_to_color_image(nifti_image *img)
/* convert to color image (but remain grayscale) */
{
  int
  i,
  size,size3;
  double
  dval,
  *dimg=NULL;
  char fcname[64]="niik_image_convert_to_color_image";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  /* nv is already 3 */
  if(img->nv==3) {
    return 1;
  }
  size3=size=img->nx*img->ny*img->nz;
  size3*=3;
  if((dimg=(double *)calloc(size3,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return 0;
  }
  for(i=0; i<size; i++) {
    dval=niik_image_get_voxel(img,i);
    dimg[i]=dimg[i+size]=dimg[i+size*2]=dval;
  }
  img->ndim=img->dim[0]=6;
  img->nt=img->dim[4]=1;
  img->nu=img->dim[5]=1;
  img->nv=img->dim[6]=3;
  img->nvox=size3;
  free(img->data);
  if((img->data=(void *)calloc(img->nvox,img->nbyper))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return 0;
  }
  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(dimg);
  return 1;
} /* make 3 bands (nv) */


nifti_image *niik_image_montage(nifti_image **imglist,int nimg,int tile,int xlo,int ylo,int zlo,int xhi,int yhi,int zhi)
/* montage image
 * -returns combined image
 * -imglist is the list of images
 *  they have to have the same dimensions
 * -nimg is the number of images
 * -tile is the #columns, #rows is determined by nimg/tile
 * -xlo|ylo|zlo|xhi|yhi|zhi are the crop index
 * -if negative, then automatically calculate
 */
{
  nifti_image
  *outimg;
  int
  ncol,nrow,
       m,n,i,j,k,t,u,v,
       xdim,ydim,zdim,
       nn,ii,jj,kk;
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return NULL;
  }
  if(xlo<0) xlo=0;
  if(ylo<0) ylo=0;
  if(zlo<0) zlo=0;
  if(xhi<0) xhi=imglist[0]->nx-1;
  if(yhi<0) yhi=imglist[0]->ny-1;
  if(zhi<0) zhi=imglist[0]->nz-1;
  if(xhi>=imglist[0]->nx-1) xhi=imglist[0]->nx-1;
  if(yhi>=imglist[0]->ny-1) yhi=imglist[0]->ny-1;
  if(zhi>=imglist[0]->nz-1) zhi=imglist[0]->nz-1;
  for(n=1; n<nimg; n++) {
    if(niik_image_cmp_dim(imglist[0],imglist[n])!=0) {
      fprintf(stderr,"ERROR: niik_image_cmp_dim %s %s\n",imglist[0]->fname,imglist[n]->fname);
      return NULL;
    }
    /*if(imglist[n]->ndim>3) {
      fprintf(stderr,"ERROR: ndim = %i (img %i)\n",imglist[n]->ndim,n);
      return NULL; }*/
    /*fprintf(stdout,"ndim  %i  [%3i %3i %3i %3i %3i %3i]\n",imglist[n]->ndim,
      imglist[n]->nx,imglist[n]->ny,imglist[n]->nz,
      imglist[n]->nt,imglist[n]->nu,imglist[n]->nv);*/
  }
  fprintf(stdout,"[niik_image_montage] tile %3i\n",tile);
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  free(outimg->data);
  ncol=tile;
  nrow=nimg/tile;
  /*fprintf(stdout,"[%3i %3i %3i %3i %3i %3i] \n",
    outimg->nx,outimg->ny,outimg->nz,outimg->nt,outimg->nu,outimg->nv);
    fprintf(stdout,"new image %i %i | %i\n",ncol,nrow,nimg);*/
  outimg->nx=xdim=outimg->dim[1]=xhi-xlo+1;
  outimg->ny=ydim=outimg->dim[2]=yhi-ylo+1;
  outimg->nz=zdim=outimg->dim[3]=zhi-zlo+1;
  outimg->nx*=ncol;
  outimg->ny*=nrow;
  outimg->nvox=outimg->nx*outimg->ny*outimg->nz*outimg->nt*outimg->nu*outimg->nv;
  outimg->data=(void *)calloc(outimg->nvox,sizeof(outimg->nbyper));
  /*  fprintf(stdout,"[%3i %3i %3i %3i %3i %3i] \n",
      outimg->nx,outimg->ny,outimg->nz,outimg->nt,outimg->nu,outimg->nv);*/
  for(m=nn=0; m<nrow; m++) {
    for(n=0; n<ncol; nn++,n++) {
      if(nn>nimg) continue;
      fprintf(stdout,"  %i [%i %i]  %s\n",nn,n,m,imglist[nn]->fname);
      for(v=0; v<imglist[nn]->nv; v++) {
        for(u=0; u<imglist[nn]->nu; u++) {
          for(t=0; t<imglist[nn]->nt; t++) {
            for(k=zlo; k<=zhi; k++) {
              kk=k-zlo;
              for(j=ylo; j<=yhi; j++) {
                jj=m*ydim+j-ylo;
                for(i=xlo; i<=xhi; i++) {
                  ii=n*xdim+i-xlo;
                  /*if(nn) fprintf(stdout,"[%3i %3i %3i %3i] -> [%3i %3i %3i]  [tuv]=[%3i %3i %3i]\n",i,j,k,nn,ii,jj,k,t,u,v);*/
                  niik_image_set_voxel(outimg,
                                       ii+
                                       jj*outimg->nx+
                                       kk*outimg->nx*outimg->ny+
                                       t *outimg->nx*outimg->ny*outimg->nz+
                                       u *outimg->nx*outimg->ny*outimg->nz*outimg->nt+
                                       v *outimg->nx*outimg->ny*outimg->nz*outimg->nt*outimg->nu,

                                       niik_image_get_voxel(imglist[nn],
                                           i+
                                           j*imglist[nn]->nx+
                                           k*imglist[nn]->nx*imglist[nn]->ny+
                                           t*imglist[nn]->nx*imglist[nn]->ny*imglist[nn]->nz+
                                           u*imglist[nn]->nx*imglist[nn]->ny*imglist[nn]->nz*imglist[nn]->nt+
                                           v*imglist[nn]->nx*imglist[nn]->ny*imglist[nn]->nz*imglist[nn]->nt*imglist[nn]->nu));
                }
              }
            }
          }
        }
      }
    }
  }
  return outimg;
} /* niik_image_montage */

int niik_image_combine_and_write_as_type(char *filename,nifti_image **imglist,int nimg,char dim,double scale_max,int datatype) {
  nifti_image *tmpimg=NULL;
  char fcname[65]="niik_image_combine_and_write_as_type";
  int dim_num=0,verbose=0;
  if(verbose>=1) niik_fc_display(fcname,1);
  switch(dim) {
  case 't':
    dim_num=4;
    break;
  case 'u':
    dim_num=5;
    break;
  case 'v':
    dim_num=6;
    break;
  case 'w':
    dim_num=7;
    break;
  default:
    fprintf(stderr,"[%s] ERROR: t-dimension is used (%i)\n",fcname,dim);
    dim_num=4;
    break;
  }
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return 0;
  }
  if((tmpimg=niik_image_combine(imglist,nimg,dim_num,scale_max))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine\n",fcname);
    return 0;
  }
  if(!niik_image_type_convert(tmpimg,datatype)) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine\n",fcname);
    return 0;
  }
  if(!niik_image_write(filename,tmpimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_write\n",fcname);
    return 0;
  }
  tmpimg=niik_image_free(tmpimg);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_combine_and_write_as_type */

int niik_image_combine_and_write(char *filename,nifti_image **imglist,int nimg,char dim,double scale_max) {
  int niik_image_combine_and_write_as_type(char *filename,nifti_image **imglist,int nimg,char dim,double scale_max,int datatype);
  char fcname[65]="niik_image_combine_and_write";
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return 0;
  }
  if(!niik_image_combine_and_write_as_type(filename,imglist,nimg,dim,scale_max,imglist[0]->datatype)) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write_as_type\n",fcname);
    return 0;
  }
  return 1;
}


/*** end of 4,5,6th dimensions ****/


nifti_image *niik_image_average_multiple(nifti_image **imglist,int nimg,double *faclist) {
  nifti_image *outimg=NULL;
  int n,i;
  double
  fac,
  *dimg;
  /* check inputs */
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return NULL;
  }
  if(nimg<=1) {
    fprintf(stderr,"ERROR: nimg is zero\n");
    return NULL;
  }
  /* check input's dimensions */
  for(n=1; n<nimg; n++) {
    for(i=0; i<imglist[0]->dim[0]; i++) {
      if(imglist[0]->dim[i]!=imglist[n]->dim[i]) {
        fprintf(stderr,"ERROR: img[%i]'s dimension is different [%i vs %i]\n",n,imglist[0]->dim[i],imglist[n]->dim[i]);
        return NULL;
      }
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: outimg=niik_image_copy(imglist[0])\n");
    return NULL;
  }
  if((dimg=(double *)calloc(outimg->nvox,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return NULL;
  }
  for(n=0; n<nimg; n++) {
    if(faclist==NULL) fac=1;
    else fac=faclist[n];
    for(i=0; i<outimg->nvox; i++) {
      dimg[i] += fac*niik_image_get_voxel(imglist[n],i);
    }
  }
  for(i=0; i<outimg->nvox; i++) {
    dimg[i]/=nimg;
  }
  if(!niik_image_set_voxels_from_double_vector(outimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector(outimg,dimg)\n");
    outimg=niik_image_free(outimg);
    return NULL;
  }
  free(dimg);
  return outimg;
}

int niik_image_multiply_2_images(nifti_image *changed,nifti_image *mult) {
  char fcname[64]="niik_image_multiply_2_images";
  int i;
  if(changed==NULL) {
    fprintf(stderr,"[%s] ERROR: change_image is null\n",fcname);
    return 0;
  }
  if(mult==NULL) {
    fprintf(stderr,"[%s] ERROR: multiply_image is null\n",fcname);
    return 0;
  }
  if(changed->nvox!=mult->nvox) {
    fprintf(stderr,"[%s] ERROR: #voxel is different %i %i\n",fcname,changed->nvox,mult->nvox);
    return 0;
  }
  for(i=0; i<changed->nvox; i++) {
    niik_image_mul_voxel(changed,i,niik_image_get_voxel(mult,i));
  }
  return 1;
}

int niik_image_add_2_images(nifti_image *changed,nifti_image *add) {
  const char *fcname=__func__;
  int i;
  if(changed==NULL) {
    fprintf(stderr,"[%s] ERROR: change_image is null\n",fcname);
    return 0;
  }
  if(add==NULL) {
    fprintf(stderr,"[%s] ERROR: add_image is null\n",fcname);
    return 0;
  }
  if(changed->nvox!=add->nvox) {
    fprintf(stderr,"[%s] ERROR: #voxel is different %i %i\n",fcname,changed->nvox,add->nvox);
    return 0;
  }
  for(i=0; i<changed->nvox; i++) {
    niik_image_add_voxel(changed,i,niik_image_get_voxel(add,i));
  }
  return 1;
}



nifti_image *niik_image_maximize_multiple(nifti_image **imglist,int nimg,double *faclist) {
  nifti_image *outimg=NULL;
  int n,i;
  double
  fac=1,
  *dimg;
  /* check inputs */
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return NULL;
  }
  if(nimg<=1) {
    fprintf(stderr,"ERROR: nimg is zero\n");
    return NULL;
  }
  /* check input's dimensions */
  for(n=1; n<nimg; n++) {
    for(i=0; i<imglist[0]->dim[0]; i++) {
      if(imglist[0]->dim[i]!=imglist[n]->dim[i]) {
        fprintf(stderr,"ERROR: img[%i]'s dimension is different [%i vs %i]\n",n,imglist[0]->dim[i],imglist[n]->dim[i]);
        return NULL;
      }
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: outimg=niik_image_copy(imglist[0])\n");
    return NULL;
  }
  if((dimg=(double *)calloc(outimg->nvox,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return NULL;
  }
  if(faclist==NULL) fac=1;
  else fac=faclist[0];
  for(i=0; i<outimg->nvox; i++) {
    dimg[i]=fac*niik_image_get_voxel(imglist[0],i);
  }
  for(n=1; n<nimg; n++) {
    if(faclist==NULL) fac=1;
    else fac=faclist[n];
    for(i=0; i<outimg->nvox; i++) {
      dimg[i] = NIIK_DMAX(dimg[i],fac*niik_image_get_voxel(imglist[n],i));
    }
  }
  if(!niik_image_set_voxels_from_double_vector(outimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector(outimg,dimg)\n");
    outimg=niik_image_free(outimg);
    return NULL;
  }
  free(dimg);
  return outimg;
}

nifti_image *niik_image_add_multiple(nifti_image **imglist,int nimg,double *faclist) {
  nifti_image *outimg=NULL;
  int n,i;
  double
  fac=1,
  *dimg;
  const char *fcname="niik_image_add_multiple";
  /* check inputs */
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return NULL;
  }
  if(nimg<=1) {
    fprintf(stderr,"ERROR: nimg is zero\n");
    return NULL;
  }
  /* check input's dimensions */
  for(n=1; n<nimg; n++) {
    for(i=0; i<imglist[0]->dim[0]; i++) {
      if(imglist[0]->dim[i]!=imglist[n]->dim[i]) {
        fprintf(stderr,"ERROR: img[%i]'s dimension is different [%i vs %i]\n",n,imglist[0]->dim[i],imglist[n]->dim[i]);
        return NULL;
      }
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: outimg=niik_image_copy(imglist[0])\n");
    return NULL;
  }
  if((dimg=(double *)calloc(outimg->nvox,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return NULL;
  }
  if(faclist==NULL) fac=1;
  else {
    fac=faclist[0];
    fprintf(stdout,"[%s] %9.4f   %s\n",fcname,fac,imglist[0]->fname);
  }
  for(i=0; i<outimg->nvox; i++) {
    dimg[i]=fac*niik_image_get_voxel(imglist[0],i);
  }
  for(n=1; n<nimg; n++) {
    if(faclist==NULL) fac=1;
    else {
      fac=faclist[n];
      fprintf(stdout,"[%s] %9.4f   %s\n",fcname,fac,imglist[n]->fname);
    }
    for(i=0; i<outimg->nvox; i++) {
      dimg[i] += fac*niik_image_get_voxel(imglist[n],i);
    }
  }
  if(!niik_image_set_voxels_from_double_vector(outimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector(outimg,dimg)\n");
    outimg=niik_image_free(outimg);
    return NULL;
  }
  free(dimg);
  return outimg;
}


nifti_image *niik_image_mul_multiple(nifti_image **imglist,int nimg) {
  nifti_image *outimg=NULL;
  int n,i;
  double
  *dimg;
  /* check inputs */
  if(imglist==NULL) {
    fprintf(stderr,"ERROR: imglist is null\n");
    return NULL;
  }
  if(nimg<=1) {
    fprintf(stderr,"ERROR: nimg is zero\n");
    return NULL;
  }
  /* check input's dimensions */
  for(n=1; n<nimg; n++) {
    for(i=0; i<imglist[0]->dim[0]; i++) {
      if(imglist[0]->dim[i]!=imglist[n]->dim[i]) {
        fprintf(stderr,"ERROR: img[%i]'s dimension is different [%i vs %i]\n",n,imglist[0]->dim[i],imglist[n]->dim[i]);
        return NULL;
      }
    }
  }
  if((outimg=niik_image_copy(imglist[0]))==NULL) {
    fprintf(stderr,"ERROR: outimg=niik_image_copy(imglist[0])\n");
    return NULL;
  }
  if((dimg=(double *)calloc(outimg->nvox,sizeof(double)))==NULL) {
    fprintf(stderr,"ERROR: calloc\n");
    return NULL;
  }
  for(i=0; i<outimg->nvox; i++) {
    dimg[i]=niik_image_get_voxel(imglist[0],i);
  }
  for(n=1; n<nimg; n++) {
    for(i=0; i<outimg->nvox; i++) {
      dimg[i] *= niik_image_get_voxel(imglist[n],i);
    }
  }
  if(!niik_image_set_voxels_from_double_vector(outimg,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector(outimg,dimg)\n");
    outimg=niik_image_free(outimg);
    return NULL;
  }
  free(dimg);
  return outimg;
} /* niik_image_mul_multiple */





/***************************************************
 *
 * ROTATE FUNCTIONS
 *
 ***************************************************/


int niik_image_rotate(nifti_image *img,char *rstr) {
  double *dimg;
  const char *fcname="niik_image_rotate";
  niikpt ct;
  int
  verbose=1,
  xydim,xdim,
  ni,nj,nk,
  n,i,j,k,
  ii,jj,kk;
  niikmat *afmat[2];
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(rstr==NULL) {
    fprintf(stderr,"[%s] ERROR: rstr is null\n",fcname);
    return 0;
  }
  if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector(img)\n");
    return 0;
  }
  xdim  = img->nx;
  xydim = img->nx*img->ny;
  afmat[0]=afmat[1]=NULL;

  if     (!strncmp(rstr,"flipxyz",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(i=0; i<img->nvox; i++) {
      dimg[img->nvox-i-1] = niik_image_get_voxel(img,i);
    }
  }

  else if(!strncmp(rstr,"swap-xy",7)) {
    if(verbose) fprintf(stdout,"[%s] %s\n",fcname,rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->nx; j++) {
        for(i=0; i<img->ny; i++,n++) {
          ii=j;
          jj=i;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + ii*xdim + jj);
        }
      }
    }
    if(verbose) fprintf(stdout,"[%s] header\n",fcname);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&   img->dim[1],   &img->dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    if(img->qform_code) {
      if(verbose) fprintf(stdout,"[%s] qform\n",fcname);
      afmat[0] =
        niikmat_multiply_free12(niikmat_affine_matrix_val(0,1,0,0,
                                1,0,0,0,
                                0,0,1,0,
                                0,0,0,1),
                                niikmat_mat44_matrix(img->qto_xyz));
      img->qto_xyz = niikmat_make_mat44(afmat[0]);
      img->qto_ijk = nifti_mat44_inverse(img->qto_xyz);
      nifti_mat44_to_quatern(img->qto_xyz,
                             &img->quatern_b,&img->quatern_c,&img->quatern_d,
                             &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                             NULL,NULL,NULL,
                             &img->qfac);
    }
    if(img->sform_code) {
      if(verbose) fprintf(stdout,"[%s] sform\n",fcname);
      afmat[1] =
        niikmat_multiply_free12(niikmat_affine_matrix_val(0,1,0,0,
                                1,0,0,0,
                                0,0,1,0,
                                0,0,0,1),
                                niikmat_mat44_matrix(img->sto_xyz));
    }
    img->sto_xyz = niikmat_make_mat44(afmat[1]);
    img->sto_ijk = nifti_mat44_inverse(img->sto_xyz);
    if(verbose>1) niikmat_display_msg("swap xy\n",afmat[0]);
  }

  else if(!strncmp(rstr,"axi2cor",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->ny; k++) {
      for(j=0; j<img->nz; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=i;
          jj=k;
          kk=j;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&   img->dim[2],   &img->dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
  }

  else if(!strncmp(rstr,"axi2sag",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nx; k++) {
      for(j=0; j<img->nz; j++) {
        for(i=0; i<img->ny; i++,n++) {
          ii=j;
          jj=k;
          kk=i;
          dimg[n]=niik_image_get_voxel(img,kk*img->nx + jj + ii*xydim);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&   img->dim[2],&   img->dim[2]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[2]);
  }

  else if(!strncmp(rstr,"cor2axi",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->ny; k++) {
      for(j=0; j<img->nz; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=i;
          jj=k;
          kk=j;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&   img->dim[2],&   img->dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
  }

  else if(!strncmp(rstr,"sag2axi",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->ny; k++) {
      for(j=0; j<img->nx; j++) {
        for(i=0; i<img->nz; i++,n++) {
          ii=j;
          jj=k;
          kk=i;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %4i %4i %4i\n",img->nx,img->ny,img->nz);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&   img->dim[2],&   img->dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %4i %4i %4i\n",img->nx,img->ny,img->nz);
  } /* sag2axi */

  else if(!strncmp(rstr,"sag2cor",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nx; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nz; i++,n++) {
          ii=k;
          jj=j;
          kk=i;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
  }

  else if(!strncmp(rstr,"cor2axi",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->ny; k++) {
      for(j=0; j<img->nz; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=i;
          jj=k;
          kk=j;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&   img->dim[2],&   img->dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
  } /* cor2axi */

  else if(!strncmp(rstr,"cor2sag",7)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nx; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nz; i++,n++) {
          ii=k;
          jj=j;
          kk=i;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    if(img->qform_code) {
      NIIK_RET0(((afmat[0]=niikmat_mat44_matrix(img->qto_xyz))==NULL),fcname,"niikmat_mat44_matrix for qform");
      NIIK_RET0(((afmat[0]=
                    niikmat_multiply_free12(niikmat_affine_matrix_val(0,0,1,0,
                                            0,1,0,0,
                                            1,0,0,0,
                                            0,0,0,1),
                                            afmat[0]))==NULL),fcname,"niikmat_multiply_free12");
      img->qto_xyz = niikmat_make_mat44(afmat[0]);
      img->qto_ijk = nifti_mat44_inverse(img->qto_xyz);
      nifti_mat44_to_quatern(img->qto_xyz,
                             &img->quatern_b,&img->quatern_c,&img->quatern_d,
                             &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                             NULL,NULL,NULL,
                             &img->qfac);
    }
    if(img->sform_code) {
      NIIK_RET0(((afmat[1]=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix for sform");
      NIIK_RET0(((afmat[1]=
                    niikmat_multiply_free12(niikmat_affine_matrix_val(0,0,1,0,
                                            0,1,0,0,
                                            1,0,0,0,
                                            0,0,0,1),
                                            afmat[1]))==NULL),fcname,"niikmat_multiply_free12");
      img->sto_xyz = niikmat_make_mat44(afmat[1]);
      img->sto_ijk = nifti_mat44_inverse(img->sto_xyz);
    }
  }

  else if(!strncmp(rstr,"flipx",5)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=img->nx-i-1;
          jj=j;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
  }

  else if(!strncmp(rstr,"flipy",5)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=i;
          jj=img->ny-j-1;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
  }

  else if(!strncmp(rstr,"flipz",5)) {
    if(verbose) fprintf(stdout,"[%s] %s\n",fcname,rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=i;
          jj=j;
          kk=img->nz-k-1;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    ct=niikpt_val((img->nx-1.0)*img->dx,
                  (img->ny-1.0)*img->dy,
                  (img->nz-1.0)*img->dz,
                  0);
    if(img->qform_code) {
      NIIK_RET0(((afmat[0]=niikmat_mat44_matrix(img->qto_xyz))==NULL),fcname,"niikmat_mat44_matrix for qform");
      NIIK_RET0(((afmat[0]=
                    niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,
                                            0,0,0,
                                            1,1,-1,
                                            0,0,0,
                                            ct.x,ct.y,ct.z),
                                            afmat[0]))==NULL),fcname,"niikmat_multiply_free12");
      img->qto_xyz = niikmat_make_mat44(afmat[0]);
      img->qto_ijk = nifti_mat44_inverse(img->qto_xyz);
      nifti_mat44_to_quatern(img->qto_xyz,
                             &img->quatern_b,&img->quatern_c,&img->quatern_d,
                             &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                             NULL,NULL,NULL,
                             &img->qfac);
    }
    if(img->sform_code) {
      NIIK_RET0(((afmat[1]=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix for sform");
      NIIK_RET0(((afmat[1]=
                    niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,
                                            0,0,0,
                                            1,1,-1,
                                            0,0,0,
                                            ct.x,ct.y,ct.z),
                                            afmat[1]))==NULL),fcname,"niikmat_multiply_free12");
      img->sto_xyz = niikmat_make_mat44(afmat[1]);
      img->sto_ijk = nifti_mat44_inverse(img->sto_xyz);
    }
  } /* flipz */

  else if(!strncmp(rstr,"right",5)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->nx; j++) {
        for(i=0; i<img->ny; i++,n++) {
          ii=j;
          jj=img->ny-i-1;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
  } else if(!strncmp(rstr,"left",4)) {
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) %s\n",rstr);
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->nx; j++) {
        for(i=0; i<img->ny; i++,n++) {
          ii=img->nx-j-1;
          jj=i;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(verbose) fprintf(stdout,"-v (niik_image_rotate) header\n");
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&   img->dim[1],&   img->dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
  }

  else if(!strncmp(rstr,"norm",4)) {
    if(verbose) fprintf(stdout,"[niik_image_rotate] %s\n",rstr);
    afmat[0]=afmat[1]=NULL;
    if(img->qform_code) {
      NIIK_RET0(((afmat[0]=niikmat_mat44_matrix(img->qto_xyz))==NULL),fcname,"niikmat_mat44_matrix for qform");
    }
    if(img->sform_code) {
      NIIK_RET0(((afmat[1]=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix for sform");
    }
    if(img->qform_code) {
      nifti_mat44_to_orientation( img->qto_xyz, &ni,&nj,&nk );
    }
    if(img->sform_code) {
      nifti_mat44_to_orientation( img->sto_xyz, &ni,&nj,&nk );
    } else {
      fprintf(stderr,"[niik_image_rotate] ERROR: missing qform or sform\n");
      return 0;
    }
    ct=niikpt_zero();
    /*ct.x = (img->nx-1.0) * img->dx * 0.5;
      ct.y = (img->ny-1.0) * img->dy * 0.5;
      ct.z = (img->nz-1.0) * img->dz * 0.5;*/
    if(verbose>1) niikmat_display_msg("initial matrix\n",afmat[0]);
    /* fprintf(stdout,"  %i %i %i\n",ni,nj,nk); */
    if((ni%2)==0) {
      if(verbose) fprintf(stdout,"[niik_image_rotate] flip in x\n");
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; i++,n++) {
            ii=img->nx-i-1;
            jj=j;
            kk=k;
            dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
          }
        }
      }
      if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
        fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
        free(dimg);
      }
      ni--;
      if(img->qform_code) {
        NIIK_RET0(((afmat[0]=
                      niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,-1,1,1,0,0,0,ct.x,ct.y,ct.z),
                                              afmat[0]))==NULL),fcname,"niikmat_multiply_free12");
      }
      if(img->sform_code) {
        NIIK_RET0(((afmat[1]=
                      niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,-1,1,1,0,0,0,ct.x,ct.y,ct.z),
                                              afmat[1]))==NULL),fcname,"niikmat_multiply_free12");
      }
      if(verbose>1) niikmat_display_msg("flip x\n",afmat[0]);
    }

    if((nj%2)==0) {
      if(verbose) fprintf(stdout,"[niik_image_rotate] flip in y\n");
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; i++,n++) {
            ii=i;
            jj=img->ny-j-1;
            kk=k;
            dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
          }
        }
      }
      if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
        fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
        free(dimg);
      }
      nj--;
      if(img->qform_code) {
        afmat[0] =
          niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,1,-1,1,0,0,0,ct.x,ct.y,ct.z),
                                  afmat[0]);
      }

      if(img->sform_code) {
        niikmat_display(afmat[1]);
        afmat[1] =
          niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,1,-1,1,0,0,0,ct.x,ct.y,ct.z),
                                  afmat[1]);
      }
      if(verbose>1) niikmat_display_msg("flip matrix\n",niikmat_affine_matrix_new(0,0,0,0,0,0,1,-1,1,0,0,0,ct.x,ct.y,ct.z));
      if(verbose>1) niikmat_display_msg("flip y\n",afmat[0]);
    }

    if((nk%2)==0) {
      if(verbose) fprintf(stdout,"[niik_image_rotate] flip in z\n");
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; i++,n++) {
            ii=i;
            jj=j;
            kk=img->nz-k-1;
            dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
          }
        }
      }
      if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
        fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
        free(dimg);
      }
      if(img->qform_code) {
        afmat[0] =
          niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,1,1,-1,0,0,0,ct.x,ct.y,ct.z),
                                  afmat[0]);
      }
      if(img->sform_code) {
        afmat[1] =
          niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,1,1,-1,0,0,0,ct.x,ct.y,ct.z),
                                  afmat[1]);
      }
      if(verbose>1) niikmat_display_msg("flip z\n",afmat[0]);
      nk--;
    }

    if(ni>2) {
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; i++,n++) {
            ii=i;
            jj=j;
            kk=k;
            if(nj<2)      {
              ii=j;
              jj=i;
            } else if(nk<2) {
              ii=k;
              kk=i;
            }
            dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
          }
        }
      }
      if(nj<2) {
        NIIK_ISWAP(&ni,&nj);
        NIIK_ISWAP(&img->nx,&img->ny);
        NIIK_ISWAP(&   img->dim[1],&   img->dim[2]);
        NIIK_FSWAP(&img->dx,&img->dy);
        NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
        if(img->qform_code) {
          afmat[0] =
            niikmat_multiply_free12(niikmat_affine_matrix_val(0,1,0,0,
                                    1,0,0,0,
                                    0,0,1,0,
                                    0,0,0,1),afmat[0]);
        }
        if(img->sform_code) {
          afmat[1] =
            niikmat_multiply_free12(niikmat_affine_matrix_val(0,1,0,0,
                                    1,0,0,0,
                                    0,0,1,0,
                                    0,0,0,1),afmat[1]);
        }
        if(verbose>1) niikmat_display_msg("swap xy\n",afmat[0]);
      }
      if(nk<2) {
        NIIK_ISWAP(&ni,&nk);
        NIIK_ISWAP(&img->nx,&img->nz);
        NIIK_ISWAP(&   img->dim[1],&   img->dim[3]);
        NIIK_FSWAP(&img->dx,&img->dz);
        NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
        if(img->qform_code) {
          NIIK_RET0(((afmat[0] =
                        niikmat_multiply_free12(niikmat_affine_matrix_val(0,0,1,0,
                                                0,1,0,0,
                                                1,0,0,0,
                                                0,0,0,1),
                                                afmat[0]))==NULL),
                    fcname,
                    "niikmat_multiply_free12, qform");
        }
        if(img->sform_code) {
          NIIK_RET0(((afmat[1] =
                        niikmat_multiply_free12(niikmat_affine_matrix_val(0,0,1,0,
                                                0,1,0,0,
                                                1,0,0,0,
                                                0,0,0,1),
                                                afmat[1]))==NULL),
                    fcname,
                    "niikmat_multiply_free12");
        }
        if(verbose>1) niikmat_display_msg("swap xz\n",afmat[0]);
      }
      if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
        fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
        free(dimg);
        return 0;
      }
    }
    if(nj>4) {
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; i++,n++) {
            ii=i;
            jj=k;
            kk=j;
            dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
          }
        }
      }
      NIIK_ISWAP(&nj,&nk);
      NIIK_ISWAP(&img->ny,&img->nz);
      NIIK_ISWAP(&   img->dim[2],&   img->dim[3]);
      NIIK_FSWAP(&img->dy,&img->dz);
      NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
      if(img->qform_code) {
        afmat[0] =
          niikmat_multiply_free12(niikmat_affine_matrix_val(1,0,0,0,
                                  0,0,1,0,
                                  0,1,0,0,
                                  0,0,0,1),afmat[0]);
      }
      if(img->sform_code) {
        afmat[1] =
          niikmat_multiply_free12(niikmat_affine_matrix_val(1,0,0,0,
                                  0,0,1,0,
                                  0,1,0,0,
                                  0,0,0,1),afmat[1]);
      }
      if(verbose>1) niikmat_display_msg("swap yz\n",afmat[0]);
      if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
        fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
        free(dimg);
        return 0;
      }
    }
    if(img->qform_code) {
      if(verbose>1) niikmat_display_msg("final qform\n",afmat[0]);
      /*nifti_disp_matrix_orient("old orientations:\n",img->qto_xyz);*/
      img->qto_xyz = niikmat_make_mat44(afmat[0]);
      /*niikmat_display(afmat[0]);
      nifti_disp_matrix_orient("new orientations:\n",img->qto_xyz);*/
      niikmat_inverse_update(afmat[0]);
      img->qto_ijk = niikmat_make_mat44(afmat[0]);
      nifti_mat44_to_quatern(img->qto_xyz,
                             &img->quatern_b,&img->quatern_c,&img->quatern_d,
                             &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                             NULL,NULL,NULL,
                             &img->qfac);
    }
    if(img->sform_code) {
      nifti_disp_matrix_orient("old orientations:\n",img->sto_xyz);
      img->sto_xyz = niikmat_make_mat44(afmat[1]);
      nifti_disp_matrix_orient("new orientations:\n",img->sto_xyz);
      niikmat_inverse_update(afmat[1]);
      img->sto_ijk = niikmat_make_mat44(afmat[1]);
    }
  } /* norm */

  else if(!strncmp(rstr,"LAS",3)) {
    if(!niik_image_rotate(img,"norm")) {
      fprintf(stderr,"ERROR: niik_image_rotate\n");
      return 0;
    }

    if(verbose) fprintf(stdout,"[niik_image_rotate] %s\n",rstr);
    afmat[0]=afmat[1]=NULL;
    if(img->qform_code) afmat[0]=niikmat_mat44_matrix(img->qto_xyz);
    if(img->sform_code) afmat[1]=niikmat_mat44_matrix(img->sto_xyz);
    if(img->qform_code) {
      nifti_mat44_to_orientation( img->qto_xyz, &ni,&nj,&nk );
    }
    if(img->sform_code) {
      nifti_mat44_to_orientation( img->sto_xyz, &ni,&nj,&nk );
    } else {
      fprintf(stderr,"[niik_image_rotate] ERROR: missing qform or sform\n");
      return 0;
    }
    ct=niikpt_zero();
    if(verbose>1) niikmat_display_msg("initial matrix\n",afmat[0]);

    if(verbose) fprintf(stdout,"[niik_image_rotate] flip in x\n");
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++,n++) {
          ii=img->nx-i-1;
          jj=j;
          kk=k;
          dimg[n]=niik_image_get_voxel(img,kk*xydim + jj*xdim + ii);
        }
      }
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
      free(dimg);
    }
    ni--;
    if(img->qform_code) {
      afmat[0]=
        niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,-1,1,1,0,0,0,ct.x,ct.y,ct.z),
                                afmat[0]);
    }
    if(img->sform_code) {
      afmat[1]=
        niikmat_multiply_free12(niikmat_affine_matrix_new(0,0,0,0,0,0,-1,1,1,0,0,0,ct.x,ct.y,ct.z),
                                afmat[1]);
    }
    if(verbose>1) niikmat_display_msg("flip x\n",afmat[0]);

    if(img->qform_code) {
      if(verbose>1) niikmat_display_msg("final qform\n",afmat[0]);
      /*nifti_disp_matrix_orient("old orientations:\n",img->qto_xyz);*/
      img->qto_xyz = niikmat_make_mat44(afmat[0]);
      /*niikmat_display(afmat[0]);
      nifti_disp_matrix_orient("new orientations:\n",img->qto_xyz);*/
      niikmat_inverse_update(afmat[0]);
      img->qto_ijk = niikmat_make_mat44(afmat[0]);
      nifti_mat44_to_quatern(img->qto_xyz,
                             &img->quatern_b,&img->quatern_c,&img->quatern_d,
                             &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                             NULL,NULL,NULL,
                             &img->qfac);
    }
    if(img->sform_code) {
      nifti_disp_matrix_orient("old orientations:\n",img->sto_xyz);
      img->sto_xyz = niikmat_make_mat44(afmat[1]);
      nifti_disp_matrix_orient("new orientations:\n",img->sto_xyz);
      niikmat_inverse_update(afmat[1]);
      img->sto_ijk = niikmat_make_mat44(afmat[1]);
    }
  } /* LAS */

  else {
    fprintf(stderr,"ERROR: unknown rotate string %s\n",rstr);
    free(dimg);
    return 0;
  }
  /* update image data */
  if(verbose) fprintf(stdout,"-v (niik_image_rotate) update image data\n");
  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERORR: nifti_k_set_voxel_values_vector(img,dimg)\n");
    free(dimg);
    return 0;
  }
  free(dimg);
  return 1;
}



int niik_image_restructure(nifti_image *img,char *rstr) {
  char fcname[32]="niik_image_restructure";
  int verbose=0;
  int
  nx,xy,
  n=0,i,j,k,ii,jj,kk;
  double *dimg=NULL;
  niikpt c,cc[2];
  niikmat *mat[2];

  if(verbose>=1)niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  nx=img->nx;
  xy=img->nx*img->ny;
  mat[0]=mat[1]=NULL;
  c.x=(img->nx-1.0)*img->dx/2.0;
  c.y=(img->ny-1.0)*img->dy/2.0;
  c.z=(img->nz-1.0)*img->dz/2.0;
  c.w=0;
  cc[0] = c; /*niikpt_affine_transform_m44(img->qto_xyz,c);*/
  cc[1] = c; /*niikpt_affine_transform_m44(img->sto_xyz,c);*/
  if(!strncmp(rstr,"+x+y+z",6)) {
    fprintf(stdout,"[%s] +x+y+z, no change\n",fcname);
    return 1;
  }
  NIIK_RET0(((dimg=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector");

  if(!strncmp(rstr,"-x+y+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=img->nx-1-i;
          jj=j;
          kk=k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    for(n=0; n<2; n++)
      mat[n] = niikmat_affine_matrix_new(0,0,0,
                                         0,0,0,
                                         -1,1,1,
                                         0,0,0,
                                         cc[n].x,cc[n].y,cc[n].z);
  }

  else if(!strncmp(rstr,"+x-y+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=img->ny-1-j;
          kk=k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    for(n=0; n<2; n++)
      mat[n] = niikmat_affine_matrix_new(0,0,0,
                                         0,0,0,
                                         1,-1,1,
                                         0,0,0,
                                         cc[n].x,cc[n].y,cc[n].z);
  }

  else if(!strncmp(rstr,"+x+y-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=j;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    for(n=0; n<2; n++)
      mat[n] = niikmat_affine_matrix_new(0,0,0,
                                         0,0,0,
                                         1,1,-1,
                                         0,0,0,
                                         cc[n].x,cc[n].y,cc[n].z);
  }

  else if(!strncmp(rstr,"-x-y+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=img->nx-1-i;
          jj=img->ny-1-j;
          kk=k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    for(n=0; n<2; n++)
      mat[n] = niikmat_affine_matrix_new(0,0,0,
                                         0,0,0,
                                         -1,-1,1,
                                         0,0,0,
                                         cc[n].x,cc[n].y,cc[n].z);
  }

  else if(!strncmp(rstr,"-x+y-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=img->nx-1-i;
          jj=j;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    for(n=0; n<2; n++)
      mat[0] = niikmat_affine_matrix_new(0,0,0,
                                         0,0,0,
                                         -1,1,-1,
                                         0,0,0,
                                         c.x,c.y,c.z);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+x-y-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=img->ny-1-j;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    mat[0] = niikmat_affine_matrix_new(0,0,0,
                                       0,0,0,
                                       1,-1,-1,
                                       0,0,0,
                                       c.x,c.y,c.z);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-x-y-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; i++) {
          ii=img->nx-1-i;
          jj=img->ny-1-j;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + jj*nx + ii);
        }
      }
    }
    mat[0] = niikmat_affine_matrix_new(0,0,0,
                                       0,0,0,
                                       -1,-1,-1,
                                       0,0,0,
                                       c.x,c.y,c.z);
    mat[1] = niikmat_copy(mat[0]);
  }


  else if(!strncmp(rstr,"+x+z+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_affine_matrix_val(1,0,0,0,
                                       0,0,1,0,
                                       0,1,0,0,
                                       0,0,0,1);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-x+z+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=img->nx-1-i;
          jj=k;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+x-z+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=img->nz-1-k;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+x+z-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-x-z+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-x+z-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+x-z-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-x-z-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(k=0; k<img->nz; k++) {
        for(i=0; i<img->nx; i++) {
          ii=i;
          jj=k;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + kk*nx + ii);
        }
      }
    }
    NIIK_ISWAP(&img->ny,&img->nz);
    NIIK_ISWAP(&img->   dim[2],&img->   dim[3]);
    NIIK_FSWAP(&img->dy,&img->dz);
    NIIK_FSWAP(&img->pixdim[2],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(1,0,0,0,
                                        0,0,1,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }




  else if(!strncmp(rstr,"+y+x+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=i;
          kk=k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_affine_matrix_val(0,1,0,0,
                                       1,0,0,0,
                                       0,0,1,0,
                                       0,0,0,1);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y+x+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y-x+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=img->nx-1-i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y+x-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y-x+z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=img->nx-1-i;
          kk=k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y+x-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y-x-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=img->nx-1-i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y-x-z",6)) {
    for(k=n=0; k<img->nz; k++) {
      for(i=0; i<img->nx; i++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=img->nx-1-i;
          kk=img->nz-1-k;
          dimg[n++]=niik_image_get_voxel(img, kk*xy + ii*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        1,0,0,0,
                                        0,0,1,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y+z+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=k;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_affine_matrix_val(0,1,0,0,
                                       0,0,1,0,
                                       1,0,0,0,
                                       0,0,0,1);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y+z+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=k;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y-z+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=img->nz-1-k;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y+z-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=k;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y-z+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=img->nz-1-k;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y+z-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=k;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+y-z-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=j;
          jj=img->nz-1-k;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-y-z-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          ii=img->ny-1-j;
          jj=img->nz-1-k;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, jj*xy + ii*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,1,0,0,
                                        0,0,1,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z+x+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=i;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_affine_matrix_val(
               0,0,1,0,
               1,0,0,0,
               0,1,0,0,
               0,0,0,1);
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z+x+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=i;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z-x+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=img->nx-1-i;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z+x-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=i;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z-x+y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=img->nx-1-i;
          kk=j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z+x-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=i;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z-x-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=img->nx-1-i;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z-x-y",6)) {
    for(j=n=0; j<img->ny; j++) {
      for(i=0; i<img->nx; i++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=img->nx-1-i;
          kk=img->ny-1-j;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + kk*nx + jj);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->ny);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[2]);
    NIIK_FSWAP(&img->dx,&img->dy);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[2]);
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        1,0,0,0,
                                        0,1,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }




  else if(!strncmp(rstr,"+z+y+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=j;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_affine_matrix_val(0,0,1,0,
                                       0,1,0,0,
                                       1,0,0,0,
                                       0,0,0,1);
    mat[1] = niikmat_copy(mat[0]);

  }

  else if(!strncmp(rstr,"-z+y+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=j;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z-y+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=img->ny-1-j;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z+y-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=j;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z-y+x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=img->ny-1-j;
          kk=i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z+y-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=j;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"+z-y-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=k;
          jj=img->ny-1-j;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else if(!strncmp(rstr,"-z-y-x",6)) {
    for(i=n=0; i<img->nx; i++) {
      for(j=0; j<img->ny; j++) {
        for(k=0; k<img->nz; k++) {
          ii=img->nz-1-k;
          jj=img->ny-1-j;
          kk=img->nx-1-i;
          dimg[n++]=niik_image_get_voxel(img, ii*xy + jj*nx + kk);
        }
      }
    }
    NIIK_ISWAP(&img->nx,&img->nz);
    NIIK_ISWAP(&img->   dim[1],&img->   dim[3]);
    NIIK_FSWAP(&img->dx,&img->dz);
    NIIK_FSWAP(&img->pixdim[1],&img->pixdim[3]);
    mat[0] = niikmat_multiply_free12
             (niikmat_affine_matrix_val(0,0,1,0,
                                        0,1,0,0,
                                        1,0,0,0,
                                        0,0,0,1),
              niikmat_affine_matrix_new(0,0,0,
                                        0,0,0,
                                        -1,-1,-1,
                                        0,0,0,
                                        c.x,c.y,c.z));
    mat[1] = niikmat_copy(mat[0]);
  }

  else {
    fprintf(stdout,"[%s] unknown restructure option, %s\n",fcname,rstr);
    return 0;
  }

  if(img->qform_code) {
    mat[0] = niikmat_multiply_free2(mat[0],
                                    niikmat_mat44_matrix(img->qto_xyz));
    img->qto_xyz = niikmat_make_mat44(mat[0]);
    img->qto_ijk = nifti_mat44_inverse(img->qto_xyz);
    nifti_mat44_to_quatern(img->qto_xyz,
                           &img->quatern_b,&img->quatern_c,&img->quatern_d,
                           &img->qoffset_x,&img->qoffset_y,&img->qoffset_z,
                           NULL,NULL,NULL,
                           &img->qfac);
  } /* qform */
  if(img->sform_code) {
    /*fprintf(stdout,"sform modification matrix:\n");
      niikmat_display(mat[1]);*/
    mat[1] = niikmat_multiply_free2(mat[1],
                                    niikmat_mat44_matrix(img->sto_xyz));
    img->sto_xyz = niikmat_make_mat44(mat[1]);
    img->sto_ijk = nifti_mat44_inverse(img->sto_xyz);
  } /* sform */
  niikmat_free(mat[0]);
  niikmat_free(mat[1]);
  NIIK_RET0((!niik_image_set_voxels_from_double_vector(img,dimg)),fcname,"niik_image_set_voxels_from_double_vector");
  free(dimg);
  if(verbose>=1)niik_fc_display(fcname,0);
  return 1;
} /* niik_image_restructure */


int niik_image_subsample(nifti_image *img,int sub_ratio) {
  nifti_image *tmpimg=NULL;
  int i,j,k,m,n;
  niikpt p;
  NIIK_RET0((img==NULL),__func__,"img is null");
  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_image_copy\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(!niik_image_resample_3d_update(img,
                                    img->dx*sub_ratio,
                                    img->dy*sub_ratio,
                                    img->dz*sub_ratio,
                                    img->nx/sub_ratio+1,
                                    img->ny/sub_ratio+1,
                                    img->nz/sub_ratio+1,
                                    NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
    return 0;
  }
  niik_image_clear(img);
  for(k=n=0; k<tmpimg->nz; k++) {
    p.z=k*tmpimg->dz;
    for(j=0; j<tmpimg->ny; j++) {
      p.y=j*tmpimg->dy;
      for(i=0; i<tmpimg->nx; i++,n++) {
        p.x=i*tmpimg->dx;
        m=niik_image_check_index_niikpt(img,p);
        niik_image_add_voxel(img,m,niik_image_get_voxel(tmpimg,n));
      }
    }
  }
  tmpimg=niik_image_free(tmpimg);
  return 1;
}

int niikmat_write_as_linear_xfm(char *fname,niikmat *mat) {
  FILE *fp=NULL;
  NIIK_RET0((mat==NULL),__func__,"mat is null");
  NIIK_RET0((fname==NULL),__func__,"fname is null");
  NIIK_RET0(((fp = fopen(fname,"w"))==NULL),__func__,"cannot open file");
  fprintf(fp,"MNI Transform File\n");
  fprintf(fp,"Transform_Type = Linear;\n");
  fprintf(fp,"Linear_Transform = \n");
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f\n",
          mat->m[0][0],mat->m[0][1],mat->m[0][2],mat->m[0][3]);
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f\n",
          mat->m[1][0],mat->m[1][1],mat->m[1][2],mat->m[1][3]);
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f ;\n",
          mat->m[2][0],mat->m[2][1],mat->m[2][2],mat->m[2][3]);
  fclose(fp);
  return 1;
}

int niikmat_convert_from_xfm(niikmat *m,nifti_image *srcimg,nifti_image *tgtimg)
/*  convert from xfm to niikmat : [target_cosines_and_origin]*[xfm]*inv(source_cosines_and_origin) */
{
  niikmat *sm=NULL;
  NIIK_RET0((m==NULL),__func__,"m is null");
  NIIK_RET0((srcimg==NULL),__func__,"srcimg is null");
  NIIK_RET0((tgtimg==NULL),__func__,"tgtimg is null");
  /* delta*inv(S2)*M*S*inv(delta) */
  NIIK_RET0(((sm = niikmat_mat44_matrix(srcimg->sto_xyz))==NULL),
            __func__,"niikkmat_mat44_matrix");
  NIIK_RET0((!niikmat_multiply_mat1_free2
             (sm,
              niikmat_scale_matrix(1.0/srcimg->dx,1.0/srcimg->dy,1.0/srcimg->dz))),__func__,"niikmat_multiply_mat1_free");
  NIIK_RET0((!niikmat_multiply_mat1_free2(m,sm)),__func__,"niikmat_multiply_mat1_free2");
  NIIK_RET0((!niikmat_multiply_mat2_free1(niikmat_mat44_matrix(tgtimg->sto_ijk),m)),__func__,"niikmat_multiply_mat2_free1");
  NIIK_RET0((!niikmat_multiply_mat2_free1(niikmat_scale_matrix(tgtimg->dx,tgtimg->dy,tgtimg->dz),m)),__func__,"niikmat_multiply_mat2_free1");
  return 1;
} /* niikmat_convert_from_xfm */

int niikmat_convert_to_xfm(niikmat *m,nifti_image *srcimg,nifti_image *tgtimg)
/*  convert from xfm to niikmat : inv([target_cosines_and_origin])*[xfm]*[source_cosines_and_origin]
 * 2013-10-10, Kunio Nakamura
 * -don't use sto_ijk --> sto_ijk suffers from floating point error due to the use of float32
 * -niikmat_inverse is more accurate with double (float64)
 */
{
  int verbose=0;
  const char *fcname="niikmat_convert_to_xfm";
  niikmat *sm=NULL;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((m==NULL),fcname,"m is null");
  NIIK_RET0((srcimg==NULL),fcname,"srcimg is null");
  NIIK_RET0((tgtimg==NULL),fcname,"tgtimg is null");
  if(verbose>0) {
    fprintf(stdout,"m[0] = ");
    niikmat_display(m);
  }

  NIIK_EXIT((niikmat_multiply_mat2_free1(niikmat_scale_matrix(1./tgtimg->dx,
                                         1./tgtimg->dy,
                                         1./tgtimg->dz),
                                         m)==0),
            fcname,"niikmat_multiply_mat2_free1 scale_matrix",9);

  if(verbose>0) {
    fprintf(stdout,"m[1] = ");
    niikmat_display(m);
  }

  NIIK_EXIT((niikmat_multiply_mat2_free1(niikmat_mat44_matrix(tgtimg->sto_xyz),
                                         m)==0),
            fcname,"niikmat_multiply_mat2_free1 tgt->sto_xyz",9);

  if(verbose>0) {
    fprintf(stdout,"m[2] = ");
    niikmat_display(m);
  }

  NIIK_EXIT((niikmat_multiply_mat1_free2(m,
                                         niikmat_scale_matrix(srcimg->dx,
                                             srcimg->dy,
                                             srcimg->dz))==0),
            fcname,"niikmat_multiply_mat2_free1 mov->scale",9);

  if(verbose>0) {
    fprintf(stdout,"m[3] = ");
    niikmat_display(m);
  }

  NIIK_RET0(((sm = niikmat_mat44_matrix(srcimg->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix");

  NIIK_RET0((!niikmat_inverse_update(sm)),fcname,"niikmat_inverse_matrix");

  NIIK_EXIT((!niikmat_multiply_mat1_free2(m,sm)),
            fcname,"niikmat_multiply_mat2_free1 mov->sto_ijk",9);

  if(verbose>0) {
    fprintf(stdout,"[%s] m = \n",fcname);
    niikmat_display(m);
  }
  return 1;
} /* niikmat_convert_to_xfm */

niikmat *niikmat_read_xfm(char *fname) {
  niikmat *m=NULL;
  const char *fcname="niikmat_read_xfm";
  char *cptr=NULL,asc[65536],c;
  FILE *fp=NULL;
  int n=0,verbose=0;
  NIIK_RET0((fname==NULL),fcname,"fname is null");
  if(verbose>=1) niik_fc_display(fcname,1);
  if(( fp = fopen(fname,"r") ) == NULL ) {
    fprintf(stderr,"[%s] ERROR: fopen %s\n",fcname,fname);
    return NULL;
  }
  for(n=0; n<65536; n++) asc[n]=0;
  n=0;
  while((c = getc(fp)) != EOF) {
    asc[n++] = c;
  }
  fclose(fp);
  if(verbose>=1) fprintf(stdout,"[[ %s ]]\n",asc);
  NIIK_RET0(((cptr = strstr(asc,"Linear_Transform ="))==NULL),fcname,"can't find \"Linear Transform =\" string");
  cptr += 19;
  if(verbose>=1) fprintf(stdout,"[[ %s ]]\n",cptr);
  m=niikmat_identity(4,4);
  sscanf(cptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &m->m[0][0],&m->m[0][1],&m->m[0][2],&m->m[0][3],
         &m->m[1][0],&m->m[1][1],&m->m[1][2],&m->m[1][3],
         &m->m[2][0],&m->m[2][1],&m->m[2][2],&m->m[2][3]);
  if(verbose>=1) niikmat_display(m);
  if(verbose>=1) niik_fc_display(fcname,0);
  return m;
}







/********************************************************************************
 *
 * get local things
 *
 ********************************************************************************/

int niik_image_get_neighbors_square_index(nifti_image *img,int *ijk,int *nei,int *nnei,int len)
/* len should be odd # like 3,5,7...
 */
{
  char fcname[64]="niik_image_get_neighbors_square_index";
  fprintf(stdout,"[%s] under development? --Can't remember (Kunio,2013-11-21)]\n",fcname);
  return 1;
} /* niik_image_get_neighbors_square_index */

int niik_image_get_neighbors_radius(nifti_image *img,nifti_image *mask,int *ijk,int *nei,int *nnei,double radius) {
  unsigned char
  *bimg=NULL;
  int
  verbose=0,
  mi,mj,mk,ni,nj,nk,
  i,j,k,n1,n2,n3;
  double r1,r2,r3;
  char fcname[32]="niik_image_get_neighbors";
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  if(mask!=NULL) {
    NIIK_RET0((mask->datatype!=NIFTI_TYPE_UINT8),fcname,"mask is not uint8");
  }
  mi=ijk[0]-(radius/img->dx)-1;
  mj=ijk[1]-(radius/img->dy)-1;
  mk=ijk[2]-(radius/img->dz)-1;
  ni=ijk[0]+(radius/img->dx)+1;
  nj=ijk[1]+(radius/img->dy)+1;
  nk=ijk[2]+(radius/img->dz)+1;
  mi=(mi<0)?0:mi;
  mj=(mj<0)?0:mj;
  mk=(mk<0)?0:mk;
  ni=(ni>=img->nx)?img->nx-1:ni;
  nj=(nj>=img->ny)?img->ny-1:nj;
  nk=(nk>=img->nz)?img->nz-1:nk;
  radius*=radius;
  *nnei=0;
  if(mask!=NULL) {
    /*NIIK_RET0(((bimg=niik_image_get_voxels_as_uint8_vector(mask))==NULL),fcname,"niik_image_get_voxels_as_uint8_vector");*/
    bimg=mask->data;
    for(k=mk; k<=nk; k++) {
      n3=k*img->nx*img->ny;
      r3=NIIK_SQ((k-ijk[2])*img->dz);
      for(j=mj; j<=nj; j++) {
        n2=j*img->nx;
        r2=NIIK_SQ((j-ijk[1])*img->dy);
        for(i=mi; i<=ni; i++) {
          n1=i+n2+n3;
          if(bimg[n1]==0) continue;
          r1=NIIK_SQ((i-ijk[0])*img->dx);
          if(r1+r2+r3>radius) continue;
          nei[*nnei]=n1;
          *nnei = *nnei + 1;
        }
      }
    }
  }  /* with mask */
  else {
    for(k=mk; k<=nk; k++) {
      n3=k*img->nx*img->ny;
      r3=NIIK_SQ((k-ijk[2])*img->dz);
      for(j=mj; j<=nj; j++) {
        n2=j*img->nx;
        r2=NIIK_SQ((j-ijk[1])*img->dy);
        for(i=mi; i<=ni; i++) {
          n1=i+n2+n3;
          r1=NIIK_SQ((i-ijk[0])*img->dx);
          if(r1+r2+r3>radius) continue;
          nei[*nnei]=n1;
          *nnei = *nnei + 1;
        }
      }
    }
  } /* withing mask */
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
}


/*
 * principal axes : eigenvectors
 */

#define NRC_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
int djacobi(double **a, int n, double d[], double **v, int *nrot) {
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  b=(double *)calloc(n+1,sizeof(double));
  z=(double *)calloc(n+1,sizeof(double));
  for (ip=1; ip<=n; ip++) {
    for (iq=1; iq<=n; iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1; ip<=n; ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1; i<=50; i++) {
    sm=0.0;
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free(z);
      free(b);
      return 1;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
            && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((double)(fabs(h)+g) == (double)fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=1; j<=ip-1; j++) {
            NRC_ROTATE(a,j,ip,j,iq);
          }
          for (j=ip+1; j<=iq-1; j++) {
            NRC_ROTATE(a,ip,j,j,iq);
          }
          for (j=iq+1; j<=n; j++) {
            NRC_ROTATE(a,ip,j,iq,j);
          }
          for (j=1; j<=n; j++) {
            NRC_ROTATE(v,j,ip,j,iq);
          }
          ++(*nrot);
        }
      }
    }
    for (ip=1; ip<=n; ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  fprintf(stderr,"[%s:%i:%s] ERROR: Too many iterations in routine jacobi\n",__FILE__,__LINE__,__func__);
  return 0;
}
#undef NRC_ROTATE

int niikmat_jacobi(niikmat *m,niikvec *eigval,niikmat *eigvec) {
  double **a=NULL;
  double *d;
  double **v;
  int i,j,n,nrot;
  NIIK_RET0((m==NULL),__func__,"m is null");
  NIIK_RET0((eigval==NULL),__func__,"eigval is null");
  NIIK_RET0((eigvec==NULL),__func__,"eigvec is null");
  NIIK_RET0((m->row!=m->col),__func__,"not a square matrix,m");
  NIIK_RET0((eigvec->row!=eigvec->col),__func__,"not a square matrix, eigvec");
  NIIK_RET0((eigvec->row!=eigval->num),__func__,"length mismatch");
  n=m->row;
  a=(double **)calloc(n+1,sizeof(double *));
  v=(double **)calloc(n+1,sizeof(double *));
  for(i=1; i<=n; i++) {
    a[i]=(double *)calloc(n+1,sizeof(double));
    v[i]=(double *)calloc(n+1,sizeof(double));
    for(j=1; j<=n; j++) {
      a[i][j]=m->m[i-1][j-1];
    }
  }
  d=(double *)calloc(n+1,sizeof(double));
  //fprintf(stdout,"running djacobi\n");
  if(!djacobi(a,n,d,v,&nrot)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: djacobi\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  //fprintf(stdout,"finished djacobi\n");
  for(i=1; i<=n; i++) {
    for(j=1; j<=n; j++) {
      eigvec->m[i-1][j-1]=v[i][j];
    }
    eigval->v[i-1]=d[i];
    free(v[i]);
    free(a[i]);
  }
  free(a);
  free(v);
  free(d);
  return 1;
}



int niik_image_principal_axes(nifti_image *img,nifti_image *maskimg,niikmat *eigvec,niikvec *eigval) {
  char fcname[32]="niik_image_principal_axes";
  int i,j,k;
  int verbose=0;
  niikpt c,p;
  niikmat *cov;
  double d,w=0;
  if(verbose>0) niik_fc_display(fcname,1);
  c=niikpt_zero();
  NIIK_RET0((!niikpt_image_get_centroid_update(img,maskimg,&c)),__func__,"niikpt_image_get_centroid_update");
  //c.w=0; niikpt_disp(c);
  cov=niikmat_init(3,3);
  for(k=0; k<img->nz; k++) {
    p.z=k*img->dz;
    for(j=0; j<img->ny; j++) {
      p.y=j*img->dy;
      for(i=0; i<img->nx; i++) {
        p.x=i*img->dx;
        if(maskimg!=NULL) if(niik_image_get_voxel(maskimg,i+j*img->nx+k*img->nx*img->ny)==0.0) continue;
        d=niik_image_get_voxel(img,i+j*img->nx+k*img->nx*img->ny);
        cov->m[0][0] += d * (p.x-c.x) * (p.x-c.x);
        cov->m[0][1] += d * (p.x-c.x) * (p.y-c.y);
        cov->m[0][2] += d * (p.x-c.x) * (p.z-c.z);
        cov->m[1][1] += d * (p.y-c.y) * (p.y-c.y);
        cov->m[1][2] += d * (p.y-c.y) * (p.z-c.z);
        cov->m[2][2] += d * (p.z-c.z) * (p.z-c.z);
        w+=d;
      }
    }
  }
  cov->m[1][0]=cov->m[0][1];
  cov->m[2][0]=cov->m[0][2];
  cov->m[2][1]=cov->m[1][2];
  NIIK_RET0((w==0.0),fcname,"div by 0 (w)");
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      cov->m[i][j]/=w;
  if(verbose>0) niikmat_display_msg("COV matrix = \n",cov);
  NIIK_RET0((!niikmat_jacobi(cov,eigval,eigvec)),__func__,"niikmat_jacobi");

  if(fabs(eigvec->m[0][0])>fabs(eigvec->m[0][1]) && fabs(eigvec->m[0][0])>fabs(eigvec->m[0][2])) {}
  else {
    for(i=0; i<3; i++) NIIK_DSWAP(&eigvec->m[0][i],&eigvec->m[1][i]);
    NIIK_DSWAP(&eigval->v[0],&eigval->v[1]);
  }
  if(fabs(eigvec->m[1][1])>fabs(eigvec->m[1][0]) && fabs(eigvec->m[1][1])>fabs(eigvec->m[1][2])) {}
  else {
    for(i=0; i<3; i++) NIIK_DSWAP(&eigvec->m[1][i],&eigvec->m[2][i]);
    NIIK_DSWAP(&eigval->v[1],&eigval->v[2]);
  }

  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
}


/* MINC-style file history*/
int niik_image_append_history(nifti_image *img,const char *history_entry) {
  size_t history_entry_len;

  if(history_entry==NULL) return 1;

  history_entry_len=strlen(history_entry);
  if(history_entry_len==0) return 1;

  /*VF: disabling to be compatible with standard nifti*/
  /*if(img->minc_history != NULL ) {
    size_t history_length=strlen(img->minc_history);
    img->minc_history=(char*)realloc(img->minc_history,history_length+history_entry_len+1);
    strcat(img->minc_history, history_entry);
  } else {
    img->minc_history=strdup(history_entry);
  }*/
  return 1;
}

char* niik_create_minc_timestamp(int argc,char *argv[]) {
  char *timestamp;
  char cur_time[200];
  time_t t;
  struct tm *tmp;
  size_t str_len;
  int i;

  t = time(NULL);
  tmp = localtime(&t);

  strftime(cur_time, sizeof(cur_time), "%a %b %d %T %Y>>>", tmp);
  /* Get the time, overwriting newline */
  str_len=strlen(cur_time);
  for (i=0; i<argc; i++)
    str_len+=strlen(argv[i])+1;

  timestamp=(char*)malloc(str_len+3);
  strcpy(timestamp,cur_time);

  /* Copy the program name and arguments */
  for (i=0; i<argc; i++) {
    strcat(timestamp,argv[i]);
    strcat(timestamp," ");
  }

  strcat(timestamp,"\n");
  return timestamp;
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/