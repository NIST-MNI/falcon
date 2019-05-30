#include <sys/types.h>
#include <ctype.h>
#include "mri.h"
#include "tags.h"
#include "matrix.h"

#define UNUSED_SPACE_SIZE 256
#define USED_SPACE_SIZE   (3*sizeof(float)+4*3*sizeof(float))

#define NO_ERROR              0
#define ERROR_NONE            NO_ERROR
#define ERROR_NO_FILE         -1
#define ERROR_NOFILE          ERROR_NO_FILE
#define ERROR_NO_MEMORY       -2
#define ERROR_NOMEMORY        ERROR_NO_MEMORY
#define ERROR_UNSUPPORTED     -3
#define ERROR_BADPARM         -4
#define ERROR_BAD_PARM        ERROR_BADPARM
#define ERROR_BADFILE         -5
#define ERROR_BAD_FILE        ERROR_BADFILE
#define ERROR_SIZE            -6
#define ERROR_BADLOOP         -7
#define ERROR_OUT_OF_BOUNDS   -8

typedef union {
  long32  l ;
  float f ;
  int   i ;
  char  buf[4] ;
  short s[2] ;
} SWAP_LONG32 ;

typedef union {
  short  s ;
  char   buf[sizeof(short)] ;
} SWAP_SHORT ;



int Gerror = NO_ERROR ;

static int (*myclose)(FILE *stream);
static int (*error_vfprintf)(FILE *fp,const char *fmt,va_list args) = vfprintf;
static void
int_local_buffer_to_image(int *buf, MRI *mri, int slice, int frame) {
  int           y, width, height ;
  int           *pslice ;
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++) {
    pslice = &MRIIseq_vox(mri, 0, y, slice, frame) ;
    memmove(pslice, buf, width*sizeof(int)) ;
    buf += width ;
  }
}
#if 0
static void
image_to_int_buffer(int *buf, MRI *mri, int slice) {
  int y, x, width, height, depth ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++) {
    if (mri->type == MRI_UCHAR) {
      for (x = 0; x < depth; x++)
        buf[x] = (int)MRIvox(mri, x, y, slice);
    } else if (mri->type == MRI_SHORT) {
      for (x = 0; x < depth; x++)
        buf[x] = (int)MRISvox(mri, x, y, slice);
    } else if (mri->type == MRI_LONG) {
      for (x = 0; x < depth; x++)
        buf[x] = (int)MRILvox(mri, x, y, slice);
    } else if (mri->type == MRI_FLOAT) {
      for (x = 0; x < depth; x++)
        buf[x] = (int)MRIFvox(mri, x, y, slice);
    } else {
      memmove(buf, mri->slices[slice][y], width*sizeof(int)) ;
    }
    buf += width ;
  }
}
static void
image_to_long_buffer(long *buf, MRI *mri, int slice) {
  int y, x, width, height, depth ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++) {
    if (mri->type == MRI_UCHAR) {
      for (x = 0; x < depth; x++)
        buf[x] = (long)MRIvox(mri, x, y, slice);
    } else if (mri->type == MRI_INT) {
      for (x = 0; x < depth; x++)
        buf[x] = (long)MRIIvox(mri, x, y, slice);
    } else if (mri->type == MRI_SHORT) {
      for (x = 0; x < depth; x++)
        buf[x] = (long)MRISvox(mri, x, y, slice);
    } else if (mri->type == MRI_FLOAT) {
      for (x = 0; x < depth; x++)
        buf[x] = (long)MRIFvox(mri, x, y, slice);
    } else {
      memmove(buf, mri->slices[slice][y], width*sizeof(long)) ;
    }
    buf += width ;
  }
}
static void
image_to_float_buffer(float *buf, MRI *mri, int slice) {
  int y, x, width, height, depth ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth;
  for (y=0; y < height ; y++) {
    if (mri->type == MRI_UCHAR) {
      for (x = 0; x < depth; x++)
        buf[x] = (float)MRIvox(mri, x, y, slice);
    } else if (mri->type == MRI_INT) {
      for (x = 0; x < depth; x++)
        buf[x] = (float)MRIIvox(mri, x, y, slice);
    } else if (mri->type == MRI_LONG) {
      for (x = 0; x < depth; x++)
        buf[x] = (float)MRILvox(mri, x, y, slice);
    } else if (mri->type == MRI_SHORT) {
      for (x = 0; x < depth; x++)
        buf[x] = (float)MRISvox(mri, x, y, slice);
    } else {
      memmove(buf, mri->slices[slice][y], width*sizeof(float)) ;
    }
    buf += width ;
  }
}
#endif
static void
long32_local_buffer_to_image(long32 *buf, MRI *mri, int slice, int frame) {
  int           y, width, height ;
  long32          *pslice ;
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++) {
    pslice = &MRILseq_vox(mri, 0, y, slice, frame) ;
    memmove(pslice, buf, width*sizeof(long)) ;
    buf += width ;
  }
}
static void
float_local_buffer_to_image(float *buf, MRI *mri, int slice, int frame) {
  int           y, width, height ;
  float         *pslice ;
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++) {
    pslice = &MRIFseq_vox(mri, 0, y, slice, frame) ;
    memmove(pslice, buf, width*sizeof(float)) ;
    buf += width ;
  }
}
static void
short_local_buffer_to_image(short *buf, MRI *mri, int slice, int frame) {
  int           y, width, height ;
  short         *pslice ;
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++) {
    pslice = &MRISseq_vox(mri, 0, y, slice, frame) ;
    memmove(pslice, buf, width*sizeof(short)) ;
    buf += width ;
  }
}
static void
local_buffer_to_image(BUFTYPE *buf, MRI *mri, int slice, int frame) {
  int           y, width, height ;
  BUFTYPE       *pslice ;
  width = mri->width ;
  height = mri->height ;
  for (y = 0 ; y < height ; y++) {
    pslice = &MRIseq_vox(mri, 0, y, slice, frame) ;
    memmove(pslice, buf, width*sizeof(BUFTYPE)) ;
    buf += width ;
  }
}


int toupper(int c) {
  if(isalpha(c)) return c-32;
  return 0;
}
char *StrUpper(char *str) {
  char *cp ;
  for (cp = str ; *cp ; cp++)   *cp = (char)toupper(*cp) ;
  return(str) ;
}

int stricmp(const char *str1,const  char *str2) {
  char buf1[STR_LEN], buf2[STR_LEN] ;
  strcpy(buf1, str1) ;
  strcpy(buf2, str2) ;
  StrUpper(buf1) ;
  StrUpper(buf2) ;
  return(strcmp(buf1, buf2)) ;
}

short
swapShort(short s) {
  SWAP_SHORT ss ;
  char       c ;
  /* first swap bytes in word */
  ss.s = s ;
  c = ss.buf[0] ;
  ss.buf[0] = ss.buf[1] ;
  ss.buf[1] = c ;
  return(ss.s) ;
}

int
swapInt(int i) {
  SWAP_LONG32  sl ;
  short      s ;
  /* first swap bytes in each word */
  sl.i = i ;
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;
  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;
  return(sl.i) ;
}

long64 swapLong64(long64 l) {
  size_t typeSize = sizeof(long64);
  char *pVar = (char *) (&l);
  char tmp;
  double w;
  size_t i;
  for (i=0; i < typeSize/2; ++i) { // typeSize must be even
    // swap front and back
    tmp = *(pVar+2*i);
    *(pVar+2*i) = *(pVar+typeSize-1-2*i);
    *(pVar+typeSize-1-2*i) = tmp;
  }
  w = *((long64 *)(pVar)); // copy
  return w;
}

long long
freadLong(FILE *fp) {
  int  nread ;
  long long i ;
  nread = fread(&i,sizeof(long long),1,fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapLong64(i) ;
#endif
  return(i) ;
}

int freadInt(FILE *fp) {
  int  i, nread ;
  nread = fread(&i,sizeof(int),1,fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapInt(i) ;
#endif
  return(i) ;
}

short freadShort(FILE *fp) {
  int   nread ;
  short s ;
  nread = fread(&s,sizeof(short),1,fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s) ;
#endif
  if (nread != 1) {
    fprintf(stderr,"ERROR: freadShort: fread failed\n") ;
    exit(0);
  }
  return(s) ;
}

int byteswapbuffloat(void *buf, long int nbufbytes) {
  register char *cbuf,c;
  register long int n, nmax;
  nmax = nbufbytes;
  cbuf = (char *)buf;
  for (n=0; n<nmax; n+=4) {
    c = *cbuf;
    *cbuf = *(cbuf+3);
    *(cbuf+3) = c;
    c = *(cbuf+1);
    *(cbuf+1) = *(cbuf+2);
    *(cbuf+2) = c;
    cbuf += 4;
  }
  return(0);
}


float
swapFloat(float f) {
  SWAP_LONG32  sl ;
  short      s ;
  /* first swap bytes in each word */
  sl.f = f ;
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;
  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;
  return(sl.f) ;
}

float freadFloat(FILE *fp) {
  char  buf[4];
  float f;
  int   ret ;
  ret = fread(buf,4,1,fp);
  //ret = fread(&f,4,1,fp); // old way
  if (ret != 1) {
    fprintf(stderr, "ERROR: freadFloat: fread failed\n") ;
    exit(0);
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf,1);
  //f = swapFloat(f);  // old way
#endif
//error: dereferencing type-punned pointer will break strict-aliasing rules:
//  f = *((float*)buf);
  memcpy(&f,&buf,sizeof(float));
  return(f) ;
}
int freadFloatEx(float *pf, FILE *fp) {
  int   ret ;
  ret = fread(pf,sizeof(float),1,fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pf = swapFloat(*pf) ;
#endif
  return ret;
}

int freadIntEx(int *pi, FILE *fp) {
  int nread ;
  nread = fread(pi,sizeof(int),1,fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pi = swapInt(*pi) ; /* swapInt(int i) */
#endif
  return(nread);
}

// String copy will allocation.
char *strcpyalloc(const char *str) {
  char *cpstr;
  cpstr = (char *) calloc(strlen(str)+1,sizeof(char));
  strcpy(cpstr,str);
  return(cpstr);
}

char *fio_dirname(const char *pathname) {
  int l,n;
  char *dirname;
  if (pathname == NULL) return(NULL);
  char *pname = strcpyalloc(pathname);
  l = strlen(pname);
  /* strip off leading forward slashes */
  while (l > 0 && pname[l-1] == '/') {
    pname[l-1] = '\0';
    l = strlen(pname);
  }
  if (l < 2) {
    /* pname is / or . or single character */
    free(pname);
    dirname = (char *) calloc(2,sizeof(char));
    if (l==0 || pname[0] == '/') dirname[0] = '/';
    else                           dirname[0] = '.';
    return(dirname);
  }
  /* Start at the end of the path name and step back
     until a forward slash is found */
  for (n=l; n >= 0; n--)if (pname[n] == '/') break;
  if (n < 0) {
    /* no forward slash found */
    dirname = (char *) calloc(2,sizeof(char));
    dirname[0] = '.';
    free(pname);
    return(dirname);
  }
  if (n == 0) {
    /* first forward slash is the first character */
    dirname = (char *) calloc(2,sizeof(char));
    dirname[0] = '/';
    free(pname);
    return(dirname);
  }
  dirname = (char *) calloc(n+1,sizeof(char));
  memmove(dirname,pname,n);
  free(pname);
  return(dirname);
}


int
FileExists(const char *fname) {
  FILE *fp ;
  int old_errno;
  old_errno = errno;
  fp = fopen(fname, "r") ;
  if (fp)
    fclose(fp) ;
  else
    errno = old_errno;
  return(fp != NULL) ;
}


int TAGreadStart(FILE *fp, long long *plen) {
  int  tag ;
  tag = freadInt(fp) ;
  if (feof(fp))
    return(0) ;
  switch (tag) {
  case TAG_OLD_MGH_XFORM:
    *plen = (long long)freadInt(fp) ;  /* sorry - backwards compatibility
                                          with Tosa's stuff */
    *plen = *plen -1 ; // doesn't include null
    break ;
  case TAG_OLD_SURF_GEOM:    // these don't take lengths at all
  case TAG_OLD_USEREALRAS:
  case TAG_OLD_COLORTABLE:
    *plen = 0 ;
    break ;
  default:
    *plen = freadLong(fp) ;
  }
  return(tag) ;
}


int stuff_four_by_four(MATRIX *m, float m11, float m12, float m13, float m14,
                       float m21, float m22, float m23, float m24,
                       float m31, float m32, float m33, float m34,
                       float m41, float m42, float m43, float m44) {
  if (m == NULL) {
    fprintf(stderr,"ERROR: stuff_four_by_four(): matrix is NULL");
    return 0;
  }
  if (m->rows != 4 || m->cols != 4) {
    fprintf(stderr,"stuff_four_by_four(): matrix is not four-by-four");
    return 0;
  }
  *MATRIX_RELT(m, 1, 1) = m11;
  *MATRIX_RELT(m, 1, 2) = m12;
  *MATRIX_RELT(m, 1, 3) = m13;
  *MATRIX_RELT(m, 1, 4) = m14;
  *MATRIX_RELT(m, 2, 1) = m21;
  *MATRIX_RELT(m, 2, 2) = m22;
  *MATRIX_RELT(m, 2, 3) = m23;
  *MATRIX_RELT(m, 2, 4) = m24;
  *MATRIX_RELT(m, 3, 1) = m31;
  *MATRIX_RELT(m, 3, 2) = m32;
  *MATRIX_RELT(m, 3, 3) = m33;
  *MATRIX_RELT(m, 3, 4) = m34;
  *MATRIX_RELT(m, 4, 1) = m41;
  *MATRIX_RELT(m, 4, 2) = m42;
  *MATRIX_RELT(m, 4, 3) = m43;
  *MATRIX_RELT(m, 4, 4) = m44;
  return(NO_ERROR);
} /* end stuff_four_by_four() */


MATRIX *extract_r_to_i(const MRI *mri) {
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras ;
  m_voxel_to_ras = extract_i_to_r(mri) ;
  m_ras_to_voxel = MatrixInverse(m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ;
  return(m_ras_to_voxel) ;
}
MATRIX *extract_i_to_r(const MRI *mri) {
  MATRIX *m;
  float m11, m12, m13, m14;
  float m21, m22, m23, m24;
  float m31, m32, m33, m34;
  float ci, cj, ck;
  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if(m == NULL) {
    fprintf(stderr,"ERROR: extract_i_to_r(): error allocating matrix");
    exit(0);
  }
  m11 = mri->xsize * mri->x_r;
  m12 = mri->ysize * mri->y_r;
  m13 = mri->zsize * mri->z_r;
  m21 = mri->xsize * mri->x_a;
  m22 = mri->ysize * mri->y_a;
  m23 = mri->zsize * mri->z_a;
  m31 = mri->xsize * mri->x_s;
  m32 = mri->ysize * mri->y_s;
  m33 = mri->zsize * mri->z_s;
  ci = (mri->width) / 2.0;
  cj = (mri->height) / 2.0;
  ck = (mri->depth) / 2.0;
  m14 = mri->c_r - (m11 * ci + m12 * cj + m13 * ck);
  m24 = mri->c_a - (m21 * ci + m22 * cj + m23 * ck);
  m34 = mri->c_s - (m31 * ci + m32 * cj + m33 * ck);
  stuff_four_by_four(m, m11, m12, m13, m14,
                     m21, m22, m23, m24,
                     m31, m32, m33, m34,
                     0.0, 0.0, 0.0, 1.0);
  return(m);
} /* end extract_i_to_r() */

MATRIX *
MatrixAlloc( const int rows, const int cols, const int type) {
  MATRIX *mat ;
  int    row, nelts;
#ifdef _POSIX_MAPPED_FILES
  int    i;
  float  f;
#endif
  mat = (MATRIX *)calloc(1, sizeof(MATRIX)) ;
  if (!mat) {
    fprintf(stderr,"MatrixAlloc(%d, %d, %d): could not allocate mat",
            rows, cols, type) ;
    exit(0);
  }
  mat->rows = rows ;
  mat->cols = cols ;
  mat->type = type ;
  /*
    allocate a single array the size of the matrix, then initialize
    row pointers to the proper locations.
  */
  nelts = rows*cols ;
  if (type == MATRIX_COMPLEX)
    nelts *= 2 ;
  /*
    because NRC is one-based, we must leave room for a few unused
    (but valid) addresses before the start of the actual data so
    that mat->rptr[0][0] is a valid address even though it wont
    be used.
  */
  mat->data = (float *)calloc(nelts+2, sizeof(float)) ;
  mat->mmapfile = NULL;
#ifdef _POSIX_MAPPED_FILES
  if (!mat->data) { /* First try to allocate a mmap'd tmpfile */
    printf("MatrixAlloc(%d, %d): Using mmap'd tmpfile\n",
           rows, cols) ;
    if ((mat->mmapfile = tmpfile())) {
      /* This maintains identical behavior with calloc */
      f = 0;
      for (i = 0; i < nelts+2; ++i) {
        if (!(fwrite(&f, sizeof(float), 1, mat->mmapfile))) {
          printf("MatrixAlloc(%d, %d): fwrite failed", rows, cols) ;
          exit(1) ;
        }
      }
      /* This seems to matter with some implementations of mmap */
      fseek(mat->mmapfile, 0, 0) ;
      /* lseek(fileno(mat->mapfile), (nelts+2) * sizeof(float), 0) ;*/
      fflush(mat->mmapfile) ;
      mat->data = (float *)mmap(0, (nelts+2) * sizeof(float),
                                PROT_READ | PROT_WRITE, MAP_SHARED,
                                fileno(mat->mmapfile), 0) ;
      if (mat->data == MAP_FAILED) {
        mat->data = 0 ;
      }
    }
  }
#endif
  if (!mat->data) { /* we _still_ couldn't get it! */
    fprintf(stderr, "MatrixAlloc(%d, %d): allocation failed\n",
            rows, cols) ;
    exit(1) ;
  }
  mat->data += 2 ;
  /*
     silly numerical recipes in C requires 1-based stuff. The full
     data array is zero based, point the first row to the zeroth
     element, and so on.
  */
  mat->rptr = (float **)calloc(rows+1, sizeof(float *)) ;
  if (!mat->rptr) {
    free(mat->data) ;
    free(mat) ;
    fprintf(stderr,"ERROR: MatrixAlloc(%d, %d): could not allocate rptr",
            rows, cols) ;
    exit(0);
  }
  for (row = 1 ; row <= rows ; row++) {
    switch (type) {
    case MATRIX_REAL:
      mat->rptr[row] = mat->data + (row-1)*cols - 1 ;
      break ;
    case MATRIX_COMPLEX:
      mat->rptr[row] = (float *)(((CPTR)mat->data) +
                                 (row-1)*cols - 1) ;
      break ;
    default:
      fprintf(stderr,"ERROR: MatrixAlloc: unknown type %d\n",type);
      exit(0);
    }
  }
  return(mat) ;
}


int
MatrixFree(MATRIX **pmat) {
  MATRIX *mat ;
  if (!pmat) {
    fprintf(stderr,"ERROR: MatrixFree: NULL pmat POINTER!\n");
    return 0;
  }
  mat = *pmat ;
  *pmat = NULL;
  if (!mat) {
    fprintf(stderr,"ERROR: MatrixFree: NULL pmat POINTER!\n");
    return 0;
  }
  /* silly numerical recipes in C requires 1-based stuff */
  mat->data -= 2 ;
  if (mat->mmapfile) {
    int nelts ;
    nelts = mat->rows*mat->cols ;
    if (mat->type == MATRIX_COMPLEX)
      nelts *= 2 ;
#ifdef _POSIX_MAPPED_FILES
    munmap((void *) mat->data, (nelts+2) * sizeof(float)) ;
#endif
    fclose(mat->mmapfile) ;
  } else {
    free(mat->data) ;
  }
  free(mat->rptr) ;
  free(mat) ;
  return(0) ;
}



/*!
  \fn MRI *MRIallocHeader(int width, int height, int depth, int type, int nframes)
  \brief allocate an MRI data structure but not space for  the image data
*/
MRI *MRIallocHeader(int width, int height, int depth, int type, int nframes) {
  MRI  *mri ;
  int  i ;
  mri = (MRI *)calloc(1, sizeof(MRI)) ;
  if (!mri) {
    fprintf(stderr, "ERROR: MRIalloc: could not allocate MRI\n") ;
    exit(0);
  }
  // Note: changes here may need to be reflected in MRISeqchangeType()
  mri->frames = (MRI_FRAME *)calloc(nframes, sizeof(MRI_FRAME)) ;
  if (!mri->frames) {
    fprintf(stderr,"MRIalloc: could not allocate %d frame\n", nframes) ;
    exit(0);
  }
  for (i = 0 ; i < mri->nframes ; i++)
    mri->frames[i].m_ras2vox = MatrixAlloc(4,4, MATRIX_REAL) ;
  mri->imnr0 = 1 ;
  mri->imnr1 = depth;
  mri->fov = width ;
  mri->thick = 1.0 ;
  mri->scale = 1 ;
  mri->roi.dx = mri->width = width ;
  mri->roi.dy = mri->height = height ;
  mri->roi.dz = mri->depth = depth ;
  mri->yinvert = 1 ;
  mri->xsize = mri->ysize = mri->zsize = 1 ;
  mri->type = type ;
  mri->nframes = nframes ;
  mri->xi = mri->yi = mri->zi = NULL ;
  mri->slices = NULL ;
  mri->ps = 1 ;
  mri->xstart = -mri->width/2.0 ;
  mri->xend = mri->width/2.0 ;
  mri->ystart = -mri->height/2.0 ;
  mri->yend = mri->height/2.0 ;
  mri->zstart = -mri->depth/2.0 ;
  mri->zend = mri->depth/2 ;
  //
  mri->x_r = -1;
  mri->x_a = 0.;
  mri->x_s = 0.;
  //
  mri->y_r = 0.;
  mri->y_a = 0.;
  mri->y_s = -1;
  //
  mri->z_r = 0.;
  mri->z_a = 1.;
  mri->z_s = 0.;
  //
  mri->c_r = mri->c_a = mri->c_s = 0.0;
  mri->ras_good_flag = 0;
  mri->brightness = 1;
  mri->register_mat = NULL;
  mri->subject_name[0] = '\0';
  mri->path_to_t1[0] = '\0';
  mri->fname_format[0] = '\0';
  mri->gdf_image_stem[0] = '\0';
  mri->tag_data = NULL;
  mri->tag_data_size = 0;
  mri->transform_fname[0] = '\0';
  mri->pedir = NULL;
  MATRIX *tmp;
  tmp = extract_i_to_r(mri);
  AffineMatrixAlloc( &(mri->i_to_r__ ) );
  SetAffineMatrix( mri->i_to_r__, tmp );
  MatrixFree( &tmp );
  mri->r_to_i__ = extract_r_to_i(mri);
  mri->AutoAlign = NULL;
  // Chunking memory management
  mri->ischunked = 0;
  mri->chunk = NULL;
  mri->bytes_per_vox   = MRIsizeof(type);
  // These things are explicitly set to 0 here because we
  // do not yet know the true number of frames, and they
  // are set in MRIallocChunk().
  mri->bytes_per_row   = 0;
  mri->bytes_per_slice = 0;
  mri->bytes_per_vol   = 0;
  mri->bytes_total     = 0;
  mri->bvals = NULL; // For DWI
  mri->bvecs = NULL;
  return(mri) ;
}

static MRI *mghRead(char *fname, int read_volume, int frame) {
  MRI  *mri ;
  FILE  *fp = 0;
  int   start_frame, end_frame, width, height, depth, nframes, type, x, y, z,
        bpv, dof, bytes, version, ival, unused_space_size, good_ras_flag, i ;
  BUFTYPE *buf ;
  char   unused_buf[UNUSED_SPACE_SIZE+1] ;
  float  fval, xsize, ysize, zsize, x_r, x_a, x_s, y_r, y_a, y_s,
         z_r, z_a, z_s, c_r, c_a, c_s ;
  short  sval ;
  //  int tag_data_size;
  char *ext;
  int gzipped=0;
  char command[STRLEN];
  int nread;
  int tag;

  ext = strrchr(fname, '.') ;
  if (ext) {
    ++ext;
    // if mgz, then it is compressed
    if (!stricmp(ext, "mgz") || strstr(fname, "mgh.gz")) {
      gzipped = 1;
      myclose = pclose;  // assign function pointer for closing
#ifdef Darwin
      // zcat on Max OS always appends and assumes a .Z extention,
      // whereas we want .m3z
      strcpy(command, "gunzip -c ");
#else
      strcpy(command, "zcat ");
#endif
      strcat(command, fname);

      errno = 0;
      fp = popen(command, "r");
      if (!fp) {
        errno = 0;
        fprintf(stderr,"ERROR: mghRead: encountered error executing: '%s',frame %d\n",command,frame) ;
        return NULL;
      }
      if (errno) {
        pclose(fp);
        errno = 0;
        fprintf(stderr,"ERROR: mghRead: encountered error executing: '%s',frame %d\n",
                command,frame);
        return NULL;
      }
    } else if (!stricmp(ext, "mgh")) {
      myclose = fclose; // assign function pointer for closing
      fp = fopen(fname, "rb") ;
      if (!fp) {
        errno = 0;
        fprintf(stderr,"ERROR: mghRead(%s, %d): could not open file\n", fname, frame) ;
        return NULL;
      }
    }
  }

  /* keep the compiler quiet */
  xsize = ysize = zsize = 0;
  x_r = x_a = x_s = 0;
  y_r = y_a = y_s = 0;
  z_r = z_a = z_s = 0;
  c_r = c_a = c_s = 0;

  nread = freadIntEx(&version, fp) ;
  if (!nread) {
    fprintf(stderr,"ERROR: mghRead(%s, %d): read error\n",
            fname, frame) ;
    return NULL;
  }

  width = freadInt(fp) ;
  height = freadInt(fp) ;
  depth =  freadInt(fp) ;
  nframes = freadInt(fp) ;
  type = freadInt(fp) ;
  dof = freadInt(fp) ;

  unused_space_size = UNUSED_SPACE_SIZE-sizeof(short) ;

  good_ras_flag = freadShort(fp) ;
  if (good_ras_flag > 0) {    /* has RAS and voxel size info */
    unused_space_size -= USED_SPACE_SIZE ;
    xsize = freadFloat(fp) ;
    ysize = freadFloat(fp) ;
    zsize = freadFloat(fp) ;

    x_r = freadFloat(fp) ;
    x_a = freadFloat(fp) ;
    x_s = freadFloat(fp) ;
    y_r = freadFloat(fp) ;
    y_a = freadFloat(fp) ;
    y_s = freadFloat(fp) ;

    z_r = freadFloat(fp) ;
    z_a = freadFloat(fp) ;
    z_s = freadFloat(fp) ;
    c_r = freadFloat(fp) ;
    c_a = freadFloat(fp) ;
    c_s = freadFloat(fp) ;
  }
  /* so stuff can be added to the header in the future */
  fread(unused_buf, sizeof(char), unused_space_size, fp) ;

  switch (type) {
  default:
  case MRI_FLOAT:
    bpv = sizeof(float) ;
    break ;
  case MRI_UCHAR:
    bpv = sizeof(char)  ;
    break ;
  case MRI_SHORT:
    bpv = sizeof(short) ;
    break ;
  case MRI_INT:
    bpv = sizeof(int) ;
    break ;
  case MRI_TENSOR:
    bpv = sizeof(float) ;
    nframes = 9 ;
    break ;
  }
  bytes = width * height * bpv ;  /* bytes per slice */
  if(!read_volume) {
    mri = MRIallocHeader(width, height, depth, type, 1) ;
    mri->dof = dof ;
    mri->nframes = nframes ;
    if(gzipped) { // pipe cannot seek
      int count;
      for (count=0; count < mri->nframes*width*height*depth*bpv; count++)
        fgetc(fp);
    } else
      fseek(fp, mri->nframes*width*height*depth*bpv, SEEK_CUR) ;
  } else {
    if (frame >= 0) {
      start_frame = end_frame = frame ;
      if (gzipped) { // pipe cannot seek
        int count;
        for (count=0; count < frame*width*height*depth*bpv; count++)
          fgetc(fp);
      } else
        fseek(fp, frame*width*height*depth*bpv, SEEK_CUR) ;
      nframes = 1 ;
    } else {
      /* hack - # of frames < -1 means to only read in that
                many frames. Otherwise I would have had to change the whole
                MRIread interface and that was too much of a pain. Sorry.
             */
      if (frame < -1)
        nframes = frame*-1 ;

      start_frame = 0 ;
      end_frame = nframes-1 ;
      /*if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "read %d frames\n", nframes);*/
    }
    buf = (BUFTYPE *)calloc(bytes, sizeof(BUFTYPE)) ;
    mri = MRIallocSequence(width, height, depth, type, nframes) ;
    mri->dof = dof ;
    for (frame = start_frame ; frame <= end_frame ; frame++) {
      for (z = 0 ; z < depth ; z++) {
        if (fread(buf, sizeof(char), bytes, fp) != bytes) {
          // fclose(fp) ;
          myclose(fp);
          free(buf) ;
          fprintf(stderr,"ERROR: mghRead(%s): could not read %d bytes at slice %d",
                  fname, bytes, z) ;
          return NULL;
        }
        switch (type) {
        case MRI_INT:
          for (i = y = 0 ; y < height ; y++) {
            for (x = 0 ; x < width ; x++, i++) {
              ival = orderIntBytes(((int *)buf)[i]) ;
              MRIIseq_vox(mri,x,y,z,frame-start_frame) = ival ;
            }
          }
          break ;
        case MRI_SHORT:
          for (i = y = 0 ; y < height ; y++) {
            for (x = 0 ; x < width ; x++, i++) {
              sval = orderShortBytes(((short *)buf)[i]) ;
              MRISseq_vox(mri,x,y,z,frame-start_frame) = sval ;
            }
          }
          break ;
        case MRI_TENSOR:
        case MRI_FLOAT:
          for (i = y = 0 ; y < height ; y++) {
            for (x = 0 ; x < width ; x++, i++) {
              fval = orderFloatBytes(((float *)buf)[i]) ;
              MRIFseq_vox(mri,x,y,z,frame-start_frame) = fval ;
            }
          }
          break ;
        case MRI_UCHAR:
          local_buffer_to_image(buf, mri, z, frame-start_frame) ;
          break ;
        default:
          errno = 0;
          fprintf(stderr,"ERROR: mghRead: unsupported type %d",
                  mri->type) ;
          return NULL;
          break ;
        }
      }
    }
    if (buf) free(buf) ;
  }

  if(good_ras_flag > 0) {
    mri->xsize =     xsize ;
    mri->ysize =     ysize ;
    mri->zsize =     zsize ;

    mri->x_r = x_r  ;
    mri->x_a = x_a  ;
    mri->x_s = x_s  ;

    mri->y_r = y_r  ;
    mri->y_a = y_a  ;
    mri->y_s = y_s  ;

    mri->z_r = z_r  ;
    mri->z_a = z_a  ;
    mri->z_s = z_s  ;

    mri->c_r = c_r  ;
    mri->c_a = c_a  ;
    mri->c_s = c_s  ;
    if (good_ras_flag > 0)
      mri->ras_good_flag = 1 ;
  } else {
    fprintf
    (stderr,
     "-----------------------------------------------------------------\n"
     "Could not find the direction cosine information.\n"
     "Will use the CORONAL orientation.\n"
     "If not suitable, please provide the information in %s.\n"
     "-----------------------------------------------------------------\n",
     fname
    );
    setDirectionCosine(mri, MRI_CORONAL);
  }
  // read TR, Flip, TE, TI, FOV
  if (freadFloatEx(&(mri->tr), fp)) {
    if (freadFloatEx(&fval, fp)) {
      mri->flip_angle = fval;
      // flip_angle is double. I cannot use the same trick.
      if (freadFloatEx(&(mri->te), fp))
        if (freadFloatEx(&(mri->ti), fp))
          if (freadFloatEx(&(mri->fov), fp))
            ;
    }
  }
  // tag reading
  {
    long long len ;
    char *fnamedir;

    while ((tag = TAGreadStart(fp, &len)) != 0) {
      switch (tag) {
      case TAG_OLD_MGH_XFORM:
      case TAG_MGH_XFORM:
        fgets(mri->transform_fname, len+1, fp);
        // if this file exists then read the transform
        if(!FileExists(mri->transform_fname)) {
          printf("  Talairach transform %s does not exist ...\n",
                 mri->transform_fname);
          fnamedir = fio_dirname(fname);
          sprintf(mri->transform_fname,"%s/transforms/talairach.xfm",fnamedir);
          printf("   ... trying %s ...",mri->transform_fname);
          if(FileExists(mri->transform_fname)) printf("which does exist ");
          else                                 printf("which does not exist ");
          printf("\n");
          free(fnamedir);
        }
        if(FileExists(mri->transform_fname)) {
          // copied from corRead()
          if(input_transform_file
              (mri->transform_fname, &(mri->transform)) == NO_ERROR) {
            mri->linear_transform = get_linear_transform_ptr(&mri->transform);
            mri->inverse_linear_transform =
              get_inverse_linear_transform_ptr(&mri->transform);
            mri->free_transform = 1;
            /*if (DIAG_VERBOSE_ON)
              fprintf
                (stderr,
            "INFO: loaded talairach xform : %s\n", mri->transform_fname);*/
          } else {
            errno = 0;
            fprintf(stderr,"ERROR: loading transform from %s",mri->transform_fname);
            exit(0);
            mri->linear_transform = NULL;
            mri->inverse_linear_transform = NULL;
            mri->free_transform = 1;
            (mri->transform_fname)[0] = '\0';
          }
        } else {
          fprintf(stderr,
                  "Can't find the talairach xform '%s'\n",
                  mri->transform_fname);
          fprintf(stderr, "Transform is not loaded into mri\n");
        }
        break ;
      case TAG_CMDLINE:
        if (mri->ncmds > MAX_CMDS) {
          fprintf(stderr,"ERROR: mghRead(%s): too many commands (%d) in file",
                  fname,mri->ncmds);
          exit(0);
        }
        mri->cmdlines[mri->ncmds] = calloc(len+1, sizeof(char)) ;
        fread(mri->cmdlines[mri->ncmds], sizeof(char), len, fp) ;
        mri->cmdlines[mri->ncmds][len] = 0 ;
        mri->ncmds++ ;
        break ;
      default:
        TAGskip(fp, tag, (long long)len) ;
        break ;
      }
    }
  }


  // fclose(fp) ;
  myclose(fp);

  // xstart, xend, ystart, yend, zstart, zend are not stored
  mri->xstart = - mri->width/2.*mri->xsize;
  mri->xend = mri->width/2. * mri->xsize;
  mri->ystart = - mri->height/2.*mri->ysize;
  mri->yend = mri->height/2.*mri->ysize;
  mri->zstart = - mri->depth/2.*mri->zsize;
  mri->zend = mri->depth/2.*mri->zsize;
  strcpy(mri->fname, fname);
  return(mri) ;
}

#undef UNUSED_SPACE_SIZE

