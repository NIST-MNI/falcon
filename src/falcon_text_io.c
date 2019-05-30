/* Filename:     falcon_text_io.c
 * Description:  text io
 * Author:       Vladmir Fonov
 * Date:         January 17, 2019
 */

#include <stdlib.h>
#include "falcon.h"
#include "csv.h"

int niikvec_write(const char *fname,niikvec *v)
{
    return niikvec_write_ex(fname,v,NULL);
}

int niikvec_write_ex(const char *fname, niikvec *v,const char *colname) {
  if(!niik_write_double_vector_ex(fname,v->v,v->num,colname)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_write_double_vector\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  return 1;
}

int niik_write_double_vector(const char *fname,double *v,int num) {
    return niik_write_double_vector_ex(fname,v,num,NULL);
}

int niik_write_double_vector_ex(const char *fname,double *v,int num, const char *colname) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  int n;
  int csv_format=0;

  if(fname==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fname is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  csv_format=(!strncmp(fname+(strlen(fname)-4),".csv",4) || !strncmp(fname+(strlen(fname)-7),".csv.gz",7));

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    if((gp = gzopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: gzopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }
    if(csv_format)
      if(colname) gzprintf(gp,"%s\n",colname);
      else gzprintf(gp,"V1\n");

    for(n=0; n<num; n++) {
      gzprintf(gp,"%15.9f\n",v[n]);
    }
    gzclose(gp);
  } else {
    if((fp=fopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }

    if(csv_format)
      if(colname) fprintf(fp,"%s\n",colname);
      else fprintf(fp,"V1\n");

    for(n=0; n<num; n++) {
      fprintf(fp,"%15.9f\n",v[n]);
    }
    fclose(fp);
  }
  return 1;
}

int niik_write_double_vectors(const char *fname,double **v,int num,int ncol)
{
    return niik_write_double_vectors_ex(fname,v,num,ncol,NULL);
}

int niik_write_double_vectors_ex(const char *fname,double **v,int num,int ncol, const char **colnames) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  int n,i;
  int csv_format;
  if(fname==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fname is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  csv_format=(!strncmp(fname+(strlen(fname)-4),".csv",4) || !strncmp(fname+(strlen(fname)-7),".csv.gz",7));


  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    if((gp = gzopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: gzopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }

    if(csv_format) 
    {
      if(colnames)
      {
        for(i=0; i<ncol; i++) {
          gzprintf(gp,"%s",colnames[i]);
          if(i<(ncol-1))
            gzprintf(gp,",");
          else
            gzprintf(gp,"\n");
        }
      } else {
        for(i=0; i<ncol; i++) {
          gzprintf(gp,"V%d",i+1);
          if(i<(ncol-1))
            gzprintf(gp,",");
          else
            gzprintf(gp,"\n");
        }
      }
    }

    for(n=0; n<num; n++) {
      for(i=0; i<ncol; i++) {
        gzprintf(gp,"%15.9f",v[i][n]);
        if(i<(ncol-1))
          gzprintf(gp,",");
        else
          gzprintf(gp,"\n");
      }
    }
    gzclose(gp);
  } else {  
    if((fp=fopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }


    if(csv_format) 
    {
      if(colnames)
      {
        for(i=0; i<ncol; i++) {
          fprintf(fp,"%s",colnames[i]);
          if(i<(ncol-1))
            fprintf(fp,",");
          else
            fprintf(fp,"\n");
        }
      } else {
        for(i=0; i<ncol; i++) {
          fprintf(fp,"V%d",i+1);
          if(i<(ncol-1))
            fprintf(fp,",");
          else
            fprintf(fp,"\n");
        }
      }
    }

    for(n=0; n<num; n++) {
      for(i=0; i<ncol; i++) {
        fprintf(fp,"%15.9f",v[i][n]);
        if(i<(ncol-1))
          fprintf(fp,",");
        else
          fprintf(fp,"\n");
      }
    }
    fclose(fp);
  }
  return 1;
}


int niik_write_double_vector_binary(const char *fname,double *v,int num) {
  FILE *fp;
  if(fname==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fname is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if((fp=fopen(fname,"w"))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
    return 0;
  }
  fwrite(v,num,sizeof(double),fp);
  fclose(fp);
  return 1;
}



int niik_write_int_vector(const char *fname,int *v,int num) {
    return niik_write_int_vector_ex(fname,v,num,NULL);
}

int niik_write_int_vector_ex(const char *fname,int *v,int num,const char *colname) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  int n;
  int csv_format;

  if(fname==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fname is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  csv_format=(!strncmp(fname+(strlen(fname)-4),".csv",4) || !strncmp(fname+(strlen(fname)-7),".csv.gz",7));

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    if((gp = gzopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: gzopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }

    if(csv_format)
      if(colname) gzprintf(gp,"%s\n",colname);
      else gzprintf(gp,"%s\n","V1");

    for(n=0; n<num; n++) {
      gzprintf(gp,"%i\n",v[n]);
    }
    gzclose(gp);
  } else {
    if((fp=fopen(fname,"w"))==NULL) {
      fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
      return 0;
    }

    if(csv_format)
      if(colname) fprintf(fp,"%s\n",colname);
      else fprintf(fp,"%s\n","V1");

    for(n=0; n<num; n++) {
      fprintf(fp,"%i\n",v[n]);
    }
    fclose(fp);
  }
  return 1;
}

int niik_display_int_vector(int *v,int num) {
  int n;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return 0;
  }
  for(n=0; n<num; n++) {
    fprintf(stdout,"%6i ",v[n]);
  }
  fprintf(stdout,"\n");
  return 1;
}

int niik_display_double_vector(double *v,int num) {
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) {
    fprintf(stdout,"%6.5lf ",v[n]);
  }
  fprintf(stdout,"\n");
  return 1;
}

int niik_display_float_vector(float *f,int num) {
  int n;
  if(f==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: f is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) {
    fprintf(stdout,"%6.5lf ",f[n]);
  }
  fprintf(stdout,"\n");
  return 1;
}


int niikmat_write(const char *fname, niikmat *mat) 
{   
    return niikmat_write_ex(fname,mat,NULL);
}

int niikmat_write_ex(const char *fname, niikmat *mat,const char **colnames) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  int i,j;
  NIIK_RET0((fname==NULL), __func__, "fname is null");
  NIIK_RET0((  mat==NULL), __func__, "mat is null");
  if(!strncmp(fname+(strlen(fname)-4),".xfm",4)) {
    NIIK_RET0((!niikmat_write_as_linear_xfm((char *)fname,mat)),__func__,"niikmat_write_as_linear_xfm");
    return 1;
  }
  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) 
  {
    if((gp=gzopen(fname,"w"))==NULL) {
        fprintf(stderr,"[%s] ERROR: gzopen %s\n",__func__,fname);
        return 0;
    }
    if(colnames)
    {
        for(j=0; j<mat->col; j++) {
            gzprintf(gp,"%s ",colnames[j]);
        }
        gzprintf(gp,"\n");
    }
    for(i=0; i<mat->row; i++) {
        for(j=0; j<mat->col; j++) 
            gzprintf(gp,"%15.9f ",mat->m[i][j]);
        
        gzprintf(gp,"\n");
    }
    gzclose(gp);
  } else {
    if((fp=fopen(fname,"w"))==NULL) {
        fprintf(stderr,"[%s] ERROR: fopen %s\n",__func__,fname);
        return 0;
    }

    if(colnames)
    {
        for(j=0; j<mat->col; j++) 
            fprintf(fp,"%s ",colnames[j]);
        
        fprintf(fp,"\n");
    }

    for(i=0; i<mat->row; i++) {
        for(j=0; j<mat->col; j++) {
            fprintf(fp,"%15.9f ",mat->m[i][j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
  }
  return 1;
}

int niikmat_write_binary(const char *fname,niikmat *mat) {
  FILE *fp;
  int i;
  NIIK_RET0((fname==NULL),__func__,"fname is null");
  NIIK_RET0((  mat==NULL),__func__,"mat is null");
  if((fp=fopen(fname,"w"))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    fwrite(mat->m[i],mat->col,sizeof(double),fp);
  }
  fclose(fp);
  return 1;
}

niikmat *niikmat_read(const char *fname) {
  FILE *fp=NULL;
  gzFile gp=NULL;
  int m,n,nn,slen;
  niikmat *mat=NULL;
  double d;
  char *cptr,str[65536];
  int verbose=0;

  if(fname==NULL) {
    fprintf(stderr,"ERROR: fname is null\n");
    return NULL;
  }

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    if((gp = gzopen(fname,"r"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return NULL;
    }
  } else {
    if((fp = fopen(fname,"r"))==NULL) {
        fprintf(stderr,"ERROR: fopen %s\n",fname);
        return 0;
    }
  }
  m=n=nn=0;

  while(  (fp && fgets(str,65536,fp)!=NULL) ||
          (gp && gzgets(gp,str,65536)!=NULL)
   ) {
    if(verbose) fprintf(stdout,"row %i |  %s ",m,str);
    n=0;
    cptr=str;
    slen=strlen(str);
    for(;;) {
      while(cptr[0]==' ') {
        cptr++;
      }
      d=atof(cptr);
      n++;
      if(verbose) fprintf(stdout,"%2i %f | %s",n,d,cptr);
      if((cptr=strchr(cptr,' '))==NULL) {
        if(verbose) fprintf(stdout,"break -- no space\n");
        break;
      }
      while(cptr[0]==' ') {
        if(cptr==str+slen) break;
        cptr++;
      }
      if(cptr==str+slen) {
        if(verbose) fprintf(stdout,"break -- max string length\n");
        break;
      }
      if(cptr[0]=='\n') {
        if(verbose) fprintf(stdout,"break -- new line byte\n");
        break;
      }
    }
    if(!nn) nn=n;
    else if(nn!=n) {
      fprintf(stderr,"ERROR: #col (%i) is different in row = %i\n",n,m);
      fprintf(stderr,"       previously #col was %i\n",nn);
      return NULL;
    }
    m++;
  }
  if(verbose) fprintf(stdout,"%i %i\n",m,n);
  mat=niikmat_init(m,n);

  if(fp) rewind(fp);
  else gzrewind(gp);

  m=0;
  while(  (fp && fgets(str,65536,fp)!=NULL) ||
          (gp && gzgets(gp,str,65536)!=NULL)) 
    {
    n=0;
    cptr=str;
    for(;;) {
      while(cptr[0]==' ') {
        cptr++;
      }
      d=atof(cptr);
      mat->m[m][n++]=d;
      if((cptr=strchr(cptr,' '))==NULL) break;
      while(cptr[0]==' ') {
        if(cptr[1]==0) break;
        cptr++;
      }
      if(cptr[0]=='\n') break;
      if(cptr[1]==0) break;
    }
    m++;
  }
  if(fp)
    fclose(fp);
  else
    gzclose(gp);
  return mat;
}

/*reading .csv files with libcsv*/
int *niik_read_int_vector_ex(const char *fname, int *num, int column_id, const char *column_name) {
  FILE *fp;
  gzFile gp=NULL;
  int *v;
  int d;
  int _num;
  char str[4096];
  int skip_header=!strncmp(fname+(strlen(fname)-7),".csv.gz",7)||
                  !strncmp(fname+(strlen(fname)-4),".csv",4);

  if(fname==NULL) {
    fprintf(stderr,"ERROR: fname is null\n");
    return NULL;
  }
  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    /*not the most effective way to work with .gz file*/
    if((gp = gzopen(fname,"r"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return NULL;
    }
    _num=0;
    while(!gzeof(gp)) {
      if(gzgets(gp,str,4096)==NULL) {
        break;
      }
      _num ++;
    }
    gzrewind(gp);
    if(!_num || _num==1&&skip_header) /*EMPTY FILE*/
    {
        gzclose(gp);
        return NULL;
    }
    if(skip_header) {
        _num--;
        gzgets(gp,str,4096);
    }

    v=(int *)calloc(_num, sizeof(int));
    *num=0;
    while(!gzeof(gp) && (*num)<_num) {
      if( (gzgets(gp,str,4096)==NULL ) ||
          (sscanf(str,"%i",&v[*num])<1 )) {
        fprintf(stderr,"ERROR: scanf\n");
        free(v);
        gzclose(gp);
        return NULL;
      }
      (*num) ++;
    }
    gzclose(gp);
  } else {
    if((fp=fopen(fname,"r"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",fname);
      return NULL;
    }
    _num=0;
    while(!feof(fp)) {
      if((fscanf(fp,"%i",&d))<1) {
        break;
      }
      _num ++;
    }
    rewind(fp);
    if(!_num || _num==1&&skip_header) /*EMPTY FILE*/
    {
        fclose(fp);
        return NULL;
    }
    if(skip_header) {
      _num--;
      if(!fgets(str,4096,fp))
      {
        fprintf(stderr,"ERROR: fscanf\n");
        free(v);
        return NULL;
      }
    }

    v=(int *)calloc(_num,sizeof(int));
    *num=0;
    while(!feof(fp) && (*num)<_num) {
      if((fscanf(fp,"%i",&v[*num]))<1) {
        fprintf(stderr,"ERROR: fscanf\n");
        free(v);
        return NULL;
      }
      (*num) ++;
    }
    fclose(fp);
  }
  return v;
}

struct falcon_csv_reader{
  int select_column;
  const char *select_column_name;
  void *buffer;
  int is_integer;
  int count;
  int alloc;
  int skip_header;
  int fld_cnt;
  int num;
};

void cb_field_identify(void *s, size_t len, void *data) {
  struct falcon_csv_reader *rdr=(struct falcon_csv_reader*)data;
  /*first row only*/
  if(!rdr->count )
  {
    if(rdr->skip_header && rdr->select_column<0 && rdr->select_column_name) /*compare column name*/
    {
      if(!strncmp((const char*)s,rdr->select_column_name,len))
        rdr->select_column=rdr->fld_cnt; /*found column*/
    }
    rdr->fld_cnt++;
  }
}

void cb_row_identify (int c, void *data) 
{
  /*called on each row*/
  struct falcon_csv_reader *rdr=(struct falcon_csv_reader*)data;
  rdr->count++;
}

void cb_field_read (void *s, size_t len, void *data) {
  struct falcon_csv_reader *rdr=(struct falcon_csv_reader*)data;
  if(rdr->skip_header && !rdr->count) return; /*skip header*/

  if(rdr->fld_cnt==rdr->select_column) /*parse this column*/
  {
    if(rdr->is_integer) {
      int *buffer=(int*)rdr->buffer;
      buffer[rdr->num]=atoi((const char*)s);
    } else {
      double *buffer=(double*)rdr->buffer;
      buffer[rdr->num]=atof((const char*)s);
    }
    rdr->num++;
  }
  rdr->fld_cnt++;
}

void cb_row_read (int c, void *data) 
{
  /*called on each row*/
  struct falcon_csv_reader *rdr=(struct falcon_csv_reader*)data;
  rdr->count++;
  rdr->fld_cnt=0;
}


void *niik_read_vector_ex(const char *fname,int is_int, int *num, int column_id, const char *column_name) {
  FILE *fp=NULL;
  gzFile gp=NULL;

  char str[4096];
  int bytes_read;
  struct csv_parser parser;
  struct falcon_csv_reader reader;

  if(fname==NULL) {
    fprintf(stderr,"ERROR: fname is null\n");
    return NULL;
  }

  memset(&reader,0,sizeof(struct falcon_csv_reader));
  reader.is_integer=is_int;
  reader.select_column=column_id;
  reader.select_column_name=column_name;

  reader.skip_header=( !strncmp(fname+(strlen(fname)-7),".csv.gz",7)||
                       !strncmp(fname+(strlen(fname)-4),".csv",4) );

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    /*not the most effective way to work with .gz file*/
    if((gp = gzopen(fname,"rb"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return NULL;
    }
  } else {
    if((fp=fopen(fname,"rb"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",fname);
      return NULL;
    }
  }

  if (csv_init(&parser, CSV_APPEND_NULL)) {
      fprintf(stderr,"ERROR: csv_init \n");
      return NULL;
  }

  /*1st pass identify column and count rows*/
  while ( (fp && (bytes_read=fread(str, 1, 4096, fp)) > 0) ||
          (gp && (bytes_read=gzfread(str, 1, 4096, gp)) > 0) )
  {

    if (csv_parse(&parser, str, bytes_read, cb_field_identify, cb_row_identify, &reader) != bytes_read) {
      fprintf(stderr, "Error while parsing file: %s\n",
          csv_strerror(csv_error(&parser)) );
      csv_free(&parser);
      if(fp)
        fclose(fp);
      if(gp)
        gzclose(gp);
      return NULL;
    }
  }
  csv_fini(&parser, cb_field_identify, cb_row_identify, &reader);

  if(column_name && reader.select_column<0) /*check if the field is found*/
  {
    csv_free(&parser);
    if(fp)
      fclose(fp);
    if(gp)
      gzclose(gp);

    fprintf(stderr, "Can't find field %s in file %s\n",column_name,fname);
    return NULL;
  }

  *num=reader.count-(reader.skip_header?1:0);
  /*2nd pass*/
  if(reader.is_integer)
    reader.buffer=calloc(*num,sizeof(int));
  else
    reader.buffer=calloc(*num,sizeof(double));
  
  if(fp) rewind(fp);
  else   gzrewind(gp);

  reader.count=0;
  reader.fld_cnt=0;

  while ( (fp && (bytes_read=fread(str, 1, 4096, fp)) > 0) ||
          (gp && (bytes_read=gzfread(str, 1, 4096, gp)) > 0) )
  {
    if (csv_parse(&parser, str, bytes_read, cb_field_read, cb_row_read, &reader) != bytes_read) {
      fprintf(stderr, "Error while parsing file: %s\n",
          csv_strerror(csv_error(&parser)) );
      free(reader.buffer);
      csv_free(&parser);
      if(fp)
        fclose(fp);
      if(gp)
        gzclose(gp);
      return NULL;
    }
  }
  csv_fini(&parser, cb_field_read, cb_row_read, &reader);

  csv_free(&parser);
  if(fp)
    fclose(fp);
  if(gp)
    gzclose(gp);

  if(reader.num!=*num) /*DEBUG*/
    fprintf(stderr,"BUG:Inconsistent number on entries\n");

  return reader.buffer;
}

double *niik_read_double_vector(const char *fname,int *num)
{
  /**/
  return (double*)niik_read_vector_ex(fname,0,num,0,NULL);
}

int *niik_read_int_vector(const char *fname,int *num)
{
  /**/
  return (int*)niik_read_vector_ex(fname,1,num,0,NULL);
}

/*table operations*/
niiktable *niiktable_init(int ncol, int nrow)
{
  int i;
  niiktable *t;

  if((t=(niiktable *)calloc(1,sizeof(niiktable)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
    return NULL;
  }

  if((t->name=(char**)calloc(ncol,sizeof(char*)))==NULL) {
    free(t);
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
    return NULL;
  }

  if((t->col=(niikvec**)calloc(ncol,sizeof(niikvec*)))==NULL) {
    free(t->name);
    free(t);

    fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  t->ncol=ncol;

  for(i=0;i<ncol;i++)
  {
    if((t->col[i]=niikvec_init(nrow))==NULL) {
      free(t);/*possible memory leak ?*/
      fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
      return NULL;
    }
  }
  return t;
}

niiktable *niiktable_free(niiktable *t)
{
  int i;
  if(!t) return NULL;
  for(i=0;i<t->ncol;i++)
  {
    if(t->col[i])  niikvec_free(t->col[i]);
    if(t->name[i]) free(t->name[i]);
  }
  free(t->name);
  free(t->col);
  free(t);
  return NULL;
}

niiktable *niiktable_copy(const niiktable *src)
{
  int i;
  niiktable * c;
  if(!src) return NULL;
  
  if((c=niiktable_init(src->ncol,src->col[0]->num))==NULL)
  {
      fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
      return NULL;
  }
  
  for(i=0;i<src->ncol;i++)
  {
    memcpy(c->col[i]->v,src->col[i]->v,sizeof(double)*src->col[0]->num);
    if(src->name[i])
      c->name[i]=strdup(src->name[i]);
  }
  return c;
}

int niiktable_set_all(niiktable *t,double val)
{
  int i;
  if(!val) return 1;

  for(i=0;i<t->ncol;i++)
  {
    if(t->col[i]) niikvec_set_all(t->col[i],val);
  }
  return 0;
}


struct falcon_csv_reader2 {
  int is_integer;
  int count;
  int alloc;
  int skip_header;
  int fld_cnt;
  char **columns;
  niiktable * t;
};


void cb_field_identify2(void *s, size_t len, void *data) {
  struct falcon_csv_reader2 *rdr=(struct falcon_csv_reader2*)data;
  /*first row only*/
  if(!rdr->count )
  {
    if(rdr->skip_header ) /*store column name*/
    {
      rdr->columns=(char **)realloc(rdr->columns,sizeof(char*)*(rdr->fld_cnt+1));
      rdr->columns[rdr->fld_cnt]=strdup((const char*)s);
    }
    rdr->fld_cnt++;
  }
}

void cb_row_identify2 (int c, void *data) 
{
  /*called on each row*/
  struct falcon_csv_reader2 *rdr=(struct falcon_csv_reader2*)data;
  rdr->count++;
}

void cb_field_read2 (void *s, size_t len, void *data) {
  struct falcon_csv_reader2 *rdr=(struct falcon_csv_reader2*)data;

  if(rdr->skip_header && !rdr->count) return; /*skip header*/

  rdr->t->col[rdr->fld_cnt]->v[ rdr->count-(rdr->skip_header?1:0) ] = atof((const char*)s);
  rdr->fld_cnt++;
}

void cb_row_read2 (int c, void *data) 
{
  /*called on each row*/
  struct falcon_csv_reader2 *rdr=(struct falcon_csv_reader2*)data;
  rdr->count++;
  rdr->fld_cnt=0;
}


niiktable *niiktable_read(const char *fname)
{
  FILE *fp=NULL;
  gzFile gp=NULL;

  char str[4096];
  int bytes_read;
  struct csv_parser parser;
  struct falcon_csv_reader2 reader;
  niiktable * t;

  if(fname==NULL) {
    fprintf(stderr,"ERROR: fname is null\n");
    return NULL;
  }

  memset(&reader,0,sizeof(struct falcon_csv_reader2));

  reader.skip_header=( !strncmp(fname+(strlen(fname)-7),".csv.gz",7)||
                       !strncmp(fname+(strlen(fname)-4),".csv",4) );

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) {
    /*not the most effective way to work with .gz file*/
    if((gp = gzopen(fname,"rb"))==NULL) {
      fprintf(stderr,"ERROR: gzopen %s\n",fname);
      return NULL;
    }
  } else {
    if((fp=fopen(fname,"rb"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",fname);
      return NULL;
    }
  }

  if (csv_init(&parser, CSV_APPEND_NULL)) {
      fprintf(stderr,"ERROR: csv_init \n");
      return NULL;
  }

  /*1st pass identify columns and count rows*/
  while ( (fp && (bytes_read=fread(str, 1, 4096, fp)) > 0) ||
          (gp && (bytes_read=gzfread(str, 1, 4096, gp)) > 0) )
  {

    if (csv_parse(&parser, str, bytes_read, cb_field_identify2, cb_row_identify2, &reader) != bytes_read) {
      fprintf(stderr, "Error while parsing file: %s\n",
          csv_strerror(csv_error(&parser)) );
      csv_free(&parser);
      if(fp)
        fclose(fp);
      if(gp)
        gzclose(gp);
      return NULL;
    }
  }
  csv_fini(&parser, cb_field_identify2, cb_row_identify2, &reader);
  t = niiktable_init(reader.fld_cnt,reader.count-(reader.skip_header?1:0));
  reader.t=t;

  /*2nd pass*/
  /*transfer ownership */
  if(reader.skip_header) {
    free(t->name);
    t->name = reader.columns;
  }

  if(fp) rewind(fp);
  else   gzrewind(gp);

  reader.count=0;
  reader.fld_cnt=0;

  while ( (fp && (bytes_read=fread(str, 1, 4096, fp)) > 0) ||
          (gp && (bytes_read=gzfread(str, 1, 4096, gp)) > 0) )
  {
    if (csv_parse(&parser, str, bytes_read, cb_field_read2, cb_row_read2, &reader) != bytes_read) {
      fprintf(stderr, "Error while parsing file: %s\n",
          csv_strerror(csv_error(&parser)) );
      csv_free(&parser);
      if(fp)
        fclose(fp);
      if(gp)
        gzclose(gp);
      return NULL;
    }
  }
  csv_fini(&parser, cb_field_read2, cb_row_read2, &reader);

  csv_free(&parser);
  if(fp)
    fclose(fp);
  if(gp)
    gzclose(gp);

  if( (reader.count-(reader.skip_header?1:0))!=t->col[0]->num) /*DEBUG*/
    fprintf(stderr,"BUG:Inconsistent number on entries\n");

  return t;
}

#define __printf(f,g,...)  if(f) fprintf(f,__VA_ARGS__); else gzprintf(g,__VA_ARGS__) 

int niiktable_write(const char *fname,const niiktable *t)
{
  /*TODO: finish*/
  FILE *fp=NULL;
  gzFile gp=NULL;
  int i,j;
  int csv_format=0;
  char sep=' ';
  NIIK_RET0((fname==NULL), __func__, "fname is null");
  NIIK_RET0((    t==NULL), __func__, "table is null");

  csv_format=( !strncmp(fname+(strlen(fname)-7),".csv.gz",7)||
               !strncmp(fname+(strlen(fname)-4),".csv",4) );
  
  if(csv_format) sep=',';

  if(!strncmp(fname+(strlen(fname)-3),".gz",3)) 
  {
    if((gp=gzopen(fname,"w"))==NULL) {
        fprintf(stderr,"[%s] ERROR: gzopen %s\n",__func__,fname);
        return 0;
    }
  } else {
    if((fp=fopen(fname,"w"))==NULL) {
        fprintf(stderr,"[%s] ERROR: fopen %s\n",__func__,fname);
        return 0;
    }
  }

  if(csv_format)
  {
      for(j=0; j<t->ncol; j++) {
        if(t->name && t->name[j])
        {
          __printf(fp,gp,"%s",t->name[j]);
        } else {
          __printf(fp,gp,"V%d",j+1);

        }
        if(j<(t->ncol-1)) __printf(fp,gp,"%c",sep);
      }
      __printf(fp,gp,"\n");
  }

  for(i=0; i<t->col[0]->num; i++) {
        for(j=0; j<t->ncol; j++)
        {
            __printf(fp,gp,"%.15g",t->col[j]->v[i]);
            if(j<(t->ncol-1)) __printf(fp,gp,"%c",sep);
        }
        __printf(fp,gp,"\n");
  }
  if(gp) gzclose(gp);
  if(fp) fclose(fp);
  return 1;
}


niiktable *niiktable_init_ex(int ncol,int nrow,const niiktable *ref)
{
  niiktable *t;
  int i;
  if(ncol!=ref->ncol)
  {
    fprintf(stderr,"Inconsistent number of columns\n");
    return NULL;
  }

  t = niiktable_init(ncol,nrow);
  if(!t) return NULL;
  /*t->name=(char**)realloc(t->name,sizeof(char*)*ncol);*/
  for(i=0;i<ncol;i++)
    t->name[i]=strdup(ref->name[i]);
  return t;
}