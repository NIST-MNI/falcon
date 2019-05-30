
#pragma once

#ifndef __FALCON_TEXT_IO_H__
#define __FALCON_TEXT_IO_H__

/******************************************************************
 *
 * falcon_text_io.c
 *
 ******************************************************************/
niikvec *niikvec_read(const char *fname);
int niikvec_write(const char *fname,niikvec *v);
int niikvec_write_ex(const char *fname,niikvec *v,const char *colnames);

niikmat *niikmat_read(const char *fname);
int niikmat_write_ex(const char *fname,niikmat *mat,const char **colnames);
int niikmat_write(const char *fname,niikmat *mat);
int niikmat_write_binary(const char *fname,niikmat *mat);


niiktable *niiktable_init(int ncol,int nrow);
niiktable *niiktable_init_ex(int ncol,int nrow,const niiktable *ref);
niiktable *niiktable_free(niiktable *v);
niiktable *niiktable_copy(const niiktable *src);
int niiktable_set_all(niiktable *v,double val);

niiktable *niiktable_read(const char *fname);
int niiktable_write(const char *fname,const niiktable *v);

#endif /*__printf(fp,gp,"\n");*/