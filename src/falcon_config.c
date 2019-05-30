#include "falcon_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

int verify_NIIKDIR(const char* d) {
  char tmp[4096];
  strcpy(tmp,d);
  /*TODO: verify something more relevant*/
  strcat(tmp,"/data/NVK/nifti1_kunio_viewer_file_bar_ui.xml");
  /*verify presence of files*/
  return !access(tmp, R_OK);
}

const char * get_NIIKDIR(void) {
  const char *env_niikdir;

  if( (env_niikdir = getenv("FALCON_DATA"))!=NULL && verify_NIIKDIR(env_niikdir) )
    return env_niikdir;

  if(verify_NIIKDIR(FALCON_DATA_DEFAULT))
    return FALCON_DATA_DEFAULT;

  if(verify_NIIKDIR(FALCON_DATA_SOURCE))
    return FALCON_DATA_SOURCE;

  /*print error?*/
  fprintf(stderr,"ERROR: please setenv or export FALCON_DATA\n");
  fprintf(stderr,"   for example, \n");
  fprintf(stderr,"   setenv FALCON_DATA /home/kunio/kproj/\n");
  fprintf(stderr,"   setenv FALCON_DATA /lab2/Kunio/kproj/\n");
  fprintf(stderr,"   export FALCON_DATA=/home/kunio/kproj/\n");
  fprintf(stderr,"   export FALCON_DATA=/lab2/Kunio/kproj/\n");

  return NULL;
}


static int __niik_verbose__=-1;

int niik_verbose( void ) {
  if(__niik_verbose__<0) {
    const char *env_niikverbose;
    if( (env_niikverbose = getenv("FALCON_VERBOSE"))!=NULL) {
      __niik_verbose__ = atoi(env_niikverbose);
    } else {
      __niik_verbose__ = FALCON_VERBOSE;
    }
  }
  return __niik_verbose__;
}

void niik_set_verbose( int v ) {
  __niik_verbose__=v;
}


const char * get_DEBUG_PREFIX(void) {
  return getenv("FALCON_DEBUG_PREFIX");
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/