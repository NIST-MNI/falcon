#ifndef __NIIK_CONFIG_H__
#define __NIIK_CONFIG_H__ 1

/*HACK*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#cmakedefine HAVE_MINC1 1 
#cmakedefine HAVE_MINC2 1
#cmakedefine HAVE_BICPL 1

#ifndef MINC2
#define MINC2 @MINC2@
#endif

#define FALCON_DATA_SOURCE  "@FALCON_DATA_SOURCE@"
#define FALCON_DATA_DEFAULT "@FALCON_DATA_DEFAULT@"
#define FALCON_VERBOSE       @FALCON_VERBOSE@

#define NIIK_MAJOR_VERSION   @NIIK_MAJOR_VERSION@
#define NIIK_MINOR_VERSION   @NIIK_MINOR_VERSION@
#define NIIK_MICRO_VERSION   @NIIK_MICRO_VERSION@
#define NIIK_VERSION        "@NIIK_MAJOR_VERSION@.@NIIK_MINOR_VERSION@.@NIIK_MICRO_VERSION@"

/*configuration accees functions*/
const char * get_NIIKDIR(void);
int niik_verbose(void);
void niik_set_verbose(int);

/*get directory for dumping debug info*/
const char * get_DEBUG_PREFIX(void);

/*file prefix for post-mortem files*/
const char * get_POSTMORTEM_PREFIX(void);


#endif /*__NIIK_CONFIG_H__*/

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
*/
