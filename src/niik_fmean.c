/* Filename:    niik_fmean.c
 * Description: floating point mean
 * Author:      Kunio Nakamura
 * Date:        March 10, 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

int main(argc, argv)
int argc;
char *argv[];
{
  int dig=10,n;
  double d=0;

  if( argc<=2 ) {
    fprintf(stdout,"  usage: num1 num2 [num3...]\n");
    exit(0);
  }

  if(!strncmp(argv[1],"-dig",4)) {
    dig=atoi(argv[2]);
    argc-=2;
    argv+=2;
  }

  for(n=1; n<argc; n++) {
    d+=atof(argv[n]);
  }
  d/=(argc-1);

  switch(dig) {
  case  9:
    printf("%-20.9f \n",d);
    break;
  case  8:
    printf("%-20.8f \n",d);
    break;
  case  7:
    printf("%-20.7f \n",d);
    break;
  case  6:
    printf("%-20.6f \n",d);
    break;
  case  5:
    printf("%-20.5f \n",d);
    break;
  case  4:
    printf("%-20.4f \n",d);
    break;
  case  3:
    printf("%-20.3f \n",d);
    break;
  case  2:
    printf("%-20.2f \n",d);
    break;
  case  1:
    printf("%-20.1f \n",d);
    break;
  case  0:
    printf("%-20.0f \n",d);
    break;
  default:
  case 10:
    printf("%-20.10f \n",d);
    break;
  }
  exit(0);
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/