/* Filename:    fcalc.c
 * Description: floating point calculation
 * Author:      Kunio Nakamura
 * Date:        March 9, 2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

typedef struct _fc {
  struct _fc *lef;
  struct _fc *rig;
  double f;
  char op;
} FC;

FC *parsetree(char **argv,int argc,int *nc);
FC *addnodefc(char op,double f);
void displaytree(FC *curr);
double calctree(FC *curr);
void displayerror(char error_text[]);

int main(argc, argv)
int argc;
char *argv[];
{
  FC *top;
  int nc=1,dig=10;

  if( argc<=3 ) {
    printf("  simple command-line calculator\n");
    printf("    e.g.  fcalc o 3 + -2 c / 5 x 2 + 1  \n");
    printf("          results in  1.400\n");
    printf("    o/c = open/closed paranthesis\n");
    exit(0);
  }

  if(!strncmp(argv[1],"-dig",4)) {
    dig=atoi(argv[2]);
    argc-=2;
    argv+=2;
  }

  top=parsetree(argv,argc,&nc);

  switch(dig) {
  case  9:
    printf("%-20.9f \n",calctree(top));
    break;
  case  8:
    printf("%-20.8f \n",calctree(top));
    break;
  case  7:
    printf("%-20.7f \n",calctree(top));
    break;
  case  6:
    printf("%-20.6f \n",calctree(top));
    break;
  case  5:
    printf("%-20.5f \n",calctree(top));
    break;
  case  4:
    printf("%-20.4f \n",calctree(top));
    break;
  case  3:
    printf("%-20.3f \n",calctree(top));
    break;
  case  2:
    printf("%-20.2f \n",calctree(top));
    break;
  case  1:
    printf("%-20.1f \n",calctree(top));
    break;
  case  0:
    printf("%-20.0f \n",calctree(top));
    break;
  default:
  case 10:
    printf("%-20.10f \n",calctree(top));
    break;
  }
  exit(0);
}



FC *parsetree(char **argv,int argc,int *nc)
/* parse the arguments and construct a linked list */
{
  FC *top,*curr,*next;
  double f;
  int nn,num,up;

  nn=*nc;

  curr=top=NULL;

  num=1; /* number expected */
  up=0;  /* current node is not a umber but treat as a number */

  while(nn<argc) {

    switch(argv[nn][0]) {
    case '[':
    case 'o':
      nn++;
      next=parsetree(argv,argc,&nn);
      if(curr==NULL) {
        up=1;
        curr=top=next;
      } else {
        curr->rig=next;
      }
      num=0;
      break;

    case ']':
    case 'c':
      nn++;
      *nc=nn;
      return top;

    case 'x':
    case '/':
      if(top==NULL) displayerror("num or [ expected");
      if(num==1) displayerror("num or [ expected");

      next=addnodefc(argv[nn++][0],0);
      if(curr->op=='v' ||  /* current node is a number */
          up==1 ) {         /* current node is not a number but treat as a number (after []) */
        next->lef=curr;
        if(top==curr) top=next;
        curr=next;
      }

      else if(curr->op=='x' || curr->op=='/') {  /* current node is multiplication or division */
        next->lef=curr;
        if(top==curr) top=next;
        else top->rig=next;
        curr=next;
      } else if(curr->op=='+' || curr->op=='-') { /* current node is addition or subtraction */
        next->lef=curr->rig;
        curr->rig=next;
        curr=next;
      }
      num=1;
      up=0;
      break;

    case '+':
    case '-':
      if( num==1 ) { /* number starting from + or - */
        f=atof(argv[nn++]);
        if(curr==NULL) {
          curr=top=addnodefc('v',f);
        } else curr->rig=addnodefc('v',f);
        num=0;
      } else { /* operation +/- */
        next=addnodefc(argv[nn++][0],0);
        if(curr->op!='v') {
          next->lef=top;
          top=curr=next;
        } else {
          next->lef=curr;
          if(top==curr) top=next;
          curr=next;
        }
        num=1;
      }
      up=0;
      break;

    default: /* numbers */
      if(num==0) displayerror("op expected");
      f=atof(argv[nn++]);
      next=addnodefc('v',f);
      if(top==NULL)
        curr=top=next; /* beginning*/
      else {
        curr->rig=next;
      }
      num=0;
      break;
    }

    /*printf("  %i: ",nn); displaytree(top);  printf("\n");*/
  }

  *nc=nn;
  return top;
}

double calctree(FC *curr) {
  if(curr==NULL) displayerror("insufficient argument");
  switch(curr->op) {
  case 'v':
    return curr->f;
  case 'x':
    return calctree(curr->lef)*calctree(curr->rig);
  case '/':
    return calctree(curr->lef)/calctree(curr->rig);
  case '+':
    return calctree(curr->lef)+calctree(curr->rig);
  case '-':
    return calctree(curr->lef)-calctree(curr->rig);
  default:
    displayerror("unknown op");
  }
  return 0;
}

FC *addnodefc(char op,double f) {
  FC *curr;
  curr=(FC *)malloc(sizeof(FC));
  curr->lef=curr->rig=NULL;
  curr->f=f;
  curr->op=op;
  return curr;
}

void displaytree(FC *curr) {
  if(curr==NULL) return;
  switch(curr->op) {
  case 'v':
    printf(" %f ",curr->f);
    break;
  default:
    printf(" %c (",curr->op);
    displaytree(curr->lef);
    displaytree(curr->rig);
    printf(" ) ");
  }
}

void displayerror(char error_text[]) {
  fprintf(stderr,"ERROR: run-time error, %s\n",error_text);
  exit(1);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/