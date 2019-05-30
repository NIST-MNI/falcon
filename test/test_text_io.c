#include <stdio.h>
#include "falcon.h"


int main(int argc,char **argv)
{
    int i;
    niiktable *t;
    if(argc<3)
    {
        fprintf(stderr,"Usage:%s <in.csv> <out.csv>\n",argv[0]);
        return 1;
    }
    if((t=niiktable_read(argv[1]))==NULL) {
        fprintf(stderr,"Error reading %s\n",argv[1]);
        abort();
    }
    printf("Read %d entries with %d columns\n",t->col[0]->num,t->ncol);
    for(i=0;i<t->ncol;i++)
        printf("%s\t",t->name[i]);
        
    printf("\n");


    if(!niiktable_write(argv[2],t))
    {
        fprintf(stderr,"Error writing %s\n",argv[2]);
        abort();
    }

    /*ADD test to verify data*/
    niiktable_free(t);

    return 0;
}