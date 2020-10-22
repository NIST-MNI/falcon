/* FILENAME:     off_split.c
 * DESCRIPTION:  Split off file into disconnected components
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include  <volume_io.h>

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input.off> <output_prefix> \n"
          "\t--clobber clobber output files\n",
          name);
}

static int compare_ints(const void *a, const void *b) {
  return (*(const int *)a - *(const int *)b);
}

struct group_cm {
  int id;
  niikpt cm;
  int c;
};

static int compare_groups(const void *a, const void *b) {
  return (((const struct group_cm *)a)->cm.x - ((const struct group_cm *)b)->cm.x)>0;
}


int split_kobj(kobj *obj, kobj ***out, int attribute, int verbose) {
  /*code is based on https://en.wikipedia.org/wiki/Connected-component_labeling#One-pass_version*/
  /*and mincmorph*/
  const size_t max_matches=1024;
  size_t max_groups=obj->nvert+1;
  int *equiv;
  int *counts;
  int *group_labels;
  int *trans;
  int group_idx;
  int num_groups;
  int _neighbours[max_matches];
  int i;
  int last_group;
  int max_obs_matches=0;
  kvert *v;
  kedge *e;
  kface *f;
  equiv=(int*)calloc(max_groups, sizeof(int));
  counts=(int*)calloc(max_groups, sizeof(int));

  group_labels = (int*)calloc(obj->nvert+1, sizeof(int)); /*vert indexes are 1-based*/
  group_idx = 1;

  /*DEBUG*/
  /**/
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    int c=0;
    int found=0;
    for(c=0; c < v->nei; ++c) {
      int k;
      kvert *n=v->neivert[c];
      for(k=0;k<n->nei;++k)
      {
        if(n->neivert[k]==v){
          ++found;
        }
      }
    }
    if(found!=v->nei)
      abort();
  }
  /*DEBUG*/


  /*first pass*/
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    int _min_label = max_groups;
    int c;
    int prev_label;
    int matches=0;

    /*DEBUG*/
    if(v->nei>max_matches)
      abort();
    /*DEBUG*/

    for(c=0; c < v->nei && matches<max_matches; c++) {
      int _lab = group_labels[ v->neivert[c]->index ];
      if( _lab>0 ) {
        _neighbours[matches] = _lab;

        if(_lab < _min_label)
          _min_label = _lab;
        matches++;
      }
    }

    switch (matches) {
    case 0:
      /* no neighbours, make a new label and increment */
      group_labels[v->index] = group_idx;
      equiv[group_idx] = group_idx;
      counts[group_idx] = 1;
      group_idx++;
      break;

    case 1:
      /* only one neighbour, no equivalences needed */
      group_labels[ v->index ] = _min_label;
      counts[ _min_label ]++;
      break;

    default:
      /* more than one neighbour */

      /* first sort the neighbours array */
      qsort(&_neighbours, (size_t)matches, sizeof(int), &compare_ints);

      /* find the minimum possible label for this vertex,    */
      /* this is done by descending through each neighbours */
      /* equivalences until an equivalence equal to itself  */
      /* is found                                           */
      prev_label = -1;
      for(c = 0; c < matches; c++) {
        int _curr_label = _neighbours[c];

        /* recurse this label if we haven't yet */
        if(_curr_label != prev_label) {
          while(equiv[ _curr_label ] != equiv[ equiv[ _curr_label]] ) {
            _curr_label = equiv[_curr_label];
          }

          /* check against the current minimum value */
          if(equiv[_curr_label] < _min_label) {
            _min_label = equiv[_curr_label];
          }
        }

        prev_label = _neighbours[c];
      }

      /* repeat, setting equivalences to the _min_label */
      prev_label = -1;
      for(c = 0; c < matches; c++) {
        int _curr_label = _neighbours[c];

        if(_curr_label != prev_label) {
          while(equiv[_curr_label] != equiv[equiv[_curr_label]]) {
            _curr_label = equiv[_curr_label];
            equiv[_curr_label] = _min_label;
          }

          _curr_label = _neighbours[c];
          /* set the label itself */
          if(equiv[_curr_label] == 0) {             
             equiv[_curr_label] = _min_label;
          } else if( equiv[_curr_label] != _min_label ) {
            equiv[equiv[_curr_label]] = _min_label;
            equiv[_curr_label] = _min_label;
          }
        }

        prev_label = _curr_label;
      }

      /* finally set the vertex in question to the minimum value */
      group_labels[v->index] = _min_label;
      counts[_min_label]++;
      break;
    }
  }

  /* reduce the equiv and counts array */
  num_groups = 0;
  for(i = 1; i < group_idx; ++i) {
    /* if this equivalence is not resolved yet */
    if( i != equiv[i] ) {
      /* find the min label value */
      int _min_label = equiv[i];
      while(_min_label != equiv[_min_label]) {
        _min_label = equiv[_min_label];
      }

      /* update the label and its counters */
      equiv[i] = _min_label;
      counts[_min_label] += counts[i];
      counts[i] = 0;
    } else {
      num_groups++;
    }
  }
  *out=(kobj**)calloc(sizeof(kobj *),num_groups);

  //rename groups
  if(num_groups>1) {
    int g=0;
    trans = (int*)calloc(sizeof(int), group_idx+1);
    last_group=0;

    for(i = 1; i < group_idx; i++) {
      if( equiv[i] != last_group ) {
        g++;
        last_group = equiv[i];
      }
      trans[i] = g;
    }
    assert(num_groups==g);

    /*relabel labels*/
    for(i=1; i<=obj->nvert ; i++) {
      group_labels[i] = trans[ group_labels[i] ];
    }

    if(verbose)
      fprintf(stdout, "Found %d unique groups from %d\n", num_groups, group_idx);

    /*sort according to the x coordiantes */
    {
      /*for sorting groups*/
      struct group_cm* groups_cm;
      groups_cm=(struct group_cm *)calloc(sizeof(struct group_cm), g-1);

      for(v=obj->vert; v!=NULL; v=v->next) {
        int gi = group_labels[v->index]-1;
        groups_cm[gi].cm = niikpt_add(groups_cm[gi].cm, v->v);
        groups_cm[gi].c++;
      }

      for(i=0; i<num_groups; i++) {
        groups_cm[i].id=i+1;
        if(groups_cm[i].c>0)
          groups_cm[i].cm=niikpt_kmul( groups_cm[i].cm, 1.0/(double)groups_cm[i].c);
      }

      qsort(groups_cm, num_groups, sizeof(groups_cm), compare_groups);
      for(i=0; i<num_groups; i++) {
        trans[groups_cm[i].id] = i+1;
      }

      for(i=1; i<=obj->nvert ; i++) {
        group_labels[i] = trans[ group_labels[i] ];
      }

      free(groups_cm);
    }

    if(!attribute)
    {
      /*split-up object into subsets*/
      kface **last_face=(kface **)calloc(sizeof(kface *),num_groups);
      kedge **last_edge=(kedge **)calloc(sizeof(kedge *),num_groups);
      kvert **last_vert=(kvert **)calloc(sizeof(kvert *),num_groups);
      //printf("Number of groups:%d\n",num_groups);

      /*going to split up objects and update links*/
      for(i=0; i<num_groups; i++) {
        int k;
        (*out)[i] = off_obj_init();
        (*out)[i]->color=obj->color;
        (*out)[i]->spherecoo=obj->spherecoo;
        
        for(k=0;k<obj->n_comments;k++)
          off_kobj_add_comment((*out)[i],obj->comment[k]);
      }


      for(f=obj->face; f!=NULL; f=f->next) {
        /*relink*/
        int g=group_labels[f->vert[0]->index]-1; /*all corners should have the same group!*/
        if((*out)[g]->nface==0) { // initialize
          (*out)[g]->face=f;
          (*out)[g]->face->prev=NULL;
        } else {
          last_face[g]->next=f;
          f->prev=last_face[g];
        }
        last_face[g]=f;
        (*out)[g]->nface++;
        f->index=(*out)[g]->nface;
      }

      for(e=obj->edge; e!=NULL; e=e->next) {
        /*relink*/
        int g=group_labels[e->endpts[0]->index]-1; /*all corners should have the same group!*/
        if((*out)[g]->nedge==0) { // initialize
          (*out)[g]->edge=e;
          (*out)[g]->edge->prev=NULL;
        } else {
          last_edge[g]->next=e;
          e->prev=last_edge[g];
        }
        last_edge[g]=e;
        (*out)[g]->nedge++;
        e->index=(*out)[g]->nedge;
      }

      /*should go last to make sure indexes are correct*/
      for(v=obj->vert; v!=NULL; v=v->next) {
        /*relink*/
        int g=group_labels[v->index]-1;

        if((*out)[g]->nvert==0) { // initialize
          (*out)[g]->vert=v;
          (*out)[g]->vert->prev=NULL;
        } else {
          last_vert[g]->next=v;
          v->prev=last_vert[g];
        }
        last_vert[g]=v;
        (*out)[g]->nvert++;
        v->index=(*out)[g]->nvert;
      }


      /*terminate lists*/
      for(i=0; i<num_groups; i++) {
        last_vert[i]->next=NULL;
        last_face[i]->next=NULL;
        last_edge[i]->next=NULL;
      }

      obj->vert=NULL;
      obj->face=NULL;
      obj->edge=NULL;

      off_kobj_free(obj);
      free(last_face);
      free(last_edge);
      free(last_vert);
      free(trans);
    } else {
      // apply this as an attribute as additional value
      (*out)[0]=obj;
      for(v=obj->vert; v!=NULL; v=v->next) {
        v->idata=group_labels[v->index]-1;
      }
      num_groups=1;
    }
  } else {
    /*just keep input information*/
    (*out)[0]=obj;
  }
  free(equiv);
  free(counts);
  free(group_labels);
  return num_groups;
}

int main(int argc, char **argv) {
  const char *fcname="falcon_split_off";
  int clobber=0;
  int verbose=0;
  int attribute=0;
  int c;
  int i;
  float fwhm=1;
  const char *in_off=NULL;
  const char *out_prefix=NULL;
  int n;
  kobj *obj=NULL;
  kobj **out_obj;
  char *timestamp = niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"verbose", no_argument, &verbose, 1},
    {"attribute", no_argument, &attribute, 1},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "i", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<2) {
    show_usage(argv[0]);
    return 1;
  }

  in_off =argv[optind];
  out_prefix=argv[optind+1];


  /*
     if (!clobber && !access (out_prefix, F_OK))
     {
       fprintf(stderr,"%s Exists!\n", out_prefix);
       return 1;
     }*/
  if(verbose)
    niik_fc_display(fcname,1);

  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"off_kobj_read_offply",1 );

  NIIK_EXIT( (n=split_kobj( obj, &out_obj, attribute, verbose))==0,      fcname,"split_kobj",1 );

  if(attribute)
  {
    const char *meas_names[] = {"group"};
    double **meas = calloc(1,sizeof(double*));
    meas[0] = calloc(out_obj[0]->nvert,sizeof(double));

    kvert *v;

    for(v=out_obj[0]->vert;v!=NULL;v=v->next)
    {
      meas[0][v->index-1]=(double)v->idata;
    }
    NIIK_EXIT((!off_kobj_add_comment(out_obj[0],timestamp)),fcname,"off_kobj_add_comment",1);
    char tmp[1024];
    sprintf(tmp,"%s_c.ply",out_prefix);
    /* convert to a variable*/
    NIIK_EXIT((!off_kobj_write_ply_ex(tmp,out_obj[0],0,1,1,1,1,meas_names,meas)),fcname,"off_kobj_write_ply_ex",1);
    free(meas[0]);free(meas);
  } else {

    if(verbose) {
      for(i=0; i<n; i++) {
        fprintf(stdout,"Group %d nver=%d nedge=%d nface=%d\n",i,
                out_obj[i]->nvert,
                out_obj[i]->nedge,
                out_obj[i]->nface);
      }
    }

    for(i=0; i<n; i++) {
      char tmp[1024];
      sprintf(tmp,"%s_%d.ply",out_prefix,i);
      NIIK_EXIT((!off_kobj_add_comment(out_obj[i],timestamp)),fcname,"off_kobj_add_comment",1);
      NIIK_EXIT((!off_kobj_write_offply(tmp,out_obj[i],0)),fcname,"off_kobj_write_offply",1);
    }

    for(i=0; i<n; i++) {
      off_kobj_free(out_obj[i]);
    }
  }

  free(out_obj);

  if(verbose) niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
