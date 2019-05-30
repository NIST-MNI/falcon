/* Filename:     nifti1_kunio_off_curvature.c
 * Description:  Curvature calculation
 * Author:       Vladimir S. FONOV
 * Date:         March 1, 2012
 *
 * Reference:
 *
 * Rusinkiewicz, Szymon.
 * "Estimating Curvatures and Their Derivatives on Triangle Meshes,"
 * Proc. 3DPVT, 2004.
 *
 * http://gfx.cs.princeton.edu/proj/trimesh2/
 * https://gfx.cs.princeton.edu/pubs/_2004_ECA/index.php
 *
 * https://gfx.cs.princeton.edu/pubs/_2004_ECA/curvpaper.pdf
 *
 */
#include "falcon.h"
#include "falcon_surfaces.h"


#define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)
#define PREV(i) ((i)>0 ? (i)-1 : (i)+2)

#define ELEM_SWAP_D(a,b) { register double t=(a);(a)=(b);(b)=t; }
#define MAX(a,b) ( (a)>=(b)?(a):(b) )

/* Rotate a coordinate system to be perpendicular to the given normal*/
void rot_coord_sys(niikpt    old_u,
                  niikpt    old_v,
                  niikpt    new_norm,
                  niikpt    *new_u,
                  niikpt    *new_v) {
  niikpt old_norm;
  niikpt perp_old;
  niikpt dperp;
  double ndot;

  *new_u = old_u;
  *new_v = old_v;
  old_norm = niikpt_cross(old_u, old_v);
  ndot = niikpt_dot(old_norm,new_norm);

  if (ndot <= -1.0) {
    *new_u = niikpt_kmul(*new_u,-1.0);
    *new_v = niikpt_kmul(*new_v,-1.0);
    return;
  }

  perp_old = niikpt_sub(new_norm, niikpt_kmul(old_norm, 1.0*ndot));
  dperp = niikpt_kmul(niikpt_add(old_norm, new_norm), 1.0 / (1.0 + ndot) );

  *new_u = niikpt_sub(*new_u, niikpt_kmul(dperp, niikpt_dot(*new_u, perp_old)));
  *new_v = niikpt_sub(*new_v, niikpt_kmul(dperp, niikpt_dot(*new_v, perp_old)));
}


/* Reproject a curvature tensor from the basis spanned by old_u and old_v
 * (which are assumed to be unit-length and perpendicular) to the
 * new_u, new_v basis.
 */
void proj_curv( niikpt old_u,
                niikpt old_v,
                double old_ku,
                double old_kuv,
                double old_kv,
                niikpt new_u,
                niikpt new_v,
                double *new_ku,
                double *new_kuv,
                double *new_kv) {
  niikpt r_new_u, r_new_v;
  rot_coord_sys(new_u, new_v, niikpt_cross(old_u, old_v), &r_new_u, &r_new_v);

  double u1 = niikpt_dot(r_new_u, old_u);
  double v1 = niikpt_dot(r_new_u, old_v);
  double u2 = niikpt_dot(r_new_v, old_u);
  double v2 = niikpt_dot(r_new_v, old_v);

  *new_ku  = old_ku * u1*u1 + old_kuv * (2.0  * u1*v1)  + old_kv * v1*v1;
  *new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
  *new_kv  = old_ku * u2*u2 + old_kuv * (2.0  * u2*v2)  + old_kv * v2*v2;
}


/* Like the above, but for dcurv*/
/*
void proj_dcurv(const vec &old_u, const vec &old_v,
		const Vec<4> old_dcurv,
		const vec &new_u, const vec &new_v,
		Vec<4> &new_dcurv)
{
	vec r_new_u, r_new_v;
	rot_coord_sys(new_u, new_v, old_u CROSS old_v, r_new_u, r_new_v);

	float u1 = r_new_u DOT old_u;
	float v1 = r_new_u DOT old_v;
	float u2 = r_new_v DOT old_u;
	float v2 = r_new_v DOT old_v;

	new_dcurv[0] = old_dcurv[0]*u1*u1*u1 +
		       old_dcurv[1]*3.0f*u1*u1*v1 +
		       old_dcurv[2]*3.0f*u1*v1*v1 +
		       old_dcurv[3]*v1*v1*v1;
	new_dcurv[1] = old_dcurv[0]*u1*u1*u2 +
		       old_dcurv[1]*(u1*u1*v2 + 2.0f*u2*u1*v1) +
		       old_dcurv[2]*(u2*v1*v1 + 2.0f*u1*v1*v2) +
		       old_dcurv[3]*v1*v1*v2;
	new_dcurv[2] = old_dcurv[0]*u1*u2*u2 +
		       old_dcurv[1]*(u2*u2*v1 + 2.0f*u1*u2*v2) +
		       old_dcurv[2]*(u1*v2*v2 + 2.0f*u2*v2*v1) +
		       old_dcurv[3]*v1*v2*v2;
	new_dcurv[3] = old_dcurv[0]*u2*u2*u2 +
		       old_dcurv[1]*3.0f*u2*u2*v2 +
		       old_dcurv[2]*3.0f*u2*v2*v2 +
		       old_dcurv[3]*v2*v2*v2;
}*/


/* Given a curvature tensor, find principal directions and curvatures
 * Makes sure that pdir1 and pdir2 are perpendicular to normal
 */
void diagonalize_curv(
  niikpt old_u, niikpt old_v,
  double ku, double kuv, double kv,
  niikpt new_norm,
  niikpt *pdir1,
  niikpt *pdir2,
  double *k1,
  double *k2) {
  niikpt r_old_u, r_old_v;
  double c = 1.0, s = 0.0, tt = 0.0;

  /*rot_coord_sys(old_u, old_v, new_norm, &r_old_u, &r_old_v);*/

  if (kuv != 0.0) {
    /* Jacobi rotation to diagonalize*/
    double h = 0.5 * (kv - ku) / kuv;
    tt = (h < 0.0) ?
         1.0 / (h - sqrt(1.0 + h*h)) :
         1.0 / (h + sqrt(1.0 + h*h));

    c = 1.0 / sqrt(1.0 + tt*tt);
    s = tt * c;
  }

  rot_coord_sys(old_u, old_v, new_norm, &r_old_u, &r_old_v);


  *k1 = ku - tt * kuv;
  *k2 = kv + tt * kuv;

  if (fabs(*k1) >= fabs(*k2)) {
    *pdir1 = niikpt_sub( niikpt_kmul(r_old_u,c), niikpt_kmul(r_old_v, s));
  } else {
    double t=*k2;
    *k2=*k1;
    *k1=t;

    /*ELEM_SWAP_D(*k1, *k2);*/

    *pdir1 = niikpt_add( niikpt_kmul(r_old_u,s), niikpt_kmul(r_old_v, c));
  }
  *pdir2 = niikpt_cross(new_norm, *pdir1);
}


int niik_off_pointareas(kobj *obj,double *pointareas,double *cornerareas[3]) {
  int i;
  kface *f;

  for(i=0;i<obj->nvert;i++)
    pointareas[i]=0.0;
  for(i=0;i<obj->nface;i++)
    cornerareas[0][i]=cornerareas[1][i]=cornerareas[2][i]=0.0;

  /*TODO: cleanup pointareas ?*/
  /*#pragma omp parallel for private(i)*/
  for (f = obj->face,i=0; f!=NULL; f=f->next,i++) {
    /* Edges*/
    niikpt e[3];
    double area;
    double l2[3];
    double ew[3];

    e[0]=niikpt_sub(f->vert[2]->v, f->vert[1]->v);
    e[1]=niikpt_sub(f->vert[0]->v, f->vert[2]->v);
    e[2]=niikpt_sub(f->vert[1]->v, f->vert[0]->v);

    // Compute corner weights
    area = 0.5 * niikpt_mag(niikpt_cross(e[0],e[1]));

    l2[0] = niikpt_mag2(e[0]);
    l2[1] = niikpt_mag2(e[1]);
    l2[2] = niikpt_mag2(e[2]);

    ew[0] = l2[0] * (l2[1] + l2[2] - l2[0]);
    ew[1] = l2[1] * (l2[2] + l2[0] - l2[1]);
    ew[2] = l2[2] * (l2[0] + l2[1] - l2[2]);

    if (ew[0] <= 0.0) {
      cornerareas[1][i] = -0.25 * l2[2] * area / niikpt_dot(e[0], e[2]);
      cornerareas[2][i] = -0.25 * l2[1] * area / niikpt_dot(e[0], e[1]);
      cornerareas[0][i] = area - cornerareas[1][i] - cornerareas[2][i];
    } else if (ew[1] <= 0.0) {
      cornerareas[2][i] = -0.25 * l2[0] * area / niikpt_dot(e[1], e[0]);
      cornerareas[0][i] = -0.25 * l2[2] * area / niikpt_dot(e[1], e[2]);
      cornerareas[1][i] = area - cornerareas[2][i] - cornerareas[0][i];
    } else if (ew[2] <= 0.0f) {
      cornerareas[0][i] = -0.25 * l2[1] * area / niikpt_dot(e[2], e[1]);
      cornerareas[1][i] = -0.25 * l2[0] * area / niikpt_dot(e[2], e[0]);
      cornerareas[2][i] = area - cornerareas[0][i] - cornerareas[1][i];
    } else {
      int j;
      double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
      for ( j = 0; j < 3; j++)
        cornerareas[j][i] = ewscale * (ew[(j+1)%3] + ew[(j+2)%3]);
    }

    /*#pragma omp atomic*/
    pointareas[f->vert[0]->index-1] += cornerareas[0][i];
    /*#pragma omp atomic*/
    pointareas[f->vert[1]->index-1] += cornerareas[1][i];
    /*#pragma omp atomic*/
    pointareas[f->vert[2]->index-1] += cornerareas[2][i];
  }
  return 1;
}

/* Perform LDL^T decomposition of a symmetric positive definite matrix.
 * Like Cholesky, but no square roots.  Overwrites lower triangle of matrix.
 */
static inline int ldltdc(double A[3][3], double rdiag[3]) {

  double v[2]= {0.0, 0.0};
  int i;
  for (i = 0; i < 3; i++) {
    int k,j;
    for (k = 0; k < i; k++)
      v[k] = A[i][k] * rdiag[k];

    for (j = i; j < 3; j++) {
      double sum = A[i][j];
      int l;
      for (k = 0; k < i; k++)
        sum -= v[k] * A[j][k];

      if (i == j) {
        if (sum <= 0.0)
          return 0;
        rdiag[i] = 1.0 / sum;
      } else {
        A[j][i] = sum;
      }
    }
  }

  return 1;
}

/* Solve Ax=B after ldltdc */
static inline void ldltsl(double A[3][3], double rdiag[3], double B[3], double x[3]) {
  int i;
  for (i = 0; i < 3; i++) {
    double sum = B[i];
    int k;
    for ( k = 0; k < i; k++)
      sum -= A[i][k] * x[k];
    x[i] = sum * rdiag[i];
  }

  for (i = 2 ; i >= 0; i--) {
    double sum = 0.0;
    int k;
    for ( k = i + 1; k < 3; k++)
      sum += A[k][i] * x[k];

    x[i] -= sum * rdiag[i];
  }
}

int niik_off_curvature_map_rusinkiewicz(kobj *obj,
                                        double *pointareas,
                                        double *cornerareas[3],
                                        niikpt *pdir1, niikpt *pdir2,
                                        double *curv1, double *curv2,double *curv12)
/* calculates curvature map
 * -obj needs updated normal for vertices
 */
{
  kvert *v;
  kface *f;
  int i;

  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }

  /*Cleanup variables*/
  for(i=0;i<obj->nvert;i++)
  {
    curv1[i]=curv2[i]=curv12[i]=0.0;
    pdir1[i]=pdir2[2]=niikpt_zero();
  }

  /* Set up an initial coordinate system per vertex*/
  for (f=obj->face,i=0; f!=NULL; f=f->next,i++) {
    pdir1[f->vert[0]->index-1] = niikpt_sub(f->vert[1]->v, f->vert[0]->v);
    pdir1[f->vert[1]->index-1] = niikpt_sub(f->vert[2]->v, f->vert[1]->v);
    pdir1[f->vert[2]->index-1] = niikpt_sub(f->vert[0]->v, f->vert[2]->v);
  }

  for (v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    pdir1[i]  = niikpt_unit(niikpt_cross(pdir1[i], v->normal));
    pdir2[i]  = niikpt_cross(v->normal, pdir1[i]);
  }

  /* Compute curvature per-face*/
  for (f=obj->face,i=0; f!=NULL; f=f->next,i++) {
    // Edges
    niikpt e[3];
    niikpt t,n,b;
    double m[3] =    {  0.0, 0.0, 0.0 };
    double w[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
    int j;
    double diag[3] = { 0.0, 0.0, 0.0 };

    e[0]= niikpt_sub(f->vert[2]->v, f->vert[1]->v);
    e[1]= niikpt_sub(f->vert[0]->v, f->vert[2]->v);
    e[2]= niikpt_sub(f->vert[1]->v, f->vert[0]->v);

    /* N-T-B coordinate system per face*/
    t = niikpt_unit(e[0]);
    n = niikpt_cross(e[0], e[1]);
    b = niikpt_unit(niikpt_cross(n, t));

    /* Estimate curvature based on variation of normals
       along edges*/
    for (  j = 0; j < 3; j++) {
      niikpt dn ;
      double u = niikpt_dot(e[j], t);
      double v = niikpt_dot(e[j], b);
      double dnu,dnv;
      w[0][0] += u*u;
      w[0][1] += u*v;
      /*w[1][1] += v*v + u*u;
        w[1][2] += u*v; */
      w[2][2] += v*v;

      dn = niikpt_sub(f->vert[PREV(j)]->normal,
                      f->vert[NEXT(j)]->normal);

      dnu = niikpt_dot(dn,t);
      dnv = niikpt_dot(dn,b);

      m[0] += dnu*u;
      m[1] += dnu*v + dnv*u;
      m[2] += dnv*v;
    }
    w[1][1] = w[0][0] + w[2][2];
    w[1][2] = w[0][1];

    /* Least squares solution*/
    if (!ldltdc(w, diag)) {
      fprintf(stdout,"solution failed %i!\n",i);
      continue;
    }
    ldltsl(w, diag, m, m);

    // Push it back out to the vertices
    for (j = 0; j < 3; j++) {
      int vj = f->vert[j]->index-1;
      double c1, c12, c2;
      double wt;

      if(pointareas[vj]>0.0) {
        proj_curv(t, b,
                  m[0], m[1], m[2],
                  pdir1[vj], pdir2[vj],
                  &c1, &c12, &c2);

        wt = cornerareas[j][i] / pointareas[vj];
        /*#pragma omp atomic*/
        curv1[vj]  += wt * c1;
        /*#pragma omp atomic*/
        curv12[vj] += wt * c12;
        /*#pragma omp atomic*/
        curv2[vj]  += wt * c2;
      } 
    }
  }

  // Compute principal directions and curvatures at each vertex
  /*#pragma omp parallel for*/
  for (v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    diagonalize_curv(pdir1[i],  pdir2[i],
                     curv1[i],  curv12[i], curv2[i],
                     v->normal,
                     &pdir1[i], &pdir2[i],
                     &curv1[i], &curv2[i]);

  }

  return 1;
}


/* Compute derivatives of curvature.*/
/*
void TriMesh::need_dcurv()
{
	if (dcurv.size() == vertices.size())
		return;
	need_curvatures();

	dprintf("Computing dcurv... ");

	// Resize the arrays we'll be using
	int nv = vertices.size(), nf = faces.size();
	dcurv.clear(); dcurv.resize(nv);

	// Compute dcurv per-face
#pragma omp parallel for
	for (int i = 0; i < nf; i++) {
		// Edges
		vec e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
			     vertices[faces[i][0]] - vertices[faces[i][2]],
			     vertices[faces[i][1]] - vertices[faces[i][0]] };

		// N-T-B coordinate system per face
		vec t = e[0];
		normalize(t);
		vec n = e[0] CROSS e[1];
		vec b = n CROSS t;
		normalize(b);

		// Project curvature tensor from each vertex into this
		// face's coordinate system
		vec fcurv[3];
		for (int j = 0; j < 3; j++) {
			int vj = faces[i][j];
			proj_curv(pdir1[vj], pdir2[vj], curv1[vj], 0, curv2[vj],
				  t, b, fcurv[j][0], fcurv[j][1], fcurv[j][2]);

		}

		// Estimate dcurv based on variation of curvature along edges
		float m[4] = { 0, 0, 0, 0 };
		float w[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };
		for (int j = 0; j < 3; j++) {
			// Variation of curvature along each edge
			vec dfcurv = fcurv[PREV(j)] - fcurv[NEXT(j)];
			float u = e[j] DOT t;
			float v = e[j] DOT b;
			float u2 = u*u, v2 = v*v, uv = u*v;
			w[0][0] += u2;
			w[0][1] += uv;
			//w[1][1] += 2.0f*u2 + v2;
			//w[1][2] += 2.0f*uv;
			//w[2][2] += u2 + 2.0f*v2;
			//w[2][3] += uv;
			w[3][3] += v2;
			m[0] += u*dfcurv[0];
			m[1] += v*dfcurv[0] + 2.0f*u*dfcurv[1];
			m[2] += 2.0f*v*dfcurv[1] + u*dfcurv[2];
			m[3] += v*dfcurv[2];
		}
		w[1][1] = 2.0f * w[0][0] + w[3][3];
		w[1][2] = 2.0f * w[0][1];
		w[2][2] = w[0][0] + 2.0f * w[3][3];
		w[2][3] = w[0][1];

		// Least squares solution
		float d[4];
		if (!ldltdc<float,4>(w, d)) {
			//dprintf("ldltdc failed!\n");
			continue;
		}
		ldltsl<float,4>(w, d, m, m);
		Vec<4> face_dcurv(m);

		// Push it back out to each vertex
		for (int j = 0; j < 3; j++) {
			int vj = faces[i][j];
			Vec<4> this_vert_dcurv;
			proj_dcurv(t, b, face_dcurv,
				   pdir1[vj], pdir2[vj], this_vert_dcurv);
			float wt = cornerareas[i][j] / pointareas[vj];
			dcurv[vj] += wt * this_vert_dcurv;
		}
	}

	dprintf("Done.\n");


}*/


int re_init_curvature(off_curvature_t *crv, kobj *obj) {
  int n;
  for(n=0;n<3;n++) {
    crv->cornerareas[n]=(double*)realloc(crv->cornerareas[n],obj->nface*sizeof(double));
  }
  crv->pointareas = (double*)realloc(crv->pointareas,obj->nvert*sizeof(double));
  crv->pdir1  = (niikpt*) realloc(crv->pdir1, obj->nvert*sizeof(niikpt));
  crv->pdir2  = (niikpt*) realloc(crv->pdir2, obj->nvert*sizeof(niikpt));
  crv->curv1  = (double*) realloc(crv->curv1, obj->nvert*sizeof(double));
  crv->curv2  = (double*) realloc(crv->curv2, obj->nvert*sizeof(double));
  crv->curv12 = (double*) realloc(crv->curv12,obj->nvert*sizeof(double));
  return 1;
}


int init_curvature(off_curvature_t *crv,kobj *obj) {
  crv->cornerareas[0]=(double*)calloc(obj->nface,sizeof(double));
  crv->cornerareas[1]=(double*)calloc(obj->nface,sizeof(double));
  crv->cornerareas[2]=(double*)calloc(obj->nface,sizeof(double));

  crv->pointareas=(double*)calloc(obj->nvert,sizeof(double));
  crv->pdir1=(niikpt*) calloc(obj->nvert,sizeof(niikpt));
  crv->pdir2=(niikpt*) calloc(obj->nvert,sizeof(niikpt));
  crv->curv1=(double*) calloc(obj->nvert,sizeof(double));
  crv->curv2=(double*) calloc(obj->nvert,sizeof(double));
  crv->curv12=(double*)calloc(obj->nvert,sizeof(double));
  return 1;
}

int free_curvature(off_curvature_t *crv) {
  free(crv->curv12);
  free(crv->curv1);
  free(crv->curv2);
  free(crv->pdir1);
  free(crv->pdir2);
  free(crv->pointareas);
  free(crv->cornerareas[0]);
  free(crv->cornerareas[1]);
  free(crv->cornerareas[2]); 
  return 1;
}

int update_curvature(off_curvature_t *crv,kobj *obj) {
  int i;
  /*DEBUG*/

  /*DEBUG*/
  niik_off_pointareas(obj, crv->pointareas, crv->cornerareas);
  niik_off_curvature_map_rusinkiewicz(obj,
             crv->pointareas,crv->cornerareas,
             crv->pdir1,crv->pdir2,
             crv->curv1,crv->curv2,
             crv->curv12);
  return 1;
}

double mean_curvature(off_curvature_t *crv,kobj *obj) {
  int i;
  double mean_crv=0.0;
//#pragma omp parallel for private(i) reduction(+:mean_crv)
  for(i=0; i<obj->nvert; i++) {
    mean_crv+=crv->curv1[i]+crv->curv2[i];
  }

  return mean_crv/(2*obj->nvert);
}



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/