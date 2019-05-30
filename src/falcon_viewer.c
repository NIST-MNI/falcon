/* Filename:     nifti1_kunio_viewer.c
 * Description:  nifti image viewer by Kunio
 * Author:       Kunio Nakamura
 * Date:         February 25, 2012
 */

#include <getopt.h>

#include <gtk/gtk.h>
#include <glib/gprintf.h>
#include <gdk/gdk.h>
#include <gdk/gdkkeysyms.h>
#include "mincio.h"

#include "falcon.h"
#include "falcon_surfaces.h"

#define NVK_MAJOR_VERSION (0)
#define NVK_MINOR_VERSION (4)
#define NVK_MICRO_VERSION (6)

#define max_window_size 4096

static char * NVK_version[] = {
  "--------------------------------------------------\n"
  "  version 0.0     Kunio 2012-02-25\n"
  "  -testing beta version\n"
  "--------------------------------------------------\n"
  "  version 0.1.0   Kunio 2012-12-18\n"
  "  -added read function for MINC-1\n"
  "--------------------------------------------------\n"
  "  version 0.2.0   Kunio 2013-01-10\n"
  "  -added read function for MINC-2\n"
  "  --but not sure if the values are correct\n"
  "--------------------------------------------------\n"
  "  version 0.3.0   January 27, 2014, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added writing minc (minc2) file for mask editing\n"
  "--------------------------------------------------\n"
  "  version 0.4.0   February 13, 2014, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added adaptive pen (uses tolerance as in seed-fill)\n"
  "  -eraser did not change\n"
  "  -added check for x,y,z-dim for mask and newly loaded image\n"
  "--------------------------------------------------\n"
  "  version 0.4.1   2014-06-20  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed t-dimension image in mask editing mode\n"
  "--------------------------------------------------\n"
  "  version 0.4.2   2014-06-23  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed t-dimension image smart-pen image\n"
  "--------------------------------------------------\n"
  "  version 0.4.3   2014-06-23  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added option -mask-autosave\n"
  "  --automatic saving\n"
  "--------------------------------------------------\n"
  "  version 0.4.4   2014-06-25  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added a menu item to toggle between automatic saving\n"
  "--------------------------------------------------\n"
  "  version 0.4.5   2014-06-26  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed a bug, undoing without editing\n"
  "--------------------------------------------------\n"
  "  version 0.4.6   2014-09-06  Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed a bug, ndim check\n"
  "--------------------------------------------------\n"
};

static char * NVK_0_0_version[] = {
  "--------------------------------------------------\n"
  "  version 0.0     Kunio 2012-02-25\n"
  "  -testing beta version\n"
  "--------------------------------------------------\n"
  "  version 0.0.1   Kunio 2012-02-25\n"
  "  -fixed niikmat function problems\n"
  "  -added icon\n"
  "--------------------------------------------------\n"
  "  version 0.0.2   Kunio 2012-04-03\n"
  "  -fixed sampling distance for drawing objects\n"
  "  --now the sampling distance depends on the zoom\n"
  "    level\n"
  "--------------------------------------------------\n"
  "  version 0.0.3   Kunio 2012-05-07\n"
  "  -added better scroll handling for updating\n"
  "--------------------------------------------------\n"
  "  version 0.0.4   Kunio 2012-05-31\n"
  "  -added rounding for surfaces\n"
  "--------------------------------------------------\n"
  "  version 0.0.5   Kunio 2012-07-24\n"
  "  -added '-cidx' for initial image index\n"
  "--------------------------------------------------\n"
  "  version 0.0.6   Kunio 2012-08-28\n"
  "  -added a menu under view to switch on/off the naviation red-line and green-dot\n"
  "  -there may be a bug when a mouse moving outside the window with <Ctrl> key pressed\n"
  "--------------------------------------------------\n"
};

static char * NVK_usage[] = {
  "  command-line optional usage:\n"
  "  -u  --help       : show this message\n"
  "  -U               : show run-time usage\n"
  "  --version        : show version info\n"
  "  --debug          : set a debug num [default=0]\n"
  "  --zoom           : set 3 zoom values for x,y,z-axis [default=2.0 2.0 2.0]\n"
  "                     for example, --zoom 3,2,2\n"
  "  --idx            : initial position (ijk-index) [default=0,0,0]\n"
  "                     -for example, -idx 120,120,120\n"
  "  --cidx           : initial position = image center (ijk-index) [default = off]\n"
  "  --flipxyz        : flips in all x,y,z [default = no flip]\n"
  "  -o <obj.off> --off <obj1.off>[,<obj2.off>...] \n"
  "                   : add off object to display\n"
  "                     can be repeated arbitrary number of times\n"
  "                    -color comes from object itself\n"
  "                    -if there's no color, then it's green\n"
  "  --no-voxel-check --no-pixel-check\n"
  "                   : skip voxel size checking [default = ON]\n"
  "  --minmax <m>,<M>\n"
  "                   : pre-sets min and max values\n"
  "                   : for example, --minmax 0 140\n"
};

static char * NVK_runtime_usage[] = {
  "  run-time usage: \n"
  "\n"
  "  MOUSE usage: \n"
  "  mouse-over:          -show voxel intensity and voxel coordinate (ijk)\n"
  "                       -if qform_code is used, then another coordinate\n"
  "                        (xyz) in blue\n"
  "  left-click           -moves the current position to that voxel\n"
  "                       -'current position' has a little green marker\n"
  "                        along red axes\n"
  "  right-up/dorn drag   -zoom out/in respectively\n"
  "  middle drag          -move viewing plane\n"
  "  shift-right          -changes min/max intensities\n"
  "  scroll down/up       -goes thru planes (increment/decrement respectively\n"
  "\n"
  "  MOUSE usage with mask image\n"
  "  CONTROL + left       -adds with a pen\n"
  "  CONTROL + right      -deletes with a pen\n"
  "  CONTROL + middle     -fill \n"
  "  SHIFT + middle       -fill with intensity range (tol) \n"
  "\n"
  "  Window\n"
  "  -numbers on the lower left show the min/max intensities for \n"
  "   interactive scaling\n"
  "  -numbers on the lower right show the current coordiante: \n"
  "   that is, x,y,z,(t,u) positions\n"
  "  -intensity at the current position is shown on left of 'x'\n"
  "  -mouse-over position and its intensity is shown on the lower\n"
  "   left, next to the min/max spin buttons\n"
  "  -if qform_code is used, then the xyz coordinate after qform\n"
  "   transformation is shown in blue.\n"
  "\n"
  "\n"
};



/***********************************************************
 *
 * GLOBAL VARIABLES
 *
 *
 ***********************************************************/


int verbose=1; /* verbose level */

static GtkWidget *g_main_window;     /* main window */
static GtkWidget *g_spin_min;
static GtkWidget *g_spin_max;        /* spin for min/max intensities */
static GtkWidget *g_spin_dim[9];
static GtkWidget *g_spin_xdim;
static GtkWidget *g_spin_ydim;
static GtkWidget *g_spin_zdim;       /* spin for 3d space */
static GtkWidget *g_spin_tdim;       /* spin for 4th dim */
static GtkWidget *g_spin_udim;       /* spin for 5th dim */
static GtkWidget *g_lbl_voxel;       /* label for voxel position */
static GtkWidget *g_lbl_dim[9];      /* label for dimensions */
static GtkWidget *g_lbl_gray;        /* label for voxel intensity */
static char g_lbl_graystr[64];       /* string for g_lbl_gray */
static int g_window_width;
static int g_window_height;          /* window width and height; can be changed */
static GdkPixbuf *g_pixbuf[3];       /* pixel buffer (where we put pixel values) -axial/coronal/sagittal */
static GtkWidget *g_gtkimg[3];       /* image widget */
static nifti_image *g_image;         /* input image */
static kobj **g_obj_list=NULL;            /* input off object */
static int g_obj_list_num;           /* input off object */
static bbox *g_off_bb;               /* input off object's bounding box */
static float *g_image_data;          /* image data in float? or double? */
static int g_idx[9];                 /* current position (index) */
static char g_axis[3];               /* axis for left[1] and right[2] images */
static char g_flip[4];               /* flip for 3 directions (x,y,z = 1,2,3) */
static double g_zoom[4];             /* zoom level */
static double g_shift[4][3];         /* image shift */
static niikmat *g_matrix[4];         /* transformation (zoom/shift) matrix */
static niikpt g_img_delta,g_img_fov; /* image min/max */
static int g_xy[3];                  /* temp window x/y position */
static gchar g_curr_filename[4096];  /* current folder name */
static int g_coordinates=1;          /* coordinate system */

int g_draw_red_line=1;               /* flag to draw a red line */
int g_draw_green_dot=1;              /* flag to draw a green dot */

static nifti_image *g_maskimg;       /* input mask image */
static unsigned char *g_mask_data;   /* mask image data in uint8 */
static GtkWidget *g_spin_pen_size;   /* mask editing pen size */
static GtkWidget *g_spin_edit_blue;  /* mask editing blue factor */
static GtkWidget *g_spin_fill_tol;   /* mask editing -filling tolerance */
double g_pen_size;                   /* mask editing pen size */
double g_edit_blue;                  /* mask editing blue factor */
double g_edit_fill_tol;              /* mask editing -filling tolerance */
int **g_edit_undo_list;              /* mask editing undo
				      * undo_list[max_undo][long vector]
				      * undo_list[max_undo][0] = size of vector
				      * undo_list[max_undo][1] =  1 for writing and
				      *                          -1 for deleting */
int g_edit_max_undo;                 /* mask editing max undo number */
int g_edit_undo_idx;                 /* mask editing undo number */
int g_automatic_mask_save=0;         /* mask is automatically saved every minute */


/***********************************************************
 *
 * FUNCTION PROTOTYPES
 *
 * -these are internal functions
 *
 ***********************************************************/

static void view_change_axis1_action ();
static void view_change_axis2_action ();
static void NVK_window_expand();
static void NVK_window_shrink();
static void NVK_view_zoom_out();
static void NVK_view_zoom_in();

static gboolean scroll_event_imagen( GtkWidget *widget,GdkEventScroll *event,int nimg );
static gboolean scroll_event_image1( GtkWidget *widget,GdkEventScroll *event );
static gboolean scroll_event_image2( GtkWidget *widget,GdkEventScroll *event );
static gboolean motion_notify_event_imagen( GtkWidget *widget,GdkEventMotion *event, int nimg );
static gboolean motion_notify_event_image1( GtkWidget *widget,GdkEventMotion *event );
static gboolean motion_notify_event_image2( GtkWidget *widget,GdkEventMotion *event );
static gboolean button_release_event_imagen( GtkWidget *widget,GdkEventMotion *event,int nimg );
static gboolean button_release_event_image1( GtkWidget *widget,GdkEventMotion *event );
static gboolean button_release_event_image2( GtkWidget *widget,GdkEventMotion *event );
static gboolean button_press_event_imagen( GtkWidget *widget,GdkEventMotion *event,int nimg );
static gboolean button_press_event_image1( GtkWidget *widget,GdkEventMotion *event );
static gboolean button_press_event_image2( GtkWidget *widget,GdkEventMotion *event );

static void NVK_update_matrix(int nimg);
static niikpt NVK_window2image(int nimg,int xi,int yi);
static int NVK_display_image();
static int NVK_display_imagen(int nimg);

static int NVK_update_zoom(int y, int y0, int nimg);
static int NVK_update_shift(int x, int x0, int y, int y0, int nimg);

static void NVK_update_lbl_gray();
static void NVK_update_lbl_voxel(niikpt q);

static int NVK_update_mask_seed_fill_tol(int x, int y, int z, int nimg, char dir, double lotol,double hitol);
static int NVK_update_mask_seed_fill(int x, int y, int z, int nimg, char dir);
static int NVK_update_mask_dilation (int x, int y, int z, int nimg, char dir, double pen_size);
static int NVK_update_mask_erosion  (int x, int y, int z, int nimg, char dir, double pen_size);

static void NVK_change_spin_min();
static void NVK_change_spin_max();

static void NVK_change_spin_xdim();
static void NVK_change_spin_ydim();
static void NVK_change_spin_zdim();
static void NVK_change_spin_tdim();
static void NVK_change_spin_udim();

static void NVK_write_mask(char *filename,nifti_image *maskimg);
static int NVK_write_minc_mask(char *filename,nifti_image *maskimg);
static gboolean NVK_automatic_save_mask_file(GtkWidget *widget);

GdkPixbuf *add_create_pixbuf(const gchar * filename);







GdkPixbuf *add_create_pixbuf(const gchar *filename) {
  GdkPixbuf *pixbuf;
  GError *error = NULL;
  pixbuf = gdk_pixbuf_new_from_file(filename,&error);
  if(!pixbuf) {
    fprintf(stderr," %s\n",error->message);
    g_error_free(error);
  }
  return pixbuf;
}



/* -calculates the index for a diretion
 * -not to be confused with nimg which can be 1 or 2
 * -this number is x,y,z -> 1,2,3 */
int dir2num(char dir) {
  switch(dir) {
  case 'x':
    return 1;
  case 'y':
    return 2;
  case 'z':
    return 3;
  }
  fprintf(stderr,"ERROR: unknown dir, %c\n",dir);
  exit(0);
  return 0;
}


/* -updates zoom for a particular direction
 * -returns nonzero if changed
 * -retruns zero if no change was made
 * -y and y0 are the vertical mouse potision in the image */
static int NVK_update_zoom(int y, int y0, int nimg) {
  double dy;
  dy = 2.0 * (y-y0) / g_window_height;
  if(fabs(dy)<0.1) return 0;
  if(verbose>2) fprintf(stdout,"-d3 updating zoom nimg=%i, axis[%i]=%c, num=%i, zoom=%8.3f\n",
                          nimg,nimg,g_axis[nimg],dir2num(g_axis[nimg]),g_zoom[dir2num(g_axis[nimg])]);
  g_zoom[dir2num(g_axis[nimg])] *= (1.0 + dy/3.0);
  return 1;
}

/* -updates shift
 * -returns zero if no change
 * -returns non-zero if shift was updated
 * */
static int NVK_update_shift(int x, int x0, int y, int y0, int nimg) {
  double dx,dy,czoom;
  int num;
  num=dir2num(g_axis[nimg]);
  czoom = 0.6*g_zoom[num] + 0.4;
  dx = (x0 - x) / czoom;
  dy = (y0 - y) / czoom;
  if(fabs(dx)<3) if(fabs(dy)<3) return 0;
  g_shift[num][1] -= dx;
  g_shift[num][2] -= dy;
  if(verbose>2) fprintf(stdout,"-d3 shift[%i][1-2] = %f %f\n",num,g_shift[num][1],g_shift[num][2]);
  return 1;
} /* NVK_update_shift */


static int NVK_update_intensity(int x, int x0, int y, int y0, int nimg) {
  double dx,dy,imin,imax,cmax;
  int
  i,ax,ay;
  ax=x-x0;
  ay=y0-y;
  if(fabs(ax)<4) if(fabs(ay)<4) return 0;  /* skip if small change */
  gtk_spin_button_get_range(GTK_SPIN_BUTTON(g_spin_min),&imin,&imax);
  if(imax<imin) {
    cmax=imax;
    imax=imin;
    imin=cmax;
  }
  cmax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_max));
  if     (cmax<imax/256) imax=(imax-imin)/512.0;
  else if(cmax<imax/128) imax=(imax-imin)/256.0;
  else if(cmax<imax/64)  imax=(imax-imin)/128.0;
  else if(cmax<imax/32)  imax=(imax-imin)/64.0;
  else if(cmax<imax/16)  imax=(imax-imin)/32.0;
  else if(cmax<imax/8)   imax=(imax-imin)/16.0;
  else if(cmax<imax/4)   imax=(imax-imin)/8.0;
  else if(cmax<imax/2)   imax=(imax-imin)/4.0;
  else                   imax=(imax-imin)/2.0;
  dx=ax*imax/g_window_width ;
  dy=ay*imax/g_window_height;
  imin = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_min));
  if(verbose>2) fprintf(stdout,"-d3 imin=%f->%f imax=%f->%f\n",imin,imin+dx,cmax,cmax+dy);
  imin += dx;
  cmax += dy;
  i=g_xy[0];
  g_xy[0]=0;
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_min),imin);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_max),cmax);
  g_xy[0]=i;
  return 1;
} /* NVK_update_intensity */



/*
 * NVK_update_mask_dilation
 *
 * -writing with a pen
 * -with undo feature
 *
 */

static int NVK_update_mask_dilation(int x, int y, int z, int nimg, char dir, double pen_size) {
  int
  nvox=0,
  xmin,xmax,ymin,
  ymax,zmin,zmax,
  i,j,k,n0,n,nn,p;
  double dd1,dd2,ps2,lotol=1e30,hitol=-1e30;
  float *image_data=NULL;

  if(g_maskimg==NULL) return 0;

  ps2 = pen_size * pen_size;
  n0 = x + y * g_maskimg->nx + z * g_maskimg->nx * g_maskimg->ny;
  image_data = g_image_data;
  if(g_image->ndim>3) { /* move to the correct memory space for 3D image */
    if(g_image->nt>1) {
      image_data += g_image->nx * g_image->ny * g_image->nz * g_idx[4];
    }
  }
  lotol = image_data[n0] - g_edit_fill_tol;
  hitol = image_data[n0] + g_edit_fill_tol;

  switch(dir) {
  case 'x':
    n = x;
    ymin = NIIK_IMINMAX(y-pen_size/g_maskimg->pixdim[2]-0.5,0,g_maskimg->dim[2]-1);
    ymax = NIIK_IMINMAX(y+pen_size/g_maskimg->pixdim[2]+0.5,0,g_maskimg->dim[2]-1);
    zmin = NIIK_IMINMAX(z-pen_size/g_maskimg->pixdim[3]-0.5,0,g_maskimg->dim[3]-1);
    zmax = NIIK_IMINMAX(z+pen_size/g_maskimg->pixdim[3]+0.5,0,g_maskimg->dim[3]-1);
    for(k=zmin; k<=zmax; k++) {
      dd2 = NIIK_SQ((k - z)*g_maskimg->pixdim[3]);
      if(dd2 > ps2) continue;
      nn = n + k * g_maskimg->dim[1] * g_maskimg->dim[2];
      for(j=ymin; j<=ymax; j++) {
        dd1 = NIIK_SQ((j - y)*g_maskimg->pixdim[2]);
        if(dd2 + dd1 > ps2) continue;
        p = nn + j*g_maskimg->dim[1];
        if(g_mask_data[p]==0) { /* update mask */
          if(!(g_xy[0] & GDK_SHIFT_MASK)) { /* SHIFT MASK IS OFF */
            g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
            g_mask_data[p] = 1;
            nvox++;
          } else { /* SHIFT MASK IS ON */
            if(image_data[p]>lotol && image_data[p]<hitol) {
              g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
              g_mask_data[p] = 1;
              nvox++;
            } /* check the tolerance */
          } /* SHIFT MASK IS ON */
        } /* updating */
      }
    }
    break;

  case 'y':
    n = y * g_maskimg->dim[1];
    xmin = NIIK_IMINMAX(x-pen_size/g_maskimg->pixdim[1]-0.5,0,g_maskimg->dim[1]-1);
    xmax = NIIK_IMINMAX(x+pen_size/g_maskimg->pixdim[1]+0.5,0,g_maskimg->dim[1]-1);
    zmin = NIIK_IMINMAX(z-pen_size/g_maskimg->pixdim[3]-0.5,0,g_maskimg->dim[3]-1);
    zmax = NIIK_IMINMAX(z+pen_size/g_maskimg->pixdim[3]+0.5,0,g_maskimg->dim[3]-1);
    for(k=zmin; k<=zmax; k++) {
      dd2 = NIIK_SQ((k - z)*g_maskimg->pixdim[3]);
      if(dd2 > ps2) continue;
      nn = n + k * g_maskimg->dim[1] * g_maskimg->dim[2];
      for(i=xmin; i<=xmax; i++) {
        dd1 = NIIK_SQ((i - x)*g_maskimg->pixdim[1]);
        if(dd2 + dd1 > ps2) continue;
        p = nn + i ;
        if(g_mask_data[p]==0) { /* update mask */
          if(!(g_xy[0] & GDK_SHIFT_MASK)) {
            g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
            g_mask_data[p] = 1;
            nvox++;
          } else {  /* SHIFT MASK IS OFF */
            if(image_data[p]>lotol && image_data[p]<hitol) {
              g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
              g_mask_data[p] = 1;
              nvox++;
            } /* check the tolerance */
          } /* SHIFT MASK IS ON */
        } /* updating */
      }
    }
    break;

  case 'z':
    n = z * g_maskimg->dim[1] * g_maskimg->dim[2];
    xmin = NIIK_IMINMAX(x-pen_size/g_maskimg->pixdim[1]-0.5,0,g_maskimg->dim[1]-1);
    xmax = NIIK_IMINMAX(x+pen_size/g_maskimg->pixdim[1]+0.5,0,g_maskimg->dim[1]-1);
    ymin = NIIK_IMINMAX(y-pen_size/g_maskimg->pixdim[2]-0.5,0,g_maskimg->dim[2]-1);
    ymax = NIIK_IMINMAX(y+pen_size/g_maskimg->pixdim[2]+0.5,0,g_maskimg->dim[2]-1);
    /*fprintf(stdout,"  zdir  %i:%i:%i  %i:%i:%i\n",xmin,x,xmax,ymin,y,ymax);
      fprintf(stdout,"        %i\n",GDK_SHIFT_MASK);*/
    for(j=ymin; j<=ymax; j++) {
      dd2 = NIIK_SQ((j - y)*g_maskimg->pixdim[2]);
      if(dd2 > ps2) continue;
      nn = n + j * g_maskimg->dim[1];
      for(i=xmin; i<=xmax; i++) {
        dd1 = NIIK_SQ((i - x)*g_maskimg->pixdim[1]);
        if(dd2 + dd1 > ps2) continue;
        /* fprintf(stdout,"  dis  %f %f > %f\n",dd1,dd2,ps2); */
        p = nn + i;
        if(g_mask_data[p]==0) { /* update mask */
          if(!(g_xy[0] & GDK_SHIFT_MASK)) {
            g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
            g_mask_data[p] = 1;
            nvox++;
          } else {  /* SHIFT MASK IS OFF */
            if(image_data[p]>lotol && image_data[p]<hitol) {
              g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
              g_mask_data[p] = 1;
              nvox++;
            } /* check the tolerance */
          } /* SHIFT MASK IS ON */
        } /* updating */
      }
    }
    break;

  }

  return nvox;
} /* NVK_update_mask_dilation */



/*
 * NVK_update_mask_erosion
 *
 * -erasing with a pen
 * -with undo feature
 *
 */

static int NVK_update_mask_erosion(int x, int y, int z, int nimg, char dir, double pen_size) {
  int
  nvox = 0,
  xmin,xmax,ymin,
  ymax,zmin,zmax,
  i,j,k,n,nn,p;
  double dd1,dd2,ps2;

  if(g_maskimg==NULL) return 0;
  ps2 = pen_size * pen_size;

  switch(dir) {
  case 'x':
    n = x;
    ymin = NIIK_IMINMAX(y-pen_size/g_maskimg->pixdim[2]-0.5,0,g_maskimg->dim[2]-1);
    ymax = NIIK_IMINMAX(y+pen_size/g_maskimg->pixdim[2]+0.5,0,g_maskimg->dim[2]-1);
    zmin = NIIK_IMINMAX(z-pen_size/g_maskimg->pixdim[3]-0.5,0,g_maskimg->dim[3]-1);
    zmax = NIIK_IMINMAX(z+pen_size/g_maskimg->pixdim[3]+0.5,0,g_maskimg->dim[3]-1);
    for(k=zmin; k<=zmax; k++) {
      dd2 = NIIK_SQ((k - z)*g_maskimg->pixdim[3]);
      if(dd2 > ps2) continue;
      nn = n + k * g_maskimg->dim[1] * g_maskimg->dim[2];
      for(j=ymin; j<=ymax; j++) {
        dd1 = NIIK_SQ((j - y)*g_maskimg->pixdim[2]);
        if(dd2 + dd1 > ps2) continue;
        p = nn + j*g_maskimg->dim[1];
        if(g_mask_data[p]>0) {
          /* update mask */
          g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
          g_mask_data[p] = 0;
          nvox++;
        }
      }
    }
    break;

  case 'y':
    n = y * g_maskimg->dim[1];
    xmin = NIIK_IMINMAX(x-pen_size/g_maskimg->pixdim[1]-0.5,0,g_maskimg->dim[1]-1);
    xmax = NIIK_IMINMAX(x+pen_size/g_maskimg->pixdim[1]+0.5,0,g_maskimg->dim[1]-1);
    zmin = NIIK_IMINMAX(z-pen_size/g_maskimg->pixdim[3]-0.5,0,g_maskimg->dim[3]-1);
    zmax = NIIK_IMINMAX(z+pen_size/g_maskimg->pixdim[3]+0.5,0,g_maskimg->dim[3]-1);
    if(verbose) fprintf(stdout,"-d3 mask_erosion y %i [x,z] [%3i %3i] -> [%3i %3i]  %i\n",y,xmin,zmin,xmax,zmax,n);
    for(k=zmin; k<=zmax; k++) {
      dd2 = NIIK_SQ((k - z)*g_maskimg->pixdim[3]);
      if(dd2 > ps2) continue;
      nn = n + k * g_maskimg->dim[1] * g_maskimg->dim[2];
      for(i=xmin; i<=xmax; i++) {
        dd1 = NIIK_SQ((i - x)*g_maskimg->pixdim[1]);
        if(dd2 + dd1 > ps2) continue;
        p = nn+i;
        if(g_mask_data[p]>0) {
          /* update mask */
          g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
          g_mask_data[p] = 0;
          nvox++;
        }
      }
    }
    break;

  case 'z':
    n = z * g_maskimg->dim[1] * g_maskimg->dim[2];
    xmin = NIIK_IMINMAX(x-pen_size/g_maskimg->pixdim[1]-0.5,0,g_maskimg->dim[1]-1);
    xmax = NIIK_IMINMAX(x+pen_size/g_maskimg->pixdim[1]+0.5,0,g_maskimg->dim[1]-1);
    ymin = NIIK_IMINMAX(y-pen_size/g_maskimg->pixdim[2]-0.5,0,g_maskimg->dim[2]-1);
    ymax = NIIK_IMINMAX(y+pen_size/g_maskimg->pixdim[2]+0.5,0,g_maskimg->dim[2]-1);
    /*fprintf(stdout,"  zdir  %i:%i:%i  %i:%i:%i\n",xmin,x,xmax,ymin,y,ymax);*/
    for(j=ymin; j<=ymax; j++) {
      dd2 = NIIK_SQ((j - y)*g_maskimg->pixdim[2]);
      if(dd2 > ps2) continue;
      nn = n + j * g_maskimg->dim[1];
      for(i=xmin; i<=xmax; i++) {
        dd1 = NIIK_SQ((i - x)*g_maskimg->pixdim[1]);
        if(dd2 + dd1 > ps2) continue;
        /* fprintf(stdout,"  dis  %f %f > %f\n",dd1,dd2,ps2); */
        p = nn + i;
        if(g_mask_data[p]>0) {
          /* update mask */
          g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = p;
          g_mask_data[p] = 0;
          nvox++;
        }
      }
    }
    break;

  }

  return nvox;
} /* NVK_update_mask_erosion */


/*
 * default seed fill
 *
 * -filling function
 * -uses the seed fill with tolerance NVK_update_mask_seed_fill_tol
 *  with really low/high tolerance
 */

static int NVK_update_mask_seed_fill(int x, int y, int z, int nimg, char dir) {
  return NVK_update_mask_seed_fill_tol(x,y,z,nimg,dir,-1e30,1e30);
}




/*
 * actual seed fill function
 *
 * -filling function
 * -faster than the previous image viewer
 * -creates 2 large arrays and upates the current positions
 * -the old one was slower because it check the whole image
 *  at each iteration
 * -this one saves the current locations
 * -supports undo
 *
 */

static int NVK_update_mask_seed_fill_tol(int x, int y, int z, int nimg, char dir, double lotol,double hitol) {
  int
  nvox = 0,
  n,num, *v1,*v2,
  idx[4],area,xdim,
  cnt,
  verbo = 0;

  if(g_maskimg==NULL) return 0;
  if(verbo) fprintf(stdout,"  seed fill function\n");

  xdim = g_image->dim[1];
  area = g_image->dim[1]*g_image->dim[2];

  switch(dir) {
  case 'x':

    if(verbo) fprintf(stdout,"    x-dir\n");
    v1 = (int *)calloc(g_image->dim[2]*g_image->dim[3],sizeof(int));
    v2 = (int *)calloc(g_image->dim[2]*g_image->dim[3],sizeof(int));

    v1[0] = x + y*xdim + z*area;
    num = cnt = 1;

    g_edit_undo_idx = (g_edit_undo_idx+1)%g_edit_max_undo;
    g_edit_undo_list[g_edit_undo_idx][0]=2;
    g_edit_undo_list[g_edit_undo_idx][1]=1;

    if(g_mask_data[v1[0]]==0) {
      g_mask_data[v1[0]]=1;
      g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[0];
    }

    if(verbo) fprintf(stdout,"  enter loop\n");
    while(cnt) {
      num = cnt;
      for(n=0,cnt=0; n<num; n++) {
        idx[2] = (int)floor((v1[n]%area)/xdim);
        idx[3] = (int)floor(v1[n]/area);

        if(idx[2]>0) {
          if(g_mask_data[v1[n]-xdim]==0) {
            if(g_image_data[v1[n]-xdim]>lotol) {
              if(g_image_data[v1[n]-xdim]<hitol) {
                g_mask_data[v1[n]-xdim]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-xdim;
                nvox++;
                v2[cnt++] = v1[n]-xdim;
              }
            }
          }
        }
        if(idx[2]<g_image->dim[2]-1) {
          if(g_mask_data[v1[n]+xdim]==0) {
            if(g_image_data[v1[n]+xdim]>lotol) {
              if(g_image_data[v1[n]+xdim]<hitol) {
                g_mask_data[v1[n]+xdim]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+xdim;
                nvox++;
                v2[cnt++] = v1[n]+xdim;
              }
            }
          }
        }

        if(idx[3]>0) {
          if(g_mask_data[v1[n]-area]==0) {
            if(g_image_data[v1[n]-area]>lotol) {
              if(g_image_data[v1[n]-area]<hitol) {
                g_mask_data[v1[n]-area]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-area;
                nvox++;
                v2[cnt++] = v1[n]-area;
              }
            }
          }
        }
        if(idx[3]<g_image->dim[3]-1) {
          if(g_mask_data[v1[n]+area]==0) {
            if(g_image_data[v1[n]+area]>lotol) {
              if(g_image_data[v1[n]+area]<hitol) {
                g_mask_data[v1[n]+area]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+area;
                nvox++;
                v2[cnt++] = v1[n]+area;
              }
            }
          }
        }
      } /* each voxels in the list */

      /* update the list */
      for(n=0; n<cnt; n++) {
        v1[n]=v2[n];
      }

      if(verbo) fprintf(stdout,"  count = %i \n",cnt);
    } /* while there's a change */

    free(v1);
    free(v2);


    break;


  case 'y':
    if(verbo) fprintf(stdout,"    y-dir\n");
    v1 = (int *)calloc(g_image->dim[1]*g_image->dim[3],sizeof(int));
    v2 = (int *)calloc(g_image->dim[1]*g_image->dim[3],sizeof(int));

    if(verbo) fprintf(stdout,"    first pos\n");
    v1[0] = x + y*xdim + z*area;
    num = cnt = 1;

    if(verbo) fprintf(stdout,"    first pos %i\n",v1[0]);
    g_edit_undo_idx = (g_edit_undo_idx+1)%g_edit_max_undo;
    g_edit_undo_list[g_edit_undo_idx][0]=2;
    g_edit_undo_list[g_edit_undo_idx][1]=1;

    if(g_mask_data[v1[0]]==0) {
      g_mask_data[v1[0]]=1;
      g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[0];
    }

    if(verbo) fprintf(stdout,"  enter loop\n");
    while(cnt) {
      num = cnt;
      for(n=0,cnt=0; n<num; n++) {
        idx[1] = v1[n]%xdim;
        idx[3] = (int)floor(v1[n]/area);

        if(idx[1]>0) {
          if(g_mask_data[v1[n]-1]==0) {
            if(g_image_data[v1[n]-1]>lotol) {
              if(g_image_data[v1[n]-1]<hitol) {
                g_mask_data[v1[n]-1]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-1;
                nvox++;
                v2[cnt++] = v1[n]-1;
              }
            }
          }
        }
        if(idx[1]<g_image->dim[1]-1) {
          if(g_mask_data[v1[n]+1]==0) {
            if(g_image_data[v1[n]+1]>lotol) {
              if(g_image_data[v1[n]+1]<hitol) {
                g_mask_data[v1[n]+1]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+1;
                nvox++;
                v2[cnt++] = v1[n]+1;
              }
            }
          }
        }

        if(idx[3]>0) {
          if(g_mask_data[v1[n]-area]==0) {
            if(g_image_data[v1[n]-area]>lotol) {
              if(g_image_data[v1[n]-area]<hitol) {
                g_mask_data[v1[n]-area]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-area;
                nvox++;
                v2[cnt++] = v1[n]-area;
              }
            }
          }
        }
        if(idx[3]<g_image->dim[3]-1) {
          if(g_mask_data[v1[n]+area]==0) {
            if(g_image_data[v1[n]+area]>lotol) {
              if(g_image_data[v1[n]+area]<hitol) {
                g_mask_data[v1[n]+area]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+area;
                nvox++;
                v2[cnt++] = v1[n]+area;
              }
            }
          }
        }
      } /* each voxels in the list */

      /* update the list */
      for(n=0; n<cnt; n++) {
        v1[n]=v2[n];
      }

      if(verbo) fprintf(stdout,"  count = %i \n",cnt);
    } /* while there's a change */

    free(v1);
    free(v2);

    break;

  case 'z':

    if(verbo) fprintf(stdout,"    z-dir\n");
    v1 = (int *)calloc(g_image->dim[1]*g_image->dim[2],sizeof(int));
    v2 = (int *)calloc(g_image->dim[1]*g_image->dim[2],sizeof(int));

    if(verbo) fprintf(stdout,"    first pos\n");
    v1[0] = x + y*xdim + z*area;
    num = cnt = 1;

    if(verbo) fprintf(stdout,"    first pos %i\n",v1[0]);
    g_edit_undo_idx = (g_edit_undo_idx+1)%g_edit_max_undo;
    g_edit_undo_list[g_edit_undo_idx][0]=2;
    g_edit_undo_list[g_edit_undo_idx][1]=1;

    if(g_mask_data[v1[0]]==0) {
      g_mask_data[v1[0]]=1;
      g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[0];
    }

    if(verbo) fprintf(stdout,"  enter loop\n");
    while(cnt) {
      num = cnt;
      for(n=0,cnt=0; n<num; n++) {
        idx[1] = v1[n]%xdim;
        idx[2] = (int)floor((v1[n]%area)/xdim);

        if(idx[1]>0) {
          if(g_mask_data[v1[n]-1]==0) {
            if(g_image_data[v1[n]-1]>lotol) {
              if(g_image_data[v1[n]-1]<hitol) {
                g_mask_data[v1[n]-1]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-1;
                nvox++;
                v2[cnt++] = v1[n]-1;
              }
            }
          }
        }
        if(idx[1]<g_image->dim[1]-1) {
          if(g_mask_data[v1[n]+1]==0) {
            if(g_image_data[v1[n]+1]>lotol) {
              if(g_image_data[v1[n]+1]<hitol) {
                g_mask_data[v1[n]+1]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+1;
                nvox++;
                v2[cnt++] = v1[n]+1;
              }
            }
          }
        }

        if(idx[2]>0) {
          if(g_mask_data[v1[n]-xdim]==0) {
            if(g_image_data[v1[n]-xdim]>lotol) {
              if(g_image_data[v1[n]-xdim]<hitol) {
                g_mask_data[v1[n]-xdim]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]-xdim;
                nvox++;
                v2[cnt++] = v1[n]-xdim;
              }
            }
          }
        }
        if(idx[2]<g_image->dim[2]-1) {
          if(g_mask_data[v1[n]+xdim]==0) {
            if(g_image_data[v1[n]+xdim]>lotol) {
              if(g_image_data[v1[n]+xdim]<hitol) {
                g_mask_data[v1[n]+xdim]=1;
                g_edit_undo_list[g_edit_undo_idx][g_edit_undo_list[g_edit_undo_idx][0]++] = v1[n]+xdim;
                nvox++;
                v2[cnt++] = v1[n]+xdim;
              }
            }
          }
        }
      } /* each voxels in the list */

      for(n=0; n<cnt; n++) {
        v1[n]=v2[n];
      }

      if(verbo) fprintf(stdout,"  count = %i \n",cnt);
    } /* while there's a change */

    free(v1);
    free(v2);

    break;

  } /* dir */

  return nvox;
} /* static int NVK_update_mask_seed_fill */


/*
 * NVK_change_spin_min
 * NVK_change_spin_max
 *
 * -change min/max intensity
 * -min can be larger than max
 *
 */

static void NVK_change_spin_min() {
  /*gdouble imin,imax;
    gtk_spin_button_get_range (GTK_SPIN_BUTTON(g_spin_max), &imin, &imax);
    imin = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_min)) + 1;
    gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_max), imin, imax ); */
  if(g_image==NULL) return;
  if(g_xy[0])
    NVK_display_image();
}

static void NVK_change_spin_max() {
  /* gdouble imin,imax;
     gtk_spin_button_get_range           (GTK_SPIN_BUTTON(g_spin_min), &imin, &imax);
     imax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_max)) - 1;
     gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_min), imin, imax ); */
  if(g_image==NULL) return;
  if(g_xy[0])
    NVK_display_image();
}




/*
 * void NVK_change_spin_xdim()
 * void NVK_change_spin_ydim()
 * void NVK_change_spin_zdim()
 * void NVK_change_spin_tdim()
 * void NVK_change_spin_udim()
 *
 *
 * -change the current location -up to 5th
 *
 *
 */


static void NVK_change_spin_xdim() {
  if(g_image==NULL) return;
  g_idx[1] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_spin_xdim));
  NVK_update_lbl_gray();
  if(g_xy[0]!=0)
    NVK_display_image();
}

static void NVK_change_spin_ydim() {
  if(g_image==NULL) return;
  g_idx[2] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_spin_ydim));
  NVK_update_lbl_gray();
  if(g_xy[0]!=0)
    NVK_display_image();
}

static void NVK_change_spin_zdim() {
  if(g_image==NULL) return;
  g_idx[3] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_spin_zdim));
  NVK_update_lbl_gray();
  if(g_xy[0]!=0)
    NVK_display_image();
}

static void NVK_change_spin_tdim() {
  if(g_image==NULL) return;
  g_idx[4] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_spin_tdim));
  NVK_update_lbl_gray();
  NVK_display_image();
}

static void NVK_change_spin_udim() {
  if(g_image==NULL) return;
  g_idx[5] = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(g_spin_udim));
  NVK_update_lbl_gray();
  NVK_display_image();
}








/*
 * NVK_update_lbl_voxel
 * NVK_update_lbl_gray
 *
 *
 * -updates the voxel label
 * -g_lbl_voxel label changes without clicking (left side)
 * -g_lbl_gray label changes with clicking (right side)
 *
 */

static void NVK_update_lbl_voxel(niikpt q) {
  char
  *markup,             /* pango markup
			  *  http://developer.gnome.org/gtk/2.24/GtkLabel.html#gtk-label-set-markup
			  *  http://developer.gnome.org/pango/stable/PangoMarkupFormat.html */
  lblstr[512];
  niikpt p;                /* coordiante after qform transformation */

  if(g_image==NULL) return;
  switch(g_image->datatype) {
  case NIFTI_TYPE_UINT8:
  case NIFTI_TYPE_UINT16:
  case NIFTI_TYPE_UINT32:
  case NIFTI_TYPE_UINT64:
  case NIFTI_TYPE_INT8:
  case NIFTI_TYPE_INT16:
  case NIFTI_TYPE_INT32:
  case NIFTI_TYPE_INT64:
    if(g_image->dim[5]>1) {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f %3i %3i] %4.0f",q.x,q.y,q.z,g_idx[4],g_idx[5],
              g_image_data[(int)(q.x+
                                 q.y     *g_image->dim[1]+
                                 q.z     *g_image->dim[1]*g_image->dim[2]+
                                 g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]+
                                 g_idx[5]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]*g_image->dim[4])]);
    } else if(g_image->dim[4]>1) {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f %3i] %4.0f",q.x,q.y,q.z,g_idx[4],
              g_image_data[(int)(q.x+
                                 q.y     *g_image->dim[1]+
                                 q.z     *g_image->dim[1]*g_image->dim[2]+
                                 g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3])]);
    } else {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f] %4.0f",q.x,q.y,q.z,
              g_image_data[(int)(q.x+
                                 q.y*g_image->dim[1]+
                                 q.z*g_image->dim[1]*g_image->dim[2])]);
    }
    break;  /* integer image */

  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
    sprintf(lblstr,"[%3.0f %3.0f %3.0f %3i %3i] %3.0f %3.0f %3.0f",q.x,q.y,q.z,g_idx[4],g_idx[5],
            g_image_data[(int)(q.x+
                               q.y     *g_image->dim[1]+
                               q.z     *g_image->dim[1]*g_image->dim[2])],
            g_image_data[(int)(q.x+
                               q.y     *g_image->dim[1]+
                               q.z     *g_image->dim[1]*g_image->dim[2] + g_image->dim[1]*g_image->dim[2]*g_image->dim[3])],
            g_image_data[(int)(q.x+
                               q.y     *g_image->dim[1]+
                               q.z     *g_image->dim[1]*g_image->dim[2] + g_image->dim[1]*g_image->dim[2]*g_image->dim[3]*2)] );
    break;  /* color image */

  case NIFTI_TYPE_FLOAT32:
  case NIFTI_TYPE_FLOAT64:
  case NIFTI_TYPE_FLOAT128:
  case NIFTI_TYPE_COMPLEX64:
  case NIFTI_TYPE_COMPLEX128:
  case NIFTI_TYPE_COMPLEX256:
    if(g_image->dim[5]>1 && g_image->dim[0]>=5) {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f %3i %3i] %1.5f",q.x,q.y,q.z,g_idx[4],g_idx[5],
              g_image_data[(int)(q.x+
                                 q.y     *g_image->dim[1]+
                                 q.z     *g_image->dim[1]*g_image->dim[2]+
                                 g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]+
                                 g_idx[5]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]*g_image->dim[4])]);
    } else if(g_image->dim[4]>1 && g_image->dim[0]>=4) {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f %3i] %1.5f",q.x,q.y,q.z,g_idx[4],
              g_image_data[(int)(q.x+
                                 q.y     *g_image->dim[1]+
                                 q.z     *g_image->dim[1]*g_image->dim[2]+
                                 g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3])]);
    } else {
      sprintf(lblstr,"[%3.0f %3.0f %3.0f] %1.5f",q.x,q.y,q.z,
              g_image_data[(int)(q.x+
                                 q.y*g_image->dim[1]+
                                 q.z*g_image->dim[1]*g_image->dim[2])]);
    }
    break; /* floating point */
  }


  /*
   * -if qform_code is positive then show the xyz coordinates using qto_xyz
   * -we put a color on that coordinate
   * -font is courier new for fixed width
   */

  if(g_image->qform_code>0) {
    /* gtk_label_set_text(GTK_LABEL(g_lbl_voxel),markup); */
    p.x =
      g_image->qto_xyz.m[0][0] * q.x +
      g_image->qto_xyz.m[0][1] * q.y +
      g_image->qto_xyz.m[0][2] * q.z +
      g_image->qto_xyz.m[0][3];
    p.y =
      g_image->qto_xyz.m[1][0] * q.x +
      g_image->qto_xyz.m[1][1] * q.y +
      g_image->qto_xyz.m[1][2] * q.z +
      g_image->qto_xyz.m[1][3];
    p.z =
      g_image->qto_xyz.m[2][0] * q.x +
      g_image->qto_xyz.m[2][1] * q.y +
      g_image->qto_xyz.m[2][2] * q.z +
      g_image->qto_xyz.m[2][3];

    markup = g_markup_printf_escaped ("<span font_desc='8.5' font-family=\"courier new\">%s | <span foreground=\"blue\">[%6.2f %6.2f %6.2f]</span></span> ",lblstr,p.x,p.y,p.z);
    gtk_label_set_markup (GTK_LABEL (g_lbl_voxel), markup);
    g_free(markup);
  }

  else {
    /* gtk_label_set_text(GTK_LABEL(g_lbl_voxel),markup);  */
    markup = g_markup_printf_escaped ("<span font-family=\"courier new\">%s</span> ",lblstr);
    gtk_label_set_markup (GTK_LABEL (g_lbl_voxel), markup);
    g_free(markup);
  }

} /* NVK_update_lbl_voxel */



static void NVK_update_lbl_gray() {
  if(g_image==NULL) return;
  switch(g_image->datatype) {
  case NIFTI_TYPE_UINT8:
  case NIFTI_TYPE_UINT16:
  case NIFTI_TYPE_UINT32:
  case NIFTI_TYPE_UINT64:
  case NIFTI_TYPE_INT8:
  case NIFTI_TYPE_INT16:
  case NIFTI_TYPE_INT32:
  case NIFTI_TYPE_INT64:
    sprintf(g_lbl_graystr,"%1.0f",g_image_data[g_idx[1]+
            g_idx[2]*g_image->dim[1]+
            g_idx[3]*g_image->dim[1]*g_image->dim[2]+
            g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]+
            g_idx[5]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]*g_image->dim[4]]);
    break;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
    sprintf(g_lbl_graystr,"%1.0f %1.0f %1.0f",
            g_image_data[g_idx[1]+
                         g_idx[2]*g_image->dim[1]+
                         g_idx[3]*g_image->dim[1]*g_image->dim[2]],
            g_image_data[g_idx[1]+
                         g_idx[2]*g_image->dim[1]+
                         g_idx[3]*g_image->dim[1]*g_image->dim[2] +
                         g_image->dim[1]*g_image->dim[2]*g_image->dim[3]],
            g_image_data[g_idx[1]+
                         g_idx[2]*g_image->dim[1]+
                         g_idx[3]*g_image->dim[1]*g_image->dim[2] +
                         2*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]]);
    break;
  case NIFTI_TYPE_FLOAT32:
  case NIFTI_TYPE_FLOAT64:
  case NIFTI_TYPE_FLOAT128:
  case NIFTI_TYPE_COMPLEX64:
  case NIFTI_TYPE_COMPLEX128:
  case NIFTI_TYPE_COMPLEX256:
    sprintf(g_lbl_graystr,"%1.3f",g_image_data[g_idx[1]+
            g_idx[2]*g_image->dim[1]+
            g_idx[3]*g_image->dim[1]*g_image->dim[2]+
            g_idx[4]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]+
            g_idx[5]*g_image->dim[1]*g_image->dim[2]*g_image->dim[3]*g_image->dim[4]]);
    break;
  }
  gtk_label_set_text(GTK_LABEL(g_lbl_gray),g_lbl_graystr);
}





/*
 *
 * EVENT BOX RELATED TO IMAGE (event_box[?])
 *   button_press_event
 *   button_release_event
 *   motion_notify_event
 *   scroll_event
 *
 * -these events use g_xy[0-3] for previous mouse/keyboard state and x,y location
 *
 *
 */

static gboolean button_press_event_image1( GtkWidget *widget,GdkEventMotion *event ) {
  return button_press_event_imagen( widget,event,1 );
} /* button_press_event_image1 */

static gboolean button_press_event_image2( GtkWidget *widget,GdkEventMotion *event ) {
  return button_press_event_imagen( widget,event,2 );
} /* button_press_event_image2 */


static gboolean button_press_event_imagen( GtkWidget *widget,GdkEventMotion *event,int nimg ) {
  int x,y,upidx=0;
  GdkModifierType state;
  niikpt p;

  if(g_image==NULL) return TRUE;
  if (event->is_hint) gdk_window_get_pointer (event->window, &x, &y, &state);
  else {
    x=event->x;
    y=event->y;
    state=event->state;
  }

  if(verbose>1) {
    fprintf(stdout,"-d2 press_event[%i]  w=[%3i %3i]  state = %4i -> %5.1f %5.1f %5.1f\n",nimg,x,y,state,p.x,p.y,p.z);
    if(verbose>3) {
      if(state &   GDK_SHIFT_MASK) fprintf(stdout,"-d4 GDK_SHIFT_MASK\n");
      else fprintf(stdout,"-d4 no GDK_SHIFT_MASK\n");
      if(state & GDK_CONTROL_MASK) fprintf(stdout,"-d4 GDK_CONTROL_MASK\n");
      else fprintf(stdout,"-d4 no GDK_CONTROL_MASK\n");
      if(state & GDK_BUTTON1_MASK) fprintf(stdout,"-d4 GDK_BUTTON1_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON1_MASK\n");
      if(state & GDK_BUTTON2_MASK) fprintf(stdout,"-d4 GDK_BUTTON2_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON2_MASK\n");
      if(state & GDK_BUTTON3_MASK) fprintf(stdout,"-d4 GDK_BUTTON3_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON3_MASK\n");
    }
  }

  g_xy[0] = state;
  g_xy[1] = x;
  g_xy[2] = y;

  /*
  if(state & GDK_CONTROL_MASK)
    fprintf(stdout,"press   %3i %3i   %5i  *\n",x,y,state);
  else
  fprintf(stdout,"press   %3i %3i   %5i\n",x,y,state); */


  if(state & GDK_BUTTON3_MASK) {
    /* button_press_event_imagen
     * right mouse button */
    p = NVK_window2image(nimg,x,y);
    if( state & GDK_CONTROL_MASK && g_maskimg !=NULL ) {
      /* EDITING  MODE */
      g_edit_undo_idx = (g_edit_undo_idx+1)%g_edit_max_undo;
      g_edit_undo_list[g_edit_undo_idx][1] = -1;
      g_edit_undo_list[g_edit_undo_idx][0] = 2;
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(!NVK_update_mask_erosion(p.x,p.y,p.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=nimg;
    } /* EDITING MODE */
  } /* RIGHT MOUSE */

  else if(state & GDK_BUTTON1_MASK) {
    /* button_press_event_imagen
     * left mouse button */
    p = NVK_window2image(nimg,x,y);
    if( state & GDK_CONTROL_MASK && g_maskimg !=NULL ) {
      g_edit_undo_idx = (g_edit_undo_idx+1)%g_edit_max_undo;
      g_edit_undo_list[g_edit_undo_idx][1] = 1;
      g_edit_undo_list[g_edit_undo_idx][0] = 2;
      /* EDITING  MODE */
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(!NVK_update_mask_dilation(p.x,p.y,p.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=nimg;
    } /* EDITING MODE */

  } /* LEFT MOUSE */



  switch(upidx) {
  case 0:
    return TRUE;
  case 1:
    NVK_display_imagen(1);
    break;
  case 2:
    NVK_display_imagen(2);
    break;
  case 3:
    NVK_display_image ( );
    break;
  }

  return TRUE;
} /* button_press_event_imagen */






/*
 * button_release_event
 *
 * -state is the current state, so it doesn't remember about CONTROL, SHIFT etc ...
 *
 *
 *
 *   CONTROL      LEFT MOUSE              dilation using a pen
 *   CONTROL      RIGHT MOUSE             erosion using a pen
 *   CONTROL      MIDDLE MOUSE            seed fill without tolerance
 *   SHIFT        MIDDLE MOUSE            seed fill with tolerance
 *   else         RIGHT MOUSE             change zoom
 *   else         MIDDLE MOUSE            moves the viewing plane
 *   else         LEFT MOUSE              update the current position
 *
 *
 *
 *
 */



static gboolean button_release_event_image1( GtkWidget *widget,GdkEventMotion *event ) {
  return button_release_event_imagen( widget,event,1 );
} /* button_release_event_image1 */

static gboolean button_release_event_image2( GtkWidget *widget,GdkEventMotion *event ) {
  return button_release_event_imagen( widget,event,2 );
} /* button_release_event_image2 */

static gboolean button_release_event_imagen( GtkWidget *widget,GdkEventMotion *event,int nimg ) {
  int x,y,aidx,upidx=0;
  GdkModifierType state;
  niikpt p;

  if(g_image==NULL) return TRUE;
  if (event->is_hint) gdk_window_get_pointer (event->window, &x, &y, &state);
  else {
    x=event->x;
    y=event->y;
    state=event->state;
  }

  if(verbose>2) {
    fprintf(stdout,"-d3 button_release_event_image[%i]\n",nimg);
    fprintf(stdout,"-d3 nimg=%i g_axis=%c num=%i\n",nimg,g_axis[nimg],dir2num(g_axis[nimg]));
    niikmat_display(g_matrix[dir2num(g_axis[nimg])]);
    p = NVK_window2image(nimg,x,y); /* not sure if this is right */
    fprintf(stdout,"-d3 w=[%3i %3i]  state = %4i -> %5.1f %5.1f %5.1f [%3i %3i %3i]\n",
            x,y,state,
            p.x,p.y,p.z,
            g_idx[1],g_idx[2],g_idx[3]);
    if(verbose>3) {
      if(state & GDK_SHIFT_MASK)   fprintf(stdout,"-d4 GDK_SHIFT_MASK\n");
      else fprintf(stdout,"-d4 no GDK_SHIFT_MASK\n");
      if(state & GDK_CONTROL_MASK) fprintf(stdout,"-d4 GDK_CONTROL_MASK\n");
      else fprintf(stdout,"-d4 no GDK_CONTROL_MASK\n");
      if(state & GDK_BUTTON1_MASK) fprintf(stdout,"-d4 GDK_BUTTON1_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON1_MASK\n");
      if(state & GDK_BUTTON2_MASK) fprintf(stdout,"-d4 GDK_BUTTON2_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON2_MASK\n");
      if(state & GDK_BUTTON3_MASK) fprintf(stdout,"-d4 GDK_BUTTON3_MASK\n");
      else fprintf(stdout,"-d4 no GDK_BUTTON3_MASK\n");
      if(state & GDK_RELEASE_MASK) fprintf(stdout,"-d4 GDK_RELEASE_MASK\n");
      else fprintf(stdout,"-d4 no GDK_RELEASE_MASK\n");
    }
  }


  if((g_xy[0] & GDK_BUTTON1_MASK) && (g_xy[0] & GDK_CONTROL_MASK)) {
    if(verbose>2) fprintf(stdout,"-d3 release %3i %3i   %5i state = %5i dilation\n",x,y,g_xy[0],state);
    /* EDITING  MODE
     * -don't need to update the undo idx, undo count, or undo type */
    if(g_maskimg!=NULL) {
      p = NVK_window2image(nimg,x,y);
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(!NVK_update_mask_dilation(p.x,p.y,p.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=3;
    } /* EDITING MODE */
  } /* control and left mouse */

  else if((g_xy[0] & GDK_BUTTON3_MASK) && (g_xy[0] & GDK_CONTROL_MASK)) {
    if(verbose>2) fprintf(stdout,"-d3 release %3i %3i   %5i state = %5i erosion\n",x,y,g_xy[0],state);
    if(g_maskimg!=NULL) {
      /* EDITING  MODE
       * -don't need to update the undo idx, undo count, or undo type */
      p = NVK_window2image(nimg,x,y);
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(!NVK_update_mask_erosion(p.x,p.y,p.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=3;
    } /* EDITING MODE */
  }

  else if((g_xy[0] & GDK_BUTTON2_MASK) && (g_xy[0] & GDK_CONTROL_MASK)) {
    if(g_maskimg!=NULL) {
      /* EDITING  MODE
       * -don't need to update the undo idx, undo count, or undo type */
      p = NVK_window2image(nimg,x,y);
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(!NVK_update_mask_seed_fill(p.x,p.y,p.z,nimg,g_axis[nimg]))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=3;
    }
  }

  else if((g_xy[0] & GDK_BUTTON2_MASK) && (g_xy[0] & GDK_SHIFT_MASK)) {
    if(g_maskimg!=NULL) {
      /* EDITING  MODE
       * -don't need to update the undo idx, undo count, or undo type */
      p = NVK_window2image(nimg,x,y);
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      aidx = floor(p.x) + floor(p.y)*g_image->dim[1] + floor(p.z)*g_image->dim[1]*g_image->dim[2];
      if(!NVK_update_mask_seed_fill_tol(p.x,p.y,p.z,nimg,g_axis[nimg],
                                        g_image_data[aidx] - g_edit_fill_tol,
                                        g_image_data[aidx] + g_edit_fill_tol))
        upidx=0; /* do not  update if there's no change */
      else
        upidx=3;
    }
  }

  else if(g_xy[0] & GDK_BUTTON3_MASK) {
    if(verbose>2) fprintf(stdout,"-d3 release %3i %3i   %5i state = %5i zoom\n",x,y,g_xy[0],state);
    /* button_release_event_imagen
     * right mouse button */
    p = NVK_window2image(nimg,x,y);
    if(NVK_update_zoom(y,g_xy[2],nimg)) {
      g_xy[2] = y;
      upidx=3;
    }
  } /* RIGHT MOUSE */

  else if(g_xy[0] & GDK_BUTTON2_MASK) {
    if(verbose>2) fprintf(stdout,"-d3 release %3i %3i   %5i state = %5i shift\n",x,y,g_xy[0],state);
    /* button_release_event_imagen
     * middle mouse button */
    if(NVK_update_shift(x,g_xy[1],y,g_xy[2],nimg)) {
      g_xy[1] = x;
      g_xy[2] = y;
      upidx=nimg;
    }
  }


  /* update position (idx) */
  else if(g_xy[0] & GDK_BUTTON1_MASK || g_xy[0] & GDK_MOD2_MASK) {
    if(verbose>2) fprintf(stdout,"-d3 release %3i %3i   %5i state = %5i update idx\n",x,y,g_xy[0],state);
    /* button_release_event_imagen
     * left mouse button */
    p = NVK_window2image(nimg,x,y);
    g_xy[0] = 0;
    if(verbose>2) fprintf(stdout,"-d3 GDK_BUTTON1_MASK release nimg=%i | %3i %3i %3i | %3.0f %3.0f %3.0f\n",nimg,g_idx[1],g_idx[2],g_idx[3],p.x,p.y,p.z);
    if(g_axis[nimg]!='x') {
      g_idx[1] = (int)floor(NIIK_DMINMAX(p.x,0,g_image->dim[1]));
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_xdim),g_idx[1]);
    }
    if(g_axis[nimg]!='y') {
      g_idx[2] = (int)floor(NIIK_DMINMAX(p.y,0,g_image->dim[2]));
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_ydim),g_idx[2]);
    }
    if(g_axis[nimg]!='z') {
      g_idx[3] = (int)floor(NIIK_DMINMAX(p.z,0,g_image->dim[3]));
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_zdim),g_idx[3]);
    }
    if(verbose>2) fprintf(stdout,"-d3 GDK_BUTTON1_MASK release nimg=%i | %3i %3i %3i | %3.0f %3.0f %3.0f\n",nimg,g_idx[1],g_idx[2],g_idx[3],p.x,p.y,p.z);
    upidx=3;
  }


  switch(upidx) {
  case 0:
    return TRUE;
  case 1:
    NVK_display_imagen(1);
    break;
  case 2:
    NVK_display_imagen(2);
    break;
  case 3:
    NVK_display_image ( );
    break;
  }

  g_xy[0] = state;
  g_xy[1] = x;
  g_xy[2] = y;

  return TRUE;
} /* button_release_event_imagen */











/*
 * motion_notify_event
 *
 *
 * keyboard + mouse descrption:
 *
 *
 * CONTROL + RIGHT MOUSE        erodes a mask using a pen (if there's a mask)
 * SHIFT   + RIGHT MOUSE        change intensity min/max
 * CONTROL + LEFT MOUSE         dilates a mask using a pen (if there's a mask)
 * else      LEFT MOUSE         moves the current position
 * any     + MIDDLE MOUSE       shifts the viewing plane
 *
 */

static gboolean motion_notify_event_image1( GtkWidget *widget,GdkEventMotion *event ) {
  return motion_notify_event_imagen( widget,event,1 );
}

static gboolean motion_notify_event_image2( GtkWidget *widget,GdkEventMotion *event ) {
  return motion_notify_event_imagen( widget,event,2 );
}


static gboolean motion_notify_event_imagen( GtkWidget *widget,GdkEventMotion *event,int nimg ) {
  int
  ivec[8],
       x,y,upidx=0;
  GdkModifierType state;
  niikpt p,q;

  if(g_image==NULL) return TRUE;
  if (event->is_hint) gdk_window_get_pointer (event->window, &x, &y, &state);
  else {
    x=event->x;
    y=event->y;
    state=event->state;
  }

  p = NVK_window2image(nimg,x,y);

  if(verbose>2) {
    if(verbose>5)
      fprintf(stdout,"-d6 motion_notify_event[%i] w=[%3i %3i]  state = %4i -> %5.1f %5.1f %5.1f\n",nimg,x,y,state,p.x,p.y,p.z);
    if(state & GDK_SHIFT_MASK) fprintf(stdout,"-d3 shift mask\n");
    if(state & GDK_LOCK_MASK)  fprintf(stdout,"-d3 lock mask\n");
    if(state & GDK_CONTROL_MASK)  fprintf(stdout,"-d3 control mask\n");
    if(state & GDK_MOD1_MASK)  fprintf(stdout,"-d3 mod1 mask\n");
    if(state & GDK_MOD2_MASK)  fprintf(stdout,"-d3 mod2 mask\n");
    if(state & GDK_MOD3_MASK)  fprintf(stdout,"-d3 mod3 mask\n");
    if(state & GDK_MOD4_MASK)  fprintf(stdout,"-d3 mod4 mask\n");
    if(state & GDK_MOD5_MASK)  fprintf(stdout,"-d3 mod5 mask\n");
    if(state & GDK_BUTTON1_MASK)  fprintf(stdout,"-d3 butotn1 mask\n");
    if(state & GDK_BUTTON2_MASK)  fprintf(stdout,"-d3 butotn2 mask\n");
    if(state & GDK_BUTTON3_MASK)  fprintf(stdout,"-d3 butotn3 mask\n");
    if(state & GDK_BUTTON4_MASK)  fprintf(stdout,"-d3 butotn4 mask\n");
    if(state & GDK_BUTTON5_MASK)  fprintf(stdout,"-d3 butotn5 mask\n");
  }

  /* 2012-03-10 Kunio
   * -voxel info was incorrect
   * -corrected from image fov to image dimension
   * -floor is in the right place
   * 2012-03-19
   * -changed from rounding to floor
   * 2012-05-31 Kunio
   * -rounding within viewing function NVK_display_imagen was changed so that
   *  matrix is unchanged */
  q=p;
  q.x = NIIK_DMINMAX(floor(q.x),0,g_image->nx-1.0);
  q.y = NIIK_DMINMAX(floor(q.y),0,g_image->ny-1.0);
  q.z = NIIK_DMINMAX(floor(q.z),0,g_image->nz-1.0);
  if(verbose>3) {
    fprintf(stdout,"-d4 X = %9i %9i\n",x,y);
    fprintf(stdout,"-d4 p = %9.3f %9.3f %9.3f\n",p.x,p.y,p.z);
    fprintf(stdout,"-d4 q = %9.3f %9.3f %9.3f *\n",q.x,q.y,q.z);
  }


  /* UPDATE THE POSITION / INTENSITY LABEL */
  NVK_update_lbl_voxel(q);

  if(state & GDK_BUTTON3_MASK) {

    if( state & GDK_CONTROL_MASK && g_maskimg !=NULL ) {
      /* EDITING  MODE */
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      if(verbose>2) fprintf(stdout,"-d3 editing mode position = %f %f %f\n",p.x,p.y,p.z);
      if(!NVK_update_mask_erosion(q.x,q.y,q.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0;
      else
        upidx=nimg;
      g_xy[0] = state;
    } /* EDITING MODE */

    else if(state & GDK_SHIFT_MASK) {
      if(verbose>2) fprintf(stdout,"-d3 updating intensity range %i,%i -> %i,%i\n",g_xy[1],g_xy[2],x,y);
      g_xy[0]=0;
      if(NVK_update_intensity(x,g_xy[1],y,g_xy[2],nimg)) {
        g_xy[1] = x;
        g_xy[2] = y;
        upidx=3;
      }
    }

    else {
      if(verbose>2) fprintf(stdout,"-d3 updating zoom %i -> %i\n",g_xy[2],y);
      if(NVK_update_zoom(y,g_xy[2],nimg)) {
        g_xy[2] = y;
        upidx=nimg;
      }
    }
  }  /* RIGHT MOUSE */


  else if(state & GDK_BUTTON1_MASK ) {

    if( state & GDK_CONTROL_MASK && g_maskimg !=NULL ) {
      /* EDITING  MODE */
      if(p.x < 0) return TRUE;
      if(p.y < 0) return TRUE;
      if(p.z < 0) return TRUE;
      if(p.x >= g_image->nx-1) return TRUE;
      if(p.y >= g_image->ny-1) return TRUE;
      if(p.z >= g_image->nz-1) return TRUE;
      /* fprintf(stdout,"  %f %f %f\n",p.x,p.y,p.z); */
      if(!NVK_update_mask_dilation(q.x,q.y,q.z,nimg,g_axis[nimg],g_pen_size))
        upidx=0;
      else
        upidx=nimg;
    } /* EDITING MODE */

    else { /* LEFT MOUSE but not CONTROL */

      if(verbose>2) fprintf(stdout,"-d3 GDK_BUTTON1_MASK motion_notify nimg=%i %c | %3i %3i %3i | %3.0f %3.0f %3.0f\n",
                              nimg,g_axis[nimg],g_idx[1],g_idx[2],g_idx[3],q.x,q.y,q.z);
      ivec[1] = (int)q.x;
      ivec[2] = (int)q.y;
      ivec[3] = (int)q.z;
      ivec[1] = NIIK_IMINMAX(ivec[1],0,g_image->dim[1]-1);
      ivec[2] = NIIK_IMINMAX(ivec[2],0,g_image->dim[2]-1);
      ivec[3] = NIIK_IMINMAX(ivec[3],0,g_image->dim[3]-1);
      if(g_axis[nimg]=='x') {
        if(ivec[2]==g_idx[2]) {
          if(ivec[3]==g_idx[3]) {
            return TRUE;
          }
        }
        g_idx[2] = ivec[2];
        g_idx[3] = ivec[3];
        g_xy[0] = 0;
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_ydim),g_idx[2]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_zdim),g_idx[3]);
      }
      if(g_axis[nimg]=='y') {
        if(ivec[1]==g_idx[1]) {
          if(ivec[3]==g_idx[3]) {
            return TRUE;
          }
        }
        g_idx[1] = ivec[1];
        g_idx[3] = ivec[3];
        g_xy[0] = 0;
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_xdim),g_idx[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_zdim),g_idx[3]);
      }
      if(g_axis[nimg]=='z') {
        if(ivec[1]==g_idx[1]) {
          if(ivec[2]==g_idx[2]) {
            return TRUE;
          }
        }
        g_idx[1] = ivec[1];
        g_idx[2] = ivec[2];
        g_xy[0] = 0;
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_xdim),g_idx[1]);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_ydim),g_idx[2]);
      }

      upidx=3;
      g_xy[0] = state;
    } /* NOT EDITING MODE */
  } /* LEFT MOUSE */


  else if(state & GDK_BUTTON2_MASK) {
    if(verbose>2) fprintf(stdout,"-d3 updating shift %i,%i -> %i,%i\n",g_xy[1],g_xy[2],x,y);
    if(NVK_update_shift(x,g_xy[1],y,g_xy[2],nimg)) {
      g_xy[1] = x;
      g_xy[2] = y;
      upidx=nimg;
    }
  } /* button2 = middle mouse */

  g_xy[0] = state;


  switch(upidx) {
  case 0:
    return TRUE;
  case 1:
    NVK_display_imagen(1);
    break;
  case 2:
    NVK_display_imagen(2);
    break;
  case 3:
    NVK_display_image ( );
    break;
  }


  return TRUE;
} /* motion_nofity_event_imagen */



/*
 * scroll_event
 *
 * -moves through the plane
 *
 */

static gboolean scroll_event_image1( GtkWidget *widget,GdkEventScroll *event ) {
  return scroll_event_imagen( widget,event,1 );
}

static gboolean scroll_event_image2( GtkWidget *widget,GdkEventScroll *event ) {
  return scroll_event_imagen( widget,event,2 );
}

static gboolean scroll_event_imagen( GtkWidget *widget,GdkEventScroll *event,int nimg ) {
  int x,y,n,txy;
  GdkModifierType state;
  niikpt p;

  if(g_image==NULL) return TRUE;

  x=event->x;
  y=event->y;
  state=event->state;

  p = NVK_window2image(nimg,x,y);
  if(verbose>2) fprintf(stdout,"-d3 w=[%3i %3i]  state = %4i -> %5.1f %5.1f %5.1f\n",x,y,state,p.x,p.y,p.z);

  /* UP/DOWN
   * -change the current location according to the image */
  n=dir2num(g_axis[nimg]); /* get index for direction */

  switch(event->direction) {
  case GDK_SCROLL_UP:
    if(g_idx[n]==0) return TRUE;
    g_idx[n]--;
    /* 2012-05-07, Kunio
     * -added temporary xy for changing g_xy[0] to non-negative */
    txy=g_xy[0];
    g_xy[0]=1;
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_dim[n]),g_idx[n]);
    g_xy[0]=txy;
    break;
  case GDK_SCROLL_DOWN:
    if(g_idx[n]==g_image->dim[n]-1) return TRUE;
    g_idx[n]++;
    /* 2012-05-07, Kunio
     * -added temporary xy for changing g_xy[0] to non-negative */
    txy=g_xy[0];
    g_xy[0]=1;
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_dim[n]),g_idx[n]);
    g_xy[0]=txy;
    break;
  case GDK_SCROLL_LEFT:
    break;
  case GDK_SCROLL_RIGHT:
    break;
  }

  return TRUE;
} /* scroll_event_image1 */





/*
 * updates the matrix
 *
 *
 * -there are 3 matrices for x,y,z directions
 * -display function assumes the image directions
 *  so we can't rotate the image
 * -we can flip though
 * -the matrix goes from the window coordinate (+plane location)
 *  to the image space
 * -we can invert and get the window coordinate from image voxel
 *  location
 *
 */

static void NVK_update_matrix(int nimg) {
  int num;
  niikmat *A=NULL;

  if(g_image==NULL) return;

  switch(g_axis[nimg]) {
  case 'z':
    num = 3;

    verbose-=2;
    if(verbose>1) fprintf(stdout,"-d2 dir=%i, axis[%i]=%c\n",num,nimg,g_axis[nimg]);
    niikmat_identity_update(g_matrix[num]);                                /* identity matrix */
    if(verbose>5) {
      fprintf(stdout,"  identity\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix(g_img_delta.x,g_img_delta.y,g_img_delta.z),
                                g_matrix[num]);                                    /* pixel spacing */
    if(verbose>5) {
      fprintf(stdout,"  scale\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(-g_img_fov.x/2.0,-g_img_fov.y/2.0,-g_idx[3]*g_img_delta.z),
                                g_matrix[num]);                                    /* image centering */
    if(verbose>5) {
      fprintf(stdout,"  centering 1\n");
      niikmat_display(g_matrix[num]);
    }

    if(g_flip[1])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('x'),
                                  g_matrix[num] );                                  /* flip in x */
    if(g_flip[2])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('y'),
                                  g_matrix[num] );                                   /* flip in y */
    if(g_flip[3])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('z'),
                                  g_matrix[num] );                                   /* flip in z */
    if(verbose>5) {
      fprintf(stdout,"  flipping\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(g_shift[num][1],g_shift[num][2],0),
                                g_matrix[num]);                                     /* shift */
    if(verbose>5) {
      fprintf(stdout,"  shift\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix(g_zoom[num],g_zoom[num],1),
                                g_matrix[num]);                                     /* zooming */
    if(verbose>5) {
      fprintf(stdout,"  zooming\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(0.5*g_window_width,
                                0.5*g_window_height,
                                0),
                                g_matrix[num]);                                   /* window centering */
    if(verbose>5) {
      fprintf(stdout,"  window centering\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_inverse_update(g_matrix[num]);                                 /* inverse */
    if(verbose>5) {
      fprintf(stdout,"  inversion\n");
      niikmat_display(g_matrix[num]);
    }

    verbose+=2;


    if(verbose>2) {
      fprintf(stdout,"-d3 z-dir matrix\n");
      niikmat_display(g_matrix[num]);
    }
    break;

  case 'y':
    num = 2;
    if(verbose>1) fprintf(stdout,"-d2 dir=%i, axis[%i]=%c\n",num,nimg,g_axis[nimg]);
    niikmat_identity_update(g_matrix[num]);                                /* identity matrix */
    A = niikmat_init(4,4);
    niikmat_identity_update(A);
    A->m[1][2]=A->m[2][1]=1;
    A->m[1][1]=A->m[2][2]=0;
    if(verbose>5) {
      fprintf(stdout,"  A\n");
      niikmat_display(A);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix    (g_img_delta.x,
                                g_img_delta.y,
                                g_img_delta.z),
                                g_matrix[num]);                                   /* pixel spacing */
    if(verbose>5) {
      fprintf(stdout,"  pixel spacing scaling\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(-g_img_fov.x/2.0,
                                -g_idx[2]*g_img_delta.y,
                                -g_img_fov.z/2.0),
                                g_matrix[num]);                                   /* image centering */
    if(verbose>5) {
      fprintf(stdout,"  image centering\n");
      niikmat_display(g_matrix[num]);
    }

    if(g_flip[1])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('x'),
                                  g_matrix[num] );                                /* flip in x */
    if(g_flip[2])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('y'),
                                  g_matrix[num] );                                /* flip in y */
    if(g_flip[3])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('z'),
                                  g_matrix[num] );                                /* flip in z */
    if(verbose>5) {
      fprintf(stdout,"  flipping\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1( A, g_matrix[num]);                               /* image orientation */
    if(verbose>5) {
      fprintf(stdout,"  image orientation\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(g_shift[num][1],
                                g_shift[num][2],
                                0.0),
                                g_matrix[num]);                                   /* window centering */
    if(verbose>5) {
      fprintf(stdout,"  centeringf 2\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix(g_zoom[num],
                                g_zoom[num],
                                1.0 ),
                                g_matrix[num]);                                   /* zooming */
    if(verbose>5) {
      fprintf(stdout,"  zooming\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(0.5*g_window_width,
                                0.5*g_window_height,
                                0.0),
                                g_matrix[num]);                                   /* window centering */
    if(verbose>5) {
      fprintf(stdout,"  window centering\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_inverse_update(g_matrix[num]);                              /* inverse */
    if(verbose>5) {
      fprintf(stdout,"  inverse\n");
      niikmat_display(g_matrix[num]);
    }

    if(verbose>2) {
      fprintf(stdout,"-d3 y-dir matrix\n");
      niikmat_display(g_matrix[num]);
    }
    break;

  case 'x':
    num=1;
    if(verbose>1) fprintf(stdout,"-d2 dir=%i, axis[%i]=%c\n",num,nimg,g_axis[nimg]);
    niikmat_identity_update(g_matrix[num]);                                /* identity matrix */
    A = niikmat_init(4,4);
    niikmat_clear(A);
    A->m[0][1]=A->m[1][2]=A->m[2][0]=A->m[3][3]=1.0;
    if(verbose>5) {
      fprintf(stdout,"  A\n");
      niikmat_display(A);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix    (g_img_delta.x,
                                g_img_delta.y,
                                g_img_delta.z),
                                g_matrix[num]);                                   /* pixel spacing */
    if(verbose>5) {
      fprintf(stdout,"  scaling (pixel spacing)\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(-g_idx[1]*g_img_delta.x,
                                -g_img_fov.y/2.0,
                                -g_img_fov.z/2.0),
                                g_matrix[num]);                                   /* image centering */
    if(verbose>5) {
      fprintf(stdout,"  centering 1\n");
      niikmat_display(g_matrix[num]);
    }

    if(g_flip[1])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('x'),
                                  g_matrix[num] );                                /* flip in x */
    if(g_flip[2])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('y'),
                                  g_matrix[num] );                                /* flip in y */
    if(g_flip[3])
      niikmat_multiply_mat2_free1(niikmat_flip_matrix('z'),
                                  g_matrix[num] );                                /* flip in z */
    if(verbose>5) {
      fprintf(stdout,"  flip\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1( A, g_matrix[num]);                               /* image orientation */
    if(verbose>5) {
      fprintf(stdout,"  image orientation\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(g_shift[num][1],
                                g_shift[num][2],
                                0.0),
                                g_matrix[num]);                                   /* translation */
    if(verbose>5) {
      fprintf(stdout,"  translation\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_scale_matrix(g_zoom[num],
                                g_zoom[num],
                                1.0 ),
                                g_matrix[num]);                                   /* zooming */
    if(verbose>5) {
      fprintf(stdout,"  zooming\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_multiply_mat2_free1(niikmat_translate_matrix(0.5*g_window_width,
                                0.5*g_window_height,
                                0.0),
                                g_matrix[num]);                                   /* window centering */
    if(verbose>5) {
      fprintf(stdout,"  window centering\n");
      niikmat_display(g_matrix[num]);
    }

    niikmat_inverse_update(g_matrix[num]);                              /* inverse */
    if(verbose>5) {
      fprintf(stdout,"  inversion\n");
      niikmat_display(g_matrix[num]);
    }

    if(verbose>2) {
      fprintf(stdout,"-d3 x-dir matrix\n");
      niikmat_display(g_matrix[num]);
    }
    break;

  default:
    fprintf(stderr,"ERROR: unknown nimg for NVK_window2image, %i\n",nimg);
    exit(0);
  }

  return;
}


/*
 * calculate the positions
 *
 * -assumes g_matrix is already updated
 *
 */

static niikpt NVK_window2image(int nimg,int xi,int yi) {
  niikpt p;
  niikmat *A;
  A = g_matrix[dir2num(g_axis[nimg])];
  /* 2012-05-31, Kunio
   * -added +0.5 for rounding */
  p.x =
    A->m[0][0] * xi +
    A->m[0][1] * yi +
    A->m[0][3] + 0.5 ;
  p.y =
    A->m[1][0] * xi +
    A->m[1][1] * yi +
    A->m[1][3] + 0.5;
  p.z =
    A->m[2][0]   * xi +
    A->m[2][1]   * yi +
    A->m[2][3] + 0.5;
  p.w =0 ;
  return p;
}


/************************************
 * image display function
 *
 *
 *   int NVK_display_image();
 *   int NVK_display_imagen(int nimg);
 *
 * -NVK_display_imagen is the main display function
 * -updates either left/right image (nimg=1/2 respectively)
 * -object is allowed and displayed on the screen
 *
 * To Do
 * -multiple objects would be nice (for cortex)
 * -faster display (?)
 *
 *
 * Bug
 * -blue mask tends to darken the image
 * -Can't rotate the image (arbitrary or 90 degrees)
 *
 *
 ************************************/

static int NVK_display_image() {
  /* int n;
    #pragma omp parallel for
  for(n=1;n<=2;n++){
  NVK_display_imagen(n); }*/
  NVK_display_imagen(1);
  NVK_display_imagen(2);
  return 1;
}

static int NVK_display_imagen(int nimg) {
  kface *f;
  int
      xo,yo,zo,
      xi,yi,vi,
      xp,yp,
      vizo,viyo,vixo,
      m,n,mm=0,nn=0,
      coplanar,
      i,j,k,
      xyzdim=0,
      isize1,
      isize2,
      isect[3][2],
      bxmin,bxmax,bymin,bymax,bzmin,bzmax,
      rowstride,rowstride2,
      xdim,area,size;
  double
    isectfrac,
    spinmin,spinmax,
    iscl,ioff;
  guchar *pixels,pix;
  niikpt
    pp[4],                /* plane triangle for g_obj */
    qq[2];                /* intersection points */
  niikmat *A, *invA, *Ao;
  char ucblue255;
  /*static int dispindex=0;*/


  if(g_image==NULL) return 1;

  pp[2].w=0;

  /*  dispindex++;
      if(verbose) fprintf(stdout,"-d1 display nimg=%i cnt=%i\n",nimg,dispindex);*/

  pixels     = gdk_pixbuf_get_pixels (g_pixbuf[nimg]);
  rowstride2 = rowstride  = gdk_pixbuf_get_rowstride  (g_pixbuf[nimg]);
  rowstride -= 3*g_window_width;       /* actual offset */
  size = rowstride2 * g_window_height; /* # bytes in the pixel buffer */

  /* intensity scaling calculation */
  spinmin = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_min));
  spinmax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_max));
  iscl = 255.0 / (spinmax - spinmin);
  ioff = -iscl * spinmin + 0.5; /* 0.5 for rounding */
  /* fprintf(stdout,"  iscl = %f  ioff = %f\n",iscl,ioff); */

  xdim = g_image->dim[1];
  area = g_image->dim[1] * g_image->dim[2];
  xyzdim=g_image->nx * g_image->ny * g_image->nz;
  isize1=g_image->nx * g_image->ny * g_image->nz;
  if(g_image->ndim>3) isize1*=g_image->nt;
  if(g_image->ndim>4) isize1*=g_image->nu;
  isize2=isize1 * 2;

  NVK_update_matrix(nimg);

  A = g_matrix[dir2num(g_axis[nimg])];

  invA = niikmat_inverse(A);

  for(n=0; n<size; n++) pixels[n]=0; /* black background */

  ucblue255 = 255 * g_edit_blue;


  switch(g_axis[nimg]) {


  case 'x':
    vizo  = g_idx[1];
    n     = area * g_image->dim[3];
    if(g_image->ndim>3) {
      vizo += n    * g_idx[4];
      n    *= g_image->dim[4];
    }
    if(g_image->ndim>4) {
      vizo += n    * g_idx[5];
    }

    if(g_image->datatype == NIFTI_TYPE_RGB24) {
      vizo  = g_idx[1];
      n = area * g_image->dim[3];
      nn = n*2;
    }

    for(yi=0,vi=0,yp=-1; yi<g_window_height; yi++) { /* sup/inf dir */
      /* 2012-05-31, Kunio
       * -+ 0.5 was added because rounding was removed from the matrix
       * -similarly for the rest of loops in this function */
      yo = (int) floor(A->m[2][1]*yi + A->m[2][3] + 0.5);
      if(yo<0) {
        vi+=rowstride2;
        continue;
      }
      if(yo>=g_image->dim[3]) {
        vi+=rowstride2;
        continue;
      }
      if(yo == yp) {
        memcpy(pixels+vi,pixels+vi-rowstride2,(size_t)rowstride2);
        vi += rowstride2;
        continue;
      }
      yp = yo;
      viyo = vizo + yo * area;
      for(xi=0,xp=-1; xi<g_window_width; xi++) { /* left-right */
        xo = floor(A->m[1][0]*xi + A->m[1][1]*yi + A->m[1][3] + 0.5);
        if(xo<0) {
          vi+=3;
          continue;
        }
        if(xo>=g_image->dim[2]) {
          vi+=3;
          continue;
        }
        if(xo==xp) {
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          continue;
        }

        /* COLOR IMAGE */
        if(g_image->datatype == NIFTI_TYPE_RGB24) {
          pixels[vi++] = (unsigned char)g_image_data[xo*xdim+viyo];
          pixels[vi++] = (unsigned char)g_image_data[xo*xdim+viyo+n];
          pixels[vi++] = (unsigned char)g_image_data[xo*xdim+viyo+nn];
        } else if(g_maskimg!=NULL) { /* mask image */
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo],0,255);
          if(g_image->nv==3) {
            if(g_mask_data[xo*xdim+viyo]>0) {
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize1],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize2],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
            } else {
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize1],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize2],0,255);
            }
          } else if(g_mask_data[(xo*xdim+viyo)%xyzdim]>0) { /* color image with a mask */
            pix *= 1.0-g_edit_blue;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
          } else {
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
          }
        } /* with mask */

        else if(g_image->nv==3) { /* COLOR IMAGE (Kunio's format) */
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize1],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo+isize2],0,255);
        }

        /* GRAY IMAGE */
        else {
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo*xdim+viyo],0,255);
          pixels[vi++] = pix;
          pixels[vi++] = pix;
          pixels[vi++] = pix;
        }
        xp = xo;
      } /* x-dir */
      vi += rowstride;
    } /* y-dir */


    /* wrote iamge for x-dir */





    /*
     *  x-direction object display
     *
     *
     * -if an object was an input, then
     *  we need to put it on the image
     * -this section was removed
     */



    if(g_obj_list_num) {
      /* get the points on the viewing plane */
      pp[0].x = pp[1].x = pp[2].x = g_idx[1];
      pp[0].y = pp[0].z = pp[1].z = pp[2].y = 0;
      pp[1].y = pp[2].z = 1e6;

      if(verbose>2) {
        fprintf(stdout,"-d3 plane triangle: %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f\n",
                pp[0].x,pp[0].y,pp[0].z,
                pp[1].x,pp[1].y,pp[1].z,
                pp[2].x,pp[2].y,pp[2].z);
        off_display_bbox_info(g_off_bb);
      }

      pp[0]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[0]);
      pp[1]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[1]);
      pp[2]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[2]);


      /* pre-calculate the matrix from image (xyz) space to window space */
      Ao = niikmat_inverse(g_matrix[dir2num(g_axis[nimg])]);
      niikmat_multiply_mat1_free2(Ao,niikmat_mat44_matrix(g_image->sto_ijk));


      /* go thru the bounding boxes -only in 2 dimensions */
      bxmin = floor( (pp[0].x - g_off_bb->origin.x) / g_off_bb->delta - 0.5);
      bxmax = floor( (pp[0].x - g_off_bb->origin.x) / g_off_bb->delta + 0.5);

      bxmin = NIIK_IMINMAX(bxmin-1,0,g_off_bb->xdim-1);
      bxmax = NIIK_IMINMAX(bxmax+1,0,g_off_bb->xdim-1);

      for(k=0; k<g_off_bb->zdim; k++) {
        nn = k * g_off_bb->area;

        for(j=0; j<g_off_bb->ydim; j++) {

          for(i=bxmin; i<=bxmax; i++) {
            n = i + j * g_off_bb->xdim + nn;

            if(verbose>7)
              fprintf(stdout,"-d8 bbox index=%i[%i,%i,%i]\n",n,i,j,k);

            /* within each bounding box, test if the surface intersect with the viewing plane */
            for(m=0; m<g_off_bb->ndata[n]; m++) {
              f = g_off_bb->data[n][m];

              if(f->pmax.x < pp[0].x) continue; /* skip if no overlap */
              if(f->pmin.x > pp[0].x) continue;

              if(!off_check_tri_tri_intersect_with_isectline(pp[0],pp[1],pp[2],
                  f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
                  &coplanar,
                  &qq[0],
                  &qq[1]) ) {
                continue;
              }

              /* if(coplanar) { continue; } */
              if(verbose>7) fprintf(stdout,"-d8 intersect: %0.2f %0.2f %0.2f -> %0.2f %0.2f %0.2f\n",
                                      qq[0].x,qq[0].y,qq[0].z,qq[1].x,qq[1].y,qq[1].z);

              /* find the window positions */
              isect[0][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[0].x +
                                  Ao->m[0][1]  * qq[0].y +
                                  Ao->m[0][2]  * qq[0].z +
                                  Ao->m[0][3]);

              isect[0][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[0].x +
                                  Ao->m[1][1]  * qq[0].y +
                                  Ao->m[1][2]  * qq[0].z +
                                  Ao->m[1][3]);

              isect[1][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[1].x +
                                  Ao->m[0][1]  * qq[1].y +
                                  Ao->m[0][2]  * qq[1].z +
                                  Ao->m[0][3]);

              isect[1][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[1].x +
                                  Ao->m[1][1]  * qq[1].y +
                                  Ao->m[1][2]  * qq[1].z +
                                  Ao->m[1][3]);

              /* linearly interpolate between the two pointer */
              /*niikmat_display(g_matrix[dir2num(g_axis[nimg])]);
              fprintf(stderr,"ERROR: %f\n",fabs(g_matrix[dir2num(g_axis[nimg])]->m[1][0])/5.0);*/

              for(isectfrac=0;
                  isectfrac<=1.0;
                  isectfrac+=fabs(g_matrix[dir2num(g_axis[nimg])]->m[1][0])/5.0) {
                isect[2][0] = isect[1][0] * isectfrac + isect[0][0] * (1.0-isectfrac);
                isect[2][1] = isect[1][1] * isectfrac + isect[0][1] * (1.0-isectfrac);

                if(isect[2][0]<0) continue;
                if(isect[2][1]<0) continue;
                if(isect[2][0]>=g_window_width)  continue;
                if(isect[2][1]>=g_window_height) continue;

                mm = isect[2][1]*rowstride2 + isect[2][0]*3;
                if(f->color!=NULL) {
                  pixels[mm  ] = 255 * f->color[0];
                  pixels[mm+1] = 255 * f->color[1];
                  pixels[mm+2] = 255 * f->color[2];
                  /* fprintf(stdout,"color %3i %3i %3i\n",pixels[mm],pixels[mm+1],pixels[mm+2]); */
                } else { /* no color -> green */
                  pixels[mm  ] = 0;
                  pixels[mm+1] = 255;
                  pixels[mm+2] = 0;
                }
              }
            } /* each triangle in the bounding box */

          }
        }
      } /* ijk of the bounding boxes */

      Ao=niikmat_free(Ao);
    } /* object for x-dir */



    break;




  case 'y':
    vizo  = g_idx[2] * xdim;
    if(g_image->ndim>3) {
      n     = area * g_image->dim[3];
      vizo += n    * g_idx[4];
    }
    if(g_image->ndim>4) {
      n    *= g_image->dim[4];
      vizo += n    * g_idx[5];
    }

    if(g_image->datatype == NIFTI_TYPE_RGB24) {
      vizo  =  g_idx[2] * xdim;
      n     = area * g_image->dim[3];
      nn    = n*2;
    }


    for(yi=0,vi=0,yp=-1; yi<g_window_height; yi++) { /* sup/inf dir */
      yo = (int) floor(A->m[2][1]*yi + A->m[2][3] + 0.5);
      if(yo<0) {
        vi+=rowstride2;
        continue;
      }
      if(yo>=g_image->dim[3]) {
        vi+=rowstride2;
        continue;
      }
      /*if(yo == yp){
      memcpy(pixels+vi,pixels+vi-rowstride2,(size_t)rowstride2);
      vi += rowstride2;
      continue; }*/
      yp = yo;
      viyo = vizo + yo * area;

      for(xi=0,xp=-1; xi<g_window_width; xi++) { /* left/right dir */
        xo = floor(A->m[0][0]*xi + A->m[0][1]*yi + A->m[0][3] + 0.5);
        if(xo<0) {
          vi+=3;
          continue;
        }
        if(xo>=g_image->dim[1]) {
          vi+=3;
          continue;
        }
        if(xo==xp) {
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          continue;
        }
        if(g_image->datatype == NIFTI_TYPE_RGB24) {
          vixo = xo+viyo;
          pixels[vi++] = (unsigned char)g_image_data[vixo];
          pixels[vi++] = (unsigned char)g_image_data[vixo+n];
          pixels[vi++] = (unsigned char)g_image_data[vixo+nn];
        } else if(g_maskimg!=NULL) {
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          if(g_image->nv==3) {
            if(g_mask_data[xo+viyo]>0) {
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
            } else {
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
            }
          } /* color image with a mask */

          else if(g_mask_data[(xo+viyo)%xyzdim]>0) {
            pix *= 1.0-g_edit_blue;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
          } else {
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
          }
        }

        /* COLOR IMAGE (Kunio's format) */
        else if(g_image->nv==3) {
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
        }

        /* GRAY IMAGE */
        else {
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          pixels[vi++] = pix;
          pixels[vi++] = pix;
          pixels[vi++] = pix;
        }
        xp = xo;
      } /* x-dir */
      vi += rowstride;
    } /* y-dir */


    /* wrote image for y-dir */






    /*
     *  y-direction object display
     *
     *
     *  -if an object was an input, then
     *  we need to put it on the image
     * -this section was removed
     */

    if(g_obj_list_num) {
      /* get the points on the viewing plane */
      pp[0].x = 0;
      pp[0].y = g_idx[2];
      pp[0].z = 0;

      pp[1].x = 1e6;
      pp[1].y = g_idx[2];
      pp[1].z = 0;

      pp[2].x = 0;
      pp[2].y = g_idx[2];
      pp[2].z = 1e6;

      pp[0]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[0]);
      pp[1]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[1]);
      pp[2]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[2]);

      if(verbose>2) {
        fprintf(stdout,"-d3 plane triangle: %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f\n",
                pp[0].x,pp[0].y,pp[0].z,
                pp[1].x,pp[1].y,pp[1].z,
                pp[2].x,pp[2].y,pp[2].z);
        off_display_bbox_info(g_off_bb);
      }


      /* pre-calculate the matrix from image (xyz) space to window space */
      Ao = niikmat_inverse(g_matrix[dir2num(g_axis[nimg])]);

      /*
      niikmat_multiply_mat1_free2(Ao,
                                  niikmat_scale_matrix(1.0/g_image->pixdim[1],
                                                       1.0/g_image->pixdim[2],
                                                       1.0/g_image->pixdim[3]));
      */
      niikmat_multiply_mat1_free2(Ao,niikmat_mat44_matrix(g_image->sto_ijk));

      /* go thru the bounding boxes -only in 2 dimensions */
      bymin = floor( (pp[0].y - g_off_bb->origin.y) / g_off_bb->delta - 0.5);
      bymax = floor( (pp[0].y - g_off_bb->origin.y) / g_off_bb->delta + 0.5);

      bymin = NIIK_IMINMAX(bymin-1, 0, g_off_bb->ydim-1);
      bymax = NIIK_IMINMAX(bymax+1, 0, g_off_bb->ydim-1);

      for(k=0; k<g_off_bb->zdim; k++) {
        nn = k * g_off_bb->area;

        for(j=bymin; j<=bymax; j++) {

          for(i=0; i<g_off_bb->xdim; i++) {
            n = i + j * g_off_bb->xdim + nn;

            if(verbose>7)
              fprintf(stdout,"-d8 bbox index=%i[%i,%i,%i]\n",n,i,j,k);

            /* within each bounding box, test if the surface intersect with the viewing plane */
            for(m=0; m<g_off_bb->ndata[n]; m++) {
              f = g_off_bb->data[n][m];

              if(f->pmax.y < pp[0].y) continue; /* skip if no overlap */
              if(f->pmin.y > pp[0].y) continue;

              if(!off_check_tri_tri_intersect_with_isectline(pp[0],pp[1],pp[2],
                  f->vert[0]->v,
                  f->vert[1]->v,
                  f->vert[2]->v,
                  &coplanar,
                  &qq[0],&qq[1]) ) {
                continue;
              }

              /* if(coplanar) { continue; } */
              if(verbose>7) fprintf(stdout,"-d8 intersect: %0.2f %0.2f %0.2f -> %0.2f %0.2f %0.2f\n",
                                      qq[0].x,qq[0].y,qq[0].z,qq[1].x,qq[1].y,qq[1].z);

              /* find the window positions */
              isect[0][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[0].x +
                                  Ao->m[0][1]  * qq[0].y +
                                  Ao->m[0][2]  * qq[0].z +
                                  Ao->m[0][3]);
              isect[0][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[0].x +
                                  Ao->m[1][1]  * qq[0].y +
                                  Ao->m[1][2]  * qq[0].z +
                                  Ao->m[1][3]);

              isect[1][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[1].x +
                                  Ao->m[0][1]  * qq[1].y +
                                  Ao->m[0][2]  * qq[1].z +
                                  Ao->m[0][3]);
              isect[1][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[1].x +
                                  Ao->m[1][1]  * qq[1].y +
                                  Ao->m[1][2]  * qq[1].z +
                                  Ao->m[1][3]);

              /* linearly interpolate between the two pointer */
              for(isectfrac=0; isectfrac<=1.0; isectfrac+=fabs(g_matrix[dir2num(g_axis[nimg])]->m[0][0])/5.0) {
                isect[2][0] = isect[1][0] * isectfrac + isect[0][0] * (1.0-isectfrac);
                isect[2][1] = isect[1][1] * isectfrac + isect[0][1] * (1.0-isectfrac);

                if(isect[2][0]<0) continue;
                if(isect[2][1]<0) continue;
                if(isect[2][0]>=g_window_width)  continue;
                if(isect[2][1]>=g_window_height) continue;

                mm = isect[2][1]*rowstride2 + isect[2][0]*3;
                if(f->color!=NULL) { /* colored object */
                  pixels[mm  ] = 255 * f->color[0];
                  pixels[mm+1] = 255 * f->color[1];
                  pixels[mm+2] = 255 * f->color[2];
                  /* fprintf(stdout,"color %3i %3i %3i\n",pixels[mm],pixels[mm+1],pixels[mm+2]); */
                } else { /* no color -> green */
                  pixels[mm  ] = 0;
                  pixels[mm+1] = 255;
                  pixels[mm+2] = 0;
                }
              }
            } /* each triangle in the bounding box */

          }
        }
      } /* ijk of the bounding boxes */

      Ao=niikmat_free(Ao);
    } /* object for y-dir */



    break;









  case 'z':
    vizo  = g_idx[3] * area;
    if(g_image->ndim>3) {
      n     = area * g_image->dim[3];
      vizo += n    * g_idx[4];
    }
    if(g_image->ndim>4) {
      n    *= g_image->dim[4];
      vizo += n    * g_idx[5];
    }
    /*fprintf(stdout,"%3i %3i %3i %3i %3i | %10i | %10i\n",g_idx[1],g_idx[2],g_idx[3],g_idx[4],g_idx[5],vizo,g_image->nvox); */

    if(g_image->datatype == NIFTI_TYPE_RGB24) {
      vizo  = g_idx[3] * area;
      n     = area * g_image->dim[3];
      nn    = n * 2;
    }

    for(yi=0,vi=0,yp=-1; yi<g_window_height; yi++) {
      yo = (int) floor(A->m[1][1]*yi + A->m[1][3] + 0.5);
      if(yo<0) {
        vi+=rowstride2;
        continue;
      }
      if(yo>=g_image->dim[2]) {
        vi+=rowstride2;
        continue;
      }
      if(yo == yp) {
        memcpy(pixels+vi,pixels+vi-rowstride2,(size_t)rowstride2);
        vi += rowstride2;
        continue;
      }
      yp = yo;
      viyo = vizo + yo * xdim;

      for(xi=0,xp=-1; xi<g_window_width; xi++) {
        xo = floor(A->m[0][0]*xi + A->m[0][1]*yi + A->m[0][3] + 0.5);
        if(xo==xp) {
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          pixels[vi] = pixels[vi-3];
          vi++;
          continue;
        }
        if(xo<0) {
          vi+=3;
          continue;
        }
        if(xo>=g_image->dim[1]) {
          vi+=3;
          continue;
        }
        if(g_image->datatype == NIFTI_TYPE_RGB24) {
          vixo = xo+viyo;
          pixels[vi++] = (unsigned char)g_image_data[vixo];
          pixels[vi++] = (unsigned char)g_image_data[vixo+n];
          pixels[vi++] = (unsigned char)g_image_data[vixo+nn];
        }

        else if(g_maskimg!=NULL) {
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          if(g_image->nv==3) {
            if(g_mask_data[xo+viyo]>0) {
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = pix;
              pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
              pix *= 1.0-g_edit_blue;
              pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
            } else {
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
              pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
            }
          } /* color image with a mask */

          else if(g_mask_data[(xo+viyo)%xyzdim]>0) {
            pix *= 1.0-g_edit_blue;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = NIIK_UCMINMAX(pix + ucblue255,0,255);
          } else {
            pixels[vi++] = pix;
            pixels[vi++] = pix;
            pixels[vi++] = pix;
          }
        }


        /* COLOR IMAGE (Kunio's format) */
        else if(g_image->nv==3) {
          /*fprintf(stdout,"    %i  color\n",g_image->nv);*/
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize1],0,255);
          pixels[vi++] = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo+isize2],0,255);
        }

        /* GRAY IMAGE */
        else {
          pix = (unsigned char)NIIK_DMINMAX(ioff + iscl * g_image_data[xo+viyo],0,255);
          /* fprintf(stdout,"%4i %4i %f %i\n",xi,yi,g_image_data[xo+viyo],pix); */
          pixels[vi++] = pix;
          pixels[vi++] = pix;
          pixels[vi++] = pix;
        }
        xp = xo;
      } /* x-dir */
      vi += rowstride;
    } /* y-dir */




    /*
     *  z-direction object
     *
     *
     * -if an object was an input, then
     *  we need to put it on the image
     * -this section was removed
     */


    if(g_obj_list_num) {

      /*fprintf(stdout,"\t\tz-dir object\n");*/
      /* get the points on the viewing plane */

      pp[0].x = pp[0].y = pp[1].y = pp[2].x = 0;
      pp[0].z = pp[1].z = pp[2].z = g_idx[3];
      pp[1].x = pp[2].y = 1e6;

      pp[0]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[0]);
      pp[1]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[1]);
      pp[2]=niikpt_affine_transform_m44(g_image->sto_xyz,pp[2]);

      if(verbose>2) {
        fprintf(stdout,"-d3 plane triangle: %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f | %0.2f %0.2f %0.2f\n",
                pp[0].x,pp[0].y,pp[0].z,
                pp[1].x,pp[1].y,pp[1].z,
                pp[2].x,pp[2].y,pp[2].z);
        off_display_bbox_info(g_off_bb);
      }

      /* pre-calculate the matrix from image (xyz) space to window space */
      Ao = niikmat_inverse(g_matrix[dir2num(g_axis[nimg])]);

      /*
      niikmat_multiply_mat1_free2(Ao,
                                  niikmat_scale_matrix(1.0/g_image->pixdim[1],
                                                       1.0/g_image->pixdim[2],
                                                       1.0/g_image->pixdim[3]));
      */
      niikmat_multiply_mat1_free2(Ao,niikmat_mat44_matrix(g_image->sto_ijk));

      /* go thru the bounding boxes -only in 2 dimensions */
      k = floor( (pp[0].z - g_off_bb->origin.z) / g_off_bb->delta + 0.5);

      bzmin = floor( ( pp[0].z  - g_off_bb->origin.z) / g_off_bb->delta - 0.5);
      bzmax = floor( ( pp[0].z  - g_off_bb->origin.z) / g_off_bb->delta + 0.5);
      bzmin = NIIK_IMINMAX(bzmin-1,0,g_off_bb->zdim-1);
      bzmax = NIIK_IMINMAX(bzmax+1,0,g_off_bb->zdim-1);
      for(k=bzmin; k<=bzmax; k++) {
        nn = k * g_off_bb->area;

        for(j=0; j<g_off_bb->ydim; j++) {
          for(i=0; i<g_off_bb->xdim; i++) {
            n = i + j * g_off_bb->xdim + nn;

            if(verbose>7)
              fprintf(stdout,"-d8 bbox index=%i[%i,%i,%i]\n",n,i,j,k);

            /* within each bounding box, test if the surface intersect with the viewing plane */
            for(m=0; m<g_off_bb->ndata[n] ; m++) {
              f = g_off_bb->data[n][m];

              if(f->pmax.z < pp[0].z) continue;
              if(f->pmin.z > pp[0].z) continue;

              if(!off_check_tri_tri_intersect_with_isectline(pp[0],pp[1],pp[2],
                  f->vert[0]->v,
                  f->vert[1]->v,
                  f->vert[2]->v,
                  &coplanar,
                  &qq[0],
                  &qq[1]) ) {
                continue;
              }

              if(coplanar) {
                continue;
              }
              if(verbose>7) fprintf(stdout,"-d8 intersect: %0.2f %0.2f %0.2f -> %0.2f %0.2f %0.2f\n",
                                      qq[0].x,qq[0].y,qq[0].z,qq[1].x,qq[1].y,qq[1].z);

              /* find the window positions */
              isect[0][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[0].x +
                                  Ao->m[0][1]  * qq[0].y +
                                  Ao->m[0][2]  * qq[0].z +
                                  Ao->m[0][3]);

              isect[0][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[0].x +
                                  Ao->m[1][1]  * qq[0].y +
                                  Ao->m[1][2]  * qq[0].z +
                                  Ao->m[1][3]);

              isect[1][0] = floor(0.5 +
                                  Ao->m[0][0]  * qq[1].x +
                                  Ao->m[0][1]  * qq[1].y +
                                  Ao->m[0][2]  * qq[1].z +
                                  Ao->m[0][3]);

              isect[1][1] = floor(0.5 +
                                  Ao->m[1][0]  * qq[1].x +
                                  Ao->m[1][1]  * qq[1].y +
                                  Ao->m[1][2]  * qq[1].z +
                                  Ao->m[1][3]);

              /* linearly interpolate between the two pointer */
              for(isectfrac=0;
                  isectfrac<=1.0;
                  isectfrac+=fabs(g_matrix[dir2num(g_axis[nimg])]->m[0][0])/5.0) {
                isect[2][0] = isect[1][0] * isectfrac + isect[0][0] * (1.0-isectfrac);
                isect[2][1] = isect[1][1] * isectfrac + isect[0][1] * (1.0-isectfrac);

                if(isect[2][0]<0) continue;
                if(isect[2][1]<0) continue;
                if(isect[2][0]>=g_window_width)  continue;
                if(isect[2][1]>=g_window_height) continue;

                mm = isect[2][1]*rowstride2 + isect[2][0]*3;
                if(f->color!=NULL) { /* colored object */
                  pixels[mm  ] = 255 * f->color[0];
                  pixels[mm+1] = 255 * f->color[1];
                  pixels[mm+2] = 255 * f->color[2];
                  /* fprintf(stdout,"color %3i %3i %3i\n",pixels[mm],pixels[mm+1],pixels[mm+2]); */
                } else { /* no color -> gree */
                  pixels[mm  ] = 0;
                  pixels[mm+1] = 255;
                  pixels[mm+2] = 0;
                }
              }
            } /* each triangle in the bounding box */

          }
        }
      } /* ijk of the bounding boxes */

      Ao=niikmat_free(Ao);

    } /* object for z-dir */

    break;

  default:
    fprintf(stderr,"ERROR: unknown direction\n");
    return 0;
    break;
  } /* switch (g_axis[nimg]) */

  /**************************************
   *  draw a red line
   ***************************************/

  xo = floor(0.5 +
             invA->m[0][0]  * g_idx[1]  +
             invA->m[0][1]  * g_idx[2]  +
             invA->m[0][2]  * g_idx[3]  +
             invA->m[0][3]);
  yo = floor(0.5 +
             invA->m[1][0]  * g_idx[1]  +
             invA->m[1][1]  * g_idx[2]  +
             invA->m[1][2]  * g_idx[3]  +
             invA->m[1][3]);
  zo = floor(0.5 +
             invA->m[2][0]  * g_idx[1]  +
             invA->m[2][1]  * g_idx[2]  +
             invA->m[2][2] * g_idx[3]  +
             invA->m[2][3]);
  invA=niikmat_free(invA);
  if(0) {
    fprintf(stdout,"  x,y,z %4i %4i %4i\n",xo,yo,zo);
  }

  xo = NIIK_IMINMAX(xo-1,0,g_window_width-1);
  yo = NIIK_IMINMAX(yo-1,0,g_window_height-1);

  if(g_draw_red_line) {
    if(g_axis[nimg]=='z') {
      if(g_axis[3-nimg]=='y') {
        /* horizontal line */
        for(xi=0,vi=yo*rowstride2; xi<g_window_width; xi++) {
          pixels[vi++] = 255;
          pixels[vi++] = 40;
          pixels[vi++] = 40;
        }
      }
      /* vertical line */
      else if(g_axis[3-nimg]=='x') {
        for(yi=0,vi=xo*3; yi<g_window_height; yi++,vi+=rowstride2) {
          pixels[vi  ] = 255;
          pixels[vi+1] = 40;
          pixels[vi+2] = 40;
        }
      } else if(g_axis[3-nimg]=='z') {}
    } /* this image is z-dir */

    else if(g_axis[nimg]=='y') {
      /* horizontal line */
      if(g_axis[3-nimg]=='z') {
        for(xi=0,vi=yo*rowstride2; xi<g_window_width; xi++) {
          pixels[vi++] = 255;
          pixels[vi++] = 40;
          pixels[vi++] = 40;
        }
      }
      /* vertical line */
      else if(g_axis[3-nimg]=='x') {
        for(yi=0,vi=xo*3; yi<g_window_height; yi++,vi+=rowstride2) {
          pixels[vi  ] = 255;
          pixels[vi+1] = 40;
          pixels[vi+2] = 40;
        }
      } else if(g_axis[3-nimg]=='z') {}
    } /* this image is y-dir */

    else if(g_axis[nimg]=='x') {
      /* horizontal line */
      if(g_axis[3-nimg]=='z') {
        for(xi=0,vi=yo*rowstride2; xi<g_window_width; xi++) {
          pixels[vi++] = 255;
          pixels[vi++] = 40;
          pixels[vi++] = 40;
        }
      }
      /* vertical line */
      else if(g_axis[3-nimg]=='y') {
        for(yi=0,vi=xo*3; yi<g_window_height; yi++,vi+=rowstride2) {
          pixels[vi  ] = 255;
          pixels[vi+1] = 40;
          pixels[vi+2] = 40;
        }
      } else if(g_axis[3-nimg]=='z') {}
    } /* this image is y-dir */

  }  /* draw a red line */


  /*********************************************
   * green marker
   **********************************************/

  vi = xo * 3 + yo * rowstride2;
  if(g_draw_green_dot) {
    pixels[vi  ] = 0;     /* current position */
    pixels[vi+1] = 255;
    pixels[vi+2] = 0;
    if(xo>0) {  /* left position */
      pixels[vi-3] = 0;
      pixels[vi-2] = 255;
      pixels[vi-1] = 0;
    }
    if(xo<g_window_width-1) { /* right position */
      pixels[vi+3] = 0;
      pixels[vi+4] = 255;
      pixels[vi+5] = 0;
    }
    if(yo>0) { /* up position */
      pixels[vi-rowstride2]   = 0;
      pixels[vi-rowstride2+1] = 255;
      pixels[vi-rowstride2+2] = 0;
    }
    if(yo<g_window_height-1) { /* down position */
      pixels[vi+rowstride2]   = 0;
      pixels[vi+rowstride2+1] = 255;
      pixels[vi+rowstride2+2] = 0;
    }
  } /* g_draw_green_dot */

  gtk_image_set_from_pixbuf(GTK_IMAGE(g_gtkimg[nimg]),g_pixbuf[nimg]);

  return 1;
}
/* NVK_display_imagen */





/***********************************************
 *
 * static int NVK_read_image(char *fname)
 *
 *
 * -reads a nifti image
 * -updates spin buttons and global variables
 * -frees memory from previous image
 * -usually called from open_action();
 *
 *
 *
 * see also
 *
 *  void open_action();
 *  void reload_action();
 *  int main(argc,argv);
 *
 *
 ***********************************************/

static int NVK_read_image(char *fname) {
  nifti_image *img;
  double imin,imax,ran;
  char
  *cptr,
  curr_folder[4097];
  int n;
  char fcname[32]="NVK_read_image";

  if(verbose) fprintf(stdout,"-d1 reading %s in NVK_read_image\n",fname);

  if((img=niik_image_read(fname))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,fname);
    return 0;
  }
  if(0) {
    nifti_image_infodump( img );
    fprintf(stdout,"  fname  = %s\n",img->fname);
    fprintf(stdout,"  iname  = %s\n",img->iname);
    fprintf(stdout,"  iname_offset  = %i\n",img->iname_offset);
  }

  /* check if dim match with mask image */
  if(g_maskimg!=NULL) {
    if(g_maskimg->nx != img->nx) {
      fprintf(stderr,"[%s] ERORR: x-dim did not match %i, %i (mask)\n",fcname,img->nx,g_maskimg->nx);
      return 0;
    } else if(g_maskimg->ny != img->ny) {
      fprintf(stderr,"[%s] ERORR: y-dim did not match %i, %i (mask)\n",fcname,img->ny,g_maskimg->ny);
      return 0;
    } else if(g_maskimg->nz != img->nz) {
      fprintf(stderr,"[%s] ERORR: z-dim did not match %i, %i (mask)\n",fcname,img->nz,g_maskimg->nz);
      return 0;
    } /* dim check */
  } /* mask exists */


  if(g_image      != NULL) nifti_image_free(g_image);
  if(g_image_data != NULL) free(  g_image_data );

  g_image = img;
  if((g_image_data = niik_image_get_voxels_as_float_vector(g_image))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_float_vector\n");
    return 0;
  }

  if (getcwd(curr_folder, sizeof(curr_folder)) == NULL) {
    fprintf(stderr,"ERROR: getcwd \n");
    exit(1);
  }

  if(fname[0]=='/')
    sprintf(g_curr_filename,"%s",fname);
  else
    sprintf(g_curr_filename,"%s/%s",curr_folder,fname);

  if(verbose) fprintf(stdout,"-d1 g_curr_filename = %s\n",g_curr_filename);

  /* update image pixel size */
  g_img_delta.x = g_image->pixdim[1];
  g_img_delta.y = g_image->pixdim[2];
  g_img_delta.z = g_image->pixdim[3];

  /* update image FOV */
  g_img_fov.x = g_image->dim[1] * g_img_delta.x;
  g_img_fov.y = g_image->dim[2] * g_img_delta.y;
  g_img_fov.z = g_image->dim[3] * g_img_delta.z;

  /* update current index */
  for(n=1; n<=5; n++) g_idx[n] = NIIK_IMINMAX(g_idx[n],0,g_image->dim[n]-1);

  /* update the min/max spin buttons */
  imin = g_image_data[niik_get_min_index_float_vector( g_image_data, g_image->nvox )];
  imax = g_image_data[niik_get_max_index_float_vector( g_image_data, g_image->nvox )];

  ran = imax - imin;
  imin = imin - ran*0.2;
  imax = imax + ran*0.2;

  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_min), imin, imax );
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_max), imin, imax );


  /* change the increment and digits for the min/max spin buttons */
  if(ran < 1) {
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_min),0.01,0.1);
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_max),0.01,0.1);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_min),2);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_max),2);
  } else if(ran < 10) {
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_min),0.1,1);
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_max),0.1,1);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_min),1);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_max),1);
  } else {
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_min),1,10);
    gtk_spin_button_set_increments ( GTK_SPIN_BUTTON(g_spin_max),1,10);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_min),0);
    gtk_spin_button_set_digits     ( GTK_SPIN_BUTTON(g_spin_max),0);
  }


  /* update the x/y/z/t spin buttons */
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_xdim), 0, img->dim[1]-1 );
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_ydim), 0, img->dim[2]-1 );
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_zdim), 0, img->dim[3]-1 );
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_tdim), 0, img->dim[4]-1 );
  gtk_spin_button_set_range (GTK_SPIN_BUTTON(g_spin_udim), 0, img->dim[5]-1 );

  gtk_spin_button_set_value (GTK_SPIN_BUTTON(g_spin_xdim), g_idx[1]);
  gtk_spin_button_set_value (GTK_SPIN_BUTTON(g_spin_ydim), g_idx[2]);
  gtk_spin_button_set_value (GTK_SPIN_BUTTON(g_spin_zdim), g_idx[3]);
  gtk_spin_button_set_value (GTK_SPIN_BUTTON(g_spin_tdim), g_idx[4]);
  gtk_spin_button_set_value (GTK_SPIN_BUTTON(g_spin_udim), g_idx[5]);


  /* hide 4th dimension (time) */
  if(g_image->ndim<=3 || g_image->nt<=1) {
    gtk_widget_hide(g_spin_dim[4]);
    gtk_widget_hide(g_lbl_dim[4]);
    g_idx[4] = 0;
  } else {
    gtk_widget_show(g_spin_dim[4]);
    gtk_widget_show(g_lbl_dim[4]);
  }


  /* hide 5th dimension (time) */
  if(g_image->ndim<=4 || g_image->nu<=1) {
    gtk_widget_hide(g_spin_dim[5]);
    gtk_widget_hide(g_lbl_dim[5]);
    g_idx[5] = 0;
  } else {
    gtk_widget_show(g_spin_dim[5]);
    gtk_widget_show(g_lbl_dim[5]);
  }


  sprintf(curr_folder,"%s",g_curr_filename);

  cptr = strrchr (curr_folder,'/');
  gtk_window_set_title(GTK_WINDOW(g_main_window),cptr+1);


  /*  free(g_image->data);
      g_image->data = NULL; */

  /* first display of this image */
  NVK_display_image();
  return 1;
} /* static int NVK_read_image(char *fname) */





/***********************************************
 *
 * spin buttons for editing
 *
 *
 * void NVK_change_edit_pen_size
 *
 * -changes the pen size
 *
 *
 *
 * void NVK_change_edit_blue
 *
 * -changes the blue strength
 * -blue function is not perfect as it tends
 *  to darken the image
 * -but I can't figure out how to linearly change
 *  so that 100  --> totally blue in the mask
 *  and normal outside the mask and 0 --> normal
 *  everywhere
 *
 *
 * void NVK_change_seed_fill_tolerance()
 *
 * -middle button + Control or shift
 * -undo features
 *
 *
 ***********************************************/


static void NVK_change_edit_pen_size() {
  g_pen_size =  gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_pen_size));
}

static void NVK_change_edit_blue() {
  g_edit_blue =  gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_edit_blue)) / 100.0;
  NVK_display_image();
}

static void NVK_change_edit_fill_tol() {
  g_edit_fill_tol =  gtk_spin_button_get_value(GTK_SPIN_BUTTON(g_spin_fill_tol));
}



/***********************************************
 *
 * timer, save file every min
 *
 ***********************************************/

static gboolean NVK_automatic_save_mask_file(GtkWidget *widget) {
  char fcname[64]="NVK_automatic_save_mask_file";
  static char *fn=NULL;
  char tmp_fn[4096];
  char *bn;
  int verbose=0;
  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  ctm=time(NULL);
  stm=localtime(&ctm);
  if(g_maskimg == NULL) return FALSE;
  if(g_automatic_mask_save == 0) return FALSE;
  if (widget->window == NULL) return FALSE;
  sprintf(tmp_fn,"%s",g_maskimg->fname);
  bn = strrchr (g_maskimg->fname,'/');
  if(bn==NULL) bn=g_maskimg->fname;
  else bn++;
  if(fn==NULL) {
    fn=(char *)calloc(4096,sizeof(char));
    strftime(tmstr,256,"%Y-%m-%d_%T",stm);
    sprintf(fn,"tmp_NVK_save_%s_%s",tmstr,bn);
  }
  if(verbose>=0) {
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    fprintf(stdout,"[%s,%s] automatic mask save = %s\n",fcname,tmstr,fn);
  }
  NVK_write_mask(fn,g_maskimg);
  sprintf(g_maskimg->fname,"%s",tmp_fn);
  // gtk_widget_queue_draw(widget);
  return TRUE;
}

static void NVK_toggle_autosave() {
  if(g_automatic_mask_save>0) g_automatic_mask_save=0;
  else g_automatic_mask_save=1;
}

static void NVK_peek_info(){
  niikpt pp;
  double min_distance=1e10;
  int closest_obj=-1;
  int closest_index=-1;
  kvert *cv=NULL;

  int n;
  pp.x=g_idx[1];
  pp.y=g_idx[2];
  pp.z=g_idx[3];
  pp=niikpt_affine_transform_m44(g_image->sto_xyz,pp);
  fprintf(stdout,"%d,%d,%d : %0.2f,%0.2f,%0.2f\n",g_idx[1],g_idx[2],g_idx[3],
                                                  pp.x,pp.y,pp.z);
  
  /*search for closest vertex*/
  for(n=0;n<g_obj_list_num;n++)
  {
      int i;
      kvert *v;
      for(v=g_obj_list[n]->vert;v!=NULL;v=v->next)
      {
        double d2=niikpt_distance2(pp,v->v);
        if(d2<min_distance)
        {
          closest_obj=n;
          closest_index=v->index;
          min_distance=d2;
          cv=v;
        }
      }
  }
  if(closest_obj>=0 && cv)
  {
    fprintf(stdout,"%s:%d %0.2f,%0.2f,%0.2f %f\n",g_obj_list[closest_obj]->fname,closest_index,cv->v.x,cv->v.y,cv->v.z,sqrt(min_distance));
  }
}

/***********************************************
 *
 * quit!
 *
 ***********************************************/
static void quit_action () {
  gtk_main_quit();
}

static int NVK_write_minc_mask(char *filename,nifti_image *maskimg) {
  static image_metadata *meta=NULL;
  int n,verbose=0;
  static float *fimg=NULL;
  char fcname[64]="NVK_write_minc_mask";
  if(meta==NULL) {
    meta = (image_metadata *)calloc( 1, sizeof(image_metadata) ) ;
    meta->start  = calloc(3,sizeof(float));
    meta->step   = calloc(3,sizeof(float));
    meta->length = calloc(3,sizeof(int));
    meta->start[0]=maskimg->sto_xyz.m[2][3];
    meta->start[1]=maskimg->sto_xyz.m[1][3];
    meta->start[2]=maskimg->sto_xyz.m[0][3];
    meta->step[0]=maskimg->dz;
    meta->step[1]=maskimg->dy;
    meta->step[2]=maskimg->dx;
    meta->length[0]=maskimg->nz;
    meta->length[1]=maskimg->ny;
    meta->length[2]=maskimg->nx;
    fimg=niik_image_get_voxels_as_float_vector(maskimg);
  } /* define meta once */
  else {
    for(n=0; n<maskimg->nvox; n++) fimg[n]=niik_image_get_voxel(maskimg,n);
  }
  if(verbose>0) fprintf(stderr,"[%s] writing .mnc file, %s\n",fcname,filename);
  n=write_minc(filename, fimg, meta, (BOOLEAN)1);
  if(verbose>0) fprintf(stderr,"[%s] finished writing minc file, %i, %s\n",fcname,n,filename);
  if(n>0) return 0;
  return 1;
}

static void NVK_write_mask(char *filename,nifti_image *maskimg)
/* writes the mask for editing mode */
{
  char fcname[32]="NVK_write_mask";
  GtkWidget *dialog;
  int verbose=0;
  if(!strncmp(filename+(strlen(filename)-7),".nii.gz",7)) {
    if(verbose>0) fprintf(stderr,"[%s] writing .nii.gz file, %s\n",fcname,filename);
    nifti_image_write(maskimg);
    if(nifti_io_get_problem_code_kunio()>0) {
      fprintf(stderr,"[%s] ERROR: nifti_image_write\n",fcname);
      dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                       GTK_DIALOG_MODAL,
                                       GTK_MESSAGE_ERROR,
                                       GTK_BUTTONS_OK,
                                       "***ERROR***\ncould not write mask:\n%s\n",filename);
      gtk_dialog_run (GTK_DIALOG (dialog));
      gtk_widget_destroy (dialog);
      return;
    }
  } else if(!strncmp(filename+(strlen(filename)-4),".nii",4)) {
    if(verbose>0) fprintf(stderr,"[%s] writing .nii file, %s\n",fcname,filename);
    nifti_image_write(maskimg);
    if(nifti_io_get_problem_code_kunio()>0) {
      fprintf(stderr,"[%s] ERROR: nifti_image_write\n",fcname);
      dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                       GTK_DIALOG_MODAL,
                                       GTK_MESSAGE_ERROR,
                                       GTK_BUTTONS_OK,
                                       "***ERROR***\ncould not write mask:\n%s\n",filename);
      gtk_dialog_run (GTK_DIALOG (dialog));
      gtk_widget_destroy (dialog);
      return;
    }
  } else if(!strncmp(filename+(strlen(filename)-4),".mnc",4)) {
    if(NVK_write_minc_mask(filename,g_maskimg)==0) {
      fprintf(stderr,"[%s] ERROR: NVK_write_minc_mask\n",fcname);
      dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                       GTK_DIALOG_MODAL,
                                       GTK_MESSAGE_ERROR,
                                       GTK_BUTTONS_OK,
                                       "***ERROR***\ncould not write mask:\n%s\n",filename);
      gtk_dialog_run (GTK_DIALOG (dialog));
      gtk_widget_destroy (dialog);
      return;
    }
  } else {
    fprintf(stderr,"[%s] ERROR: writing file, %s\n",fcname,filename);
    dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                     GTK_DIALOG_MODAL,
                                     GTK_MESSAGE_ERROR,
                                     GTK_BUTTONS_OK,
                                     "***ERROR***\ncould not write mask:\n%s\n",filename);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
  }
} /* NVK_write_mask */




/***********************************************
 *
 * open
 *
 * -opens a dialog to choose a file
 * -checks if the new image is actuallly loaded
 *  then removes the old image from memory
 * -updates the min/max intensity spin buttons
 *  (range)
 * -updates the x/y/z/t/u spin buttons (range)
 * -but keeps the current positions and
 *  intensity range (as long as they within the
 *  range)
 * -currently, it does not remember about the
 *  last condition about the file extension
 *  filters
 *
 ***********************************************/

static void open_action () {
  GtkWidget
  *  file_chooser, *dialog;
  GtkRecentManager *  manager = NULL;
  GtkFileFilter *filter[6];
  char
  curr_folder[4096],
              *filename = NULL, *cptr;
  gchar *uri;

  if(verbose) fprintf(stdout,"-d1 file_chooser diaglog\n");
  file_chooser = gtk_file_chooser_dialog_new("Open File",
                 GTK_WINDOW(g_main_window),
                 GTK_FILE_CHOOSER_ACTION_OPEN,
                 GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                 GTK_STOCK_OPEN,   GTK_RESPONSE_ACCEPT,
                 NULL);

  if(verbose) fprintf(stdout,"-d1 filters for file_chooser \n");
  filter[0]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[0]),"*.nii.gz");
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[0]),"*.nii");
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[0]),"*.mnc.gz");
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[0]),"*.mnc");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[0]),"readable files");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[0]));

  filter[1]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[1]),"*.nii");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[1])," NIfTI");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[1]));

  filter[2]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[2]),"*.nii.gz");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[2]),"gz NIfTI");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[2]));

  filter[3]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[3]),"*.mnc");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[3]),"MINC files");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[3]));

  filter[4]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[4]),"*.mnc.gz");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[4]),"MINC gz files");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[4]));

  filter[5]=gtk_file_filter_new();
  gtk_file_filter_add_pattern (GTK_FILE_FILTER(filter[5]),"*");
  gtk_file_filter_set_name    (GTK_FILE_FILTER(filter[5]),"all files");
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(file_chooser), GTK_FILE_FILTER(filter[5]));

  if(verbose) fprintf(stdout,"-d1 current folder\n");
  sprintf(curr_folder,"%s",g_curr_filename);
  cptr = strrchr (curr_folder,'/');
  cptr[0] = '\0';

  if(verbose) fprintf(stdout,"-d1 set current_folder = %s\n",curr_folder);
  gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(file_chooser),curr_folder);

  if(verbose) fprintf(stdout,"-d1 run dialog\n");
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT) {
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    uri =      gtk_file_chooser_get_uri   (GTK_FILE_CHOOSER (file_chooser));
    manager = gtk_recent_manager_get_default ();
    gtk_recent_manager_add_item    (manager,uri);

    if(verbose) fprintf(stdout,"open file %s!\n",filename);

    if(!NVK_read_image(filename)) {
      fprintf(stdout,"ERROR: reading image\n");

      dialog = gtk_message_dialog_new( GTK_WINDOW(file_chooser),
                                       GTK_DIALOG_MODAL,
                                       GTK_MESSAGE_ERROR,
                                       GTK_BUTTONS_OK,
                                       "Could not read file:\n%s\n",filename);

      gtk_dialog_run (GTK_DIALOG (dialog));
      gtk_widget_destroy (dialog);
    }

    g_free (filename);
    g_free (uri);
  }

  gtk_widget_destroy (file_chooser);

  return;
} /* open action */





/***********************************************
 *
 * ask_quit
 *
 * -asks the user if the window should be closed
 * -just to make sure before really closing
 *
 ***********************************************/

static void NVK_ask_quit() {
  GtkWidget *dialog;
  if(g_maskimg!=NULL) {
    dialog=gtk_message_dialog_new ( GTK_WINDOW(g_main_window),
                                    GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                    GTK_MESSAGE_QUESTION,
                                    GTK_BUTTONS_YES_NO,
                                    "Do you want to save before exit [%s]?",g_maskimg->fname);
    switch(gtk_dialog_run (GTK_DIALOG (dialog))) {
    case GTK_RESPONSE_YES:
    case GTK_RESPONSE_ACCEPT:
      /* fprintf(stdout,"  writing an image %s\n",g_maskimg->fname); */
      if(!niik_image_set_voxels_from_uint8_vector(g_maskimg,g_mask_data)) {
        fprintf(stderr,"ERROR: couldn't write %s\n",g_maskimg->fname);
        exit(0);
      }
      NVK_write_mask(g_maskimg->fname,g_maskimg);
      // nifti_image_write(g_maskimg);

      quit_action(); /* gtk_main_quit(); (same things) */
      return;

    case GTK_RESPONSE_NO:
    case GTK_RESPONSE_REJECT:
      gtk_widget_destroy(dialog);
      break;
    } /* responce */
  } /* mask image exists */

  dialog=gtk_message_dialog_new ( GTK_WINDOW(g_main_window),
                                  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                  GTK_MESSAGE_QUESTION,
                                  GTK_BUTTONS_YES_NO,
                                  "Are you sure you want to quit?");
  switch(gtk_dialog_run (GTK_DIALOG (dialog))) {
  case GTK_RESPONSE_YES:
  case GTK_RESPONSE_ACCEPT:
    gtk_main_quit();
  case GTK_RESPONSE_NO:
  case GTK_RESPONSE_REJECT:
    gtk_widget_destroy(dialog);
  }
  return;
}





/***********************************************
 *
 * NVK_save_action
 *
 * -saves the mask image
 *
 ***********************************************/

static void NVK_save_action() {
  GtkWidget *dialog;
  if(g_image==NULL) return;
  if(g_maskimg!=NULL) {
    dialog=gtk_message_dialog_new ( GTK_WINDOW(g_main_window),
                                    GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                    GTK_MESSAGE_QUESTION,
                                    GTK_BUTTONS_YES_NO,
                                    "Are you sure you want to save\n [%s]?",g_maskimg->fname);
    switch(gtk_dialog_run (GTK_DIALOG (dialog))) {
    case GTK_RESPONSE_YES:
    case GTK_RESPONSE_ACCEPT:
      /* fprintf(stdout,"  writing an image %s\n",g_maskimg->fname); */
      /* we use the same pointer so that g_mask_data = g_maskimg->data
       * therefore, we don't need to update the g_maskimg anymore
       if(!nifti_k_set_voxel_values_uint8_vector(g_maskimg,g_mask_data)){
       fprintf(stderr,"ERROR: couldn't write %s\n",g_maskimg->fname);
       exit(0); } */
      NVK_write_mask(g_maskimg->fname,g_maskimg);
    // nifti_image_write(g_maskimg);

    case GTK_RESPONSE_NO:
    case GTK_RESPONSE_REJECT:
      gtk_widget_destroy(dialog);
      break;
    } /* responce */
  } /* mask image exists */
  else {
    dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                     GTK_DIALOG_MODAL,
                                     GTK_MESSAGE_ERROR,
                                     GTK_BUTTONS_OK,
                                     "There is no mask image to save");
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
  } /* no mask image  */
  return;
}




/***********************************************
 *
 * view_info
 *
 * -shows the image info as in the function in
 *  nifti1_io.c
 * -uses nifti_image_to_ascii
 * -opens a new message dialog
 *
 ***********************************************/

static void view_info_action () {
  GtkWidget *dialog;
  dialog=gtk_message_dialog_new ( GTK_WINDOW(g_main_window),
                                  GTK_DIALOG_MODAL | GTK_DIALOG_DESTROY_WITH_PARENT,
                                  GTK_MESSAGE_INFO,
                                  GTK_BUTTONS_OK,
                                  "Image info: \n%s",nifti_image_to_ascii( g_image ));
  gtk_dialog_run (GTK_DIALOG (dialog) );
  gtk_widget_destroy(dialog);
  return;
}




/***********************************************
 *
 * view reset
 *
 * -this function is linked to menu bar so the
 *  shortcut is Control+C
 * -resets the zoom and shifts
 * -Should this reset the flip?
 *
 ***********************************************/

static void view_reset_action () {
  g_zoom[1]=g_zoom[2]=g_zoom[3]=2.0;    /* zoom for x/y/z-directions */
  g_shift[1][1] = g_shift[1][2] =       /* shift within window */
                    g_shift[2][1] = g_shift[2][2] =
                                      g_shift[3][1] = g_shift[3][2] = 0;
  NVK_display_image();
}



/***********************************************
 *
 * simple re-loading function
 *
 * -this function is linked to menu bar so the
 *  shortcut is Control+R
 * -opens the same image again
 * -the filename is saved as g_curr_filename
 * -it's useful when processing the same image
 *  multiple times
 *
 ***********************************************/
static void reload_action () {
  GtkWidget * dialog;
  if(!NVK_read_image(g_curr_filename)) {
    fprintf(stderr,"ERROR: reading image\n");
    dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                     GTK_DIALOG_MODAL,
                                     GTK_MESSAGE_ERROR,
                                     GTK_BUTTONS_OK,
                                     "Could not read file:\n%s\n",g_curr_filename);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
  }
}




/***********************************************
 *
 * flip functions
 *
 * -these functions are linked to the menu bar
 * -there are no shortcut keys for them
 * -these change the display matrix in
 *  static void NVK_update_matrix(int nimg);
 * -working ok
 *
 ***********************************************/
static void NVK_view_flip_x () {
  if(g_image==NULL) return;
  g_flip[1] = (g_flip[1]==0);
  NVK_display_image();
}

static void NVK_view_flip_y () {
  if(g_image==NULL) return;
  g_flip[2] = (g_flip[2]==0);
  NVK_display_image();
}

static void NVK_view_flip_z () {
  if(g_image==NULL) return;
  g_flip[3] = (g_flip[3]==0);
  NVK_display_image();
}









/**********************************************
 *
 * axis direction
 *
 * -axis1 changes the left image
 * -axis2 changes the right image
 * -makes sure that each one is different
 * -these functions are linked to the menu bar
 * -the key is Control+A for axis2
 *  and Shift+Control+A for axis1
 *
 **********************************************/

static void view_change_axis1_action () {
  if(g_image==NULL) return;
  switch(g_axis[1]) {
  case 'x':
    g_axis[1]='y';
    if(g_axis[2]=='y') g_axis[1]='z';
    break;
  case 'y':
    g_axis[1]='z';
    if(g_axis[2]=='z') g_axis[1]='x';
    break;
  case 'z':
    g_axis[1]='x';
    if(g_axis[2]=='x') g_axis[1]='y';
    break;
  }
  NVK_display_image();
}

static void view_change_axis2_action () {
  if(g_image==NULL) return;
  switch(g_axis[2]) {
  case 'x':
    g_axis[2]='y';
    if(g_axis[1]=='y') g_axis[2]='z';
    break;
  case 'y':
    g_axis[2]='z';
    if(g_axis[1]=='z') g_axis[2]='x';
    break;
  case 'z':
    g_axis[2]='x';
    if(g_axis[1]=='x') g_axis[2]='y';
    break;
  }
  NVK_display_image();
}



/**********************************************
 *
 * resets the intensity
 *
 * -basically sets the spin buttons to min/max
 *  intensity in the image
 * -the function is linked to the menu bar so the
 *  shortcut is Control+L
 *
 **********************************************/

static void view_reset_intensity () {
  if(g_image==NULL) return;
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_min),g_image_data[niik_get_min_index_float_vector( g_image_data, g_image->nvox )]);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_max),g_image_data[niik_get_max_index_float_vector( g_image_data, g_image->nvox )]);
  NVK_display_image();
}



/**********************************************
 *
 * zoom in/out functions
 *
 * -added Kunio 2012-05-07
 * -on my mac, zoom didn't work so i'm adding this
 *
 **********************************************/

static void NVK_view_zoom_in() {
  g_zoom[1]*=1.1;
  g_zoom[2]*=1.1;
  g_zoom[3]*=1.1;
  NVK_display_image();
  return;
}

static void NVK_view_zoom_out() {
  g_zoom[1]/=1.1;
  g_zoom[2]/=1.1;
  g_zoom[3]/=1.1;
  NVK_display_image();
  return;
}

static void NVK_switch_draw_red_line() {
  g_draw_red_line = (g_draw_red_line == 0);
  NVK_display_image();
  return;
}

static void NVK_switch_draw_green_dot() {
  g_draw_green_dot = (g_draw_green_dot == 0);
  NVK_display_image();
  return;
}


/**********************************************
 *
 * window expand/shrink functions
 *
 * -because I don't want the user changing the window size,
 *  these functions can expand/shrink the window by
 *  in(de)crement of 64 pixels
 * -frees memory properly (?)
 * -these functions are linked to the menu bar so the
 *  shortcuts are Control+1 and Control+2 for shrink/expand
 *  respectively
 *
 **********************************************/

static void NVK_window_expand() {
  int n;
  if(g_image==NULL) return;

  g_window_width  = NIIK_IMINMAX(g_window_width+64, 256,max_window_size);
  g_window_height = NIIK_IMINMAX(g_window_height+64,256,max_window_size);

  for(n=1; n<3; n++) {
    if(g_pixbuf[n]!=NULL) {
      g_object_unref(g_pixbuf[n]);
    }
    g_pixbuf[n]=gdk_pixbuf_new(GDK_COLORSPACE_RGB,FALSE,8,g_window_width,g_window_height);
    gtk_image_set_from_pixbuf(GTK_IMAGE(g_gtkimg[n]),g_pixbuf[n]);
  }
  gtk_window_resize(GTK_WINDOW(g_main_window),g_window_width,g_window_height);
  NVK_display_image();
  return;
} /* NVK_window_expand */

static void NVK_window_shrink() {
  int n;
  if(g_image==NULL) return;
  g_window_width  = NIIK_IMINMAX(g_window_width-64, 256,max_window_size);
  g_window_height = NIIK_IMINMAX(g_window_height-64,256,max_window_size);
  for(n=1; n<3; n++) {
    if(g_pixbuf[n]!=NULL) {
      g_object_unref(g_pixbuf[n]);
    }
    g_pixbuf[n]=gdk_pixbuf_new(GDK_COLORSPACE_RGB,FALSE,8,g_window_width,g_window_height);
    gtk_image_set_from_pixbuf(GTK_IMAGE(g_gtkimg[n]),g_pixbuf[n]);
  }
  gtk_window_resize(GTK_WINDOW(g_main_window),g_window_width,g_window_height);
  NVK_display_image();
  return;
}/* NVK_window_shrink */



/******************************************
 *
 * coordinates function
 *
 * -I wanted to be able to use qform
 *  but it isn't working well (2012-02-25)
 * -maybe we should just assume the image directions
 *  and convert images to this viewer's convention
 *  which is
 *     x  =  Right-to-Left
 *     y  =  Anterior-to-Posterior
 *     z  =  Superior-to-Inferior
 *
 *****************************************
static void NVK_view_coordinate_xyz()
{
  return;
  fprintf(stdout,"  view coordinate xyz\n");
  g_coordinates = 2;
  NVK_display_image();
  return;
}*/


/***************************************
 *
 * undo function
 *
 * -undo the processing
 * -can't redo (undo the undo)
 * -voxel coordinates are saved in g_edit_undo_list[idx][2-]
 * -current idx is saved in g_edit_undo_idx
 * -if g_edit_undo_list[idx][1] == 0   --> undefined
 *                              == 1   --> voxels were added
 *                              == -1  --> voxels were removed
 * -so this function reverses these changes
 *
 ***************************************/

static void NVK_edit_undo_action() {
  int n;
  GtkWidget * dialog;

  /* SHOW ERROR MESSAGE */
  if(g_maskimg==NULL) {
    dialog = gtk_message_dialog_new( GTK_WINDOW(g_main_window),
                                     GTK_DIALOG_MODAL,
                                     GTK_MESSAGE_ERROR,
                                     GTK_BUTTONS_OK,
                                     "Undo is enabled only in editing mode\n");
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
    return;
  }

  if(g_edit_undo_idx<0) return;

  if(g_edit_undo_list[g_edit_undo_idx][1] == 0 ) return;
  if(verbose>2) {
    if(g_edit_undo_list[g_edit_undo_idx][1] > 0) {
      fprintf(stdout,"-d3 undo write  [%2i] %i\n",g_edit_undo_idx,g_edit_undo_list[g_edit_undo_idx][0]-2);
    } else {
      fprintf(stdout,"-d3 undo delete [%2i] %i\n",g_edit_undo_idx,g_edit_undo_list[g_edit_undo_idx][0]-2);
    }
  }
  /* DO UNDO HERE */
  if(g_edit_undo_list[g_edit_undo_idx][1] > 0) { /* we added, so remove these */
    for(n=2; n<g_edit_undo_list[g_edit_undo_idx][0]; n++) {
      g_mask_data[g_edit_undo_list[g_edit_undo_idx][n]]=0;
    }
  } else {                     /* we deleted, so add these */
    for(n=2; n<g_edit_undo_list[g_edit_undo_idx][0]; n++) {
      g_mask_data[g_edit_undo_list[g_edit_undo_idx][n]]=1;
    }
  }
  g_edit_undo_list[g_edit_undo_idx][0] = 2;
  g_edit_undo_list[g_edit_undo_idx][1] = 0;
  g_edit_undo_idx--;
  if(g_edit_undo_idx<0) g_edit_undo_idx=g_edit_max_undo-1;
  NVK_display_image();
}

/* Create a list of entries which are passed to the Action constructor.
 * This is a huge convenience over building Actions by hand. */
static GtkActionEntry entries[] = {
  { "FileMenuAction",   NULL, "_File" },  /* name, stock id, label */
  { "EditMenuAction",   NULL, "_Edit" },
  { "ViewMenuAction",   NULL, "_View" },
  { "WindowMenuAction", NULL, "_Window" },

  {
    "EditUndoAction", GTK_STOCK_UNDO,
    "_undo", "<control>z",
    "Undo",
    G_CALLBACK( NVK_edit_undo_action)
  },
  {
    "EditPeekAction", NULL,
    "_Peek", NULL,
    "Peek",
    G_CALLBACK( NVK_peek_info)
  },
  {
    "EditAutoSaveToggleAction", NULL,
    "_Autosave", NULL,
    "Autosave",
    G_CALLBACK( NVK_toggle_autosave)
  },


  {
    "ViewResetAction", GTK_STOCK_HOME,
    "_Reset", "<control>c",
    "Reset",
    G_CALLBACK( view_reset_action)
  },
  {
    "ViewInfoAction", GTK_STOCK_INFO,
    "_Info", "<control>i",
    "Info",
    G_CALLBACK( view_info_action)
  },

  {
    "ViewChangeAxis1", NULL,
    "_Axis1", "<control><shift>a",
    "Change axis 1",
    G_CALLBACK( view_change_axis1_action)
  },
  {
    "ViewChangeAxis2", NULL,
    "_Axis2", "<control>a",
    "Change axis 2",
    G_CALLBACK( view_change_axis2_action)
  },
  {
    "ViewResetIntensity", NULL,
    "_Reset Intensity", "<control>l",
    "Reset Intensity",
    G_CALLBACK( view_reset_intensity)
  },

  { "ViewFlipMenuAction",   NULL, "_Flip" },
  {
    "ViewFlipX", NULL,
    "flip _x", "",
    "flip x",
    G_CALLBACK( NVK_view_flip_x)
  },
  {
    "ViewFlipY", NULL,
    "flip _y", "",
    "flip y",
    G_CALLBACK( NVK_view_flip_y)
  },
  {
    "ViewFlipZ", NULL,
    "flip _z", "",
    "flip z",
    G_CALLBACK( NVK_view_flip_z)
  },

  { "ViewZoomMenuAction",   NULL, "_Zoom" },
  {
    "ViewZoomIn", NULL,
    "zoom _in", "<Ctrl>k",
    "zoom in",
    G_CALLBACK( NVK_view_zoom_in)
  },
  {
    "ViewZoomOut", NULL,
    "zoom _out", "<Shift><Ctrl>k",
    "zoom out",
    G_CALLBACK( NVK_view_zoom_out)
  },

  /*  { "ViewCoordinateMenuAction",   NULL, "_Coordinate" },
      { "ViewCoordinateXYZ", NULL,
      "coordinate _xyz", "",
      "Coordinate xyz",
      G_CALLBACK( NVK_view_coordinate_xyz) },*/

  { "ViewDrawLineDot",   NULL, "_Navigation Help" },
  {
    "ViewDrawGreenDot", NULL,
    "Draw green dot", NULL,
    "Draw _green dot",
    G_CALLBACK( NVK_switch_draw_green_dot)
  },
  {
    "ViewDrawRedLine", NULL,
    "Draw red line", NULL,
    "Draw _red line",
    G_CALLBACK( NVK_switch_draw_red_line)
  },

  {
    "WindowExpand", NULL,
    "_Expand", "<control>2",
    "Expand window",
    G_CALLBACK( NVK_window_expand)
  },
  {
    "WindowShrink", NULL,
    "_Shrink", "<control>1",
    "Shrink window",
    G_CALLBACK( NVK_window_shrink)
  },

  {
    "FileOpenAction",GTK_STOCK_OPEN,
    "_Open", "<control>v",
    "Open",
    G_CALLBACK( open_action)
  },

  {
    "FileReloadAction",GTK_STOCK_REFRESH,
    "_Reload", "<control>r",
    "Reload",
    G_CALLBACK( reload_action)
  },

  {
    "FileSaveAction",GTK_STOCK_SAVE,
    "_Save", "<control>s",
    "Save",
    G_CALLBACK( NVK_save_action)
  },

  {
    "FileQuitAction", GTK_STOCK_QUIT,
    "_Quit", "<control>Q",
    "Quit",
    G_CALLBACK (NVK_ask_quit)
  }
};

static guint n_entries = G_N_ELEMENTS (entries);




/***********************************************
 *                  MAIN                       *
 ***********************************************/
void usage() {
  printf("  Kunio's image viewer\n");
  printf("  version %i.%i.%i\n",NVK_MAJOR_VERSION,NVK_MINOR_VERSION,NVK_MICRO_VERSION);
}


int main(int argc,char *argv[],char *envp[]) {
  GtkWidget
  *mainv,               /* main windows' vertical box */
  *imagh,               /* images on horizontal boxes */
  *labelh,              /* labels on horizontal boxes */
  *edith,               /* edit things on horizontal box */
  *lbledit[4],          /* labels for editing */
  *event_box[3],        /* event box for image; mouse things */
  *img_align[3];        /* image alignment */
  GtkUIManager *uiman;
  /* GdkGeometry geomhints; */
  char
  **CSlist,
  tmpstr[4096],
  curr_folder[4096];
  double *dlist;
  int    *ilist;
  int    list_num;

  GtkActionGroup      *action_group;          /* Packing group for our Actions */
  GError              *error;                 /* For reporting exceptions or errors */
  GtkWidget           *menubar;               /* The actual menubar */

  char
  iconfile[4096],
           uixmlfile[4096];
  const char  *NIIKDIR;
  FILE *fp;

  double
  this_spin_pct[2],
                this_spin_max = -2,
                this_spin_min = -1;
  int
  flag_check_voxelsize = 1,
  n,m,
  nc,sc;   /* argv counters */
  const char *fcname="NVK";


  struct option long_options[] = {
    {"version", no_argument,    0, 'v'},
    {"help", no_argument, 0, 'h'},
    {"no-pixel-check",no_argument,    &flag_check_voxelsize, 0},
    {"no-voxel-check",no_argument,    &flag_check_voxelsize, 0},
    {"mask-autosave",no_argument,&g_automatic_mask_save,1},
    {"flipxyz",no_argument,0,'f'},
    {"minmax",required_argument,0,'m'},
    {"shift",required_argument,0,'s'},
    {"debug",required_argument,0,'d'},
    {"zoom",required_argument,0,'z'},
    {"cidx",required_argument,0,'c'},
    {"off",required_argument,0,'o'},
    {"obj",required_argument,0,'o'},
    {"pct",required_argument,0,'p'},
    {"idx",required_argument,0,'i'},
    {"usage",no_argument,0,'U'},
    {0, 0, 0, 0}
  };



  if((NIIKDIR = get_NIIKDIR())==NULL) {
    exit(1);
  }

  g_automatic_mask_save = 0;

  /*
   * no usage is ok
   if(argc==1){
   usage();
   exit(9); }
  */

  g_window_width = g_window_height = 384;
  g_image = NULL;
  g_idx[1] = g_idx[2] = g_idx[3] = 0;   /* initial position, x/y/z */
  g_idx[4] = g_idx[5] = 0;              /* intiial time position */
  g_zoom[1]=g_zoom[2]=g_zoom[3]=2.0;    /* zoom for x/y/z-directions */
  g_shift[1][1] = g_shift[1][2] =       /* shift within window */
                    g_shift[2][1] = g_shift[2][2] =
                                      g_shift[3][1] = g_shift[3][2] = 0;

  g_matrix[1] = niikmat_init(4,4);
  g_matrix[2] = niikmat_init(4,4);
  g_matrix[3] = niikmat_init(4,4);

  g_axis[1] = 'z';
  g_axis[2] = 'y';

  g_flip[1] = g_flip[2] = g_flip[3] = 1;

  g_spin_min = g_spin_max = NULL;
  g_edit_blue = 0.5;
  g_pen_size = 1;
  g_edit_fill_tol = 25;

  g_edit_max_undo = 20;
  g_edit_undo_idx = -1;

  g_maskimg = NULL;
  g_obj_list = NULL;
  g_obj_list_num = 0;

  g_idx[0] = 0;
  g_xy[0] = -1;
  g_coordinates = 1;

  this_spin_pct[0]=1; /* disabled by default ([0]>[1]) */
  this_spin_pct[1]=0;

  /* read options */

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "i:a:o:f:m:s:z:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"  NVK \n  current version = %i.%i.%i\n",NVK_MAJOR_VERSION,NVK_MINOR_VERSION,NVK_MICRO_VERSION);
      fprintf(stdout,"\nversion history\n");
      fprintf(stdout,"%s",*NVK_version);
      return 0;
    case 'f':
      g_flip[1] = g_flip[2] = g_flip[3] = 0;
      break;
    case 'm':
      if((dlist=niik_csstring_to_double_list(optarg,&list_num))==NULL) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }
      if(list_num!=2) {
        fprintf(stderr,"ERROR: Expected 2 arguments:%s\n",optarg);
        exit(9);
      }
      this_spin_min=dlist[0];
      this_spin_max=dlist[1];
      free(dlist);
      break;
    case 's':
      if((dlist=niik_csstring_to_double_list(optarg,&list_num))==NULL) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }
      if(list_num!=6) {
        fprintf(stderr,"ERROR: Expected 6 arguments:%s\n",optarg);
        exit(9);
      }
      g_shift[1][1] = dlist[0];
      g_shift[1][2] = dlist[1];
      g_shift[2][1] = dlist[2];
      g_shift[2][2] = dlist[3];
      g_shift[3][1] = dlist[4];
      g_shift[3][2] = dlist[5];
      free(dlist);
      break;
    case 'c':
      g_idx[0] = 1;
      break;
    case 'o':
      if((CSlist=niik_csstring_to_list(optarg,&list_num))==NULL) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }

      g_obj_list=(kobj **)realloc(g_obj_list,(g_obj_list_num+list_num)*sizeof(kobj *));
      if(verbose) fprintf(stdout,"-d1  reading obj(s)   # %i\n",list_num);

      for(n=0; n<list_num; n++) {
        if(verbose) fprintf(stdout,"-d1  reading object    %s \n",CSlist[n]);

        if((g_obj_list[g_obj_list_num+n]=off_kobj_read_offply(CSlist[n]))==NULL) {
          fprintf(stderr,"ERROR: off_kobj_read_off %s\n",CSlist[n]);
          exit(9);
        }

        if(verbose) fprintf(stdout,"-d1  object vfe %i %i %i\n",
                              g_obj_list[g_obj_list_num+n]->nvert,
                              g_obj_list[g_obj_list_num+n]->nface,
                              g_obj_list[g_obj_list_num+n]->nedge);

        if(g_obj_list[g_obj_list_num+n]->color==0) {
          if(verbose) fprintf(stdout,"-d1 add color (green)\n");
          if(!off_kobj_add_green_color(g_obj_list[g_obj_list_num+n])) {
            fprintf(stderr,"ERROR: off_kobj_add_green_color\n");
            exit(9);
          }

        } else {
          if(verbose) fprintf(stdout,"-d1 do not add color\n");
        }
        off_update_kobj_kface_pminmax(g_obj_list[g_obj_list_num+n]);
      } /* each obj */
      niik_free_list(CSlist,list_num);
      g_obj_list_num+=list_num;
      break;
    case 'p':
      this_spin_pct[0] = atof(optarg);
      if(this_spin_pct[0]>1) this_spin_pct[0]/=100.0;
      this_spin_pct[1] = 1.0 - this_spin_pct[0];
      break;
    case 'i':
      if((ilist=niik_csstring_to_int_list(optarg,&list_num))==NULL) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }
      if(list_num!=3) {
        fprintf(stderr,"ERROR: Expected 3 arguments:%s\n",optarg);
        exit(9);
      }

      g_idx[1] = ilist[0];
      g_idx[2] = ilist[1];
      g_idx[3] = ilist[2];
      free(ilist);
      break;
    case 'z':
      if((dlist=niik_csstring_to_double_list(optarg,&list_num))==NULL) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }
      if(list_num!=3) {
        fprintf(stderr,"ERROR: Expected 3 arguments:%s\n",optarg);
        exit(9);
      }
      g_zoom[1] = dlist[0];
      g_zoom[2] = dlist[1];
      g_zoom[3] = dlist[2];
      free(dlist);
      break;
    case 'U':
      fprintf(stdout,"%s\n",*NVK_runtime_usage);
      return 1;
    case '?':
    default:
      usage();
      fprintf(stdout,"%s\n",*NVK_usage);
      return 1;
    }
  }

  if((argc - optind)<1) {
    usage();
    fprintf(stdout,"%s\n",*NVK_usage);
    return 1;
  }

  if (getcwd(curr_folder, sizeof(curr_folder)) == NULL) {
    fprintf(stderr,"ERROR: getcwd \n");
    exit(9);
  }

  /* overlay object */
  if(g_obj_list_num>0) {
    if(verbose) fprintf(stdout,"-d1 bounding box for %i object(s)\n",g_obj_list_num);
    g_off_bb = off_bbox_init(7,260);
    if(!off_create_bbox_from_multiple_kobj(g_off_bb, g_obj_list, g_obj_list_num)) {
      fprintf(stderr,"ERROR: off_create_bbox_from_multiple_kobj\n");
      exit(9);
    }
  } /* bounding box */


  /******************************
   *  initialize and check gtk
   ******************************/
  if(verbose>1) fprintf(stdout,"-d2 initialize and check gtk\n");
  gtk_init(&argc,&argv);
  if((gtk_init_check(&argc,&argv))==FALSE) {
    fprintf(stderr,"ERORR: gtk_init_check failure\n");
    fprintf(stderr,"ERORR: windowing system has NOT been successfully initialized\n");
    exit(9);
  }

  /***************************
   * create the main window
   ***************************/

  if(verbose>1) fprintf(stdout,"-d2 create the main window (g_main_window)\n");
  g_main_window = gtk_window_new (GTK_WINDOW_TOPLEVEL);


  sprintf(iconfile,"%s/data/NVK/icontest1.png",NIIKDIR);
  gtk_window_set_icon(GTK_WINDOW(g_main_window), add_create_pixbuf(iconfile));


  /* gtk_window_set_title(GTK_WINDOW(g_main_window),argv[1]);*/

  gtk_window_set_resizable(GTK_WINDOW(g_main_window),FALSE);


  /* not working very well
   * -would like to fix this so that user cannot resize the window
   geomhints.min_width  = geomhints.min_height = 0;
   gtk_window_set_geometry_hints(GTK_WINDOW(g_main_window),g_main_window,
   &geomhints,GDK_HINT_MIN_SIZE );
  */

  /***************************************
   * create the image space
   **************************************/

  if(verbose>1) fprintf(stdout,"-d2 image space x2 (g_pixbuf, g_gtkimg, event_box)\n");
  for(n=1; n<3; n++) {
    if(verbose>2) fprintf(stdout,"-d3 individual image space %i\n",n);
    g_pixbuf[n] = gdk_pixbuf_new(GDK_COLORSPACE_RGB,FALSE,8,g_window_width,g_window_height);
    g_gtkimg[n] = gtk_image_new_from_pixbuf(g_pixbuf[n]);
    img_align[n]   =gtk_alignment_new (0,0,0,0);
    event_box[n]=gtk_event_box_new();
    gtk_container_add (GTK_CONTAINER (event_box[n]),img_align[n]);
    gtk_container_add (GTK_CONTAINER (img_align[n]),g_gtkimg[n]);
  }




  /**************************
   *   MENU BAR
   **************************/


  if(verbose>1) fprintf(stdout,"-d2 menu lists\n");
  action_group = gtk_action_group_new ("TestActions");
  gtk_action_group_set_translation_domain (action_group, "blah");

  uiman = gtk_ui_manager_new ();

  if(verbose>1) fprintf(stdout,"-d2 group actions\n");
  gtk_action_group_add_actions (action_group, entries, n_entries, NULL);
  gtk_ui_manager_insert_action_group (uiman, action_group, 0);



  /* Read in the UI from our XML file */
  if(verbose>1) fprintf(stdout,"-d2 reading UI from xml file\n");
  error = NULL;
  sprintf(uixmlfile,"%s/data/NVK/nifti1_kunio_viewer_file_bar_ui.xml",NIIKDIR);
  gtk_ui_manager_add_ui_from_file (uiman, uixmlfile, &error);

  if (error) {
    g_message ("building menus failed: %s", error->message);
    g_error_free (error);
  }



  /*********************************
   * set up the labels
   **********************************/

  g_lbl_voxel = gtk_label_new("     ");
  g_lbl_gray  = gtk_label_new("     ");

  lbledit[0] = gtk_label_new("Blue % = ");
  lbledit[1] = gtk_label_new("    Pen Size = ");
  lbledit[2] = gtk_label_new("    Fill Tol = ");



  /*********************************
   * set up the spin buttons
   **********************************/


  if(verbose) fprintf(stdout,"-d1 adding a spin button 1\n");
  g_spin_min = gtk_spin_button_new_with_range ( 0,1,1 );
  g_spin_max = gtk_spin_button_new_with_range ( 0,1,1 );
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_min),0);
  gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_max),1);

  g_spin_xdim = g_spin_dim[1] = gtk_spin_button_new_with_range ( 0,1,1 );
  g_spin_ydim = g_spin_dim[2] = gtk_spin_button_new_with_range ( 0,1,1 );
  g_spin_zdim = g_spin_dim[3] = gtk_spin_button_new_with_range ( 0,1,1 );
  g_spin_tdim = g_spin_dim[4] = gtk_spin_button_new_with_range ( 0,1,1 );
  g_spin_udim = g_spin_dim[5] = gtk_spin_button_new_with_range ( 0,1,1 );

  gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(g_spin_min),GTK_UPDATE_IF_VALID);
  gtk_spin_button_set_update_policy(GTK_SPIN_BUTTON(g_spin_max),GTK_UPDATE_IF_VALID);


  /*
   * EDIT SPIN BUTTONS
   *
   *   g_spin_pen_size   g_spin_edit_blue    g_spin_fill_tol
   *
   */
  g_spin_pen_size = gtk_spin_button_new_with_range ( 0,20,0.1 );
  gtk_spin_button_set_digits ( GTK_SPIN_BUTTON(g_spin_pen_size), 1);
  gtk_spin_button_set_value( GTK_SPIN_BUTTON(g_spin_pen_size), g_pen_size);

  g_spin_edit_blue = gtk_spin_button_new_with_range ( 0,100,10 );
  gtk_spin_button_set_digits ( GTK_SPIN_BUTTON(g_spin_edit_blue), 0);
  gtk_spin_button_set_value( GTK_SPIN_BUTTON(g_spin_edit_blue), g_edit_blue*100);

  g_spin_fill_tol = gtk_spin_button_new_with_range ( 0,1000,1 );
  gtk_spin_button_set_digits ( GTK_SPIN_BUTTON(g_spin_fill_tol), 0);
  gtk_spin_button_set_value( GTK_SPIN_BUTTON(g_spin_fill_tol), g_edit_fill_tol);


  g_lbl_dim[1] = gtk_label_new("  |  x");
  g_lbl_dim[2] = gtk_label_new("y");
  g_lbl_dim[3] = gtk_label_new("z");
  g_lbl_dim[4] = gtk_label_new("t");
  g_lbl_dim[5] = gtk_label_new("u");



  /*********************************
   * Pack up our objects:
   * actions -> action_group
   * action_group -> menu_manager
   **********************************/

  mainv=gtk_vbox_new(FALSE,0);
  imagh=gtk_hbox_new(FALSE,0);
  labelh=gtk_hbox_new(FALSE,0);
  edith =gtk_hbox_new(FALSE,0);


  /* Get the menubar and put it in the vertical packing box */

  if(verbose>1) fprintf(stdout,"-d2 ui manager\n");
  if((menubar = gtk_ui_manager_get_widget (uiman, "/MainMenu"))==NULL) {
    fprintf(stderr,"ERROR: gtk_ui_manager_get_widget \n");
    exit(9);
  }


  /* Make sure that the accelerators work */

  if(verbose>1) fprintf(stdout,"-d2 accelerator work\n");
  gtk_window_add_accel_group (GTK_WINDOW (g_main_window),
                              gtk_ui_manager_get_accel_group (uiman));


  /* put things on the main window */

  if(verbose>1) fprintf(stdout,"-d2 pack things in the main window\n");
  gtk_box_pack_start (GTK_BOX (mainv), menubar, FALSE, FALSE, 0);

  if(verbose>2) fprintf(stdout,"-d3 container add (main window) main vbox\n");
  gtk_container_add (GTK_CONTAINER(g_main_window),mainv);

  if(verbose>2) fprintf(stdout,"-d3 eventbox to imgh\n");
  gtk_box_pack_start(GTK_BOX(imagh), event_box[1], FALSE,TRUE,1);
  gtk_box_pack_start(GTK_BOX(imagh), event_box[2], FALSE,TRUE,1);

  if(verbose>2) fprintf(stdout,"-d3 imagh to mainv\n");
  gtk_box_pack_start(GTK_BOX(mainv),  imagh,      FALSE,TRUE,1);  /* horizontal box LEFT=img, RIGHT=options */
  gtk_box_pack_start(GTK_BOX(labelh), g_spin_min, FALSE,TRUE,1);  /* spin button for min intensity */
  gtk_box_pack_start(GTK_BOX(labelh), g_spin_max, FALSE,TRUE,1);  /* spin button for max intensity */
  if(verbose>2) fprintf(stdout,"-d3 labels to labelh\n");
  gtk_box_pack_start(GTK_BOX(labelh), g_lbl_voxel, FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_spin_udim,  FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_dim[5], FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_spin_tdim,  FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_dim[4], FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_spin_zdim,  FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_dim[3], FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_spin_ydim,  FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_dim[2], FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_spin_xdim,  FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_dim[1], FALSE,TRUE,1);
  gtk_box_pack_end(GTK_BOX(labelh), g_lbl_gray,   FALSE,TRUE,1);
  gtk_box_pack_start(GTK_BOX(mainv),labelh,       FALSE,TRUE,1);  /* horizontal box for labels */

  if(argc==3) {
    gtk_box_pack_start(GTK_BOX(mainv), edith, FALSE,TRUE,1);  /* editing tools */
    gtk_box_pack_start(GTK_BOX(edith),lbledit[0], FALSE,TRUE,1);
    gtk_box_pack_start(GTK_BOX(edith),g_spin_edit_blue, FALSE,TRUE,1);
    gtk_box_pack_start(GTK_BOX(edith),lbledit[1], FALSE,TRUE,1);
    gtk_box_pack_start(GTK_BOX(edith),g_spin_pen_size, FALSE,TRUE,1);
    gtk_box_pack_start(GTK_BOX(edith),lbledit[2], FALSE,TRUE,1);
    gtk_box_pack_start(GTK_BOX(edith),g_spin_fill_tol, FALSE,TRUE,1);
  }


  /* Connect up important signals */

  if(verbose>1) fprintf(stdout,"-d2 connect signals\n");

  gtk_widget_set_events    (event_box[1],
                            GDK_BUTTON_PRESS_MASK | GDK_POINTER_MOTION_MASK | GDK_BUTTON_RELEASE_MASK | GDK_SCROLL_MASK |
                            GDK_KEY_PRESS_MASK | GDK_KEY_RELEASE_MASK );
  gtk_widget_set_events    (event_box[2],
                            GDK_BUTTON_PRESS_MASK | GDK_POINTER_MOTION_MASK | GDK_BUTTON_RELEASE_MASK | GDK_SCROLL_MASK |
                            GDK_KEY_PRESS_MASK | GDK_KEY_RELEASE_MASK);

  g_signal_connect(G_OBJECT(event_box[1]), "motion_notify_event", G_CALLBACK (motion_notify_event_image1), NULL);
  g_signal_connect(G_OBJECT(event_box[2]), "motion_notify_event", G_CALLBACK (motion_notify_event_image2), NULL);
  g_signal_connect(G_OBJECT(event_box[1]), "scroll_event",        G_CALLBACK (scroll_event_image1),        NULL);
  g_signal_connect(G_OBJECT(event_box[2]), "scroll_event",        G_CALLBACK (scroll_event_image2),        NULL);
  g_signal_connect(G_OBJECT(event_box[1]), "button_press_event",  G_CALLBACK (button_press_event_image1),  NULL);
  g_signal_connect(G_OBJECT(event_box[2]), "button_press_event",  G_CALLBACK (button_press_event_image2),  NULL);
  g_signal_connect(G_OBJECT(event_box[1]), "button_release_event",G_CALLBACK (button_release_event_image1),NULL);
  g_signal_connect(G_OBJECT(event_box[2]), "button_release_event",G_CALLBACK (button_release_event_image2),NULL);

  g_signal_connect(G_OBJECT(g_spin_min),"value_changed",G_CALLBACK(NVK_change_spin_min),NULL);
  g_signal_connect(G_OBJECT(g_spin_max),"value_changed",G_CALLBACK(NVK_change_spin_max),NULL);

  g_signal_connect(G_OBJECT(g_spin_xdim),"value_changed",G_CALLBACK(NVK_change_spin_xdim),NULL);
  g_signal_connect(G_OBJECT(g_spin_ydim),"value_changed",G_CALLBACK(NVK_change_spin_ydim),NULL);
  g_signal_connect(G_OBJECT(g_spin_zdim),"value_changed",G_CALLBACK(NVK_change_spin_zdim),NULL);
  g_signal_connect(G_OBJECT(g_spin_tdim),"value_changed",G_CALLBACK(NVK_change_spin_tdim),NULL);
  g_signal_connect(G_OBJECT(g_spin_udim),"value_changed",G_CALLBACK(NVK_change_spin_udim),NULL);

  g_signal_connect(G_OBJECT(g_spin_pen_size ),"value_changed",G_CALLBACK(NVK_change_edit_pen_size),NULL);
  g_signal_connect(G_OBJECT(g_spin_edit_blue),"value_changed",G_CALLBACK(NVK_change_edit_blue),NULL);
  g_signal_connect(G_OBJECT(g_spin_fill_tol ),"value_changed",G_CALLBACK(NVK_change_edit_fill_tol),NULL);

  g_signal_connect (g_main_window, "delete-event", G_CALLBACK (NVK_ask_quit), NULL);
  // g_signal_connect (g_main_window, "destroy",      G_CALLBACK (NVK_ask_quit), NULL);

  g_timeout_add(60000*3, (GSourceFunc) NVK_automatic_save_mask_file, (gpointer) g_main_window);
  // g_timeout_add(10000, (GSourceFunc) NVK_automatic_save_mask_file, (gpointer) g_main_window);

  gtk_widget_show_all (g_main_window);



  /******************************
   *  read input image
   ******************************/

  if(verbose) fprintf(stdout,"-d1 current folder %s\n",curr_folder);

  if((argc - optind)==1 || (argc - optind)==2) {
    if(verbose) fprintf(stdout,"-d1 reading %s\n",argv[optind]);

    if(NVK_read_image(argv[optind])==0) {
      fprintf(stderr,"ERROR: could not open file %s\n",argv[optind]);
      exit(9);
    }

    if(g_image==NULL) {
      fprintf(stderr,"ERROR: could not read %s\n",argv[optind]);
      exit(9);
    }


    if(g_idx[0]>0) {
      g_idx[1]=g_image->nx/2;
      g_idx[2]=g_image->ny/2;
      g_idx[3]=g_image->nz/2;
      g_idx[0]=0;
    }


    if (this_spin_pct[0]<this_spin_pct[1]) {
      this_spin_min = niik_image_get_percentile(g_image,g_maskimg,this_spin_pct[0]);
      this_spin_max = niik_image_get_percentile(g_image,g_maskimg,this_spin_pct[1]);
      if(verbose>=2) fprintf(stdout,"[%s] using percentile = %9.2f -> %12.3f %12.1f\n",fcname,this_spin_pct[0],this_spin_min,this_spin_max);
    }
    if(this_spin_min < 0 && this_spin_max < 0) {
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_min),g_image_data[niik_get_min_index_float_vector(g_image_data,g_image->nvox)]);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_max),g_image_data[niik_get_max_index_float_vector(g_image_data,g_image->nvox)]);
    } else {
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_min),this_spin_min);
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(g_spin_max),this_spin_max);
    }


    if((argc - optind)==2) {
      /* try to read it */
      if(verbose) fprintf(stdout,"-d1 try reading a file %s\n",argv[optind+1]);
      if((fp=fopen(argv[optind+1],"rb"))==NULL) {
        if(verbose) fprintf(stdout,"-d1 could not read %s, copy image as maskimg\n",argv[optind+1]);

        /* copy image from g_image
         * but it's missing data (as it is stored only g_image_data
         * for saving memory */
        if((g_maskimg=niik_image_copy(g_image))==NULL) {
          fprintf(stderr,"ERROR: nifti_copy_nim_info \n");
          exit(9);
        }

        /* we have to allocate a new memory for g_maskimg
         * this can be a unsigned char image */
        g_maskimg->datatype = NIFTI_TYPE_UINT8;
        nifti_datatype_sizes( NIFTI_TYPE_UINT8, &g_maskimg->nbyper, &g_maskimg->swapsize );
        if(verbose) fprintf(stdout,"-d1 clear image\n");
        g_mask_data = (unsigned char *)calloc(g_maskimg->nvox,sizeof(char));
        g_maskimg -> data = g_mask_data;
        free(g_maskimg -> fname);
        /* update the filename otherwise, the input image is overwritten */
        g_maskimg -> fname = (char *)calloc(strlen(argv[optind+1]),sizeof(char));
        g_maskimg -> iname = (char *)calloc(strlen(argv[optind+1]),sizeof(char));
        sprintf(g_maskimg -> fname,"%s",argv[optind+1]);
        sprintf(g_maskimg -> iname,"%s",argv[optind+1]);
        g_maskimg ->nifti_type=1;

        fprintf(stdout,"\t  making a file: %s -> %s, %i\n",argv[optind+1],g_maskimg->fname,g_maskimg->nifti_type);
      } /* file does not exist -> make a mask image */

      else { /* file exists, read it */
        fclose(fp);
        if((g_maskimg=niik_image_read(argv[optind+1]))==NULL) {
          fprintf(stderr,"ERROR: niik_image_read %s\n",argv[optind+1]);
          exit(9);
        }

        /* convert to uint8 b/c it's a mask */
        if(verbose) fprintf(stdout,"-d1 convert maskimg to byte\n");
        if(!niik_image_type_convert(g_maskimg,NIFTI_TYPE_UINT8)) {
          fprintf(stderr,"ERROR: niik_image_type_convert\n");
          exit(9);
        }

        if(verbose) fprintf(stdout,"-d1 get the mask data\n");
        g_mask_data = g_maskimg -> data;
      }

      if(verbose) fprintf(stdout,"-d1 cmp image dimension\n");
      /*if(niik_image_cmp_dim(g_image,g_maskimg)!=0){
      fprintf(stderr,"ERROR: niik_iamge_cmp_dim\n");
      exit(9); }*/
      if(g_image->nx != g_maskimg->nx) {
        fprintf(stderr,"ERROR: different xdim %i %i\n",g_image->nx,g_maskimg->nx);
        exit(9);
      }
      if(g_image->ny != g_maskimg->ny) {
        fprintf(stderr,"ERROR: different ydim %i %i\n",g_image->ny,g_maskimg->ny);
        exit(9);
      }
      if(g_image->nz != g_maskimg->nz) {
        fprintf(stderr,"ERROR: different zdim %i %i\n",g_image->nz,g_maskimg->nz);
        exit(9);
      }


      if(flag_check_voxelsize) {
        if(verbose) fprintf(stdout,"-d1 check pixdim\n");
        if(floor(g_image->pixdim[1]*1e5) != floor(g_maskimg->pixdim[1]*1e5)) {
          fprintf(stderr,"ERROR: wrong x pixdim %12.7f %12.7f\n",g_image->pixdim[1],g_maskimg->pixdim[1]);
          fprintf(stderr,"       use -no-voxel-check to skip this test\n");
          exit(9);
        }
        if(floor(g_image->pixdim[2]*1e5) != floor(g_maskimg->pixdim[2]*1e5)) {
          fprintf(stderr,"ERROR: wrong y pixdim %12.7f %12.7f\n",g_image->pixdim[2],g_maskimg->pixdim[2]);
          fprintf(stderr,"       use -no-voxel-check to skip this test\n");
          exit(9);
        }
        if(floor(g_image->pixdim[3]*1e5) != floor(g_maskimg->pixdim[3]*1e5)) {
          fprintf(stderr,"ERROR: wrong z pixdim %12.7f %12.7f\n",g_image->pixdim[3],g_maskimg->pixdim[3]);
          fprintf(stderr,"       use -no-voxel-check to skip this test\n");
          exit(9);
        }
      } else {
        if(verbose) fprintf(stdout,"-d1 don't check pixdim\n");
      }

      if(verbose) fprintf(stdout,"-d1 create undo list\n");
      m = NIIK_IMAX(g_maskimg->dim[1],NIIK_IMAX(g_maskimg->dim[2],g_maskimg->dim[3])) + 5;
      g_edit_undo_list=(int **)calloc(g_edit_max_undo,sizeof(int *));
      for(n=0; n<g_edit_max_undo; n++) {
        g_edit_undo_list[n] = (int *)calloc(m*m,sizeof(int));
        g_edit_undo_list[n][0] = 2;
        g_edit_undo_list[n][1] = 0;
      }

      if(verbose) fprintf(stdout,"-d1 set the title\n");
      sprintf(tmpstr,"%s, %s",argv[optind],argv[optind+1]);
      gtk_window_set_title(GTK_WINDOW(g_main_window),tmpstr);

      if(g_automatic_mask_save>0) {
        NVK_automatic_save_mask_file(g_main_window);
      } /* automatic save */

    } /* using g_maskimg */


    if(verbose) fprintf(stdout,"-d1 finish preparing the image(s)\n");
  } /* commnad-line arguments */

  else {
    sprintf(g_curr_filename,"%s/",curr_folder);
  }

  if(verbose>1) fprintf(stdout,"-d2 gtk main starts here\n");
  gtk_main ();

  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/