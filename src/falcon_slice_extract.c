/* FILENAME:    flacon_slice_extract.c
* DESCRIPTION:  Extract slice for visualization
* AUTHOR:       Vladimir S. Fonov
* DATE:         Oct 12, 2018
*/

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include "falcon.h"
#include "falcon_cortex.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (1)


#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


void prog_version() {
    fprintf(stdout,"  falcon_slice_extract history\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  0.0.1  Sep 14, 2016, Vladimir S. FONOV - moved from niikmath\n");
    fprintf(stdout,"\n");
}

void prog_usage() {
    fprintf(stdout,"falcon_slice_extract\n");
    fprintf(stdout,"  usage: [options] <input.mnc/.nii> <out pattern have to have %%s for direction and %%d for slice no> \n");
    fprintf(stdout,"\n");
    fprintf(stdout,"\n  optional usage:\n");
    fprintf(stdout,"  -d/--debug            Added debugging output\n");
    fprintf(stdout,"  -o/--obj  <o>[,o,o]   Add object(s) to render\n");
    fprintf(stdout,"  -x/-y/-z <n1>[,n2,n2] Slices list, -1 to extract all slices for that axis\n");
    fprintf(stdout,"  -X/--XYZ <vx,vy,vz>   Resample to this step sizes\n");
    fprintf(stdout,"  -S/--scale <d>        Scaling for images, default (1)\n");
    fprintf(stdout,"  --min <f> image minimum default auto\n");
    fprintf(stdout,"  --max <f> image maximum default auto\n");
    fprintf(stdout,"  -P/--pct <f> use image percentile as max\n");
}

int main(int argc,char *argv[],char *envp[]) {
    int clobber=0;
    int verbose=niik_verbose();
    nifti_image *img=NULL;
    kobj **obj=NULL;
    const char *in_scan=NULL;
    const char *out_prefix=NULL;
    const char *fcname="flacon_slice_extract";
    int i,n;
    int flag_debug=0;
    char **CSlist=NULL;
    char **objlist=NULL;
    int n_obj=0;
    int *xslices=NULL,*yslices=NULL,*zslices=NULL;
    char fname[4096];
    double xyz[4]= {0,0,0,0};
    double scale=1.0;

    int interp=NIIK_INTERP_LINEAR;
    double    imin=NIIKMAX,
              imax=NIIKMAX,
              percent=NIIKMAX;
    static bbox *off_bb=NULL;               /* input off object's bounding box */

    extract_slice_info slices_x;
    extract_slice_info slices_y;
    extract_slice_info slices_z;

    struct option long_options[] = {
        {"clobber", no_argument,  &clobber, 1},
        {"debug",   no_argument,  &flag_debug, 1},
        {"version", no_argument, 0, 'v'},
        {"help",    no_argument, 0, 'h'},
        {"usage",   no_argument, 0, 'U'},
        {"XYZ",     required_argument, 0, 'X'},
        {"obj",     required_argument, 0, 'o'},
        {"scale",   required_argument, 0, 'S'},
        {"min",     required_argument, 0, 'I'},
        {"max",     required_argument, 0, 'A'},
        {"pct",     required_argument, 0, 'P'},
        {0, 0, 0, 0}
    };

    memset(&slices_x,0,sizeof(extract_slice_info));
    memset(&slices_y,0,sizeof(extract_slice_info));
    memset(&slices_z,0,sizeof(extract_slice_info));

    for (;;) {
        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long (argc, argv, "o:x:y:z:X:I:A:S:P:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 0:
            break;
        case 'S':
            scale=atof(optarg);
            break;
        case 'o':
            NIIK_EXIT(((objlist=niik_csstring_to_list(optarg,&n_obj))==NULL),fcname,"niik_csstring_to_list",9);
            break;
        case 'P':
            percent=atof(optarg);
            break;
        case 'X':
            NIIK_EXIT(((CSlist=niik_csstring_to_list(optarg,&n))==NULL),fcname,"niik_csstring_to_list",9);
            if(n!=3) {
                fprintf(stderr,"[niikmath] ERROR: xyz should have 3 numnbers\n");
                exit(1);
            }
            xyz[0]=1.0;
            for(i=0; i<3; i++) {
                xyz[i+1]=atof(CSlist[i]);
            }
            niik_free_list(CSlist,n);
            break;
        case 'x':
            NIIK_EXIT(((CSlist=niik_csstring_to_list(optarg,&n))==NULL),fcname,"niik_csstring_to_list",9);
            slices_x.slice=(int *)calloc(n,sizeof(int));
            slices_x.slice_num=n;
            for(i=0; i<n; i++) {
                slices_x.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
            break;
        case 'y':
            NIIK_EXIT(((CSlist=niik_csstring_to_list(optarg,&n))==NULL),fcname,"niik_csstring_to_list",9);
            slices_y.slice=(int *)calloc(n,sizeof(int));
            slices_y.slice_num=n;
            for(i=0; i<n; i++) {
                slices_y.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
            break;
        case 'z':
            NIIK_EXIT(((CSlist=niik_csstring_to_list(optarg,&n))==NULL),fcname,"niik_csstring_to_list",9);
            slices_z.slice=(int *)calloc(n,sizeof(int));
            slices_z.slice_num=n;
            for(i=0; i<n; i++) {
                slices_z.slice[i]=atoi(CSlist[i]);
                free(CSlist[i]);
            }
            free(CSlist);
            break;
        case 'I':
            imin=atof(optarg);
            break;
        case 'A':
            imax=atof(optarg);
            break;
        case 'v':
            prog_version();
            return 0;
        case 'h':
        case 'U':
        case '?':
        default:
            prog_usage ();
            return 1;
        }
    }

    if((argc - optind)<2) {
        prog_usage();
        return 1;
    }

    in_scan    = argv[optind];
    out_prefix = argv[optind+1];
    slices_x.fpattern=slices_y.fpattern=slices_z.fpattern=out_prefix;
    slices_x.compression=slices_y.compression=slices_z.compression=1;

    if(verbose>0) niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
    if(verbose>0) niik_fc_display(fcname,1);
    if( (img=niik_image_read(in_scan))==NULL) {
        fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,in_scan);
        exit(1);
    }

#ifdef HAVE_OPENMP
    if(verbose>0) fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif
    /*from niikmath*/
    if(n_obj>0) {
        obj=(kobj **)calloc(sizeof(kobj *), n_obj);
        for(i=0;i<n_obj;i++) {
            if((obj[i]=off_kobj_read_offply(objlist[i]))==NULL) {
                fprintf(stderr,"[%s] ERROR: off_kobj_read_off\n",fcname);
                exit(1);
            }
            off_update_kobj_kface_pminmax(obj[i]);
        }
    }

    if(xyz[0]>0.0) { /* if resampling */
        if(verbose>0) fprintf(stdout,"[niikmath] resampling        %7.3f %7.3f %7.3f\n",xyz[1],xyz[2],xyz[3]);
        NIIK_EXIT((!niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],-1,-1,-1,interp)),fcname,"niik_image_resample_3d_update",9);
    } /* resampling */

    if(n_obj>0) {
        off_bb = off_bbox_init(7,260);
        
        if(!off_create_bbox_from_multiple_kobj(off_bb, obj, n_obj)) {
           fprintf(stderr,"ERROR: off_create_bbox_from_multiple_kobj\n");
           exit(9);
        }
    }

    if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,NULL);
    if(niik_check_double_problem(imax)) {
        if(niik_check_double_problem(percent)) {
            imax=niik_image_get_max(img,NULL);
        } else {
            imax=niik_image_get_percentile(img,NULL,percent);
        }
    }

    slices_x.dmin =slices_y.dmin =slices_z.dmin =imin;
    slices_x.dmax =slices_y.dmax =slices_z.dmax =imax;
    slices_x.scale=slices_y.scale=slices_z.scale=scale;

    slices_x.slice_dir=0;
    slices_y.slice_dir=1;
    slices_z.slice_dir=2;

    /*extract all slices*/
    if(slices_x.slice_num) {
        NIIK_EXIT((!niik_image_write_slices_obj(&slices_x,img,off_bb)),fcname,"niik_image_write_slices_obj",9);
    } /* x-slices */
    if(slices_y.slice_num) {
        NIIK_EXIT((!niik_image_write_slices_obj(&slices_y,img,off_bb)),fcname,"niik_image_write_slices_obj",9);
    } /* y-slices */
    if(slices_z.slice_num) {
        NIIK_EXIT((!niik_image_write_slices_obj(&slices_z,img,off_bb)),fcname,"niik_image_write_slices_obj",9);
    } /* z-slices */
    if(objlist) niik_free_list(objlist,n_obj);
    
    if(n_obj>0) {
        for(i=0;i<n_obj;i++)
            obj[i]=off_kobj_free(obj[i]);
        free(obj);
    }
    
    img=niik_image_free(img);
    if(xslices) free(xslices);
    if(yslices) free(yslices);
    if(zslices) free(zslices);
    if(off_bb)  off_bb=off_bbox_free(off_bb);
    if(verbose>0) niik_fc_display(fcname,0);
    return 0;
} /* main */

/*
kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
