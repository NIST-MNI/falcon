/* Filename:     falcon_slice.c
 * Description:  slice extraction routine
 * Author:       Vladimir S. FONOV
 * Date:         October 16, 2018
 */

#include "falcon.h"
#include <stdlib.h>
#include <unistd.h>
#include "stb_image_write.h"

#include "falcon_cortex.h"


int niik_image_write_slices_obj(extract_slice_info *info,
        nifti_image *img, bbox *off_bb)
{
    //TIFF *tifimg;
    int   s;
    double d,dran;
    unsigned char *bout;
    int verbose=0;
    int stride_i,stride_j,stride_vol,stride_k;
    int ni,nj,nk;
    const char *fcname=__func__;
    int round_scale=ceil(info->scale);
    const char *dirs[] = { "x", "y", "z" };
    int slice_num;

    NIIK_RET0((img==NULL),fcname,"img is null");

    if(info->slice_dir==0) { /*X*/
            stride_i=img->nx; /*Y*/
            stride_j=img->nx*img->ny; /*Z*/
            stride_k=1;
            ni=img->ny;
            nj=img->nz;
            nk=img->nx;
    } else if(info->slice_dir==1) { /*Y*/
            stride_i=1; /*X*/
            stride_j=img->nx*img->ny; /*Z*/
            stride_k=img->nx;
            ni=img->nx;
            nj=img->nz;
            nk=img->ny;
    } else if(info->slice_dir==2) { /*Z*/
            stride_i=1; /*X*/
            stride_j=img->nx; /*Y*/
            stride_k=img->nx*img->ny;
            ni=img->nx;
            nj=img->ny;
            nk=img->nz;
    } else {
            fprintf(stderr,"%s: Error unsupported slice direction:%d\n",fcname,info->slice_dir);
            return 0;
    }

    /*volume*/
    stride_vol = img->nx * img->ny * img->nz * img->nt * img->nu; /*to sckip oper vector dimensions*/

    dran  = info->dmax - info->dmin;
    bout  = (unsigned char *) calloc(3*ni*nj*round_scale*round_scale, sizeof(unsigned char)); /*RGB*size*/

    slice_num=info->slice_num;
    if(info->slice[0]<0)
    {
        if(info->slice_dir==0)       slice_num=img->nx;
        else if (info->slice_dir==1) slice_num=img->ny;
        else if (info->slice_dir==2) slice_num=img->nz;
    }
    
    for(s=0; s<slice_num;s++) 
    {
        int i,j;
        char fname[4096];
        int slice;
        
        if(info->slice[0]>=0)
            slice=info->slice[s];
        else
            slice=s;

        if(slice>=nk) slice=nk-1;
        
        printf("Processing dir:%d slice:%d\n",info->slice_dir, slice);

        sprintf(fname, info->fpattern, dirs[info->slice_dir], slice);

        /*create tiff file*/
        #if 0
        if((tifimg = TIFFOpen(fname, "w")) == NULL) {
            fprintf(stderr,"%s: ERROR: could not open for writing: %s\n",fcname, fname);
            return 0;
        }
        TIFFSetField( tifimg, TIFFTAG_IMAGEWIDTH,  ni*round_scale );
        TIFFSetField( tifimg, TIFFTAG_IMAGELENGTH, nj*round_scale );
        TIFFSetField( tifimg, TIFFTAG_BITSPERSAMPLE,   8);
        TIFFSetField( tifimg, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField( tifimg, TIFFTAG_ROWSPERSTRIP, nj*round_scale);

        if(info->compression)
            TIFFSetField( tifimg, TIFFTAG_COMPRESSION,  COMPRESSION_DEFLATE);
        else
            TIFFSetField( tifimg, TIFFTAG_COMPRESSION,  COMPRESSION_NONE);

        TIFFSetField(tifimg, TIFFTAG_PHOTOMETRIC,  PHOTOMETRIC_RGB);

        TIFFSetField( tifimg, TIFFTAG_FILLORDER,    FILLORDER_LSB2MSB);
        TIFFSetField( tifimg, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField( tifimg, TIFFTAG_XRESOLUTION, 150.0);
        TIFFSetField( tifimg, TIFFTAG_YRESOLUTION, 150.0);
        TIFFSetField( tifimg, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
        TIFFSetField( tifimg, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        #endif 

        if(img->nv==3) {
            for(j=0; j<nj*round_scale; j++) {
                int _j=floor((double)j/info->scale);
                int m = (nj*round_scale-j-1)*ni*round_scale*3;

                for(i=0; i<ni*round_scale; i++) {
                    int _i=floor((double)i/info->scale);
                    
                    d = niik_image_get_voxel(img, slice*stride_k + _i*stride_i + _j*stride_j );
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                    d = niik_image_get_voxel(img, slice*stride_k + _i*stride_i + _j*stride_j + stride_vol);
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                    d = niik_image_get_voxel(img, slice*stride_k + _i*stride_i + _j*stride_j + stride_vol*2);
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                }
            } /* each voxel on a slice */
        } else {
            for(j=0; j<nj*round_scale; j++) {
                int _j=floor((double)j/info->scale);
                int m = (nj*round_scale-j-1)*ni*round_scale*3;
                for(i=0; i<ni*round_scale; i++) {
                    int _i=floor((double)i/info->scale);

                    d = niik_image_get_voxel(img, slice*stride_k + _i*stride_i + _j*stride_j);
                    /*gray-scale*/
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                    bout[m++] = (unsigned char)NIIK_DMINMAX(255.0*(d-info->dmin)/dran,0,255);
                }
            } /* each voxel on a slice */
        }

        if(off_bb)
        {
            int bxmin,bxmax,bymin,bymax,bzmin,bzmax;
            int i,j,k;

            niikpt
                ppi[3],              /* plane triangle in image space */
                ppw[3],              /* plane triangle in world space */
                qq[2];               /* intersection points  in image space*/

            int stride_x=1,stride_y=img->nx,stride_z=img->nx*img->ny;

            ppi[0]=ppi[1]=ppi[2]=niikpt_zero();

            if(info->slice_dir==0) { /*X*/
                    ppi[0].x = ppi[1].x = ppi[2].x = slice;
                    ppi[1].y = ppi[2].z = 1e6;
            } else if(info->slice_dir==1) { /*Y*/
                    ppi[0].y = ppi[1].y = ppi[2].y = slice;
                    ppi[1].x = ppi[2].z = 1e6;
            } else { /*Z*/
                    ppi[0].z = ppi[1].z = ppi[2].z = slice;
                    ppi[1].x = ppi[2].y = 1e6;
            }

            niik_ijk_to_world(img,&ppi[0],&ppw[0]);
            niik_ijk_to_world(img,&ppi[1],&ppw[1]);
            niik_ijk_to_world(img,&ppi[2],&ppw[2]);

            bxmin=bymin=bzmin=0;
            bxmax=off_bb->xdim;
            bymax=off_bb->ydim;
            bzmax=off_bb->zdim;

            /* go thru the bounding boxes -only in 2 dimensions */
            if(info->slice_dir==0) { /*X*/
                    bxmin = floor( (ppw[0].x - off_bb->origin.x) / off_bb->delta - 0.5);
                    bxmax = floor( (ppw[0].x - off_bb->origin.x) / off_bb->delta + 0.5);
            } else if(info->slice_dir==1) { /*Y*/
                    bymin = floor( (ppw[0].y - off_bb->origin.y) / off_bb->delta - 0.5);
                    bymax = floor( (ppw[0].y - off_bb->origin.y) / off_bb->delta + 0.5);
            } else { /*Z*/
                    bzmin = floor( (ppw[0].z - off_bb->origin.z) / off_bb->delta - 0.5);
                    bzmax = floor( (ppw[0].z - off_bb->origin.z) / off_bb->delta + 0.5);
            }

            bxmin = NIIK_IMINMAX(bxmin-2,0,off_bb->xdim-1);
            bxmax = NIIK_IMINMAX(bxmax+2,0,off_bb->xdim-1);

            bymin = NIIK_IMINMAX(bymin-2,0,off_bb->ydim-1);
            bymax = NIIK_IMINMAX(bymax+2,0,off_bb->ydim-1);

            bzmin = NIIK_IMINMAX(bzmin-2,0,off_bb->zdim-1);
            bzmax = NIIK_IMINMAX(bzmax+2,0,off_bb->zdim-1);

            for(k=bzmin; k<=bzmax; k++) {
                for(j=bymin; j<=bymax; j++) {
                    for(i=bxmin; i<=bxmax; i++) {
                        int n = i + j * off_bb->xdim + k * off_bb->area;
                        int m;
                        /* within each bounding box, test if the surface intersect with the viewing plane */
                        for(m=0; m<off_bb->ndata[n]; m++) {

                            kface * f = off_bb->data[n][m];
                            int coplanar;
                            niikpt isect[2];
                            double isectfrac;

                            if(info->slice_dir==0) { /*X*/
                                    if(f->pmax.x < ppw[0].x || /* skip if no overlap */
                                       f->pmin.x > ppw[0].x) continue;
                            } else if(info->slice_dir==1) { /*Y*/
                                    if(f->pmax.y < ppw[0].y || /* skip if no overlap */
                                       f->pmin.y > ppw[0].y) continue;
                            } else { /*Z*/
                                    if(f->pmax.z < ppw[0].z || /* skip if no overlap */
                                       f->pmin.z > ppw[0].z) continue;
                            }

                            if(!off_check_tri_tri_intersect_with_isectline(ppw[0],ppw[1],ppw[2],
                                    f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
                                    &coplanar, &qq[0], &qq[1]) ) {
                                continue;
                            }

                            niik_world_to_ijk(img, &qq[0], &isect[0]);
                            niik_world_to_ijk(img, &qq[1], &isect[1]);

                            /*TODO: plot a line between qq[0] and qq[1] intelligently*/
                            for(isectfrac=0.0; isectfrac<=1.0; isectfrac+=1.0/5.0) {
                                int mm;
                                niikpt p=niikpt_wavg(isect[0],isect[1],isectfrac);
                                /*move to the center of voxel and take integer*/
                                p.x=floor((p.x+0.5)*info->scale);
                                p.y=floor((p.y+0.5)*info->scale);
                                p.z=floor((p.z+0.5)*info->scale);
                                
                                if(!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z) ||
                                   p.x<0.0 ||
                                   p.y<0.0 || 
                                   p.z<0.0 ||  
                                   p.x >= img->nx*info->scale || 
                                   p.y >= img->ny*info->scale || 
                                   p.z >= img->nz*info->scale ) continue;

                                if(info->slice_dir==0) { /*X*/
                                    mm = p.y*3 + (nj*round_scale-p.z-1)*ni*3*round_scale;
                                } else if(info->slice_dir==1) { /*Y*/
                                    mm = p.x*3 + (nj*round_scale-p.z-1)*ni*3*round_scale;
                                } else { /*Z*/
                                    mm = p.x*3 + (nj*round_scale-p.y-1)*ni*3*round_scale;
                                }
        
                                if(f->color!=NULL) {
                                    bout[mm++] = (unsigned char)NIIK_DMINMAX(255 * f->color[0],0,255);
                                    bout[mm++] = (unsigned char)NIIK_DMINMAX(255 * f->color[1],0,255);
                                    bout[mm++] = (unsigned char)NIIK_DMINMAX(255 * f->color[2],0,255);
                                    /* fprintf(stdout,"color %3i %3i %3i\n",pixels[mm],pixels[mm+1],pixels[mm+2]); */
                                } else { /* no color -> green */
                                    bout[mm++] = 0;
                                    bout[mm++] = 255;
                                    bout[mm++] = 0;
                                }
                            }
                        } /* each triangle in the bounding box */
                    }
                }
            } /* ijk of the bounding boxes */
        } /* object for slice-dir */

        #if 0
        TIFFSetField(tifimg, TIFFTAG_BITSPERSAMPLE,  8);
        TIFFWriteEncodedStrip(tifimg, 0, bout, nj * ni * round_scale * round_scale * 3);
        TIFFClose(tifimg);
        #endif
        // TODO: flip Y ?
        stbi_write_png(fname, ni*round_scale, nj*round_scale,3, bout, ni*round_scale*3);

        fprintf(stdout,"    writing %s\n",fname);
    } /*next slice*/

    if(bout!=NULL) {
        free(bout);
        bout=NULL;
    }
    return 1;
} 


