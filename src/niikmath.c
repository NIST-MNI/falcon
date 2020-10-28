/* Filename:     niikmath.c
 * Description:  general nifti1 program
 * Author:       Kunio Nakamura
 * Date:         February 23, 2012
 */
#include "falcon.h"
#include "falcon_cortex.h"
#include "falcon_morph.h"
#include "falcon_seg.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#define NIIKMATH_MAJOR_VERSION (0)
#define NIIKMATH_MINOR_VERSION (14)
#define NIIKMATH_MICRO_VERSION (0)

#define NIIKMATH_CHECK_VERSION(major,minor,micro)           \
  (NIIKMATH_MAJOR_VERSION > (major) ||                  \
   (NIIKMATH_MAJOR_VERSION == (major) && NIIKMATH_MINOR_VERSION > (minor)) || \
   (NIIKMATH_MAJOR_VERSION == (major) && NIIKMATH_MINOR_VERSION == (minor) && \
    NIIKMATH_MICRO_VERSION >= (micro)))

static char * niikmath_version[] = {
  "  niikmath version history\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0     Kunio 2012-02-23\n"
  "  -initial testing version\n"
  "  -helpful to use export COMP_WORDBREAKS=\"$COMP_WORDBREAKS,\"\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.1   Kunio 2012-02-25\n"
  "  -matrix functions\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.2   Kunio 2012-02-29\n"
  "  -added aregmni for affine 12-dof MNI registration\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.3   Kunio 2012-02-29\n"
  "  -fixed symmetry calculation bug for aregmni\n"
  "  -added tiffimg for creating tiff slice image\n"
  "  -added swaput swapuv swapvt for creating switching\n"
  "   4th, 5th, and 6th dimensions but swapuv assumes\n"
  "   there's not u-dimension\n"
  "  -as defined by nifti1_io.h dimensions are x,y,z,t,u,v,w\n"
  "  -as defined by nifti1_kunio.h dimensions dimension v is for\n"
  "   colors\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.3   Kunio 2012-03-02\n"
  "  -fixed aregmni(-sform)'s memory free processes\n"
  "  -fixed imglist list size display\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.4   Kunio 2012-03-04\n"
  "  -added sobel filters (sobel, sobelx, sobely, sobelz, sobelvoxel)\n"
  "  -added laplacemap\n"
  "  -added otsu and ridler threshold methods\n"
  "  -added bounds bounds-inc bounds-color for boundary images\n"
  "  -added calcvol for volume calculation\n"
  "  -added xdirimg ydirimg zdirimg adirimg for directional image\n"
  "  -added diffinfo infodiff absdiffinfo infoabsdiff\n"
  "  -added median for median filter\n"
  "  -added closeholes for closing holes\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.5   Kunio 2012-03-05\n"
  "  -added decompose-affine for affine matrix decomposition\n"
  "  -added NBCSR for non-brain constrained symmetric registration; Ref = Chen 2008 ISMRM\n"
  "  -added rotate=* for rotating image\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.6   Kunio 2012-03-07\n"
  "  -added resample\n"
  "  -fixed NBCSR's '-outcheck'\n"
  "  -'interp' supports images with more than 3D, but interpolation is done in 3D. that is, the interpolating point must be 3d (cannot interpolate a 4D point)\n"
  "  -updated operation list\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.7   Kunio 2012-03-08\n"
  "  -added fsasc2off (for now)\n"
  "  -added shrinkwrap (for now)\n"
  "  -added obj2img (for now)\n"
  "  -added make-offe (for now)\n"
  "  -added applyaffine-obj (for now)\n"
  "  -added add-off-color (for now)\n"
  "  -added offinfo (for now)\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.8   Kunio 2012-03-09\n"
  "  -added cropmask\n"
  "  -fixed invmatrix\n"
  "  -added maskimg\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.0.9   Kunio 2012-03-11\n"
  "  -added applywarp\n"
  /*  "  -added combinewarp\n"*/
  "  ----------------------------------------------------------------\n"
  "  version 0.1.0  Kunio 2012-03-12\n"
  "  -added off2pts and obj2pts\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.1  Kunio 2012-03-15\n"
  "  -modified the usage for '-min' and '-max' for iscale\n"
  "  --there are '-imin' '-imax' '-omin' and '-omax' now for iscale\n"
  "  --for tiffimg command use '-imin' and '-imax'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.2  Kunio 2012-03-16\n"
  "  -added 'histo' for calculating histogram\n"
  "  --to be plotted using octave, excel, etc\n"
  "  -fixed bug in affine transformation using multi-thread processing\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.3  Kunio 2012-03-22\n"
  "  -added more stats for 'offinfo'\n"
  "  -added 'kmul' for constant multiplication\n"
  "  -added 'addimgs' for adding multiple images\n"
  "  -modifying 'bseg' with spatially varying erosion filter size (ongoing)\n"
  "  -modified 'aregmni-sform' so that aregmni-sform with '-out=<out.matrix>'\n"
  "   will write both into sform and out.matrix\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.4  Kunio 2012-03-25\n"
  "  -modifyied 'bseg' significantly\n"
  "  --works quite well\n"
  "  --changing the parameters like radius, thresh, ithresh, uthresh,\n"
  "    and uthresh2 should make a fairly good segmentation.\n"
  "  -added '-remove-obj-color' or '-remove-off-color' for\n"
  "   removing color from off file\n"
  "  -added a category of operations for 'niikcortex' including\n"
  "   'niikcortex-initics-expand'\n"
  "  -to add others such as 'niikcortex-initics-shrink'\n"
  "  -added 'closebrain'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.5  Kunio 2012-03-29\n"
  "  -added 'halfway-matrix'\n"
  "  -added 'calc_mtr'\n"
  "  -fixed 'aregimg' -> default histogram size is 32\n"
  "  -added optional usage for 'aregimg'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.6  Kunio 2012-04-02\n"
  "  -adding subroutines for cortex detection -> to be minor update as\n"
  "   version 0.2.0\n"
  "  -added 'distancemap'\n"
  "  -added 'vertexinfo' for showing off object's vertex\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.1.7  Kunio 2012-04-25\n"
  "  -added 'gaussnoise' for creating a gaussian noise image\n"
  "  -added 'maskthresh' for thresholding and masking at the same time\n"
  "  -added 'dbc' for differential bias correction\n"
  "     Ref = Lewis and Fox 2004 NeuroImage\n"
  "  -added 'kadd' 'kdiv' 'ksub' like 'kmul' for constant addition,\n"
  "   divsion, and subtraction\n"
  "  -changed rotate=* to permute=*\n"
  "  -added 'permute=norm' and 'permute=LAS' for (Left,Anterior,Superior\n"
  "   direction\n"
  "  -added 'avgimgs' for averaging multiple images (same dimension)\n"
  "  -added 'avgimgs-fov' for averaging multiple images with field-of-\n"
  "   view correction; require a list of matrices\n"
  "  -added 'get3dimg' for getting a 3d image from more than 3d image\n"
  "  -added 'fill-lesion' for lesion filling\n"
  "  -changed 'niikcortex-gwi' with check flag and debug flag\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.0  Kunio 2012-05-01\n"
  "  -added 'niikcortex' routines for clada like processing\n"
  "  -added voxel-wise dilation and erosion\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.1  Kunio 2012-05-24\n"
  "  -added 'montage' for combining images\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.2  Kunio 2012-06-01\n"
  "  -changed threshold calculation in 'vseg' using\n"
  "   niik_image_bseg_basic_thresh_with_mask in nifti1_kunio_bseg.c\n"
  "  -fixed 'nregimg' with affine matrix\n"
  "  -added to understand IDENTITY for '-matrix' and '-invmatrix' in\n"
  "   addition to '-matrixlist'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.3  Kunio 2012-06-07\n"
  "  -changed 'avgimg-fov' so that the bspline is turned off for\n"
  "   voxels outside the field-of-view.\n"
  "  -added an option for 'determinant' to calculate matrix\n"
  "   determinants\n"
  "  -changed non-brain registration (R4) in 'nbcsr' so that global\n"
  "   scaling is unchanged.\n"
  "  -added 'bseg0' for the simplest brain segmentation function\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.4  Kunio 2012-06-15\n"
  "  -added 'featuremap1'\n"
  "  -fixed a bug in nearest neighbor interpolation\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.5  Kunio 2012-06-23\n"
  "  -added 'setone' for set image to one (like clear)\n"
  "  -added 'padimg' for padding a 3d image\n"
  "  -modified 'bseg' for brain segmentation\n"
  "  -modified 'invmatrix' section with many similar operation names\n"
  "  -modified 'halfway-matrix' using Ref[Alexa 2002 ACM]\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.6  Kunio 2012-06-28\n"
  "  -modified 'niikcortex-initocs' for cortical surface model\n"
  "  -added 'bseg1' for brain segmentation in low-resolution (3-5mm)\n"
  "   T1-weighed MRIs\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.7  Kunio 2012-07-30\n"
  "  -added 'add-off-color-curvature' or 'add-obj-color-curvature' to\n"
  "   calculate the curvature and map on a surface object\n"
  "  -added option '-delta=<D>' for 'histo' for histogram bin size\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.8  Kunio 2012-08-21\n"
  "  -added 'bimodalfit' to fit to bimodal histogram\n"
  "  -added automatic brain ROI masking for 'vseg' or ventricle\n"
  "   segmentation\n"
  "  -added 'cuberille' to create a surface model\n"
  "  -added 'edgeinfo' to show information about an edge in object\n"
  "  -fixed 'bseg' to run without '-matrix=<M.matrix>'\n"
  "  -added options for 'niikcortex-gwi'\n"
  "  -added 'gaussPDF' for gaussian PDF\n"
  "  -added '-FWHM=<fwhm>' for sobel filters ('sobel?')\n"
  "  -added 'thinedge' for non-maximum suppression to find the binary edge mask\n"
  "  -added '-outcheck=<check.nii>' option for 'maskimg' and 'maskout\n"
  "  -modified the bimodal fitting (niik_bimodal_fit) function so that it\n"
  "   uses multi-seed / multi-level optimization\n"
  "  -corrected usage for '-seedfill-grad='\n"
  "  -added 'cuberille-prep' and 'cuberille-prep-topo'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.2.9  Kunio 2012-09-29\n"
  "  -fixed 'aregimg2' for '-nmi' option\n"
  "  -fixed '-img2=' reading\n"
  "  -added 'idiv' for image division\n"
  "  -added 'voxel-gt' 'voxel-ge' 'voxel-lt' 'voxel-le' 'voxel-eq'\n"
  "   'voxel-ne' for voxel-wise comparison\n"
  "  -added a check for 'niik_image_get_mode' for 'num'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.3.0  Kunio 2012-10-10\n"
  "  -added 'siena' for niik implementaiton of SIENA; but not fully validated\n"
  "  -added 'NBCR' for non-symmetric version of NBCSR\n"
  "  -modified 'nregimg'\n"
  "  -added 'ocvobj' for outer contour volume\n"
  "  -added 'intratemplate' for subject-specific nonlinear registration\n"
  "   and template creation\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.3.1  Kunio 2012-11-15\n"
  "  -added 'interpacket-correction' for correction inter-packet or\n"
  "   interleave patient motion (validated thru simulation)\n"
  "  -added 'interpacket-split' and interpacket-merge for splitting \n"
  "   and merging interpacket (interleave) images\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.3.2  Kunio 2012-11-19\n"
  "  -modified 'nregimg' for nonlinear registration across time points\n"
  "  -modified function 'niik_image_get_upper_threshold' to avoid greater than max\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.4.0  Kunio 2012-12-18\n"
  "  -added functions for minc-1 reading\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.5.1  Kunio 2013-03-25\n"
  "  -type conversion for dilate erode open close\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.5.2  Kunio 2013-04-16\n"
  "  -lesion filling erodes WM mask to compute mean and stdv\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.6.0  Kunio 2013-05-26\n"
  "  -fixed sobel filter (z-dir)\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.7.0  Kunio 2013-05-26\n"
  "  -updating intensity normalization (inorm3)\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.8.0  Kunio 2013-07-15\n"
  "  -lost track of previous changes...\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.9.0  Kunio 2013-07-20\n"
  "  -changed tiff so that it can write multiple slices\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.10.0  Kunio 2013-08-02\n"
  "  -added gray-dilate for grayscale dilation\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.0  2013-11-21, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed Rician noise function 'add-rice' and 'add-rice-random'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.1  2014-01-03, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -trying to port this to my new laptop\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.2  2014-01-27, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added an option for kernels\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.3  2014-01-29, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed some warning at compiling\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.4  2014-02-04, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -stats option has threshold for mask\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.5  2014-02-12, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added comment for heaviside function in --list option\n"
  "  -added skewness for 'stats'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.6  2014-02-16, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added min/max in physical coordinates for 'stats'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.7  2014-02-16 (b), Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added 'map-spectral'\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.8  2014-03-11, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added 'labelvol' for measuring volumes of label file\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.9  2014-03-20, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added 'modify-header-scl_slope' for changing the intensity scaling slope\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.10 2014-03-27, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -working on cortical pipeline\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.11 2014-03-28, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -removed niikcortex-deform and writing new program for niikcortex_deform\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.12 2014-05-04, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added the operation for avgobjs to average multiple objects\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.13 2014-05-09, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -removing niikcortex-initics[-*] and writing niikcortex_initics.c\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.14 2014-06-02, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -minor edit on absdiffinfo\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.15 2014-06-02, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -more digit for calcvol\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.16 2014-06-16, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -adding spherical coordinate for objects\n"
  "  -added add-off-color-sphere-theta and add-off-color-sphere-psi\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.17 2014-06-16, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added spherical coordinate for objects\n"
  "  -added sphere-off to create a sphere\n"
  "  -added optional usage -phase=<phase>\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.18 2014-06-21, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed output data type for removesqform and copysqform\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.12.19 2014-07-01, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added modify-header-dz to change dz\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.0 2014-08-14, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed bug in niik_image_resample_3d_update\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.1 2014-08-22, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -removed general check for FSLDIR\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.2 2014-08-22, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed ndim issues in mergetdim\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.3 2014-09-30, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added absdiffimg to calculate absolute difference images\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.4 2014-10-13, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -fixed niik_image_add_multiple error\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.13.5 2014-11-17, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -checks applyaffine\n"
  "  -added modify-header-dx and modify-header-dy\n"
  "  ----------------------------------------------------------------\n"
  "  version 0.14.0 2014-12-28, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>\n"
  "  -added principal component analysis\n"
  "  -added detailed error message\n"
  "  ----------------------------------------------------------------\n"
};

static char * niikmath_op_list[] = {
  "  ----------------------------------------------------------------\n"
  "    list of operations\n"
  "  ----------------------------------------------------------------\n"
  "    very simple operations\n"
  "  typeconvert  -converts the image data type\n"
  "  info         -displays image stats\n"
  "  infodiff     -displays difference image stats\n"
  "  infoabsdiff  -displays absolute difference image stats\n"
  "  clear        -clears an image\n"
  "  setone       -sets one for an image\n"
  "  setval       -sets voxel value\n"
  "  calcvol      -calculates volumes\n"
  "  labelvol     -calculates volumes of label file\n"
  "  ----------------------------------------------------------------\n"
  "    voxel-wise comparisons\n"
  "  voxel-gt     -greater than comparison\n"
  "  voxel-ge     -greater than or equal to comparison\n"
  "  voxel-lt     -less than comparison\n"
  "  voxel-le     -less than or equal to comparison\n"
  "  voxel-eq     -equal to comparison\n"
  "  voxel-ne     -not equal to comparison\n"
  "  ----------------------------------------------------------------\n"
  "    more simple operations\n"
  "  maskimg       -mask image with a mask\n"
  "  maskout       -mask-out image with a mask\n"
  "  resample      -resamples image with new voxel spacing or image size\n"
  "  avgimgs       -average multiple images (voxel-by-voxel)\n"
  "  mulimgs       -multiply multiple images (voxel-by-voxel)\n"
  "  maximgs       -maximize multiple images (voxel-by-voxel)\n"
  "  addimgs       -add multiple images (voxel-by-voxel)\n"
  "  subimg        -subtract two images (voxel-by-voxel)\n"
  "  absdiffimg    -absolute difference of two images (voxel-by-voxel)\n"
  "  divimg        -subtract two images (voxel-by-voxel)\n"
  "  ksub          -subtract image voxel value from a constant\n"
  "  heavisideF    -heaviside function\n"
  "  map-spectral  -put spectral colors to an image\n"
  "  map-jacobian  -put jacobian colors to an image\n"
  "  log           -natural log\n"
  "  avgimgs-fov   -average multiple images with field-of-view correction\n"
  "  kmul          -multiply a constant (-val=<K>)\n"
  "  iscale        -intensity scaling\n"
  "  gaussPDF      -gaussian normal PDF\n"
  "  histo         -creates histogram\n"
  "  binimg        -binarize image\n"
  "  gaussnoise    -image of Gaussian random\n"
  "  fill-lesion   -fill lesions as NAWM mean/stdv\n"
  "  tiffimg       -creates tiff of a slice\n"
  "  bimodalfit    -fits the intensity histogram to bimodal Gaussian distributions\n"
  "  brainhistofit -fits intensity histogram to brain tissues (GM,WM,CSF)\n"
  "  ----------------------------------------------------------------\n"
  "  pseudoT2     -calculates T2 map from dual echo images\n"
  "  ----------------------------------------------------------------\n"
  "    changing image size\n"
  "  montage      -image montage\n"
  "  cropimg      -crop image\n"
  "  cropmask     -crop image according to mask\n"
  "  padimg       -pads image with zeros\n"
  "  ----------------------------------------------------------------\n"
  "    threshold operations\n"
  "  thresh       -threshold image\n"
  "  otsu         -threshold image using Otsu's algorithm\n"
  "  ridler       -threshold image using Ridler's algorithm\n"
  "  maskthresh   -threshold image and mask image (output is an image)\n"
  "  ----------------------------------------------------------------\n"
  "    segmentation operations\n"
  "  bseg         -brain segmentation\n"
  "  vseg         -ventricle segmentation\n"
  "  segcolor     -write colored segmentation results\n"
  "  ----------------------------------------------------------------\n"
  "    intensity normalization (and related functions)\n"
  "  iscale       -intensity scaling\n"
  "  inorm        -intensity normalization\n"
  "  inorm2       -intensity normalization using histograms\n"
  "               -registration is not required\n"
  "  inorm3       -intensity normalization using linear regression\n"
  "               -registration is required\n"
  "  dbc          -differential bias correction\n"
  "                [Lewis and Fox 2004 NeuroImage]\n"
  "  ----------------------------------------------------------------\n"
  "    boundary images\n"
  "  bounds       -creates a boundary image\n"
  "  bounds-inc   -creates a boundary image (boundary within mask)\n"
  "  bounds-color -creates a boundary image with colors\n"
  "  ----------------------------------------------------------------\n"
  "    create maps\n"
  "  laplacemap   -creates laplace map\n"
  "  distancemap  -creates signed distance map\n"
  "  xdirimg      -x-directional map (voxel intensity for voxel position)\n"
  "  ydirimg      -y-directional map (voxel intensity for voxel position)\n"
  "  zdirimg      -z-directional map (voxel intensity for voxel position)\n"
  "  adirimg      -directional maps (voxel intensity for voxel position, 3D)\n"
  "  ----------------------------------------------------------------\n"
  "    image permutation\n"
  "  permute-norm         -converts to more or less RAS directions\n"
  "  permute-left         -turns left  [x'=-y; y'=x;  z'=z]\n"
  "  permute-right        -turns right [x'=y;  y'=-x; z'=z]\n"
  "  permute-cor2sag      -rotates     [x'=z;  y'=y;  z'=x]\n"
  "  permute-cor2axi      -rotates     [x'=x;  y'=z;  z'=y]\n"
  "  permute-sag2cor      -rotates     [x'=z;  y'=z;  z'=x]\n"
  "  permute-sag2axi      -rotates     [x'=y;  y'=z;  z'=x]\n"
  "  permute-cor2axi      -rotates     [x'=x;  y'=z;  z'=y]\n"
  "  permute-axi2sag      -rotates     [x'=y;  y'=z;  z'=x]\n"
  "  permute-axi2cor      -rotates     [x'=x;  y'=z;  z'=y]\n"
  "  permute-flipxyz      -flips in all x,y,z directions [x'=-x; y'=-y; z'=-z]\n"
  "  ----------------------------------------------------------------\n"
  "    morphologic operations\n"
  "               -uses -radius=<R> to set the kernel radius\n"
  "  open         -morphologic opening\n"
  "  close        -morphologic closing\n"
  "  dilate       -morphologic dilation\n"
  "  erode        -morphologic erosion\n"
  "  closeholes   -morphologic closing by seed-filling background from\n"
  "                image edge\n"
  "  closebrain   -series of morphologic operations to find brain ROI:\n"
  "                (1) dilation, (2) close holes, and (3) erosion\n"
  "  ----------------------------------------------------------------\n"
  "    simple filters\n"
  "  gauss        -gaussian filter\n"
  "  median       -median filter\n"
  "  sobelx       -x-directional sobel filter\n"
  "  sobely       -y-directional sobel filter\n"
  "  sobelz       -z-directional sobel filter\n"
  "  sobelm       -sobel filter magnitude sqrt(x^2+y^2+z^2)\n"
  "  sobel        -sobel filter for all directions\n"
  "  sobelvoxel   -sobel filter for a voxel\n"
  "  ----------------------------------------------------------------\n"
  "    more than 3 dimensions\n"
  "    -v (6th)-dimension is used as color in this set of programs\n"
  "  mergetdim         -merges images in t(4th)-dim\n"
  "  mergeudim         -merges images in u(5th)-dim\n"
  "    swap routines assume the largest dimension is unused\n"
  "  swaput            -swaps u and t dimensions\n"
  "  swapuv            -swaps u and v dimensions\n"
  "  swapvt            -swaps v and t dimensions\n"
  "  ----------------------------------------------------------------\n"
  "    changing header info\n"
  "  qformradio        -puts arbitrary qform with center in middle\n"
  "  removesqform      -removes sform and qform\n"
  "  copysqform        -copies sform and qform to output image\n"
  "  addsform          -adds sform to an image\n"
  "  modify-header-scl_slope\n"
  "                    -changes scl_slope\n"
  "  modify-header-dz  -changes dz\n"
  "  ----------------------------------------------------------------\n"
  "    matrix operations\n"
  "  invmatrix         -inverts matrix\n"
  "  multiplymatrix    -multiplies (multiple) matrix(matrices)\n"
  "  affine-decompose decompose-affine\n"
  "                    -decomposes matrix into affine parameters\n"
  "  affinematrix      -creates affine matrix\n"
  "  halfway-matrix    -calculates halfway matrix\n"
  "  avgmatrix         -calculates average of matrices (approximate)\n"
  "  determinant       -calculates determinant of a matrix\n"
  "  ----------------------------------------------------------------\n"
  "    affine registration / transformation\n"
  "  aregmni           -affine registration to MNI152 in FSL\n"
  "  aregmni-sform     -affine registration to MNI152 in FSL\n"
  "                    -output matrix is put into image header sform\n"
  "  applysform        -applies affine transformation using sform\n"
  "  writesform        -writes sform affine transformation\n"
  "  NBCSR             -non-brain constrained symmetric registration\n"
  "  NBCR              -non-brain constrained registration (non-symmetric)\n"
  "  NBCSR-rigid       -non-brain constrained symmetric registration\n"
  "                     rigid (6-dof) version\n"
  "  applyaffine       -applies affine transformation using input matrix\n"
  "  aregimg           -affine registration\n"
  "  aregimg2          -multi-seed/multi-resolution affine registration\n"
  "  ----------------------------------------------------------------\n"
  "    nonlinear registration / transformation\n"
  "  applywarp         -nonlinearly transforms images\n"
  "  nregimg           -nonlinearly registers images (under construction)\n"
  "  intratemplate     -nonlinearly registers images from the same subject and create a template\n"
  "  ----------------------------------------------------------------\n"
  "    OFF object processing\n"
  "  cuberille         -creates a surface model from a binary image\n"
  "  applyaffine-obj   -applies affine transformation using input\n"
  "                     matrix on an off object\n"
  "  fsasc2off         -converts FreeSurfer's ASCII format to off object\n"
  "  shrinkwrap        -shrink-wraps a mask image (not well validated)\n"
  "  sphere-off        -creates a sphere\n"
  "  obj2img           -creates boundary at the object surface\n"
  "  avgobjs           -averages objects, same object config required\n"
  "  obj2pts off2pts   -creates point file with a list of 3d coordinates\n"
  "                     for each vertex\n"
  "  make-offe         -write complete off file with offe file\n"
  "  combine-obj       -combine multiple off files\n"
  "  add-off-color\n"
  "  add-obj-color\n"
  "                    -adds one color to an object\n"
  "  add-off-color-curvature\n"
  "  add-obj-color-curvature\n"
  "                    -adds colors to an object based on curvatures\n"
  "  add-off-color-sphere-theta\n"
  "                    -creates a color-map of spherical coordinates (theta)\n"
  "  add-off-color-sphere-psi\n"
  "                    -creates a color-map of spherical coordinates (psi)\n"
  "  offinfo           -displays off information\n"
  "  vertexinfo        -displays vertex information\n"
  "  edgeinfo          -displays edge information\n"
  "  faceinfo          -displays face information\n"
  "  ----------------------------------------------------------------\n"
  "    Cortical surface reconstruction\n"
  "  niikcortex-gwi     -gray matter - white matter interface mask\n"
  "  niikcortex-wm      -white matter segmentation for niikcortex\n"
  "  niikcortex-mimg    -modifying image (obsolete)\n"
  "  niikcortex-initocs -creates initial pial surface\n"
  "  niikcortex-nonctx  -creates (non-)cortex label\n"
  "  niikcortex-deform  -deform cortical surface model\n"
  "  niikcortex-thick-color  - apply colour according to cortical thickness\n"
  "  niikcortex-edge    - edge-based cortical atrophy method\n"
  "  ----------------------------------------------------------------\n"
  "  inter-packet (interleave) processing\n"
  "  interpacket-correction     -corrects interpacket patient motion by registering\n"
  "  interpacket-split          -splits interpacket (interleave) images by specified number\n"
  "  interpacket-merge          -merges interpacket (interleave) images\n"
  "  ----------------------------------------------------------------\n"
  "  To see the specific usage and options, type operation: e.g.,\n"
  "  ----------------------------------------------------------------\n"
  "  Other (less common)operations:\n"
  "  grid-lines            -creates image with grid-lines\n"
  "                         e.g. to show nonlinear deformation fields\n"
  "  > niikmath offinfo\n"
  "  ----------------------------------------------------------------\n"
};

static char * niikmath_op_niikcortex_list[] = {
  "  ----------------------------------------------------------------\n"
  "    Cortical surface reconstruction\n"
  "  ----------------------------------------------------------------\n"
  "  niikcortex-gwi     -gray matter - white matter interface mask\n"
  "  niikcortex-wm      -white matter segmentation for niikcortex\n"
  "  niikcortex-mimg    -modifying image (obsolete)\n"
  "  niikcortex-initocs -creates initial pial surface\n"
  "  niikcortex-nonctx  -creates (non-)cortex label\n"
  "  niikcortex-deform  -deform cortical surface model\n"
  "  niikcortex-thick-color  - apply colour according to cortical thickness\n"
  "  niikcortex-edge    - edge-based cortical atrophy method\n"
  "  ----------------------------------------------------------------\n"
};

static char * niikmath_option_list[] = {
  "  ----------------------------------------------------------------\n"
  "    general optional usage\n"
  "  for more detailed usage, write operation (e.g. 'niikmath bseg')\n"
  "  ----------------------------------------------------------------\n"
  "  -u -help          -help and general usage\n"
  "  --version         -version info\n"
  "  -list             -list of operations\n"
  "  ----------------------------------------------------------------\n"
  "    common options\n"
  "  -in=<img.nii>     -set input image\n"
  "  -out=<out>        -set output filename\n"
  "  -obj=<obj.off>    -set input object filename\n"
  "  -imglist=<in1.nii>,<in2.nii>,...\n"
  "                    -list of input images\n"
  "  ----------------------------------------------------------------\n"
  "    matrix options\n"
  "  -The matrix format is same as FSL\n"
  "  -The matrix is in ascii format and should not have empty rows\n"
  "  -matrix=<mat>     -set input matrix\n"
  "  -matrixlist=<mat1>,<mat2>,...\n"
  "                    -list of of input matrices\n"
  "                    -separate by comma\n"
  "  -invmatrix=<mat>  -set and use inverse of input matrix\n"
  "  ----------------------------------------------------------------\n"
  "    XFM options\n"
  "  -The transformation format is linear MINC XFM with single transform\n"
  "  -xfm=<trans.xfm>     -set input matrix\n"
  "  -invxfm=<trans.xfm>  -set and use inverse of input matrix\n"
  "  ----------------------------------------------------------------\n"
  "    interpolation options\n"
  "  -nn                -nearest neighbor\n"
  "  -linear            -trilinear\n"
  "  -bspline           -b-spline\n"
  "  ----------------------------------------------------------------\n"
  "    data types\n"
  "  -uint8             -unsigned 8-bit integer\n"
  "  -uint16            -unsigned 16-bit integer\n"
  "  -uint32            -unsigned 32-bit integer\n"
  "  -uint64            -unsigned 64-bit integer\n"
  "  -int8              -8-bit signed integer\n"
  "  -int16             -16-bit signed integer\n"
  "  -int32             -32-bit signed integer\n"
  "  -int64             -64-bit signed integer\n"
  "  -float32           -32-bit floating point\n"
  "  -float64           -64-bit floating point\n"
  /*  "  ----------------------------------------------------------------\n"
  "    Parameters\n"
  "  -thresh=<T>        -set threshold to T\n"
  "  -uthresh=<T2>      -set upper threshold to T2\n"
  "  -ithresh=<IT>      -set initial threshold to IT\n"
  "  -radius=<R>        -set radius to T\n"
  "  -FWHM=<fwhm>       -set FWHM to <fwhm>\n"
  "  -mask=<mask.nii>   -set mask to <mask.nii>\n"
  "  -iter=<iter>       -set iteration\n"
  "  -min=<min>         -set min value\n"
  "  -max=<max>         -set max value\n"
  "  -ijk=<i>,<j>,<k>   -set ijk\n"
  "  -xyz=<i>,<j>,<k>   -set xyz\n"
  "  -rgb=<R>,<G>,<B>   -set RGB\n"*/
  "  ----------------------------------------------------------------\n"
};

static char * niikmath_matrix_descrip[] = {
  "  matrix format\n"
  "  -ascii text file\n"
  "  -numbers are separated by space(s)\n"
  "  -rows are read once by fgets\n"
  "  -number of columns must be same for each row"
  "  -affine transformation is 4-by-4\n"
  "  -but support for N-by-M matrix not just 4x4\n"
  "   e.g. invmatrix and multiplymatrix\n"
  "  -example of identity matrix\n"
  "     1 0 0 0\n"
  "     0 1 0 0\n"
  "     0 0 1 0\n"
  "     0 0 0 1\n"
  "  -example of translation matrix\n"
  "     1 0 0 5.0\n"
  "     0 1 0 3.2\n"
  "     0 0 1 1.12\n"
  "     0 0 0 1\n"
};


void usage() {
  fprintf(stdout,"niikmath general tool\n");
  fprintf(stdout,"  usage: niikmath operation [options]\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *mni_ven_list[4],
  *tmpimglist[512],
  *dist_map=NULL, /* distance map */
   *mni_img=NULL,  /* for mni registration (aregmni), mni image */
    *mni_seg=NULL,  /* for mni registration (aregmni), mni brain mask */
     *brainmask=NULL,
      **imglist=NULL,
        **maskimglist=NULL,
          *img=NULL,
           *img2=NULL,
            *gradimg=NULL,
             *refimg=NULL,
              *outimg=NULL,
               *warpimg=NULL,
                **warpimglist=NULL,
                  *maskroi=NULL,
                   *maskref=NULL,
                    *maskimg=NULL;
  kobj
  **objlist=NULL,
    *obj=NULL,
     *obj2=NULL;
  kvert
  **vertlist=NULL,
    *vert,
    *v2,
    *v=NULL;
  kface
  *f,*face;
  kedge *edge;
  FILE *fp=NULL;

  niikpt
  ctr,
  pt;

  /* MRI *mgzimg; */
  niikvec
  *FWHM_list=NULL,
   *vec=NULL;

  niikmat
  *tstats=NULL,
   *mat=NULL,
    *histomat=NULL,
     *premat=NULL,
      *postmat=NULL,
       **matlist=NULL,
         **xfmlist=NULL,
           *cmap=NULL,
            *atv=NULL,
             *afmat2=NULL,
              *afmat=NULL;
  float
  *fimg=NULL;

  double
  *MRI_TE=NULL,*MRI_TR=NULL,
   xsum,ysum,wsum,xssq,yssq,xysum,d1,d2,
   bimodal_fit[6],
   conmat_tp,conmat_fp,conmat_fn,conmat_tn,
   affpar[20],
   daffpar[20],
   delta=NIIKMAX,
   inval=NIIKMAX,
   imin=NIIKMAX,
   imax=NIIKMAX,
   omin=NIIKMAX,
   omax=NIIKMAX,
   thresh=NIIKMAX,
   uthresh=NIIKMAX,
   uthresh2=NIIKMAX,
   ithresh=NIIKMAX,
   radius=NIIKMAX,
   radius2=NIIKMAX,
   FWHM=NIIKMAX,
   *vlist=NULL,
    dmean=NIIKMAX,
    dstdv=NIIKMAX,
    dval=NIIKMAX,
    dtmp=NIIKMAX,
    percent=NIIKMAX,
    phase=NIIKMAX,
    dmin,
    dvol,
    elen=NIIKMAX,
    dd,
    ijk[4],
    xyz[4],
    rgb[4],
    dmton,dmtoff,
    niikcortex_tissue_val[12],
    *dimg=NULL;

  unsigned char *bimg=NULL;
  long *limg=NULL;

  const char  *NIIKDIR=NULL;
  char  fcname[12]="niikmath",
        **CSlist=NULL,
        **outnamelist,
        fname[4096],
        *fsasc_name=NULL,
        *FSLDIR=NULL,
        *inname=NULL,
        *outname=NULL,
        *outcheckname=NULL,
        *cptr=NULL;

  nmi_obj *nmiobj=NULL;
  int
    xfmflag=0,
    num_MRI_TE=0,num_MRI_TR=0,
    morph_kernel_shape=NIIK_MORPH_3D_KERNEL_SQUARE,
    num_method = NIIK_NUM_METHOD_RUNGE_KUTTA, /*numerical method */
    debug=0,
    step=0,
    kernel=1,
    bbox_num=-8,
    histo_num=-128,
    histo_anum=0,
    tiff_xslice=-1,
    tiff_yslice=-1,
    tiff_zslice=-1,
    *tiff_xslices=NULL,
    *tiff_yslices=NULL,
    *tiff_zslices=NULL,
     flag_grad=0,
     flag_remesh=1,
     flag_bbox=1,
     iter=-1,
     ndegree=-1,
     dof=6,
     nmi_num[2],
     areg_cost=NIIK_REGISTER_UNKNOWN,
     nummatlist=0,
     numxfmlist=0,
     numimglist=0,
     nummaskimglist=0,
     numvlist=0,
     numobjlist=0,
     i=0,j=0,k=0,m=0,n=0,
     size=0,
     interp=NIIK_INTERP_LINEAR,
     datatype=-1,
     idx[9],
     img_dim_min[9],img_dim_max[9],
     warp_map_type=NIIK_WARP_MAP_UNKNOWN,
     *niikcortex_ctx_label=NULL,niikcortex_ctx_label_size=0,
      niikcortex_ctx_id=3,
      num=0,
      sc,nc;

  struct tm *stm;
  time_t ctm;
  char tmstr[256];
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  if(argc==1) {
    fprintf(stdout,"\n  niikmath : Kunio's nifti-related general program\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"  --help -u         : help and general usage\n");
    fprintf(stdout,"  --list            : list operations\n");
    fprintf(stdout,"  --version         : show version info   %i.%i.%i\n",
            NIIKMATH_MAJOR_VERSION,NIIKMATH_MINOR_VERSION,NIIKMATH_MICRO_VERSION);
    fprintf(stdout,"\n");
    exit(1);
  }

#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  srand(time(NULL));
  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  /*niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);*/
  niik_disp_exec_info(fcname,NIIKMATH_MAJOR_VERSION,NIIKMATH_MINOR_VERSION,NIIKMATH_MICRO_VERSION);

  /* check for the environmental variables: NIIKDIR */
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    exit(1);
  }

  /* affine parameters */
  for(n=0; n<20; n++) affpar[n]=0;
  for(n=7; n<=10; n++) affpar[n]=1;
  for(n=0; n<20; n++) daffpar[n]=0;

  /* niikcortex parameters */
  for(n=0; n<12; n++)
    niikcortex_tissue_val[n]=-1;

  nc = sc = 1;
  ijk[0] = 0;
  xyz[0] = 0;
  idx[0] = 0;
  rgb[0] = 0;

  nmi_num[0] = nmi_num[1] = 0;

  while(nc<argc) {
    /* fprintf(stdout,"  argv[%i] = %s\n",nc,argv[nc]); */

    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9) ||
          !strncmp(argv[nc],"-version",8)  ) {
        fprintf(stdout,"%s",*niikmath_version);
        exit(1);
      }

      else if(!strncmp(argv[nc],"--help-matrix",13)) {
        fprintf(stdout,"%s",*niikmath_matrix_descrip);
        exit(1);
      }

      else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*niikmath_option_list);
        exit(1);
      }

      else if(!strncmp(argv[nc],"--list",5)) {
        fprintf(stdout,"%s",*niikmath_op_list);
        exit(1);
      }

      else if(!strncmp(argv[nc],"-warp-disp-bspline=",19)) {
        cptr=argv[nc]+19;
        if(warpimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading warp      %s\n",cptr);
        if((warpimg=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
        if(warpimg->nu!=3) {
          fprintf(stderr,"ERROR: nu must be 3\n");
          exit(1);
        }
        warp_map_type=NIIK_WARP_MAP_DISP_BSPLINE;
      } /* -warp-disp-bspline=<warpimg.nii.gz> */

      else if(!strncmp(argv[nc],"-matrix=IDENTITY",16)) {
        if(afmat!=NULL) {
          fprintf(stderr,"ERROR: afmat is already defined\n");
          exit(1);
        }
        if((afmat = niikmat_identity(4,4))==NULL) {
          fprintf(stderr,"ERROR: niikmat_identity\n");
          exit(1);
        }
      } /* -matrix=IDENTITY */

      else if(!strncmp(argv[nc],"-ref=MNI152_1mm",15)) {
        if(refimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        NIIK_EXIT((niik_check_fsldir_exists()>0),fcname,"fsldir is not defined",1);
        sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
        fprintf(stdout,"[niikmath] reading refimg    %s\n",fname);
        if((refimg=niik_image_read(fname))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
      } /* reference image */

      else if(!strncmp(argv[nc],"-ref=MNI152_2mm",15)) {
        if(refimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        NIIK_EXIT((niik_check_fsldir_exists()>0),fcname,"fsldir is not defined",1);
        sprintf(fname,"%s/data/standard/MNI152_T1_2mm.nii.gz",FSLDIR);
        fprintf(stdout,"[niikmath] reading refimg    %s\n",fname);
        if((refimg=niik_image_read(fname))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
      } /* reference image */

      else if(!strncmp(argv[nc],"-warp-disp-map=",15)) {
        cptr=argv[nc]+15;
        if(warpimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading warp      %s\n",cptr);
        if((warpimg=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
        warp_map_type=NIIK_WARP_MAP_DISP;
      } /* -warp-disp-map=<warpimg.nii.gz> */

      else if(!strncmp(argv[nc],"-backward-euler",15)) {
        num_method = NIIK_NUM_METHOD_BACKWARD_EULER;
      }

      else if(!strncmp(argv[nc],"-seedfill-grad=",15)) {
        flag_grad=atoi(argv[nc]+15);
      } /* seedfill-grad for seed-filling subroutines */

      else if(!strncmp(argv[nc],"-warp-loc-map=",14)) {
        cptr=argv[nc]+14;
        if(warpimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading warp      %s\n",cptr);
        if((warpimg=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
        warp_map_type=NIIK_WARP_MAP_LOC;
      } /* -warp-loc-map=<warpimg.nii.gz> */


      else if(!strncmp(argv[nc],"-tiff-xslices=",14)) {
        cptr=argv[nc]+14;
        NIIK_EXIT(((CSlist=niik_csstring_to_list(cptr,&n))==NULL),fcname,"niik_csstring_to_list",9);
        tiff_xslices=(int *)calloc(n+1,sizeof(int));
        tiff_xslices[0]=n+1;
        for(i=0; i<n; i++) {
          tiff_xslices[i+1]=atoi(CSlist[i]);
          free(CSlist[i]);
        }
        free(CSlist);
      } /* tiff-xslices for tiffimg */
      else if(!strncmp(argv[nc],"-tiff-yslices=",14)) {
        cptr=argv[nc]+14;
        NIIK_EXIT(((CSlist=niik_csstring_to_list(cptr,&n))==NULL),fcname,"niik_csstring_to_list",9);
        tiff_yslices=(int *)calloc(n+1,sizeof(int));
        tiff_yslices[0]=n+1;
        for(i=0; i<n; i++) {
          tiff_yslices[i+1]=atoi(CSlist[i]);
          free(CSlist[i]);
        }
        free(CSlist);
      } /* tiff-yslice for tiffimg */
      else if(!strncmp(argv[nc],"-tiff-zslices=",14)) {
        cptr=argv[nc]+14;
        NIIK_EXIT(((CSlist=niik_csstring_to_list(cptr,&n))==NULL),fcname,"niik_csstring_to_list",9);
        tiff_zslices=(int *)calloc(n+1,sizeof(int));
        tiff_zslices[0]=n+1;
        for(i=0; i<n; i++) {
          tiff_zslices[i+1]=atoi(CSlist[i]);
          free(CSlist[i]);
        }
        free(CSlist);
      } /* tiff-zslice for tiffimg */

      else if(!strncmp(argv[nc],"-forward-euler",14)) {
        num_method = NIIK_NUM_METHOD_FORWARD_EULER;
      } else if(!strncmp(argv[nc],"-num-midpoint",13)) {
        num_method = NIIK_NUM_METHOD_FORWARD_EULER;
      }

      else if(!strncmp(argv[nc],"-tiff-xslice=",13)) {
        cptr=argv[nc]+13;
        tiff_xslice=atoi(cptr);
      } /* tiff-xslice for tiffimg */
      else if(!strncmp(argv[nc],"-tiff-yslice=",13)) {
        cptr=argv[nc]+13;
        tiff_yslice=atoi(cptr);
      } /* tiff-yslice for tiffimg */
      else if(!strncmp(argv[nc],"-tiff-zslice=",13)) {
        cptr=argv[nc]+13;
        tiff_zslice=atoi(cptr);
      } /* tiff-zslice for tiffimg */

      else if(!strncmp(argv[nc],"-distancemap=",13)) {
        cptr=argv[nc]+13;
        if(dist_map!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: distance map is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading           %s\n",cptr);
        if((dist_map=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -distance=<map.nii.gz> */

      else if(!strncmp(argv[nc],"-runge-kutta",12)) {
        num_method = NIIK_NUM_METHOD_RUNGE_KUTTA;
      }

      else if(!strncmp(argv[nc],"-niikcortex-",12)) {
        if     (!strncmp(argv[nc],"-niikcortex-edge-tissue=",24)) {
          if((CSlist=niik_csstring_to_list(argv[nc]+24,&n))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",argv[nc]);
            exit(1);
          }
          if(n!=6) {
            fprintf(stderr,"[niikmath] ERROR: -niikcortex-edge-tissue=<CSF_mean>,<CSF_stdv>,<GM_mean>,<GM_stdv>,<WM_mean>,<WM_stdv>\n");
            exit(1);
          }
          if(tstats!=NULL) {
            fprintf(stderr,"[niikmath] ERROR: tstats is already defined\n");
            exit(1);
          }
          tstats=niikmat_init(3,2);
          for(i=k=0; i<3; i++) {
            for(j=0; j<2; j++,k++) {
              tstats->m[i][j]=atof(CSlist[k]);
              free(CSlist[k]);
            }
          }
          free(CSlist);
          niikmat_display(tstats);
        } /* niikcortex-edge-tissue */
        else if(!strncmp(argv[nc],"-niikcortex-white-only",22)) {
          niikcortex_ctx_id = 1;
        } else if(!strncmp(argv[nc],"-niikcortex-pial-only",21)) {
          niikcortex_ctx_id = 2;
        } else if(!strncmp(argv[nc],"-niikcortex-ics-ran=",20)) {
          niikcortex_tissue_val[6] = atof(argv[nc]+20);
        } else if(!strncmp(argv[nc],"-niikcortex-ocs-ran=",20)) {
          niikcortex_tissue_val[7] = atof(argv[nc]+20);
        }

        else if(!strncmp(argv[nc],"-niikcortex-tissue=",19)) {
          if((CSlist=niik_csstring_to_list(argv[nc]+19,&n))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",argv[nc]);
            exit(1);
          }
          if(n!=8) {
            fprintf(stderr,"ERROR: need '-niikcortex-tissue=<GM>,<WM>,<CSF>,<Brain>,<WM_Surf>,<Pial_Surf>,<WM_Range>,<Pial_Range>'\n");
            exit(1);
          }
          for(n=0; n<8; n++) {
            niikcortex_tissue_val[n]=atof(CSlist[n]);
            free(CSlist[n]);
          }
          free(CSlist);
        } /* tissue intensities */

        else if(!strncmp(argv[nc],"-niikcortex-brain=",18)) {
          niikcortex_tissue_val[3] = atof(argv[nc]+18);
        } else if(!strncmp(argv[nc],"-niikcortex-label=",18)) {
          if((niikcortex_ctx_label = niik_read_int_vector(argv[nc]+18,&niikcortex_ctx_label_size))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niik_read_int_vector %s\n",argv[nc]+18);
            exit(1);
          }
        } else if(!strncmp(argv[nc],"-niikcortex-Rdgm=",17)) {
          radius2 = atof(argv[nc]+17);
        } else if(!strncmp(argv[nc],"-niikcortex-step=",17)) {
          delta = atof(argv[nc]+17);
        } else if(!strncmp(argv[nc],"-niikcortex-ics=",16)) {
          niikcortex_tissue_val[4] = atof(argv[nc]+16);
        } else if(!strncmp(argv[nc],"-niikcortex-ocs=",16)) {
          niikcortex_tissue_val[5] = atof(argv[nc]+16);
        } else if(!strncmp(argv[nc],"-niikcortex-csf=",16)) {
          niikcortex_tissue_val[2] = atof(argv[nc]+16);
        } else if(!strncmp(argv[nc],"-niikcortex-wm=",15)) {
          niikcortex_tissue_val[1] = atof(argv[nc]+15);
        } else if(!strncmp(argv[nc],"-niikcortex-gm=",15)) {
          niikcortex_tissue_val[0] = atof(argv[nc]+15);
        }
      } /* niikcortex */

      else if(!strncmp(argv[nc],"-register-sad",13)) {
        areg_cost=NIIK_REGISTER_SAD;
      } /* set affine registration's cost method */
      else if(!strncmp(argv[nc],"-register-ssq",13)) {
        areg_cost=NIIK_REGISTER_SSQ;
      } /* set affine registration's cost method */

      else if(!strncmp(argv[nc],"-maskimglist=",13)) {
        if(maskimglist!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: -maskimglist is already defined\n");
          exit(1);
        }
        if((CSlist=niik_csstring_to_list(argv[nc]+13,&nummaskimglist))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list\n");
          exit(1);
        }
        fprintf(stdout,"[niikmath] mask imglist size = %i\n",nummaskimglist);
        maskimglist = (nifti_image **)calloc(nummaskimglist, sizeof(nifti_image *));
        for(n=0; n<nummaskimglist; n++) {
          fprintf(stdout,"[niikmath] reading image:    %s\n",CSlist[n]);
          if((maskimglist[n]=niik_image_read(CSlist[n]))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",CSlist[n]);
            exit(1);
          }
        }
        free(CSlist);
      } /* -maskimglist=<mask1.nii>,<mask2.nii>... */

      else if(!strncmp(argv[nc],"-matrixlist=",12)) {
        if(matlist!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: -matrixlist is already defined\n");
          exit(1);
        }
        if((CSlist=niik_csstring_to_list(argv[nc]+12,&nummatlist))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",argv[nc]);
          exit(1);
        }
        fprintf(stdout,"[niikmath] matrixlist size = %i\n",nummatlist);
        matlist = (niikmat **)calloc(nummatlist,sizeof(niikmat *));
        for(n=0; n<nummatlist; n++) {
          fprintf(stdout,"[niikmath] reading matrix:   %s\n",CSlist[n]);
          if(!strncmp(CSlist[n],"IDENTITY",8)) {
            if((matlist[n] = niikmat_identity(4,4))==NULL) {
              fprintf(stderr,"ERROR: niikmat_identity(4,4)\n");
              exit(1);
            }
          } else if((matlist[n]=niikmat_read(CSlist[n]))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niikmat_read %s\n",cptr);
            exit(1);
          }
          free(CSlist[n]);
        }
        free(CSlist);
      } /* -matrixlist=<reg1>,<reg2>... */

      else if(!strncmp(argv[nc],"-brainmask=",11)) {
        cptr=argv[nc]+11;
        if(brainmask!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: maskimg is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading brain     %s\n",cptr);
        if((brainmask=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -brainmask=<brainmask.nii.gz> */

      else if(!strncmp(argv[nc],"-invmatrix=",11)) {
        if(afmat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: matrix is already defined\n");
          exit(1);
        }
        cptr=argv[nc]+11;
        if(!strncmp(cptr,"IDENTITY",8)) {
          if((afmat = niikmat_identity(4,4))==NULL) {
            fprintf(stderr,"ERROR: niikmat_identity\n");
            exit(1);
          }
        } else {
          fprintf(stdout,"[niikmath] reading matrix:   %s\n",cptr);
          if((afmat=niikmat_read(cptr))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niikmat_read\n");
            exit(1);
          }
          if(!niikmat_inverse_update(afmat)) {
            fprintf(stderr,"[niikmath] ERROR: niik_inverse_update(afmat)\n");
            exit(1);
          }
        }
      } /* -invmatrix=<matrix.txt> */

      else if(!strncmp(argv[nc],"-invxfm=",8)) {
        if(afmat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: matrix is already defined\n");
          exit(1);
        }
        cptr=argv[nc]+8;
        if(!strncmp(cptr,"IDENTITY",8)) {
          if((afmat = niikmat_identity(4,4))==NULL) {
            fprintf(stderr,"ERROR: niikmat_identity\n");
            exit(1);
          }
        } else {
          fprintf(stdout,"[niikmath] reading xfm:   %s\n",cptr);
          if((afmat=niikmat_read_xfm(cptr))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niikmat_read_xfm\n");
            exit(1);
          }
          if(!niikmat_inverse_update(afmat)) {
            fprintf(stderr,"[niikmath] ERROR: niik_inverse_update(afmat)\n");
            exit(1);
          }
        }
      } /* -invmatrix=<matrix.txt> */

      else if(!strncmp(argv[nc],"-histo-num=",11)) {
        histo_num=atoi(argv[nc]+11);
      } /* histo-num */

      else if(!strncmp(argv[nc],"-histo-avg=",11)) {
        histo_anum=atoi(argv[nc]+11);
      } /* histo-avg */

      else if(!strncmp(argv[nc],"-remesh=off",11)) {
        flag_remesh = 0;
      }

      else if(!strncmp(argv[nc],"-FWHM-list=",11)) {
        if((FWHM_list = niikvec_init_from_ascii(argv[nc]+11))==NULL) {
          fprintf(stderr,"[%s] ERROR: niikmat_init_from_ascii\n",fcname);
          exit(1);
        }
      } /* -FWHM-list=<...> */

      else if(!strncmp(argv[nc],"-outcheck=",10)) {
        outcheckname=argv[nc]+10;
      } /* -outcheck=<outimg.nii.gz> */

      else if(!strncmp(argv[nc],"-uthresh2=",10)) {
        uthresh2=atof(argv[nc]+10);
      } /* -uthresh2=<UT2> */

      else if(!strncmp(argv[nc],"-bbox-num=",10)) {
        bbox_num=atoi(argv[nc]+10);
      } /* bbox num */

      else if(!strncmp(argv[nc],"-xfmlist=",9)) {
        if(xfmlist!=NULL) {
          fprintf(stderr,"[%s] ERROR: -xfmlist is already defined\n",fcname);
          exit(1);
        }
        if((CSlist=niik_csstring_to_list(argv[nc]+9,&numxfmlist))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_csstring_to_list %s\n",fcname,argv[nc]);
          exit(1);
        }
        fprintf(stdout,"[%s] xfmlist size = %i\n",fcname,numxfmlist);
        xfmlist = (niikmat **)calloc(numxfmlist,sizeof(niikmat *));
        for(n=0; n<numxfmlist; n++) {
          fprintf(stdout,"[%s] reading matrix:   %s\n",fcname,CSlist[n]);
          if(!strncmp(CSlist[n],"IDENTITY",8)) {
            if((xfmlist[n] = niikmat_identity(4,4))==NULL) {
              fprintf(stderr,"ERROR: niikmat_identity(4,4)\n");
              exit(1);
            }
          } else if((xfmlist[n]=niikmat_read_xfm(CSlist[n]))==NULL) {
            fprintf(stderr,"[%s] ERROR: niikmat_read %s\n",fcname,cptr);
            exit(1);
          }
          free(CSlist[n]);
        }
        free(CSlist);
      } /* -xfmlist=<reg1>,<reg2>... */

      else if(!strncmp(argv[nc],"-bbox=off",9)) {
        flag_bbox = 0;
      }

      else if(!strncmp(argv[nc],"-initocs-",9)) {
        if     (!strncmp(argv[nc],"-initocs-smooth-iter=",21)) {
          iter = atoi(argv[nc]+21);
        } else if(!strncmp(argv[nc],"-initocs-init-thick=",20)) {
          delta = atof(argv[nc]+20);
        } else if(!strncmp(argv[nc],"-initocs-max-thick=",19)) {
          omax = atof(argv[nc]+19);
        } else if(!strncmp(argv[nc],"-initocs-gthresh=",17)) {
          thresh = atof(argv[nc]+17);
        } else if(!strncmp(argv[nc],"-initocs-brain=",15)) {
          niikcortex_tissue_val[3] = atof(argv[nc]+15);
        } else if(!strncmp(argv[nc],"-initocs-step=",14)) {
          delta = atof(argv[nc]+14);
        } else {
          fprintf(stderr,"[niikmath] ERROR: unknown option for initocs\n");
          fprintf(stderr,"                  %s\n",argv[nc]);
          exit(1);
        }
      } /* -initiocs- */

      else if(!strncmp(argv[nc],"-imglist=",9)) {
        if(imglist!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: -imglist is already defined\n");
          exit(1);
        }
        if((CSlist=niik_csstring_to_list(argv[nc]+9,&numimglist))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+9,&numimglist)\n");
          exit(1);
        }
        /*for(n=0;n<numimglist;n++){
          fprintf(stdout,"CS %i: %s\n",n,CSlist[n]); }*/
        fprintf(stdout,"[niikmath] imglist size = %i\n",numimglist);
        imglist = (nifti_image **)calloc(numimglist,sizeof(nifti_image *));
        for(n=0; n<numimglist; n++) {
          fprintf(stdout,"[niikmath] reading image:    %s\n",CSlist[n]);
          if((imglist[n]=niik_image_read(CSlist[n]))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",CSlist[n]);
            exit(1);
          }
        }
        niik_free_list(CSlist,numimglist);
      } /* -imglist=<reg1>,<reg2>... */

      else if(!strncmp(argv[nc],"-uthresh=",9)) {
        uthresh=atof(argv[nc]+9);
      } /* -uthresh=<uthresh> */
      else if(!strncmp(argv[nc],"-ithresh=",9)) {
        ithresh=atof(argv[nc]+9);
      } /* -ithresh=<ithresh> */

      else if(!strncmp(argv[nc],"-nmi-num=",9) ||
              !strncmp(argv[nc],"-nmi-dim=",9)) {
        if((CSlist=niik_csstring_to_list(argv[nc]+9,&n))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+9,&n)\n");
          exit(1);
        }
        if(n!=2) {
          fprintf(stderr,"[niikmath] ERROR: 2 need numbers\n");
          exit(1);
        }
        nmi_num[0]=atoi(CSlist[0]);
        nmi_num[1]=atoi(CSlist[1]);
        niik_free_list(CSlist,n);
      } /* nmi-num */

      else if(!strncmp(argv[nc],"-radius2=",9)) {
        radius2=atof(argv[nc]+9);
      } /* -radius=<radius.nii.gz> */

      else if(!strncmp(argv[nc],"-objlist=",9)) {
        if(objlist!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: '-objlist' is already defined\n");
          exit(1);
        }
        if((CSlist=niik_csstring_to_list(argv[nc]+9,&numobjlist))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n");
          exit(1);
        }
        objlist = (kobj **)calloc(numobjlist,sizeof(kobj *));
        for(n=0; n<numobjlist; n++) {
          fprintf(stdout,"[niikmath] reading object:   %s\n",CSlist[n]);
          if((objlist[n]=off_kobj_read_offply(CSlist[n]))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: off_kobj_read_off\n");
            exit(1);
          }
        }
        niik_free_list(CSlist,numobjlist);
      } /* -vlist=<v1>,<v2>... */

      else if(!strncmp(argv[nc],"-postmat=",9)) {
        if(postmat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: post-matrix is already defined\n");
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading matrix:   %s\n",argv[nc]+9);
        if((postmat=niikmat_read(argv[nc]+9))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niikmat_read\n");
          exit(1);
        }
      } /* -postmat=<postmat.matrix> */
      else if(!strncmp(argv[nc],"-premat=",8)) {
        if(premat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: pre-matrix is already defined\n");
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading matrix:   %s\n",argv[nc]+8);
        if((premat=niikmat_read(argv[nc]+8))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niikmat_read\n");
          exit(1);
        }
      } /* -premat=<premat.matrix> */

      else if(!strncmp(argv[nc],"-matrix=I",9)) {
        NIIK_EXIT((afmat!=NULL),fcname,"afmat is already defined",9);
        NIIK_EXIT(((afmat = niikmat_identity(4,4))==NULL),fcname,"niikmat_identity",9);
      } /* -matrix=I */

      else if(!strncmp(argv[nc],"-degree=",8)) {
        ndegree = atoi(argv[nc]+8);
      } /* -degree=<nd> */

      else if(!strncmp(argv[nc],"-diamond",8)) {
        morph_kernel_shape = NIIK_MORPH_3D_KERNEL_DIAMOND;
      }

      else if(!strncmp(argv[nc],"-matrix=",8)) {
        if(afmat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: matrix is already defined\n");
          exit(1);
        }
        cptr=argv[nc]+8;
        if(!strncmp(cptr,"IDENTITY",8)) {
          if((matlist[n] = niikmat_identity(4,4))==NULL) {
            fprintf(stderr,"ERROR: niikmat_identity(4,4)\n");
            exit(1);
          }
        } else {
          fprintf(stdout,"[niikmath] reading matrix:   %s\n",cptr);
          if((afmat=niikmat_read(cptr))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niikmat_read\n");
            exit(1);
          }
        } /* not identity, actually read the matrix */
      } /* -matrix=<matrix.txt> */

      else if(!strncmp(argv[nc],"-xfm=",5)) {
        if(afmat!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: matrix is already defined\n");
          exit(1);
        }
        cptr=argv[nc]+5;
        if(!strncmp(cptr,"IDENTITY",8)) {
          if((matlist[n] = niikmat_identity(4,4))==NULL) {
            fprintf(stderr,"ERROR: niikmat_identity(4,4)\n");
            exit(1);
          }
        } else {
          fprintf(stdout,"[niikmath] reading xfm matrix:   %s\n",cptr);
          if((afmat=niikmat_read_xfm(cptr))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: niikmat_read\n");
            exit(1);
          }
        } /* not identity, actually read the matrix */
      } /* -xfm=<xfm.xfm> */

      else if(!strncmp(argv[nc],"-affpar=",8)) {
        if((CSlist=niik_csstring_to_list(argv[nc]+8,&n))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list\n");
          exit(1);
        }
        if(n==6 || n==9 || n==12 || n==15) {
          for(m=0; m<n; m++) {
            if(m>=7)
              affpar[m+2] = atof(CSlist[m]);
            else
              affpar[m+1] = atof(CSlist[m]);
          }
        } else {
          fprintf(stderr,"ERROR: affpar is invalid\n");
          fprintf(stderr,"       -affpar requires 6,9,12, or 15 parameters\n");
          fprintf(stderr,"        for 3 rotations, 3 translation, 3 scaling,\n");
          fprintf(stderr,"        3 shears, and 3 centers\n");
          exit(1);
        }
        niik_free_list(CSlist,n);
      } /* affpar */

      else if(!strncmp(argv[nc],"-thresh=",8)) {
        thresh=atof(argv[nc]+8);
      } /* -thresh=<thresh> */

      else if(!strncmp(argv[nc],"-radius=",8)) {
        radius=atof(argv[nc]+8);
      } /* -radius=<radius.nii.gz> */

      else if(!strncmp(argv[nc],"-percent=",9)) {
        percent=atof(argv[nc]+9);
      } /* output max */

      else if(!strncmp(argv[nc],"-refmask=",9)) {
        cptr=argv[nc]+9;
        if(maskref!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: maskref is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading ref mask  %s\n",cptr);
        if((maskref=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -refmask=<maskref.nii.gz> */

      else if(!strncmp(argv[nc],"-inname=",8)) {
        inname=argv[nc]+8;
      } /* -inname=<in_name> */

      else if(!strncmp(argv[nc],"-kernel=",8)) {
        kernel = atoi(argv[nc]+8);
      } /* kernel */

      else if(!strncmp(argv[nc],"-invxfm=",8)) {
        xfmflag=1;
        fprintf(stdout,"[%s] reading xfm (inv) %s\n",fcname,argv[nc]+8);
        NIIK_EXIT(((afmat = niikmat_read_xfm(argv[nc]+8))==NULL),
                  fcname,"niikmat_read_xfm",9);
        NIIK_EXIT((!niikmat_inverse_update(afmat)),
                  fcname,"niikmat_inverse_update",9);
      } /* -invxfm=<xfm.xfm> */

      else if(!strncmp(argv[nc],"-bspline",8)) {
        interp=NIIK_INTERP_BSPLINE;
      }

      else if(!strncmp(argv[nc],"-float16",8)) {
        datatype=NIFTI_TYPE_FLOAT16;
      } else if(!strncmp(argv[nc],"-float32",8)) {
        datatype=NIFTI_TYPE_FLOAT32;
      } else if(!strncmp(argv[nc],"-float64",8)) {
        datatype=NIFTI_TYPE_FLOAT64;
      } else if(!strncmp(argv[nc],"-double",7)) {
        datatype=NIFTI_TYPE_FLOAT64;
      } else if(!strncmp(argv[nc],"-uint64",7)) {
        datatype=NIFTI_TYPE_UINT64;
      } else if(!strncmp(argv[nc],"-uint32",7)) {
        datatype=NIFTI_TYPE_UINT32;
      } else if(!strncmp(argv[nc],"-uint16",7)) {
        datatype=NIFTI_TYPE_UINT16;
      }

      else if(!strncmp(argv[nc],"-fsasc=",7)) {
        fsasc_name=argv[nc]+7;
      } /* fsasc */

      else if(!strncmp(argv[nc],"-delta=",7)) {
        delta=atof(argv[nc]+7);
      } /* fsasc */

      else if(!strncmp(argv[nc],"-linear",7)) {
        interp=NIIK_INTERP_LINEAR;
      }

      else if(!strncmp(argv[nc],"-vlist=",7)) {
        if(vlist!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: '-vlist' is already defined\n");
          exit(1);
        }
        cptr=argv[nc]+7;
        numvlist=1;
        while((cptr=strchr(cptr,','))!=NULL) {
          numvlist++;
          cptr[0] = '\0';
          cptr++;
        }
        fprintf(stdout,"[niikmath] vlist size = %i\n",numvlist);
        vlist = (double *)calloc(numvlist,sizeof(double));
        cptr = argv[nc]+7;
        for(n=0; n<numvlist; n++) {
          vlist[n]=atof(cptr);
          if((cptr = strchr(cptr,'\0'))==NULL) {
            fprintf(stderr,"[niikmath] ERROR: missing next value\n");
            exit(1);
          }
          cptr++;
          /*fprintf(stdout,"vlist %i  %f\n",n,vlist[n]);*/
        }
      } /* -vlist=<v1>,<v2>... */

      else if(!strncmp(argv[nc],"-phase=",7)) {
        phase=atof(argv[nc]+7);
      } /* -phase=<value> */

      /*else if(!strncmp(argv[nc],"-inmgz=",7)){
      cptr=argv[nc]+7;
      if(img!=NULL){
      fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
      exit(1); }
      fprintf(stdout,"\treading           %s\n",cptr);
      if((mgzimg=MRIread(cptr))==NULL){
      fprintf(stderr,"[niikmath] ERROR: mri_read\n");
      exit(1); }
      fprintf(stdout,"\t    read the file!\n");
      } *//* -inmgz=<img.nii.gz> */

      else if(!strncmp(argv[nc],"-mean=",6)) {
        dmean=atof(argv[nc]+6);
      } /* -dmean=<dmean> */
      else if(!strncmp(argv[nc],"-stdv=",6)) {
        dstdv=atof(argv[nc]+6);
      } /* -dstdv=<dstdv> */

      else if(!strncmp(argv[nc],"-debug",6)) {
        debug=1;
      } /* -debug */

      else if(!strncmp(argv[nc],"-step=",6)) {
        step=atoi(argv[nc]+6);
      } /* -debug */

      else if(!strncmp(argv[nc],"-iter=",6)) {
        iter=atoi(argv[nc]+6);
      } /* -iter=<N> */

      else if(!strncmp(argv[nc],"-FWHM=",6)) {
        FWHM=atof(argv[nc]+6);
      } /* -FWHM=<fwhm> */
      else if(!strncmp(argv[nc],"-fwhm=",6)) {
        FWHM=atof(argv[nc]+6);
      } /* -FWHM=<fwhm> */

      else if(!strncmp(argv[nc],"-mask=",6)) {
        cptr=argv[nc]+6;
        if(maskimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: maskimg is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading mask      %s\n",cptr);
        if((maskimg=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -mask=<mask.nii.gz> */

      else if(!strncmp(argv[nc],"-uint8",6)) {
        datatype=NIFTI_TYPE_UINT8;
      } else if(!strncmp(argv[nc],"-int64",6)) {
        datatype=NIFTI_TYPE_INT64;
      } else if(!strncmp(argv[nc],"-int32",6)) {
        datatype=NIFTI_TYPE_INT32;
      } else if(!strncmp(argv[nc],"-int16",6)) {
        datatype=NIFTI_TYPE_INT16;
      }

      else if(!strncmp(argv[nc],"-imin=",6)) {
        imin=atof(argv[nc]+6);
      } /* input min */
      else if(!strncmp(argv[nc],"-imax=",6)) {
        imax=atof(argv[nc]+6);
      } /* input max */

      else if(!strncmp(argv[nc],"-omin=",6)) {
        omin=atof(argv[nc]+6);
      } /* output min */
      else if(!strncmp(argv[nc],"-omax=",6)) {
        omax=atof(argv[nc]+6);
      } /* output max */

      else if(!strncmp(argv[nc],"-elen=",6)) {
        elen=atof(argv[nc]+6);
      } /* -elen=<edge_length> */

      else if(!strncmp(argv[nc],"-img2=",6)) {
        cptr=argv[nc]+6;
        if(img2!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img2 is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading           %s\n",cptr);
        if((img2=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -img=<img.nii.gz> */

      else if(!strncmp(argv[nc],"-img=",5)) {
        cptr=argv[nc]+5;
        if(img!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading           %s\n",cptr);
        if((img=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",cptr);
          exit(1);
        }
      } /* -img=<img.nii.gz> */

      else if(!strncmp(argv[nc],"-xfm=",5)) {
        xfmflag=1;
        fprintf(stdout,"[%s] reading xfm       %s\n",fcname,argv[nc]+5);
        NIIK_EXIT(((afmat = niikmat_read_xfm(argv[nc]+5))==NULL),fcname,"niikmat_read_xfm",9);
      } /* -xfm=<xfm.xfm> */

      else if(!strncmp(argv[nc],"-dof=",5)) {
        dof=atof(argv[nc]+5);
      } /* -dof=<dof> */

      else if(!strncmp(argv[nc],"-num=",5)) {
        num=atoi(argv[nc]+5);
      } /* -num=<N> */

      else if(!strncmp(argv[nc],"-ref=",5)) {
        cptr=argv[nc]+5;
        if(refimg!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: img is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading refimg    %s\n",cptr);
        if((refimg=niik_image_read(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: nifti_image_read\n");
          exit(1);
        }
      } /* reference image */

      else if(!strncmp(argv[nc],"-obj=",5)) {
        cptr=argv[nc]+5;
        if(obj!=NULL) {
          fprintf(stderr,"[niikmath] ERROR: obj is already defined, %s\n",cptr);
          exit(1);
        }
        fprintf(stdout,"[niikmath] reading object    %s\n",cptr);
        if((obj=off_kobj_read_offply(cptr))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: off_kobj_read_off\n");
          exit(1);
        }
      } /* object */

      else if(!strncmp(argv[nc],"-rgb=",5)) {
        if(rgb[0]>0) {
          fprintf(stderr,"[niikmath] ERROR: -rgb is already defined\n");
          exit(1);
        }
        rgb[0]=1;
        if((CSlist=niik_csstring_to_list(argv[nc]+5,&m))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n");
          exit(1);
        }
        if(m!=3) {
          fprintf(stderr,"[niikmath] ERROR: rgb should have 3 numnbers\n");
          exit(1);
        }
        for(n=0; n<m; n++) {
          rgb[n+1]=atof(CSlist[n]);
          /*fprintf(stdout,"%i %f\n",n,rgb[n+1]);*/
        }
        niik_free_list(CSlist,m);
      } /* -rgb=<R>,<G>,<B> */

      else if(!strncmp(argv[nc],"-ijk=",5)) {
        /* get voxel location */
        if(ijk[0]>0) {
          fprintf(stderr,"[niikmath] ERROR: -ijk is already defined\n");
          exit(1);
        }
        ijk[0]=1;
        if((CSlist=niik_csstring_to_list(argv[nc]+5,&m))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n");
          exit(1);
        }
        if(m!=3) {
          fprintf(stderr,"[niikmath] ERROR: ijk should have 3 numnbers\n");
          exit(1);
        }
        for(n=0; n<m; n++) {
          ijk[n+1]=atof(CSlist[n]);
        }
        niik_free_list(CSlist,m);
      } /* -ijk=<i>,<j>,<k> */

      else if(!strncmp(argv[nc],"-xyz=",5)) {
        /* get voxel location */
        if(xyz[0]>0) {
          fprintf(stderr,"[niikmath] ERROR: -xyz is already defined\n");
          exit(1);
        }
        xyz[0]=1;
        if((CSlist=niik_csstring_to_list(argv[nc]+5,&m))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n");
          exit(1);
        }
        if(m!=3) {
          fprintf(stderr,"[niikmath] ERROR: xyz should have 3 numnbers\n");
          exit(1);
        }
        for(n=0; n<m; n++) {
          xyz[n+1]=atof(CSlist[n]);
        }
        niik_free_list(CSlist,m);
      } /* -xyz=<x>,<y>,<z> */

      else if(!strncmp(argv[nc],"-int8",5)) {
        datatype=NIFTI_TYPE_INT8;
      }

      else if(!strncmp(argv[nc],"-out=",5)) {
        outname=argv[nc]+5;
      } /* -out=<out.nii.gz> */

      else if(!strncmp(argv[nc],"-val=",5)) {
        inval=atof(argv[nc]+5);
      } /* -val=<value> */

      else if(!strncmp(argv[nc],"-roi=",5)) {
        if((CSlist=niik_csstring_to_list(argv[nc]+5,&m))==NULL) {
          fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n");
          exit(1);
        }
        if(m!=6) {
          fprintf(stderr,"[niikmath] ERROR: roi should have 6 numnbers\n");
          fprintf(stderr,"                  xmin ymin zmin xmax ymax zmax\n");
          exit(1);
        }
        img_dim_min[0]=img_dim_max[0]=1;
        img_dim_min[1]=atoi(CSlist[0]);
        img_dim_min[2]=atoi(CSlist[1]);
        img_dim_min[3]=atoi(CSlist[2]);
        img_dim_max[1]=atoi(CSlist[3]);
        img_dim_max[2]=atoi(CSlist[4]);
        img_dim_max[3]=atoi(CSlist[5]);
        niik_free_list(CSlist,m);
      } /* ROI */

      else if(!strncmp(argv[nc],"-TE=",4)) {
        if((CSlist=niik_csstring_to_list(argv[nc]+4,&m))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n",fcname);
          exit(1);
        }
        num_MRI_TE=m;
        MRI_TE=(double *)calloc(num_MRI_TE,sizeof(double));
        for(m=0; m<num_MRI_TE; m++) {
          MRI_TE[m] = atof(CSlist[m]);
        }
        fprintf(stdout,"[niikmath] TE: ");
        for(m=0; m<num_MRI_TE; m++) {
          fprintf(stdout,"%12.3f ",MRI_TE[m]);
        }
        fprintf(stdout,"\n");
        niik_free_list(CSlist,num_MRI_TE);
      } /* TE */

      else if(!strncmp(argv[nc],"-TR=",4)) {
        if((CSlist=niik_csstring_to_list(argv[nc]+4,&m))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_csstring_to_list(argv[nc]+5,&m)\n",fcname);
          exit(1);
        }
        num_MRI_TR=m;
        MRI_TR=(double *)calloc(num_MRI_TR,sizeof(double));
        for(m=0; m<num_MRI_TR; m++) {
          MRI_TR[m] = atof(CSlist[m]);
        }
        fprintf(stdout,"[niikmath] TR: ");
        for(m=0; m<num_MRI_TR; m++) {
          fprintf(stdout,"%12.3f ",MRI_TR[m]);
        }
        fprintf(stdout,"\n");
        niik_free_list(CSlist,num_MRI_TR);
      } /* TR */

      else if(!strncmp(argv[nc],"-in=",4)) {
        cptr=argv[nc]+4;
        fprintf(stdout,"[%s] reading           %s\n",fcname,cptr);
        NIIK_EXIT((img!=NULL),fcname,"img is already opened",9);
        NIIK_EXIT(((img=niik_image_read(cptr))==NULL),fcname,"niik_image_read",9);
      } /* -in=<img.nii.gz> */

      else if(!strncmp(argv[nc],"-rx=",4)) {
        affpar[1]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-ry=",4)) {
        affpar[2]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-rz=",4)) {
        affpar[3]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-tx=",4)) {
        affpar[4]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-ty=",4)) {
        affpar[5]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-tz=",4)) {
        affpar[6]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-sg=",4)) {
        affpar[7]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-sx=",4)) {
        affpar[8]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-sy=",4)) {
        affpar[9]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-sz=",4)) {
        affpar[10]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-kx=",4)) {
        affpar[11]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-ky=",4)) {
        affpar[12]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-kz=",4)) {
        affpar[13]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-cx=",4)) {
        affpar[14]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-cy=",4)) {
        affpar[15]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-cz=",4)) {
        affpar[16]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-px=",4)) {
        affpar[14]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-py=",4)) {
        affpar[15]=atof(argv[nc]+4);
      } else if(!strncmp(argv[nc],"-pz=",4)) {
        affpar[16]=atof(argv[nc]+4);
      }

      else if(!strncmp(argv[nc],"-nmi",4) ||
              !strncmp(argv[nc],"-NMI",4)) {
        areg_cost=NIIK_REGISTER_NMI;
      } /* set affine registration's cost method */

      else if(!strncmp(argv[nc],"-cc",3) ||
              !strncmp(argv[nc],"-CC",3)) {
        areg_cost=NIIK_REGISTER_CC;
      } /* set affine registration's cost method */

      else if(!strncmp(argv[nc],"-nn",3)) {
        interp=NIIK_INTERP_NN;
      }

      else if(!strncmp(argv[nc],"-h",2)) {
        fprintf(stdout,"%s",*niikmath_option_list);
        exit(1);
      }

      else if(!strncmp(argv[nc],"-u",6)) {
        fprintf(stdout,"%s",*niikmath_option_list);
        exit(1);
      }

      else {
        fprintf(stderr,"[niikmath] ERROR: unknown option %s\n",argv[nc]);
        exit(1);
      }
      nc++;
    }

    else {
      argv[sc++]=argv[nc++];
    }
  }
  argc=sc;

  /*img=niik_read_minc(argv[1]);
    exit(0);*/


  if(argc!=2) {
    fprintf(stderr,"[niikmath] ERROR: wrong usage for %s\n",argv[0]);
    fprintf(stderr,"           possibly missing operations\n");
    fprintf(stderr,"           please type --list for possible operations\n");
    exit(1);
  }


  /**********************************************
   *
   * main processing
   *
   **********************************************/

  if(!strncmp(argv[1],".",1)) {

  } /* nothing */

  else if(!strncmp(argv[1],"add-off-color-sphere-theta",25)) {
    fprintf(stdout,"[%s] add-off-color-sphere-theta\n",fcname);
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: add-off-color-sphere-theta -obj=<obj.off> -out=<out.off>\n");
      exit(1);
    }
    vec=niikvec_init(obj->nvert);
    for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
      vec->v[n]=v->sph.the;
    }
    NIIK_EXIT((!off_kobj_apply_color_map(obj,vec->v,0,NIIK_PI,NIIK_COLORMAP_SPECTRAL)),fcname,"off_kobj_apply_color_map",1);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,obj,0)),fcname,"off_kobj_write_off",1);
    obj=off_kobj_free(obj);
    exit(0);
  } /* add-off-color-sphere-theta */

  else if(!strncmp(argv[1],"add-off-color-sphere-psi",23)) {
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: add-off-color-sphere-psi -obj=<obj.off> -out=<out.off>\n");
      exit(1);
    }
    vec=niikvec_init(obj->nvert);
    for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
      vec->v[n]=v->sph.psi;
    }
    NIIK_EXIT((!off_kobj_apply_color_map(obj,vec->v,0,NIIK_PI2,NIIK_COLORMAP_SPECTRAL)),fcname,"off_kobj_apply_color_map",1);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,obj,0)),fcname,"off_kobj_write_off",1);
    obj=off_kobj_free(obj);
    exit(0);
  } /* add-off-color-sphere-psi */

  else if(!strncmp(argv[1],"add-off-color-curvature",23) ||
          !strncmp(argv[1],"add-obj-color-curvature",23)) {
    fprintf(stdout,"[niikmath] curvature color map\n");
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -out=<out.off>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optoinal usage:\n");
      fprintf(stdout,"  -iter=<i>     : number of local neighbors [default=3]\n");
      fprintf(stdout,"  -val=<v>      : range of values for coloring [default=10]\n");
      exit(1);
    }
    vec=niikvec_init(obj->nvert);
    if(iter<0) iter=3;
    if(niik_check_double_problem(inval)) inval=10;
    fprintf(stdout,"[niikmath]   local nei = %i\n",iter);
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    if(!niik_off_curvature_map_update(obj,vec,iter)) {
      fprintf(stderr,"ERROR: niik_off_curvature_map_update\n");
      exit(1);
    }
    fprintf(stdout,"             min,avg,max = %.4f %.4f %.4f\n",
            niik_get_min_from_double_vector(vec->v,vec->num),
            niik_get_mean_from_double_vector(vec->v,vec->num),
            niik_get_max_from_double_vector(vec->v,vec->num));
    fprintf(stdout,"[niikmath]   range = %.4f\n",inval);
    if(!off_kobj_apply_color_map(obj,vec->v,-inval,inval,NIIK_COLORMAP_SPECTRAL)) {
      fprintf(stderr,"ERROR: off_kobj_apply_color_map\n");
      exit(1);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[%s] ERROR: off_kobj_write_offply(outname,obj,0)\n",fcname);
      exit(1);
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* add-off-color-curvature */

  else if(!strncmp(argv[1],"modify-header-scl_slope",23)) {
    fprintf(stdout,"[%s] modify-header-scl_slope\n",fcname);
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: modify-header-scl_slope -in=<img.nii> -out=<out.nii> -val=<new_scl_slope>\n");
      exit(1);
    }
    img->scl_slope=inval;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_iamge_write",9);
    exit(0);
  } /* OP = modify-header-scl_slope */

  else if(!strncmp(argv[1],"map-jacobian-overlay",20)) {
    if(imglist==NULL || vlist==NULL || outname==NULL) {
      fprintf(stderr,"[%s] usage: niikmath map-jacobian-overlay -imglist=<jdet.mnc>,<img.mnc> -out=<out.nii> -vlist=<min>,<max>\n",fcname);
      fprintf(stderr,"[%s] optional usage:\n",fcname);
      fprintf(stderr,"  -mask=<mask.nii>          mask image\n");
      fprintf(stderr,"  -imin=<min>               min intensity value [default=image's min]\n");
      fprintf(stderr,"  -imax=<max>               max intensity value [default=image's max]\n");
      fprintf(stderr,"  -val=<bg_factor>          factor for background [default = 0.2, range 0 - 1]\n");
      fprintf(stderr,"                            0 = black, 1 = image\n");
      fprintf(stderr,"  -thresh=<thresh>          threshold for mask image [default = none]\n");
      exit(1);
    }
    cmap=niik_colormap_get_jacobian(NIIK_VAL_USMAX);
    if(niik_check_double_problem(imin)) imin=niik_image_get_min(imglist[1],NULL);
    if(niik_check_double_problem(imax)) imax=niik_image_get_max(imglist[1],NULL);
    if(niik_check_double_problem(inval)) inval=0.2;
    if(!niik_check_double_problem(thresh)) {
      if(maskimg!=NULL) {
        fprintf(stdout,"[%s] threshold image %12.6f\n",fcname,thresh);
        NIIK_EXIT((!niik_image_threshold(maskimg,thresh)),fcname,"niik_image_threshold",1);
      }
    }
    NIIK_EXIT((!niik_image_type_convert(imglist[0],NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);
    NIIK_EXIT(((outimg=niik_image_map_color_overlay_mask(imglist[1],imin,imax,imglist[0],vlist[0],vlist[1],inval,maskimg,cmap))==NULL),fcname,"niik_image_map_color",9);
    for(n=0; n<2; n++) imglist[n]=niik_image_free(imglist[n]);
    free(imglist);
    if(datatype>=0) NIIK_EXIT((!niik_image_type_convert(outimg,datatype)),fcname,"niik_image_type_convert",1);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_iamge_write",9);
    outimg=niik_image_free(outimg);
    exit(0);
  } /* map-jacobian-overlay */

  else if(!strncmp(argv[1],"map-spectral",12)) {
    if(img==NULL || vlist==NULL || outname==NULL) {
      fprintf(stderr,"[%s] usage: niikmath map-spectral -in=<img.nii> -out=<out.nii> -vlist=<min>,<max>\n",fcname);
      exit(1);
    }
    cmap=niik_colormap_get_spectral(NIIK_VAL_USMAX);
    NIIK_EXIT((!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);
    NIIK_EXIT(((outimg=niik_image_map_color(img,vlist[0],vlist[1],cmap))==NULL),fcname,"niik_image_map_color",9);
    img=niik_image_free(img);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_iamge_write",9);
    outimg=niik_image_free(outimg);
    exit(0);
  } /* create spectral image */

  else if(!strncmp(argv[1],"map-jacobian",12)) {
    if(img==NULL || vlist==NULL || outname==NULL) {
      fprintf(stderr,"[%s] usage: niikmath map-jacobian -in=<img.nii> -out=<out.nii> -vlist=<min>,<max>\n",fcname);
      exit(1);
    }
    cmap=niik_colormap_get_jacobian(NIIK_VAL_USMAX);
    NIIK_EXIT((!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);
    NIIK_EXIT(((outimg=niik_image_map_color(img,vlist[0],vlist[1],cmap))==NULL),fcname,"niik_image_map_color",9);
    img=niik_image_free(img);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_iamge_write",9);
    outimg=niik_image_free(outimg);
    exit(0);
  } /* create jacobian image */

  else if(!strncmp(argv[1],"interpacket-correction",23)) {
    if(img==NULL || num==0 || outname==NULL) {
      fprintf(stdout,"  usage: %s interpacket-correction -in=<img.nii> -num=<#packets> -out=<out.nii>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -dof=<DOF>    : registration degrees of freedom [default=6]\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  required args:\n");
      fprintf(stdout,"  <img.nii>     : input image\n");
      fprintf(stdout,"  <#packets>    : # of packet\n");
      fprintf(stdout,"  <out.nii>     : output image\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  example:\n");
      fprintf(stdout,"  %s interpacket-correction -in=ASSERT_104-HMR-1_f-p_910402_m0_t1p_IPM01.nii.gz -num=2 -out=ASSERT_104_HMR-1_f-p_910402_m0_t1p_IPM03.nii.gz\n",argv[0]);
      exit(1);
    }
    if((outimg=niik_image_interpacket_motion_correction(img,num,dof))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_interpacket_motion_correction\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    exit(0);
  } /* interpacket-correction */

  /* THIS FUNCTION DOES NOT WORK YET
  else if(!strncmp(argv[1],"cuberille-prep-topo",19)){
    if(maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -mask=<mask.nii> -out=<out.off>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:");
      fprintf(stdout,"  -val=0       : remove voxels instead of adding");
      exit(1); }
    if(!niik_image_correct_topology_for_cuberille(maskimg)){
      fprintf(stderr,"[niikmath] ERROR: niik_image_correct_topology_for_cuberille\n");
      exit(1); }
    if(niik_check_double_problem(inval)) inval=1;
    if(!niik_image_correct_for_cuberille(maskimg,(int)inval)){
      fprintf(stderr,"[niikmath] ERROR: niik_image_correct_for_cuberille\n");
      exit(1); }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niik_image_write(outname,maskimg)){
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1); }
    exit(0);
  }  OP = cuberille-prep-topo */

  else if(!strncmp(argv[1],"interpacket-split",18)) {
    if(img==NULL || num==0 || outname==NULL) {
      fprintf(stdout,"  usage: %s interpacket-split -in=<img.nii> -num=<#packets> -out=<out>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  required args:\n");
      fprintf(stdout,"  <img.nii>     : input image\n");
      fprintf(stdout,"  <#packets>    : # of packet\n");
      fprintf(stdout,"  <out>         : output prefix\n");
      exit(1);
    }
    if((imglist=niik_image_split_interpacket(img,num))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_split_interpacket\n",fcname);
      exit(1);
    }
    for(n=0; n<num; n++) {
      sprintf(fname,"%s-%i.nii.gz",outname,n+1);
      fprintf(stdout,"[%s] writing output  %s\n",fcname,fname);
      niik_image_append_history(imglist[n],timestamp);
      if(!niik_image_write(fname,imglist[n])) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,fname);
        exit(1);
      }
    } /* num */
    exit(0);
  } /* OP = interpacket-split */

  else if(!strncmp(argv[1],"interpacket-merge",18)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s interpacket-merge -imglist=<img1.nii>,<img2.nii>[,<img3.nii>...] -out=<out.nii>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  required args:\n");
      fprintf(stdout,"  <img1.nii>...   : input images \n");
      fprintf(stdout,"  <out.nii>       : output image\n");
      exit(1);
    }
    if((outimg=niik_image_merge_interpacket(imglist,numimglist))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_merge_interpacket\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    exit(0);
  } /* OP = interpacket-merge */

  else if(!strncmp(argv[1],"modify-header-dx",17)) {
    fprintf(stdout,"[%s] modify-header-dx\n",fcname);
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: modify-header-scl_slope -in=<img.nii> -out=<out.nii> -val=<dx>\n");
      exit(1);
    }
    img->dx=img->pixdim[1]=inval;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_iamge_write",9);
    exit(0);
  } /* OP = modify-header-dx */

  else if(!strncmp(argv[1],"modify-header-dy",17)) {
    fprintf(stdout,"[%s] modify-header-dy\n",fcname);
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: modify-header-scl_slope -in=<img.nii> -out=<out.nii> -val=<dy>\n");
      exit(1);
    }
    img->dy=img->pixdim[2]=inval;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_iamge_write",9);
    exit(0);
  } /* OP = modify-header-dy */

  else if(!strncmp(argv[1],"modify-header-dz",17)) {
    fprintf(stdout,"[%s] modify-header-dz\n",fcname);
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: modify-header-scl_slope -in=<img.nii> -out=<out.nii> -val=<dz>\n");
      exit(1);
    }
    img->dz=img->pixdim[3]=inval;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_iamge_write",9);
    exit(0);
  } /* OP = modify-header-dz */


  else if(!strncmp(argv[1],"applyaffine-xfm",17)) {
    if(img==NULL || refimg==NULL || outname==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: applyaffine-xfm -in=<img.mnc> -ref=<ref.mnc> -out=<out.mnc> -xfm=<transform.xfm>\n");
      exit(1);
    }
    NIIK_EXIT((xfmflag==0),fcname,"-xfm=<mat.xfm> should be used",9);
    NIIK_EXIT((!niikmat_convert_from_xfm(afmat,img,refimg)),fcname,"niikmat_convert_from_xfm",9);
    NIIK_EXIT((!niik_image_affine_transform_3d_update(img,refimg,afmat,interp)),fcname,"niik_image_affine_transform_3d_update",9);
    fprintf(stdout,"[%s] writing image:    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
    refimg=niik_image_free(refimg);
    afmat=niikmat_free(afmat);
    exit(0);
  } /* applyaffine-xfm */

  else if(!strncmp(argv[1],"applyaffine-inv",17)) {
    if(img==NULL || outname==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<in.matrix> -out=<out.nii> [-ref=<ref.nii>]\n",argv[0],argv[1]);
      exit(1);
    }
    if(refimg==NULL) {
      if((refimg=niik_image_copy(img))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_copy(img)\n");
        exit(1);
      }
    }
    niikmat_inverse_update(afmat);
    fprintf(stdout,"[niikmath] affine transformation\n");
    niikmat_display(afmat);
    if(!niik_image_affine_transform_3d_update(img,refimg,afmat,interp)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update(img,img,afmat,interp)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing image:    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    if(img==refimg)
      img=niik_image_free(img);
    else {
      img=niik_image_free(img);
      refimg=niik_image_free(refimg);
    }
  } /* OP = applyaffine */

  else if(!strncmp(argv[1],"applyaffine-obj",17)) {
    if(obj==NULL || outname==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -matrix=<in.matrix> / -xfm=<in.xfm> -out=<out.off> [-ref=<ref.nii>]\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] affine transformation\n");
    niikmat_display(afmat);
    if(!niik_affine_transform_off(obj,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niik_affine_transform_off(obj,afmat)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(1);
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* OP = applyaffine-obj */

  else if(!strncmp(argv[1],"remove-obj-color",16) ||
          !strncmp(argv[1],"remove-off-color",16)) {
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<obj.off> -out=<out.off>\n",argv[0],argv[1]);
      exit(1);
    }
    if(!off_kobj_remove_color(obj)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_remove_color(obj)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(1);
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* remove-off-color */

  else if(!strncmp(argv[1],"world2voxel-int",15)) {
    if(img==NULL || !(xyz[0]>0)) {
      fprintf(stdout,"  usage: [-options] -in=<img.nii> -xyz=<x>,<y>,<z>\n");
      exit(1);
    }
    pt = niikpt_val(xyz[1],xyz[2],xyz[3],0);
    if(afmat!=NULL) {
      pt = niikpt_affine_transform(afmat,pt);
    }
    pt = niikpt_affine_transform_m44(img->sto_ijk,pt);
    fprintf(stdout,"  voxel: %19i %19i %19i\n",(int)floor(0.5+pt.x),(int)floor(0.5+pt.y),(int)floor(0.5+pt.z));
    img=niik_image_free(img);
    exit(0);
  } /* OP = voxel2world */

  else if(!strncmp(argv[1],"write-ascii-img",15)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.txt>\n",argv[0],argv[1]);
      exit(1);
    }
    if((fp=fopen(outname,"w"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",outname);
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      fprintf(fp,"%12.9f ",niik_image_get_voxel(img,i));
    }
    fclose(fp);
    fprintf(stdout,"%3i %3i %3i %3i %3i %3i\n",img->nx,img->ny,img->nz,img->nt,img->nu,img->nv);
    exit(0);
  }  /* OP = write-ascii-img */

  else if(!strncmp(argv[1],"removeiscaling",14)) {
    fprintf(stdout,"[%s] remove intensity scaling\n",fcname);
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"[%s] usage: removeiscaling -in=<in.nii> -out=<out.nii>\n",fcname);
      exit(1);
    }
    img->scl_slope=img->scl_inter=0.0;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = removeiscaling */

  else if(!strncmp(argv[1],"decompose-affine",14) ||
          !strncmp(argv[1],"affine-decompose",14)) {
    if(afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -matrix=<in.matrix>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"    niikmat_decompose_affine  dof=%i\n",dof);
    niikmat_display(afmat);
    if(!niikmat_decompose_affine(afmat,affpar,dof)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_decompose_affine\n");
      exit(1);
    }
    if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update(afmat,affpar)\n");
      exit(1);
    }
    fprintf(stdout,"      estimated matrix\n");
    niikmat_display(afmat);
    fprintf(stdout,"      r  = %15.9f %15.9f %15.9f\n",affpar[1],affpar[2],affpar[3]); /* in degrees */
    fprintf(stdout,"      t  = %15.9f %15.9f %15.9f\n",affpar[4],affpar[5],affpar[6]);
    fprintf(stdout,"      s  = %15.9f %15.9f %15.9f   includes global scaling of %-15.9f\n",affpar[7]*affpar[8],affpar[7]*affpar[ 9],affpar[7]*affpar[10],affpar[7]);
    fprintf(stdout,"      sh = %15.9f %15.9f %15.9f\n",affpar[11],affpar[12],affpar[13]);
    fprintf(stdout,"      c  = %15.9f %15.9f %15.9f\n",affpar[14],affpar[15],affpar[16]);
    afmat=niikmat_free(afmat);
  } /* OP = affine-decompose */

  else if(!strncmp(argv[1],"cuberille-prep",14)) {
    fprintf(stdout,"[niikmath] preparation for cuberille\n");
    fprintf(stdout,"           checks and removes (or adds) shared corners\n");
    if(maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -mask=<mask.nii> -out=<out.off>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:");
      fprintf(stdout,"  -val=0       : remove voxels instead of adding");
      exit(1);
    }
    if(niik_check_double_problem(inval)) inval=1;
    if(!niik_image_correct_for_cuberille(maskimg,(int)inval)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_correct_for_cuberille\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(maskimg,timestamp);
    if(!niik_image_write(outname,maskimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    exit(0);
  } /* OP = cuberille-prep */

  else if(!strncmp(argv[1],"multiplymatrix",14)) {
    if(matlist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -matrixlist=<in1.matrix>,<in2.matrix>... -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -matrixlist is a list of matrices to multiply\n");
      fprintf(stdout,"  -output = [inN] * [inN-1] ... [in1]\n");
      fprintf(stdout,"  -for more info on matrix, please '--help-matrix'\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    niikmat_display(matlist[0]);
    for(n=1; n<nummatlist; n++) {
      fprintf(stdout,"[niikmath] multiplying\n");
      niikmat_display(matlist[n]);
      if(!niikmat_multiply_mat2_free1(matlist[n],matlist[0])) {
        fprintf(stderr,"[niikmath] ERROR: niikmat_multiply_mat2_free1\n");
        exit(1);
      }
      matlist[n]=NULL;
    }
    fprintf(stdout,"[niikmath] final matrix\n");
    niikmat_display(matlist[0]);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niikmat_write(outname,matlist[0])) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write %s\n",outname);
      exit(1);
    }
    matlist[0]=niikmat_free(matlist[0]);
    free(matlist);
    matlist=NULL;
    nummatlist=0;
  } /* OP = multiplymatrix */

  else if(!strncmp(argv[1],"add-rice-WRONG",14)) {
    fprintf(stdout,"[%s] add-rice: adding Rician noise (this was the wrong version)\n",fcname);
    if(img==NULL || outname==NULL || vlist==NULL ) {
      fprintf(stdout,"[%s] usage: -in=<in.nii> -out=<out.nii> -vlist=<v>,<s>\n",fcname);
      exit(1);
    }
    NIIK_EXIT(((dimg=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector",9);
    if(!strncmp(argv[1],"add-rice-random",15)) {
      fprintf(stdout,"[%s] real pseudo-random generator     %.6f\n",fcname,inval);
      for(i=0; i<img->nvox; i++) {
        dimg[i]+=NIIK_RiceRnd2(vlist[0],vlist[1]);
      }
    } else if(!strncmp(argv[1],"add-rice",8)) {
      fprintf(stdout,"[%s] reproducible pseudo-random generator    %.6f %.6f\n",fcname,vlist[0],vlist[1]);
      for(i=0; i<img->nvox; i++) {
        dimg[i]+=NIIK_RiceRnd(vlist[0],vlist[1]);
      }
    }
    NIIK_EXIT((!niik_image_set_voxels_from_double_vector(img,dimg)),fcname,"niik_image_set_voxels_from_double_vector",9);
    free(dimg);
    if(datatype>=0) {
      NIIK_EXIT((!niik_image_type_convert(outimg,datatype)),fcname,"niik_image_type_convert",9);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    outimg=niik_image_free(outimg);

    exit(0);
  } /* add-rice-WRONG */

  else if(!strncmp(argv[1],"halfway-matrix",14) ||
          !strncmp(argv[1],"matrix-halfway",14)) {
    if(afmat==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -matrix=<input.matrix> -out=<fwd.matrix>[,<inv.matrix>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -<input.matrix> is a required input matrix\n");
      fprintf(stdout,"  -<fwd.matrix> is required, but <inv.matrix> is optional\n");
      fprintf(stdout,"  -halfway matrix is defined as :\n");
      fprintf(stdout,"     <input.matrix> = <fwd.matrix> x <fwd.matrix>\n");
      fprintf(stdout,"     <inv.matrix> = INVERSE(<fwd.matrix>)\n");
      fprintf(stdout,"     thus,\n");
      fprintf(stdout,"     <input.matrix> = INVERSE(<inv.matrix>) x <fwd.matrix>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    /* calculate halfway matrix */
    afmat2=niikmat_halfway_matrix2(afmat);
    niikmat_display(afmat2);
    if((CSlist=niik_csstring_to_list(outname,&nummatlist))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: CSlist=niik_csstring_to_list(outname,&nummatlist)\n");
      exit(1);
    }
    switch(nummatlist) {
    case 2:
      afmat=niikmat_inverse(afmat2);
      fprintf(stdout,"[niikmath] writing inverse matrix    %s\n",CSlist[1]);
      if(!niikmat_write(CSlist[1],afmat)) {
        fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
        exit(1);
      }
    case 1:
      fprintf(stdout,"[niikmath] writing forward matrix    %s\n",CSlist[0]);
      if(!niikmat_write(CSlist[0],afmat2)) {
        fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
        exit(1);
      }
      break;
    default:
      fprintf(stderr,"[niikmath] ERROR: need 1 or 2 output matrix filenames\n");
      exit(1);
    }
    exit(0);
    /*if(nummatlist>2) {
      fprintf(stderr,"[niikmath] ERROR: need 1 or 2 output matrix filenames\n");
      exit(1); }
      fprintf(stdout,"[niikmath] output fwd matrix     %s\n",CSlist[0]);
    if(nummatlist==2) fprintf(stdout,"[niikmath] output inv matrix     %s\n",CSlist[1]);
    fprintf(stdout,"  calculate halfway parameters\n");
    if(!niikmat_halfway_matrix(afmat,affpar,dof)){
      fprintf(stderr,"[niikmath] ERROR: niikmat_halfway_matrix(afmat,affpar,dof)\n");
      exit(1); }
    fprintf(stdout,"[niikmath] calculate forward halfway matrix\n");
    if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)){
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update(afmat,affpar)\n");
      exit(1); }
    fprintf(stdout,"[niikmath] writing forward matrix    %s\n",CSlist[0]);
    if(!niikmat_write(CSlist[0],afmat)){
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1); }
    if(nummatlist==2) {
      fprintf(stdout,"[niikmath] calculate inverse halfway matrix\n");
      niikmat_inverse_update(afmat);
      fprintf(stdout,"[niikmath] writing inverse matrix    %s\n",CSlist[1]);
      if(!niikmat_write(CSlist[1],afmat)){
    fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
    exit(1); }
        } */
    nummatlist=0;
    exit(1);
  } /* OP = halfway-matrix */

  else if(!strncmp(argv[1],"correct-noise",13)) {
    fprintf(stdout,"  removes zeros and shift intensities\n");
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    imin=niik_image_get_max(img,NULL);
    if(!niik_check_double_problem(inval))
      imin=inval;
    for(i=0; i<img->nvox; i++) {
      dval=niik_image_get_voxel(img,i);
      if(dval>1e-6)
        imin=NIIK_DMIN(imin,dval);
    }
    fprintf(stdout,"[%s] bg-noise subtraction:  %12.9g\n",fcname,imin);
    for(i=0; i<img->nvox; i++) {
      dval = niik_image_get_voxel(img,i);
      if(dval>imin)
        niik_image_add_voxel(img,i,-imin);
      else
        niik_image_set_voxel(img,i,0);
    }
    fprintf(stdout,"[%s] writing image     %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = correct-noise */

  else if(!strncmp(argv[1],"intratemplate",13)) {
    if(imglist==NULL || outname==NULL || refimg==NULL || maskimg==NULL) {
      fprintf(stdout,"  usage: %s %s -ref=<refimg.nii> -imglist=<img1.nii>,<img2.nii>[,<img3.nii>...] -mask=<refmask.nii> -out=<warp1.nii>,<warp2.nii>[,...]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -matrixlist=<img1.matrix>,<img2.matrix>[,...]  : list of affine registration matrix\n");
      fprintf(stdout,"  -outcheck=<check.nii>           check image\n");
      exit(1);
    }
    /*
     * check and prepare matrix list
     */
    if(nummatlist>0) {
      if(numimglist!=nummatlist) {
        fprintf(stderr,"[%s] ERROR: #image list [%i] and #matrix list [%i] are different\n",fcname,numimglist,nummatlist);
        exit(1);
      }
    } else {
      nummatlist=numimglist;
      matlist=(niikmat **)calloc(nummatlist,sizeof(niikmat *));
      for(n=0; n<nummatlist; n++) {
        matlist[n]=niikmat_identity(4,4);
      }
    }
    if((CSlist=niik_csstring_to_list(outname,&n))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_csstring_to_list %s\n",fcname,outname);
      exit(1);
    }
    if(numimglist!=n) {
      fprintf(stderr,"[%s] ERROR: #image list [%i] and #output list [%i] are different\n",fcname,numimglist,n);
      exit(1);
    }
    for(n=0; n<numimglist; n++) {
      free(CSlist[n]);
    }
    free(CSlist);

    /*
     * prepare output warp image
     */
    warpimglist=(nifti_image **)calloc(numimglist,sizeof(nifti_image *));
    fprintf(stdout,"[%s] prepare output warp image\n",fcname);
    if((warpimglist[0]=niik_image_copy_as_type(refimg,NIFTI_TYPE_FLOAT32))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
      exit(1);
    }
    warpimglist[0]->ndim=warpimglist[0]->dim[0]=5;
    warpimglist[0]->nt=warpimglist[0]->dim[4]=1;
    warpimglist[0]->nu=warpimglist[0]->dim[5]=3;
    warpimglist[0]->nvox=warpimglist[0]->nx*warpimglist[0]->ny*warpimglist[0]->nz*warpimglist[0]->nt*warpimglist[0]->nu;
    free(warpimglist[0]->data);
    warpimglist[0]->data=(float *)calloc(warpimglist[0]->nvox,warpimglist[0]->nbyper);
    for(n=1; n<nummatlist; n++) {
      if((warpimglist[n]=niik_image_copy(warpimglist[0]))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
        exit(1);
      }
    }

    /*
     * template creationg and image warps
     */
    afmat=niikmat_identity(4,4);
    if(iter<0) iter=2;
    /*if(!niik_image_nregister_intrasubject_template_test1(refimg,maskimg,
                                                         warpimglist,imglist,matlist,numimglist)){
      fprintf(stderr,"[%s] ERROR: niik_image_nregister_intrasubject_template_test1\n",fcname);
      exit(1); }*/
    if(!niik_image_nregister_intrasubject_template_test1(refimg,maskimg,
        warpimglist,imglist,matlist,numimglist)) {
      fprintf(stderr,"[%s] ERROR: niik_image_nregister_intrasubject_template_test1\n",fcname);
      exit(1);
    }

    /*
     * post-processing
     */
    if((CSlist=niik_csstring_to_list(outname,&n))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_csstring_to_list %s\n",fcname,outname);
      exit(1);
    }
    for(n=0; n<nummatlist; n++) {
      fprintf(stdout,"[%s] writing image     %s\n",fcname,CSlist[n]);
      niik_image_append_history(warpimglist[n],timestamp);
      if(!niik_image_write(CSlist[n],warpimglist[n])) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,fname);
        exit(1);
      }
      free(CSlist[n]);
    }
    free(CSlist);

    if(outcheckname!=NULL) {
      if(!niik_image_nregister_average_with_fov(refimg,warpimglist,imglist,matlist,numimglist,NIIK_INTERP_BSPLINE)) {
        fprintf(stderr,"[%s] ERROR: niik_image_nregister_average_with_fov\n",fcname);
        exit(1);
      }
      fprintf(stdout,"[%s] writing image     %s\n",fcname,outcheckname);
      niik_image_append_history(refimg,timestamp);
      if(!niik_image_write(outcheckname,refimg)) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outcheckname);
        exit(1);
      }
    } /* outcheckname */
    exit(0);
  } /* OP = intratemplate */

  else if(!strncmp(argv[1],"add-off-color",13) ||
          !strncmp(argv[1],"add-obj-color",13)) {
    if(obj==NULL || outname==NULL || rgb[0]==0) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -out=<out.off> -rgb=<R>,<G>,<B>\n",argv[0],argv[1]);
      /*fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -out=<out.off>       change the output filename\n");*/
      fprintf(stdout,"\n  R G B should have range of [0,1]\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] adding color %8.5f %8.5f %8.5f\n",rgb[1],rgb[2],rgb[3]);
    if(!off_kobj_add_one_color(obj,rgb[1],rgb[2],rgb[3])) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_add_one_color(obj,rgb[1],rgb[2],rgb[3])\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(1);
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* OP = add-off-color */


  else if(!strncmp(argv[1],"remove-faces",12)) {
    fprintf(stdout,"[niikmath] remove-faces\n");
    if(obj==NULL || outname==NULL || !(xyz[0]!=0 || ijk[0]!=0)) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -out=<out.off> [-xyz=<x>,<y>,<z> | -ijk=<i>,<j>,<k>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -val=<D>           : distance from xyz or ijk point in which\n                       the faces are kept [default=3mm]\n");
      fprintf(stdout,"  -removes faces so that the faces around the point [xyz|ijk] is kept\n");
      fprintf(stdout,"  -this is useful to look at small area of surface on geomview.\n");
      exit(1);
    }
    if(xyz[0]==0 && ijk[0]>0 && img!=NULL) {
      xyz[0]=1;
      xyz[1]=ijk[1]*img->dx;
      xyz[2]=ijk[2]*img->dy;
      xyz[3]=ijk[3]*img->dz;
    }
    pt=niikpt_val(xyz[1],xyz[2],xyz[3],0);
    if(niik_check_double_problem(inval)) inval=3.0;
    fprintf(stdout,"[niikmath] start removal %6.2f %6.2f %6.2f\n",pt.x,pt.y,pt.z);
    for(face=obj->face; face!=NULL; face=face->next) {
      if(niikpt_distance(face->vert[0]->v,pt)>inval)
        if(niikpt_distance(face->vert[1]->v,pt)>inval)
          if(niikpt_distance(face->vert[2]->v,pt)>inval) {
            if(debug) {
              fprintf(stdout,"removing %i\n",face->index);
              off_display_face_info(face,20);
            }
            off_kface_remove(face,obj);
            face=obj->face;
          }
    }
    face=obj->face;
    if(niikpt_distance(face->vert[0]->v,pt)>inval)
      if(niikpt_distance(face->vert[1]->v,pt)>inval)
        if(niikpt_distance(face->vert[2]->v,pt)>inval) {
          if(debug) {
            fprintf(stdout,"removing %i\n",face->index);
            off_display_face_info(face,20);
          }
          off_kface_remove(face,obj);
        }
    off_kobj_update_all_index(obj);
    fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
    for(vert=obj->vert; vert!=NULL; vert=vert->next) vert->v.w=0;
    for(face=obj->face; face!=NULL; face=face->next) {
      face->vert[0]->v.w=face->vert[1]->v.w=face->vert[2]->v.w=1;
    }
    for(vert=obj->vert; vert!=NULL; vert=vert->next)
      if(vert->v.w==0)
        off_kvert_remove(vert,obj);
    off_kobj_update_all_index(obj);
    fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
    obj->edge=NULL;
    off_kobj_update_all_index(obj);
    fprintf(stdout,"        vfe %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(1);
    }
    img=niik_image_free(img);
    obj=off_kobj_free(obj);
    exit(0);
  } /* remove-faces */

  else if(!strncmp(argv[1],"removesqform",12)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -removes img's sform and qform (\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] removing qform and sform\n");
    if(!img->sform_code) {
      fprintf(stdout,"sto_xyz: %s\n",niik_nifti1_xform_string(img->sform_code));
      for(i=0; i<4; i++) {
        fprintf(stdout,"  %15.8f %15.8f %15.8f %15.8f\n",img->sto_xyz.m[i][0],img->sto_xyz.m[i][1],img->sto_xyz.m[i][2],img->sto_xyz.m[i][3]);
      }
    }
    if(!img->qform_code) {
      fprintf(stdout,"qto_xyz: %s\n",niik_nifti1_xform_string(img->qform_code));
      for(i=0; i<4; i++) {
        fprintf(stdout,"  %15.8f %15.8f %15.8f %15.8f\n",img->qto_xyz.m[i][0],img->qto_xyz.m[i][1],img->qto_xyz.m[i][2],img->qto_xyz.m[i][3]);
      }
    }
    img->qform_code = 0;
    img->sform_code = 0;
    niik_image_append_history(img,timestamp);

    if(outname!=NULL) {
      fprintf(stdout,"[%s] writing image:    %s\n",fcname,outname);
      NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",1);
    } else {
      sprintf(fname,"%s",img->fname);
      NIIK_RET0((!niik_image_write(fname,img)),fcname,"niik_image_write");
      fprintf(stdout,"[%s] writing image:    %s\n",fcname,fname);
      NIIK_EXIT((!niik_image_write(fname,img)),fcname,"niik_image_write",1);
    }
    img=niik_image_free(img);
  } /* OP = removesqform */

  else if(!strncmp(argv[1],"affinematrix",12)) {
    if(outname==NULL) {
      /*fprintf(stdout,"  usage: %s %s -vlist=<Rx>,<Ry>,<Rz>,<Tx>,<Ty>,<Tz>,<Sx>,<Sy>,<Sz>,<Kx>,<Ky>,<Kz>,<Px>,<Py>,<Pz>\n",argv[0],argv[1]);*/
      fprintf(stdout,"  usage: %s %s -affpar=<Rx>,<Ry>,<Rz>,<Tx>,<Ty>,<Tz>,<Sx>,<Sy>,<Sz>,<Kx>,<Ky>,<Kz>,<Px>,<Py>,<Pz>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  alternative optional usage:\n");
      fprintf(stdout,"  -rx -ry -rz         x,y,z rotation in degrees respectively\n");
      fprintf(stdout,"  -tx -ty -tz         x,y,z translation in mm respectively\n");
      fprintf(stdout,"  -sx -sy -sz -sg     x,y,z,global scaling respectively\n");
      fprintf(stdout,"  -kx -ky -kz         x,y,z skewing respectively\n");
      fprintf(stdout,"  -px -py -pz         x,y,z pivot point in mm from [0,0,0] image space\n");
      exit(1);
    }
    if(vlist!=NULL) {
      for(n=1,m=0; n<=3; n++) affpar[n]=vlist[m++];
      for(n=4; n<=6; n++) affpar[n]=vlist[m++];
      affpar[7]=1;
      for(n=8; n<=16; n++) affpar[n]=vlist[m++];
    }
    niik_aregister_display_affine(affpar);
    if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar\n");
      exit(1);
    }
    niikmat_display(afmat);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1);
    }
    afmat=niikmat_free(afmat);
  } /* OP = affinematrix */

  else if(!strncmp(argv[1],"affine2rigid",12)) {
    if(outname==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -matrix=<in.matrix> -out=<rigid_out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  alternative optional usage:\n");
      fprintf(stdout,"  -px=<PX> -py=<PY> -pz=<PZ>        : x,y,z pivot point in mm from [0,0,0] image space\n");
      fprintf(stdout,"  -dof=<DOF>                        : DOF\n");
      exit(1);
    }
    fprintf(stdout,"[%s] niikmat_decompose_affine dof=%i\n",fcname,dof);
    niikmat_display(afmat);
    if(!niikmat_decompose_affine(afmat,affpar,dof)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_decompose_affine\n");
      exit(1);
    }
    if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update(afmat,affpar)\n");
      exit(1);
    }
    if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update(afmat,affpar)\n");
      exit(1);
    }
    fprintf(stdout,"      estimated matrix\n");
    niikmat_display(afmat);
    fprintf(stdout,"      r  = %15.9f %15.9f %15.9f\n",affpar[1],affpar[2],affpar[3]); /* in degrees */
    fprintf(stdout,"      t  = %15.9f %15.9f %15.9f\n",affpar[4],affpar[5],affpar[6]);
    fprintf(stdout,"      s  = %15.9f %15.9f %15.9f   includes global scaling of %-15.9f\n",affpar[7]*affpar[8],affpar[7]*affpar[ 9],affpar[7]*affpar[10],affpar[7]);
    fprintf(stdout,"      sh = %15.9f %15.9f %15.9f\n",affpar[11],affpar[12],affpar[13]);
    fprintf(stdout,"      c  = %15.9f %15.9f %15.9f\n",affpar[14],affpar[15],affpar[16]);
    affpar[7]=affpar[8]=affpar[9]=affpar[10]=1.0;
    affpar[11]=affpar[12]=affpar[13]=0.0;
    if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update(afmat,affpar)\n");
      exit(1);
    }
    fprintf(stdout,"     Reconstructed matrix:\n");
    niikmat_display(afmat);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1);
    }
    afmat=niikmat_free(afmat);
    exit(0);
  } /* OP = affine2rigid */

  else if(!strncmp(argv[1],"merge-slices",12)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: merge-slices -imglist=<img1.nii>[,<img2.nii>...] -out=<out.nii>\n");
      exit(1);
    }
    fprintf(stdout,"[%s] merging slices\n",fcname);
    NIIK_EXIT(((outimg=niik_image_copy(imglist[0]))==NULL),fcname,"niik_image_copy",9);
    outimg->dim[0]=outimg->ndim=3;
    outimg->dim[3]=outimg->nz=numimglist;
    outimg->nvox=outimg->nx*outimg->ny*outimg->nz;
    num=imglist[0]->nvox;
    free(outimg->data);
    outimg->data=(void *)calloc(outimg->nvox,sizeof(outimg->nbyper));
    for(k=j=0; k<outimg->nz; k++) {
      NIIK_EXIT((num!=imglist[k]->nvox),fcname,"2d image size did not match",9);
      for(i=0; i<imglist[k]->nvox; j++,i++) {
        niik_image_set_voxel(outimg,j,niik_image_get_voxel(imglist[k],i));
      }
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = merge-slices */

  else if(!strncmp(argv[1],"featuremap",10)) {
    if(outname==NULL || (img==NULL && imglist==NULL)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<outname.nii> -val=<kernel_size>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -kernel=<kernel_size>         kernel_size is 1 for 3x3x3 square matrix [default=1]\n");
      exit(1);
    }
    if(ijk[0]>0) {
      ijk[0]=ijk[1]+ijk[2]*img->nx+ijk[3]*img->nx*img->ny;
      fprintf(stdout,"[%s] voxel %3i %3i %3i  %9i\n",fcname,(int)ijk[1],(int)ijk[2],(int)ijk[3],(int)ijk[0]);
      switch(argv[1][10]) {
      case '1':
        vec=niikvec_init(10);
        if(!niik_image_feature_voxel_type1(img,(int)ijk[0],kernel,vec)) {
          fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type1\n",fcname);
          exit(1);
        }
        niikvec_display(vec);
        exit(0);
      case '2':
        vec=niikvec_init(24);
        if(!niik_image_feature_voxel_type2(img,(int)ijk[0],vec)) {
          fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type2\n",fcname);
          exit(1);
        }
        niikvec_display(vec);
        exit(0);
      }
      exit(1);
    }
    switch(argv[1][10]) {
    case '1':
      NIIK_EXIT(((outimg=niik_image_feature_map_type1(img,(int)kernel))==NULL),fcname,"niik_image_feature_map_type1",1);
      break;
    case '2':
      NIIK_EXIT(((outimg=niik_image_feature_map_type2(img))==NULL),fcname,"niik_image_feature_map_type2",1);
      break;
    case '4':
      if(img!=NULL) {
        NIIK_EXIT(((outimg=niik_image_feature_map_type4(img))==NULL),fcname,"niik_image_feature_map_type4",1);
      } else if(imglist!=NULL) {
        NIIK_EXIT(((outimg=niik_image_feature_map_type4_multi(imglist,numimglist))==NULL),fcname,"niik_image_feature_map_type4_multi",1);
        for(n=0; n<numimglist; n++) imglist[n]=niik_image_free(imglist[n]);
      }
      break;
    default:
      fprintf(stderr,"[%s] ERROR: unknown option featuremap %c\n",fcname,argv[1][10]);
      exit(1);
    } /* switch */
    img=niik_image_free(img);
    if(datatype>=0) NIIK_EXIT((!niik_image_type_convert(outimg,datatype)),fcname,"niik_image_type_convert",1);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = featuremap1 */

  /* used this to test minc things
   * knakamura@mrs.mni.mcgill.ca
   *
   else if(!strncmp(argv[1],"test-qsform",11)){
   niik_mat44_to_cosines_start_step(img->sto_xyz,(double)img->dx,(double)img->dy,(double)img->dz,NULL,NULL,NULL,NULL,NULL);
   exit(0);
   }
  */

  else if(!strncmp(argv[1],"typeconvert",11)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s typeconvert -in=<img.nii> -out=<out.nii> [-val=<multiply_constant>]\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -converts the image type and writes <out.nii>\n");
      /* fprintf(stdout,"  -list of types = -uint8 -uint16 -uint32 -uint64 -int8 -int16 -int32 -int64 -float32 -float64\n");*/
      fprintf(stdout,"\n");
      exit(1);
    }
    mat44_display(img->qto_xyz);
    if(datatype<0) {
      fprintf(stderr,"[%s] ERROR: please select the new datatype\n",fcname);
      fprintf(stderr,"           -uint8  -int8      (unsigned) byte\n");
      fprintf(stderr,"           -uint16 -int16     (unsigned) short\n");
      fprintf(stderr,"           -uint32 -int32     (unsigned) int\n");
      fprintf(stderr,"           -uint64 -int64     (unsigned) long\n");
      fprintf(stderr,"           -float32 -float64  (double) floating point\n");
      exit(1);
    }
    fprintf(stdout,"[%s]   convert image:  %s to %s\n",fcname,nifti_datatype_string(img->datatype),nifti_datatype_string(datatype));
    fprintf(stdout,"[%s] qform:\n",fcname);
    mat44_display(img->qto_xyz);
    fprintf(stdout,"[%s] sform:\n",fcname);
    mat44_display(img->sto_xyz);
    if(niik_check_double_problem(inval)) inval=1.0;
    for(i=0; i<img->nvox; i++) niik_image_mul_voxel(img,i,inval);
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
      exit(1);
    }
    mat44_display(img->sto_xyz);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
  } /* op = typeconvert */

  else if(!strncmp(argv[1],"absdiffinfo",11) ||
          !strncmp(argv[1],"infoabsdiff",11)) {
    fprintf(stdout,"[%s] absdiffinfo\n",fcname);
    if(imglist==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii> [-mask=<mask.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"  -statistical info about the difference image = abs(in1 - in2)\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(numimglist!=2) {
      fprintf(stderr,"[%s] ERROR: require 2 images\n",fcname);
      exit(1);
    }
    if(outname!=NULL) {
      fprintf(stderr,"[%s] ERROR: there's no output\n",fcname);
      exit(1);
    }
    if(niik_image_cmp_dim(imglist[0],imglist[1])!=0) {
      fprintf(stderr,"[%s] ERROR: niik_image_cmp_dim\n",fcname);
      exit(1);
    }
    if(maskimg!=NULL) {
      if(niik_image_cmp_dim(imglist[0],maskimg)!=0) {
        fprintf(stderr,"[%s] ERROR: niik_image_cmp_dim\n",fcname);
        exit(1);
      }
    }
    for(i=0; i<imglist[0]->nvox; i++) {
      niik_image_set_voxel(imglist[0],i,
                           fabs(niik_image_get_voxel(imglist[0],i) - niik_image_get_voxel(imglist[1],i)));
    }
    nifti_image_free(imglist[1]);
    if(!niik_image_display_stats(imglist[0],maskimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_display_stats\n",fcname);
      exit(1);
    }
    nifti_image_free(imglist[0]);
    nifti_image_free(maskimg);
    imglist[0]=maskimg=NULL;
    free(imglist);
    imglist=NULL;
  } /* OP = absdiffinfo infoabsdiff */

  else if(!strncmp(argv[1],"histothresh",11)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> [-val=<frac>]\n",argv[0],argv[1]);
      fprintf(stdout,"  -calculates upper threshold\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(dval)) dval=0.01;
    fprintf(stdout,"upper threshold %15.9f\n",niik_image_get_upper_threshold(img,maskimg,dval));
    exit(0);
  } /* OP = histothresh */

  else if(!strncmp(argv[1],"applyaffine",11)) {
    if(img==NULL || outname==NULL || afmat==NULL) {
      fprintf(stdout,"  transform image using affine matrix\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<in.matrix> -out=<out.nii> [-ref=<ref.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -ref=<ref.nii>                : image dimension and pixel spacing of <ref.nii> is used\n");
      fprintf(stdout,"  -nn | -linear | -bspline      : choice of interpolation [default=linear]\n");
      fprintf(stdout,"  -invmatrix=<in.matrix>        : matrix is inverted before transformation\n");
      fprintf(stdout,"                                : this replaces -matrix=<in.matrix>\n");
      exit(1);
    }
    if(refimg==NULL) {
      if((refimg=niik_image_copy(img))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy(img)\n",fcname);
        exit(1);
      }
    }
    if(xfmflag>0) {
      fprintf(stdout,"[%s] xfm file:\n",fcname);
      niikmat_display(afmat);
      NIIK_EXIT((!niikmat_convert_from_xfm(afmat,img,refimg)),
                fcname,"niikmat_convert_from_xfm",1);
      fprintf(stdout,"[%s] xfm file converted:\n",fcname);
      niikmat_display(afmat);
    }
    fprintf(stdout,"[%s] affine transformation: %s\n",fcname,niik_interpolate_string(interp));
    niikmat_display(afmat);
    if(!niik_image_affine_transform_3d_update(img,refimg,afmat,interp)) {
      fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update(img,img,afmat,interp)\n",fcname);
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
        exit(1);
      }
    }
    fprintf(stdout,"[%s] writing image:    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write\n",fcname);
      exit(1);
    }
    if(img==refimg)
      img=niik_image_free(img);
    else {
      img=niik_image_free(img);
      refimg=niik_image_free(refimg);
    }
  } /* OP = applyaffine */

  else if(!strncmp(argv[1],"avgimgs-fov",11)) {
    if(imglist==NULL || outname==NULL || !(matlist!=NULL || xfmlist!=NULL) || refimg==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>,... -matrixlist=<mat1.matrix>,<mat2.matrix>... -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[%s] fov-adjusted image averaging\n",fcname);
    if(matlist==NULL) {
      nummatlist=numxfmlist;
    }
    if(nummatlist != numimglist) {
      fprintf(stderr,"[%s] ERROR: number of images and matrices did not match %i %i\n",fcname,numimglist,nummatlist);
      exit(1);
    }
    if(matlist==NULL) {
      matlist=xfmlist;
      for(n=0; n<nummatlist; n++) {
        NIIK_EXIT((!niikmat_convert_from_xfm(matlist[n],imglist[n],refimg)),
                  fcname,"niikmat_convert_from_xfm",1);
        /* fprintf(stdout,"[%s] matrix %i\n",fcname,n+1); niikmat_display(matlist[n]); */
      }
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(refimg,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    }
    if(!niik_image_average_with_fov(refimg,imglist,matlist,numimglist,interp)) {
      fprintf(stderr,"ERROR: niik_image_average_with_fov\n");
      exit(1);
    }
    /*fprintf(stdout,"image info %s\n%s\n\n",refimg->fname,nifti_image_to_ascii(refimg));*/
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(refimg,timestamp);
    if(!niik_image_write(outname,refimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    refimg=niik_image_free(refimg);
    exit(0);
  } /* OP = avgimgs-fov */

  else if(!strncmp(argv[1],"determinant",11)) {
    if(afmat!=NULL) {
      dval = niikpt_det4(afmat->m[0][0],afmat->m[0][1],afmat->m[0][2],afmat->m[0][3],
                         afmat->m[1][0],afmat->m[1][1],afmat->m[1][2],afmat->m[1][3],
                         afmat->m[2][0],afmat->m[2][1],afmat->m[2][2],afmat->m[2][3],
                         afmat->m[3][0],afmat->m[3][1],afmat->m[3][2],afmat->m[3][3]);
      fprintf(stdout,"determinant = %-19.9f\n",dval);
      afmat=niikmat_free(afmat);
      exit(0);
    } else if(matlist!=NULL) {
      for(n=0; n<nummatlist; n++) {
        afmat = matlist[n];
        dval = niikpt_det4(afmat->m[0][0],afmat->m[0][1],afmat->m[0][2],afmat->m[0][3],
                           afmat->m[1][0],afmat->m[1][1],afmat->m[1][2],afmat->m[1][3],
                           afmat->m[2][0],afmat->m[2][1],afmat->m[2][2],afmat->m[2][3],
                           afmat->m[3][0],afmat->m[3][1],afmat->m[3][2],afmat->m[3][3]);
        fprintf(stdout,"determinant = %-19.9f\n",dval);
        afmat=niikmat_free(afmat);
      }
      free(matlist);
    } else {
      fprintf(stdout,"  usage: %s %s -matrix=<bspline_img.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    exit(0);
  } /* OP = determinant */

  else if(!strncmp(argv[1],"bspline-inv",11)) {
    fprintf(stdout,"  estiamte from bspline coefficients\n");
    if(img==NULL || refimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<bspline_img.nii> -ref=<ref_img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(!niik_image_interpolate_inverse_3d_bspline_coeff(img,refimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_inverse_3d_bspline_coeff\n");
      exit(1);
    }
    /* fprintf(stdout,"   %3i %3i %3i   %3i %3i %3i\n",refimg->nx,refimg->ny,refimg->nz,refimg->nt,refimg->nu,refimg->nv); */
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(refimg,timestamp);
    if(!niik_image_write(outname,refimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    refimg=niik_image_free(refimg);
    img=niik_image_free(img);
    exit(0);
  } /* OP = bspline */

  else if(!strncmp(argv[1],"fill-lesion",11)) {
    if(imglist==NULL || numimglist!=2 || maskimg==NULL || outname==NULL ) {
      fprintf(stdout,"  usage: %s %s -imglist=<img.nii>,<wm_mask.nii> -mask=<lesion.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -matrixlist=<les-reg.matrix>,<wm-reg.matrix>\n");
      fprintf(stdout,"                   lesreg.matrix is the matrix to transform lesion mask to <img.nii>\n");
      fprintf(stdout,"                   and <wm-reg.matrix> is the matrix to transform WM mask to <img.nii>\n");
      exit(1);
    }
    if(!strncmp(argv[1],"fill-lesion-feature",19)) {
      fprintf(stderr,"ERROR: UNDER CONSTRUCTION\n");
      exit(1);
      if(!niik_image_fill_lesion_with_feature(imglist[0],maskimg,NULL,imglist[1],NULL)) {
        fprintf(stderr,"ERROR: niik_image_fill_lesion_with_feature\n");
        exit(1);
      }
      exit(0);
    }
    if(nummatlist==0) {
      if(!niik_image_fill_lesion(imglist[0],maskimg,imglist[1])) {
        fprintf(stderr,"ERROR: niik_image_fill_lesion\n");
        exit(1);
      }
    } else if(nummatlist==1) {
      if(!niik_image_fill_lesion_with_matrix(imglist[0],maskimg,matlist[0],imglist[1],NULL)) {
        fprintf(stderr,"ERROR: niik_image_fill_lesion_with_matrix\n");
        exit(1);
      }
    } else if(nummatlist==2) {
      if(!niik_image_fill_lesion_with_matrix(imglist[0],maskimg,matlist[0],imglist[1],matlist[1])) {
        fprintf(stderr,"ERROR: niik_image_fill_lesion_with_matrix\n");
        exit(1);
      }
    } else {
      fprintf(stderr,"ERROR: niik_image_fill_lesion_with_matrix\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(imglist[0],timestamp);
    if(!niik_image_write(outname,imglist[0])) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(1);
    }
    exit(0);
  } /* OP = fill-lesion */

  else if(!strncmp(argv[1],"combinewarp",11)) {
    if(img==NULL || refimg==NULL || outname==NULL || warpimg==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<movimg.nii> -ref=<refimg.nii> -out=<out.nii> -warp-disp-bspline=<warp.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] combine warp\n");
    if(!niik_image_combine_warp(img,refimg,warpimg,premat,postmat,warp_map_type)) {
      fprintf(stderr,"ERROR: niik_image_combine_warp\n");
      exit(1);
    }
    exit(0);
  } /* OP = combinewarp */

  else if(!strncmp(argv[1],"voxel2world",11)) {
    if(img==NULL || !(ijk[0]>0)) {
      fprintf(stdout,"  usage: [-options] -in=<img.nii> -ijk=<x>,<y>,<z>\n");
      exit(1);
    }
    pt = niikpt_val(ijk[1],ijk[2],ijk[3],0);
    pt = niikpt_affine_transform_m44(img->sto_xyz,pt);
    if(afmat!=NULL) {
      pt = niikpt_affine_transform(afmat,pt);
    }
    fprintf(stdout,"  world: %19.12f %19.12f %19.12f\n",pt.x,pt.y,pt.z);
    img=niik_image_free(img);
    exit(0);
  } /* OP = voxel2world */

  else if(!strncmp(argv[1],"world2voxel",11)) {
    if(img==NULL || !(xyz[0]>0)) {
      fprintf(stdout,"  usage: [-options] -in=<img.nii> -xyz=<x>,<y>,<z>\n");
      exit(1);
    }
    pt = niikpt_val(xyz[1],xyz[2],xyz[3],0);
    if(afmat!=NULL) {
      pt = niikpt_affine_transform(afmat,pt);
    }
    pt = niikpt_affine_transform_m44(img->sto_ijk,pt);
    fprintf(stdout,"  voxel: %19.12f %19.12f %19.12f\n",pt.x,pt.y,pt.z);
    img=niik_image_free(img);
    exit(0);
  } /* OP = voxel2world */

  else if(!strncmp(argv[1],"jacobianmap",11)) {
    fprintf(stdout,"[niikmath] Jacobian map calculation *** not validated!!!\n");
    if(warpimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -warp-disp-map=<warp_disp_map.nii.gz> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if((outimg = niik_image_jacobian_map(warpimg,maskimg,warp_map_type))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_jacobian_map\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    warpimg=niik_image_free(warpimg);
    outimg=niik_image_free(outimg);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = jacobianmap */

  else if(!strncmp(argv[1],"gray-dilate",6)) {
    fprintf(stdout,"[%s] grayscale dilation\n",fcname);
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: gray-dilate -in=<img.nii> -out=<out.nii> -radius=<R> -val=<V>\n");
      exit(1);
    }
    NIIK_EXIT((niik_check_double_problem(inval)),fcname,"missing -val=<V>",9);
    NIIK_EXIT((niik_check_double_problem(radius)),fcname,"missing -radius=<R>",9);
    fprintf(stdout,"niik_image_morph_gray_dilate\n");
    NIIK_EXIT((!niik_image_morph_gray_dilate(img,radius,inval)),fcname,"niik_image_morph_gray_dilate",9);
    fprintf(stdout,"[%s] writing output image:   %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = gray-dilate */

  else if(!strncmp(argv[1],"stats-label",11)) {
    if(img==NULL || maskimg==NULL) {
      fprintf(stdout,"  usage: %s stats-label -in=<in.nii> -mask=<label.nii> -val=<label>\n",argv[0]);
      fprintf(stdout,"\n");
      exit(1);
    }
    NIIK_EXIT((outname!=NULL),fcname,"There's no output",9);
    NIIK_EXIT((niik_check_double_problem(inval)),fcname,"ERROR: please use -val=<label>",9);
    NIIK_EXIT((!niik_image_label_display_stats(img,maskimg,inval)),fcname,"niik_image_label_display_stats",9);
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
  } /* stats-label */

  else if(!strncmp(argv[1],"writesform2",11)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"  -writes img's sform only, no reference\n");
      exit(1);
    }
    fprintf(stdout,"     sform\n");
    NIIK_EXIT(((afmat=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix",9);
    fprintf(stdout,"[%s] writing image:    %s\n",fcname,outname);
    NIIK_EXIT((!niikmat_write(outname,afmat)),fcname,"niikmat_write",9);
    afmat=niikmat_free(afmat);
  } /* OP = writesform2 */

  else if(!strncmp(argv[1],"writesform3",11)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"  -writes img's sform only with voxel size correction\n");
      exit(1);
    }
    fprintf(stdout,"     sform\n");
    NIIK_EXIT(((afmat=niikmat_mat44_matrix(img->sto_xyz))==NULL),fcname,"niikmat_mat44_matrix",9);
    NIIK_EXIT((!niikmat_multiply_mat1_free2
               (afmat,
                niikmat_scale_matrix(1.0/img->dx,1.0/img->dy,1.0/img->dz))),
              fcname,"niikmat_multiply_mat1_free",9);
    fprintf(stdout,"[%s] writing image:    %s\n",fcname,outname);
    NIIK_EXIT((!niikmat_write(outname,afmat)),fcname,"niikmat_write",9);
    afmat=niikmat_free(afmat);
  } /* OP = writesform3 */

  else if(!strncmp(argv[1],"unmerge-slice",13)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: unmerge-slice -in=<img.nii> -out=<outname>\n");
      exit(1);
    }
    fprintf(stdout,"[%s] unmerging slices\n",fcname);
    NIIK_EXIT((img->ndim!=3),fcname,"ndim is assumed to be 3",9);
    NIIK_EXIT(((outimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy",9);
    outimg->dim[0]=outimg->ndim=2;
    outimg->dim[3]=outimg->nz=1;
    outimg->nvox=outimg->nx*outimg->ny;
    for(k=j=0; k<img->nz; k++) {
      for(i=0; i<outimg->nvox; j++,i++) {
        niik_image_set_voxel(outimg,i,niik_image_get_voxel(img,j));
      }
      sprintf(fname,"%s-slice%i.nii.gz",outname,k);
      fprintf(stdout,"[%s] writing output  %s\n",fcname,fname);
      niik_image_append_history(outimg,timestamp);
      NIIK_EXIT((!niik_image_write(fname,outimg)),fcname,"niik_image_write",1);
    }
    exit(0);
  } /* OP = unmerge-slice */

  else if(!strncmp(argv[1],"brainhistofit",13)) {
    if(img==NULL || brainmask==NULL) {
      fprintf(stdout,"  usage: brainhistofit -img=<t1.nii> -brainmask=<mask.nii>\n");
      exit(1);
    }
    imin=niik_image_get_min(img,NULL);
    imax=0.8*niik_image_get_max(img,brainmask) + 0.2*niik_image_get_max(img,NULL);
    histomat=niikmat_init(3,5);
    if(histo_num<0) histo_num*=-1;
    fprintf(stdout,"[%s] min/max %9.5f %9.5f\n",fcname,imin,imax);
    if(!niik_image_fit_gaussian_tissues(img,brainmask,imin,imax,histo_num,0,histomat->m[0],histomat->m[1],histomat->m[2],&dval)) {
      fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian_tissues\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] GM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][0],histomat->m[0][0],histomat->m[1][0]);
    fprintf(stdout,"[%s] WM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][1],histomat->m[0][1],histomat->m[1][1]);
    fprintf(stdout,"[%s] CSF   %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][2],histomat->m[0][2],histomat->m[1][2]);
    exit(0);
  } /* brainhistofit */

  else if(!strncmp(argv[1],"absdiffimg",10)) {
    if(img==NULL || refimg==NULL || outname==NULL) {
      fprintf(stdout,"absolute difference image\n");
      fprintf(stdout,"  usage: niikmath absdiffimg -in=<in.nii> -ref=<ref.nii> -out=<out.nii>\n");
      exit(1);
    }
    NIIK_EXIT((img->nvox!=refimg->nvox),fcname,"different nvox",1);
    for(i=0; i<img->nvox; i++) {
      dval = fabs(niik_image_get_voxel(img,i) - niik_image_get_voxel(refimg,i));
      niik_image_set_voxel(img,i,dval);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",1);
    img=niik_image_free(img);
    refimg=niik_image_free(refimg);
  } /* OP = absdiffimg */

  else if(!strncmp(argv[1],"bimodalfit",10)) {
    if(img==NULL || maskimg==NULL) {
      fprintf(stdout,"bimodal fitting\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -mask=<mask.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -imin=<min>       : minimum value for histogram [default=auto]\n");
      fprintf(stdout,"  -imax=<max>       : maximum value for histogram [default=auto]\n");
      fprintf(stdout,"  -histo-num=<num>  : histogram bin number [default=200?]\n");
      exit(1);
    }
    if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,maskimg);
    if(niik_check_double_problem(imax)) {
      /*imax=(niik_image_get_max(img,maskimg) +
        niik_image_get_percentile(img,maskimg,0.99)) * 0.5;*/
      imax=niik_image_get_percentile(img,maskimg,0.999999);
    }
    if(histo_num < 0) histo_num=500;
    if(niik_check_double_problem(delta)) delta = (imax-imin)/(histo_num-1);
    if(!niik_image_bimodal_fit(img,maskimg,imin,delta,imax,histo_num,bimodal_fit+2,bimodal_fit+4,bimodal_fit,&dd)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_bimodal_fit\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] bimodal fit: distr1 %12.6f %12.6f %12.6f  [peak,mean,stdv]\n",
            bimodal_fit[0],bimodal_fit[2],bimodal_fit[4]);
    fprintf(stdout,"[niikmath] bimodal fit: distr2 %12.6f %12.6f %12.6f  [peak,mean,stdv]\n",
            bimodal_fit[1],bimodal_fit[3],bimodal_fit[5]);
    fprintf(stdout,"[niikmath] model error %12.6f\n",dd);
    omin=niik_image_get_min(img,maskimg)-5;
    omax=niik_image_get_max(img,maskimg)+5;
    for(dval=omin,dvol=0; dval<omax; dval+=0.1) {
      dvol += 0.1 * NIIK_DMIN(bimodal_fit[0]*NIIK_GaussPDF(dval-bimodal_fit[2],bimodal_fit[4]),
                              bimodal_fit[1]*NIIK_GaussPDF(dval-bimodal_fit[3],bimodal_fit[5]));
    }
    fprintf(stdout,"[niikmath] overlap %12.6f\n",dvol);
    if(debug) {
      fprintf(stdout,"*** TO PLOT ON MATLAB/OCTAVE: first create a text file of histogram ***\n");
      fprintf(stdout,"# create a histogram\n");
      fprintf(stdout,"  niikmath histogram -in=%s -mask=%s -out=tmp_histogram.txt -imin=%.3f -imax=%.3f -histo-num=%i\n",
              img->fname,maskimg->fname,imin,imax,histo_num);
      fprintf(stdout,"*** RUN THE FOLLOWING ON MATLAB/OCTAVE ***\n");
      fprintf(stdout,"  h=load('tmp_histogram.txt');\n");
      fprintf(stdout,"  x=h(1,:);\n");
      fprintf(stdout,"  y=h(2,:);\n");
      fprintf(stdout,"  y1=%.5f * exp(-(x-%.5f).^2/2/%.5f^2);\n",bimodal_fit[0],bimodal_fit[2],bimodal_fit[4]);
      fprintf(stdout,"  y2=%.5f * exp(-(x-%.5f).^2/2/%.5f^2);\n",bimodal_fit[1],bimodal_fit[3],bimodal_fit[5]);
      fprintf(stdout,"  plot(x,y,x,y1,x,y2,x,y1+y2);\n");
    }
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = bimodalfit */

  else if(!strncmp(argv[1],"applysform",10)) {
    if(img==NULL || outname==NULL || refimg==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -apply affine transformation usign sform in <img.nii>\n");
      fprintf(stdout,"  -img is the input transformed (moving) image\n");
      fprintf(stdout,"  -ref is the reference (target) image\n");
      fprintf(stdout,"   probably MNI152_T1_1mm.nii.gz or MNI152_T1_2mm.nii.gz\n");
      fprintf(stdout,"  -out.nii is the output image\n");
      fprintf(stdout,"  -instead of '-ref=<ref.nii>', '-ref=MNI152_1mm' or '-ref=MNI152_2mm'\n");
      fprintf(stdout,"   can be used\n");
      exit(1);
    }
    fprintf(stdout,"     sform\n");
    niikmat_display(niikmat_mat44_matrix(img->sto_xyz));
    afmat=niikmat_mat44_matrix(img->sto_xyz);
    niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(1.0/img->dx,1.0/img->dy,1.0/img->dz));
    niikmat_multiply_mat2_free1(niikmat_mat44_matrix(refimg->sto_ijk),afmat);
    niikmat_multiply_mat2_free1(niikmat_scale_matrix(refimg->dx,refimg->dy,refimg->dz),afmat);
    fprintf(stdout,"     using matrix\n");
    niikmat_display(afmat);
    fprintf(stdout,"     using %s\n",niik_interpolate_string(interp));
    if(!niik_image_affine_transform_3d_update(img,refimg,afmat,interp)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update\n");
      exit(1);
    }
    if(nifti_set_filenames(img,outname,0,img->byteorder)) {
      fprintf(stderr,"[niikmath] ERROR: nifti_set_filenames %s\n",outname);
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing image:    %s\n",img->fname);
    nifti_image_write(img);
    img=niik_image_free(img);
    refimg=niik_image_free(refimg);
  } /* OP = applysform */

  else if(!strncmp(argv[1],"vertexinfo",10)) {
    if(obj==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off>\n",argv[0],argv[1]);
      fprintf(stdout,"  -object's vertex info\n");
      fprintf(stdout,"  -val=<index>         : vertex index\n");
      fprintf(stdout,"  -xyz=<x>,<y>,<z>     : xyz coordinate in mm\n");
      fprintf(stdout,"  -ijk=<i>,<j>,<k>     : ijk coordinate in voxel space\n");
      fprintf(stdout,"  -in=<img.nii>        : image when using '-ijk'\n");
      exit(1);
    }
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    if(!niik_check_double_problem(inval)) {
      for(v=obj->vert; v!=NULL; v=v->next) {
        if(v->index-1==(int)inval) {
          off_display_vert_info(v,20);
        }
      }
      exit(0);
    } else if(xyz[0]>0) {
      pt.x = xyz[1];
      pt.y = xyz[2];
      pt.z = xyz[3];
    } else if(ijk[0]>0) {
      if(img!=NULL) {
        pt = niikpt_val(ijk[1] * img->dx,ijk[2] * img->dy,ijk[3] * img->dz,0);
      } else if(!niik_check_double_problem(delta)) {
        pt = niikpt_val(ijk[1]*delta,ijk[2]*delta,ijk[3]*delta,0);
      } else {
        fprintf(stderr,"[niikmath] ERROR: please use -xyz=<x>,<y>,<z>\n");
        fprintf(stderr,"           or -ijk=<i>,<j>,<k> with -delta or -in=<img.nii>\n");
        exit(1);
      }
    } else {
      fprintf(stderr,"[niikmath] ERROR: please use -xyz=<x>,<y>,<z>\n");
      fprintf(stderr,"           or -ijk=<i>,<j>,<k> with -delta or -in=<img.nii>\n");
      fprintf(stderr,"           or -val=<index>\n");
      exit(1);
    }
    dmin=1e9;
    for(v=v2=obj->vert; v!=NULL; v=v->next) {
      dval = niikpt_distance(v->v,pt);
      if(dmin>dval) {
        dmin = dval;
        v2=v;
      }
    }
    fprintf(stdout,"[niikmath]   distance %9.4f\n",dmin);
    off_display_vert_info(v2,20);
    exit(0);
  } /* vertex info */

  else if(!strncmp(argv[1],"writesform",10)) {
    if(img==NULL || outname==NULL || refimg==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"  -writes img's sform\n");
      fprintf(stdout,"  -'-ref=ref.nii' can be replaced with '-ref=MNI152_1mm' or '-ref=MNI152_2mm'\n");
      exit(1);
    }
    fprintf(stdout,"     sform\n");
    niikmat_display(niikmat_mat44_matrix(img->sto_xyz));
    afmat=niikmat_mat44_matrix(img->sto_xyz);
    niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(1.0/img->dx,1.0/img->dy,1.0/img->dz));
    niikmat_multiply_mat2_free1(niikmat_mat44_matrix(refimg->sto_ijk),afmat);
    niikmat_multiply_mat2_free1(niikmat_scale_matrix(refimg->dx,refimg->dy,refimg->dz),afmat);
    img=niik_image_free(img);
    refimg=niik_image_free(refimg);
    fprintf(stdout,"     using matrix\n");
    niikmat_display(afmat);
    fprintf(stdout,"[niikmath] writing image:    %s\n",outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1);
    }
    afmat=niikmat_free(afmat);
  } /* OP = writesform */

  else if(!strncmp(argv[1],"sphere-off",10)) {
    if(outname==NULL) {
      fprintf(stdout,"  usage: niikmath sphere-off -out=<out.nii>\n");
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -elen=<E>         target edge length [default=2mm]\n");
      fprintf(stdout,"  -radius=<R>       initial object's radius [default=image dimension]\n");
      fprintf(stdout,"  -phase=<phase>    angular phase shift for the spherical coordinate [default = 0]\n");
      fprintf(stdout,"  -xyz=<center.x>,<center.y>,<center.z>\n");
      fprintf(stdout,"                    center of sphere [default = 0,0,0]\n");
      exit(1);
    }
    if(niik_check_double_problem(phase)) {
      phase=0;
    }
    if(niik_check_double_problem(radius)) {
      radius = NIIK_DMAX(maskimg->nx*maskimg->dx,NIIK_DMAX(maskimg->ny*maskimg->dy,maskimg->nz*maskimg->dz));
    }
    if(obj==NULL) {
      fprintf(stdout,"[%s] creating a new spherical object\n",fcname);
      /*ctr.x = ctr.y = ctr.z = ctr.w = 0; */
      /*ctr=niikpt_image_get_centroid(maskimg,NULL); */
      ctr = niikpt_val(xyz[1],xyz[2],xyz[3],0);
      fprintf(stdout,"[%s]   centroid = %.3f %.3f %.3f\n",fcname,ctr.x,ctr.y,ctr.z);
      NIIK_EXIT(((obj=off_make_sphere_from_icosahedron(3,radius,ctr))==NULL),fcname,"off_make_sphere_from_icosahedron",1);
      /* phase shift */
      for(v=obj->vert; v!=NULL; v=v->next) {
        v->sph.psi+=NIIK_DEGREE2RAD(phase);
        if(v->sph.psi>NIIK_PI2) v->sph.psi-=NIIK_PI2;
        if(v->sph.psi<0) v->sph.psi+=NIIK_PI2;
      }
      off_kobj_remove_color(obj);
    }
    if(iter<0) iter=10;
    if(niik_check_double_problem(elen)) elen=2;
    NIIK_EXIT((!off_remesh_kobj(obj,elen,iter,0)),fcname,"off_remesh_kobj",1);
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->sph=niikxyz_sph(v->v);
    }
    obj->spherecoo=1;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,obj,0)),fcname,"off_kobj_write_off",1);
    obj=off_kobj_free(obj);
  } /* OP = sphere-off */
  /*moved to niikcortex_shrinkwrap.c */
#if 0
  else if(!strncmp(argv[1],"shrinkwrap",10)) {
    if(maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -mask=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -elen=<E>         target edge length [default=2mm]\n");
      fprintf(stdout,"  -obj=<obj.off>    initial object [default=auto]\n");
      fprintf(stdout,"  -radius=<R>       initial object's radius [default=image dimension]\n");
      fprintf(stdout,"  -iter=<I>         iteration [default=20]\n");
      fprintf(stdout,"  -val=<max_step_size>\n");
      fprintf(stdout,"                    maximum shrink size for each iteration [default=2mm]\n");
      fprintf(stdout,"  -phase=<phase>    change phase (angle) for spherical coordinate [default=0]\n");
      exit(1);
    }
    if(niik_check_double_problem(radius)) {
      radius = NIIK_DMAX(maskimg->nx*maskimg->dx,NIIK_DMAX(maskimg->ny*maskimg->dy,maskimg->nz*maskimg->dz));
    }
    if(niik_check_double_problem(inval)) {
      inval = 2;
    }
    if(niik_check_double_problem(phase)) {
      phase=0;
    }
    if(obj==NULL) {
      fprintf(stdout,"[%s] creating a new spherical object\n",fcname);
      ctr=niikpt_image_get_centroid(maskimg,NULL);
      fprintf(stdout,"[%s]   centroid = %.3f %.3f %.3f\n",fcname,ctr.x,ctr.y,ctr.z);
      NIIK_EXIT(((obj=off_make_sphere_from_icosahedron(3,radius,ctr))==NULL),fcname,"off_make_sphere_from_icosahedron",1);
      /* phase shift */
      for(v=obj->vert; v!=NULL; v=v->next) {
        v->sph.psi+=NIIK_DEGREE2RAD(phase);
        v->sph=niiksph_norm(v->sph);
      }
      off_kobj_remove_color(obj);
    }
    if(iter<0) iter=20;
    if(niik_check_double_problem(elen)) elen=2;
    /*
     * run shrink wrap
     */
    fprintf(stdout,"[%s] shrink wrap\n",fcname);
    NIIK_EXIT((!off_shrinkwrap_kobj_bbox_remesh(maskimg,obj,elen,iter,inval,flag_bbox,flag_remesh,debug)),fcname,"off_shrinkwrap_kobj_simple",1);
    /*
     * write output
     */
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,obj,0)),fcname,"off_kobj_write_off",1);
    obj=off_kobj_free(obj);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = shrinkwrap */
#endif
  else if(!strncmp(argv[1],"closeholes",10)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -closes holes by seed filling from the edge\n");
      fprintf(stdout,"  -out.nii is the output image\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]  mask count = %i\n",niik_image_count_mask(img));
    if(!niik_image_close_holes(img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_close_holes\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]  mask count = %i\n",niik_image_count_mask(img));
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"niik_image_write");
    img=niik_image_free(img);
  } /* OP = closeholes */

  else if(!strncmp(argv[1],"laplacemap",10)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<Small_img.nii>,<Large_img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -create laplacian map between <Small_img.nii> and <Large_img>\n");
      fprintf(stdout,"  -out.nii is the output image\n");
      exit(1);
    }
    if(numimglist!=2) {
      fprintf(stderr,"[%s] ERROR: 2 images are required in the '-imglist=' for laplacemap\n",fcname);
      exit(1);
    }
    for(n=0; n<2; n++)
      NIIK_RET0((!niik_image_type_convert(imglist[n],NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert");

    if(xyz[0]>0) {
      fprintf(stdout,"[%s] resampling to %.5f %.5f %.5f\n",fcname,xyz[1],xyz[2],xyz[3]);
      for(n=0; n<2; n++) {
        fprintf(stdout,"[%s] img %i processing\n",fcname,n+1);
        /*fprintf(stdout,"[%s] img %i   convert to byte image from [%s]\n",fcname,n+1,nifti_datatype_string(imglist[n]->datatype));*/
        /*NIIK_EXIT((!niik_image_type_convert(imglist[n],NIFTI_TYPE_UINT8)),fcname,"niik_image_convert",1);*/
        fprintf(stdout,"[%s] img %i   resampling \n",fcname,n+1);
        NIIK_EXIT((!niik_image_resample_3d_update(imglist[n],xyz[1],xyz[2],xyz[3],-1,-1,-1,NIIK_INTERP_NN)),fcname,"niik_image_resample_3d_update",1);
      }
      /*NIIK_EXIT((niik_image_write("niik_image_resample_3d_update_0.mnc",imglist[0])==0),fcname,"niik_image_write",9);
      NIIK_EXIT((niik_image_write("niik_image_resample_3d_update_1.mnc",imglist[1])==0),fcname,"niik_image_write",9);*/
    }
    fprintf(stdout,"[%s] start laplace map calculation\n",fcname);
    NIIK_RET0(((outimg=niik_image_laplace_map(imglist[0],imglist[1]))==NULL),fcname,"niik_image_laplace_map");
    for(n=0; n<2; n++) imglist[n]=niik_image_free(imglist[n]);
    free(imglist);
    imglist=NULL;
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((niik_image_write(outname,outimg)==0),fcname,"niik_image_write",9);
    outimg=niik_image_free(outimg);
  } /* OP = laplacemap */

  else if(!strncmp(argv[1],"similarity-img",14)) {
    if(img==NULL || refimg==NULL) {
      fprintf(stdout,"  usage: similarity-img -in=<img> -ref=<ref>\n");
      exit(1);
    }
    if(!niik_check_double_problem(FWHM)) {
      fprintf(stdout,"[%s:%i:%s] smoothing, %.3f\n",__FILE__,__LINE__,__func__,FWHM);
      niik_image_filter_gaussian_update(   img,(int)floor(FWHM*3.5),FWHM);
      niik_image_filter_gaussian_update(refimg,(int)floor(FWHM*3.5),FWHM);
    }
    if(xyz[0]>0) {
      fprintf(stdout,"[%s:%i:%s] resampling, %8.3f %8.3f %8.3f\n",__FILE__,__LINE__,__func__,xyz[1],xyz[2],xyz[3]);
      NIIK_EXIT((!niik_image_resample_3d_update(refimg,xyz[1],xyz[2],xyz[3],-1,-1,-1,interp)),fcname,"niik_image_resample_3d_update, refimg",1);
      NIIK_EXIT((!niik_image_affine_transform_3d_update(img,refimg,NULL,interp)),fcname,"niik_image_affine_transform_3d_update",1);
      //NIIK_EXIT((!niik_image_resample_3d_update(refimg,xyz[1],xyz[2],xyz[3],img->nx,img->ny,img->nz,interp)),fcname,"niik_image_resample_3d_update, ref",1);
      if(maskimg!=NULL) NIIK_EXIT((!niik_image_affine_transform_3d_update(maskimg,img,NULL,NIIK_INTERP_NN)),fcname,"niik_image_inverse_affine_transform_3d_update",1);
      if(maskref!=NULL) NIIK_EXIT((!niik_image_affine_transform_3d_update(maskref,img,NULL,NIIK_INTERP_NN)),fcname,"niik_image_inverse_affine_transform_3d_update",1);
      NIIK_EXIT((img->nx!=refimg->nx),fcname,"wrong nx",1);
      NIIK_EXIT((img->ny!=refimg->ny),fcname,"wrong nz",1);
      NIIK_EXIT((img->nz!=refimg->nz),fcname,"wrong ny",1);
      NIIK_EXIT((img->nvox!=refimg->nvox),fcname,"wrong nvox",1);
    }
    NIIK_EXIT((img->nvox!=refimg->nvox),fcname,"wrong nvox",1);
    NIIK_EXIT((!niik_image_aregister2_cost_cc(refimg,maskimg,img,&inval)),fcname,"niik_image_aregister2_cost_cc",1);
    fprintf(stdout,"[%s] CC   %-18.12f\n",fcname,inval);
    nmiobj=niik_aregister2_nmi_obj_init(niik_image_get_min(img,maskimg),
                                        niik_image_get_max(img,maskimg),
                                        niik_image_get_min(refimg,maskref),
                                        niik_image_get_max(refimg,maskref),
                                        niik_image_get_max(img,maskimg),
                                        niik_image_get_max(refimg,maskref));
    if(!niik_image_aregister2_cost_nmi(refimg,maskimg,img,nmiobj)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: niik_image_aregister2_cost_nmi\n",__FILE__,__LINE__,__func__);
      exit(1);
    }
    fprintf(stdout,"[%s] NMI  %-18.12f\n",fcname,nmiobj->nmi);

    // difference image
    NIIK_EXIT(((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL),fcname,"niik_image_copy_as_type",1);
    for(i=0; i<img->nvox; i++) niik_image_set_voxel(outimg,i,fabs(niik_image_get_voxel(img,i)-niik_image_get_voxel(refimg,i)));
    NIIK_EXIT((!niik_image_display_stats(outimg,maskimg)),fcname,"niik_image_display_stats",1);
    outimg=niik_image_free(outimg);

    // histogram comparison
    matlist=(niikmat **)calloc(2,sizeof(niikmat *));
    imax=floor(1+NIIK_DMAX(niik_image_get_max(img,maskimg),niik_image_get_max(refimg,maskref)));
    NIIK_EXIT(((matlist[0]=niik_image_histogram_auto(   img,maskimg,imax))==NULL),fcname,"niik_image_histogram_auto, img",1);
    NIIK_EXIT(((matlist[1]=niik_image_histogram_auto(refimg,maskref,imax))==NULL),fcname,"niik_image_histogram_auto, ref",1);
    for(i=0,dval=0; i<imax; i++) {
      dval+=fabs(matlist[0]->m[1][i]-matlist[1]->m[1][i]);
    }
    fprintf(stdout,"[%s] average absolute histogram diff   %16.0f\n",fcname,dval/imax);
    for(i=0,dval=0; i<imax; i++) {
      dval+=NIIK_SQ(matlist[0]->m[1][i]-matlist[1]->m[1][i]);
    }
    fprintf(stdout,"[%s] average squared histogram diff    %16.0f\n",fcname,dval/imax);
    exit(0);
  } /* OP similarity-img */

  else if(!strncmp(argv[1],"principalaxes",13)) {
    vec=niikvec_init(3);
    mat=niikmat_init(3,3);
    niik_image_principal_axes(img,maskimg,mat,vec);
    fprintf(stdout,"[%s]: eigen values:   %12.6f %12.6f %12.6f\n",fcname,vec->v[0],vec->v[1],vec->v[2]);
    for(i=0; i<3; i++) {
      fprintf(stdout,"[%s]: eigen vector %i: %12.6f %12.6f %12.6f\n",fcname,i,mat->m[i][0],mat->m[i][1],mat->m[i][2]);
    }
    exit(0);
  }

  else if(!strncmp(argv[1],"similarity",10)) {
    if(maskimg!=NULL && maskref!=NULL) {
      conmat_tp=conmat_fp=conmat_tn=conmat_fn=0;
      j=k=0;
      for(i=0; i<maskimg->nvox; i++) {
        dval = niik_image_get_voxel(maskimg,i);
        dtmp = niik_image_get_voxel(maskref,i);
        if(dval>0) j++;
        if(dtmp>0) k++;
        if(dval>0) {
          if(dtmp>0) conmat_tp+=1.0;
          else       conmat_fp+=1.0;
        } else {
          if(dtmp>0) conmat_fn+=1.0;
          else       conmat_tn+=1.0;
        }
      }
      fprintf(stdout,"  size mask                       %9.0f\n",(float)j);
      fprintf(stdout,"  size ref mask                   %9.0f\n",(float)k);
      fprintf(stdout,"  # true positive                 %9.0f\n",conmat_tp);
      fprintf(stdout,"  # true negative                 %9.0f\n",conmat_tn);
      fprintf(stdout,"  # false positive                %9.0f\n",conmat_fp);
      fprintf(stdout,"  # false negative                %9.0f\n",conmat_fn);
      fprintf(stdout,"  kappa                           %9.6f\n",2*conmat_tp / (j+k));
      fprintf(stdout,"  true positive rate              %9.6f\n",conmat_tp/(conmat_tp+conmat_fn));
      fprintf(stdout,"  specificity                     %9.6f\n",conmat_tn/(conmat_fp+conmat_tn));
      fprintf(stdout,"  positive predictive value       %9.6f\n",conmat_tp/(conmat_tp+conmat_fp));
      fprintf(stdout,"  negative predictive value       %9.6f\n",conmat_tn/(conmat_tn+conmat_fn));
      fprintf(stdout,"  accuracy                        %9.6f\n",(conmat_tp+conmat_tn)/maskimg->nvox);
      fprintf(stdout,"  false discovery rate            %9.6f\n",1.0-(conmat_tp/(conmat_tp+conmat_fp)));
      exit(0);
    } else {
      fprintf(stderr,"[%s] usage: niikmath similarity -mask=<mask1.nii> -refmask=<mask2.nii>\n",fcname);
      exit(1);
    }
  } /* OP = similarity */


  else if(!strncmp(argv[1],"copysqform",10)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s copysqform -in=<src.nii> -out=<dest.nii>\n\n",argv[0]);
      fprintf(stdout,"  -copies sqform from <src.nii> to <dest.nii>\n\n");
      exit(1);
    }
    fprintf(stdout,"[%s] reading outname   %s\n",fcname,outname);
    NIIK_EXIT(((outimg=niik_image_read(outname))==NULL),fcname,"niik_image_read",9);
    fprintf(stdout,"[%s] copying sform and qform from %s to %s\n",fcname,img->fname,outname);
    outimg->quatern_b  = img->quatern_b;
    outimg->quatern_c  = img->quatern_c;
    outimg->quatern_d  = img->quatern_d;
    outimg->qoffset_x  = img->qoffset_x;
    outimg->qoffset_y  = img->qoffset_y;
    outimg->qoffset_z  = img->qoffset_z;
    outimg->qfac       = img->qfac;
    outimg->qform_code = img->qform_code;
    outimg->qto_xyz = nifti_quatern_to_mat44(outimg->quatern_b,outimg->quatern_c,outimg->quatern_d,
                      outimg->qoffset_x,outimg->qoffset_y,outimg->qoffset_z,
                      outimg->dx,outimg->dy,outimg->dz,
                      outimg->qfac);
    fprintf(stdout,"  quatrn bcd = %f %f %f\n",outimg->quatern_b,outimg->quatern_c,outimg->quatern_d);
    fprintf(stdout,"  offset xyz = %f %f %f\n",outimg->qoffset_x,outimg->qoffset_y,outimg->qoffset_z);
    fprintf(stdout,"  pixdim xyz = %f %f %f\n",outimg->dx,outimg->dy,outimg->dz);
    fprintf(stdout,"  qfac       = %f\n",outimg->qfac);
    nifti_disp_matrix_orient("",outimg->qto_xyz);
    outimg->sform_code = img->sform_code;
    outimg->sto_xyz    = img->sto_xyz;
    outimg->sto_ijk    = img->sto_ijk;
    fprintf(stdout,"  sto_xyz = \n");
    mat44_display(outimg->sto_xyz);
    fprintf(stdout,"[%s] writing image:    %s\n",fcname,outimg->fname);
    sprintf(fname,"%s",outimg->fname);
    niik_image_append_history(outimg,timestamp);
    NIIK_RET0((!niik_image_write(fname,outimg)),fcname,"niik_image_write");
    img=niik_image_free(img);
    outimg=niik_image_free(outimg);
  } /* OP = copysqform */

  else if(!strncmp(argv[1],"percentile",10)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -percent=<percent> [-mask=<mask.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(inval)) {
      NIIK_EXIT((niik_check_double_problem(percent)),fcname,"please use -percent=<percent>",9);
    } else {
      percent=inval;
    }
    NIIK_EXIT((percent>1),fcname,"val must be [0,1]",9);
    fprintf(stdout,"%2i%% = %15.9f\n",(int)floor(percent*100),
            niik_image_get_percentile(img,maskimg,percent));
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = percentile */

  else if(!strncmp(argv[1],"qformradio",10)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -puts arbitrary qform (NIFTI_XFORM_ALIGNED_ANAT\n");
      fprintf(stdout,"   with center in the middle of image\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    img->quatern_b = 0;
    img->quatern_c = 1;
    img->quatern_d = 0;
    img->qfac      = -1;
    img->qform_code = NIFTI_XFORM_ALIGNED_ANAT;
    img->qoffset_x =  img->dim[1] * img->pixdim[1] * 0.5;
    img->qoffset_y = -img->dim[2] * img->pixdim[2] * 0.5;
    img->qoffset_z = -img->dim[3] * img->pixdim[3] * 0.5;
    /*
     * i'm not sure if offsets are set properly...
     * -hopefully it's initialized to zeros
     * --> changed to image center
     */
    img->qto_xyz = nifti_quatern_to_mat44(img->quatern_b,img->quatern_c,img->quatern_d,
                                          img->qoffset_x,img->qoffset_y,img->qoffset_z,
                                          img->dx,img->dy,img->dz,
                                          img->qfac);
    fprintf(stdout,"  quatrn bcd = %f %f %f\n",img->quatern_b,img->quatern_c,img->quatern_d);
    fprintf(stdout,"  offset xyz = %f %f %f\n",img->qoffset_x,img->qoffset_y,img->qoffset_z);
    fprintf(stdout,"  pixdim xyz = %f %f %f\n",img->dx,img->dy,img->dz);
    fprintf(stdout,"  qfac       = %f\n",img->qfac);
    nifti_disp_matrix_orient("",img->qto_xyz);
    fprintf(stdout,"  %f %f %f %f\n",img->qto_xyz.m[0][0],img->qto_xyz.m[0][1],img->qto_xyz.m[0][2],img->qto_xyz.m[0][3]);
    fprintf(stdout,"  %f %f %f %f\n",img->qto_xyz.m[1][0],img->qto_xyz.m[1][1],img->qto_xyz.m[1][2],img->qto_xyz.m[1][3]);
    fprintf(stdout,"  %f %f %f %f\n",img->qto_xyz.m[2][0],img->qto_xyz.m[2][1],img->qto_xyz.m[2][2],img->qto_xyz.m[2][3]);
    fprintf(stdout,"  %f %f %f %f\n",img->qto_xyz.m[3][0],img->qto_xyz.m[3][1],img->qto_xyz.m[3][2],img->qto_xyz.m[3][3]);
    if(!(datatype<0)) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert %s\n",nifti_datatype_string(datatype));
        exit(1);
      }
    }
    if(outname!=NULL) {
      niik_image_append_history(img,timestamp);
      NIIK_EXIT((niik_image_write(outname,img)==0),fcname,"niik_image_write",1);
    } else {
      if(nifti_set_filenames(img,outname,0,img->byteorder)) {
        fprintf(stderr,"[niikmath] ERROR: nifti_set_filenames %s\n",outname);
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing image:    %s\n",img->fname);
      nifti_image_write(img);
    }
    img=niik_image_free(img);
  } /* qformradio */

  else if(!strncmp(argv[1],"distancemap",10)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -out=<out_map.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -val=<max_size>         maximum allowable size [default=5mm]\n");
      exit(1);
    }
    if(niik_check_double_problem(inval)) inval = 5;
    if((outimg=niik_image_distance_map(img,inval))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_distance_map\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niik_image_write(outname,outimg)) {
      niik_image_append_history(outimg,timestamp);
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    outimg = niik_image_free(outimg);
    img = niik_image_free(img);
    exit(0);
  } /* OP = distancemap */

  else if(!strncmp(argv[1],"closebrain",10)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s closebrain -in=<in.nii> -out=<out.nii> -radius=<R>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -closes (dilates, closes holes, and erodes)\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(radius)) {
      fprintf(stderr,"[niikmath] ERROR: please set radius by -radius=<R>\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   vol   %8.3f pre-dilation\n",niik_image_get_mask_vol(img));
    /* dilation */
    fprintf(stdout,"[niikmath]   dilating R %8.3f\n",radius);
    if(!niik_image_morph_3d_radius_mask(img,NULL,NIIK_MORPH_DILATE,radius)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_morph3d_radius_mask\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   vol   %8.3f post-dilation\n",niik_image_get_mask_vol(img));
    /* close holes */
    fprintf(stdout,"[niikmath]   closing holes\n");
    if(!niik_image_close_holes(img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_close_holes\n");
      exit(1);
    }
    /* erosion */
    fprintf(stdout,"[niikmath]   vol   %8.3f pre-erosion\n",niik_image_get_mask_vol(img));
    fprintf(stdout,"[niikmath]   eroding R %8.3f\n",radius);
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_ERODE,radius)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_morph3d_radius_mask\n");
      exit(1);
    }
    /* writing output */
    fprintf(stdout,"[niikmath]   vol   %8.3f final\n",niik_image_get_mask_vol(img));
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=NULL;
  } /* OP = closebrain */

  else if(!strncmp(argv[1],"gaussnoise",10)) {
    fprintf(stdout,"[niikmath] gaussnoise (Gaussian noise)\n");
    if(img==NULL || outname==NULL || niik_check_double_problem(FWHM)) {
      fprintf(stdout,"  usage: %s %s -img=<ref.nii> -out=<out.nii> [-FWHM=<FWHM> -val=<Mean>]\n",argv[0],argv[1]);
      exit(1);
    }
    if(niik_check_double_problem(inval)) inval=0;
    if(!niik_image_gaussian_noise(img,FWHM,inval)) {
      fprintf(stderr,"ERROR: niik_image_gaussian_noise\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(1);
    }
    img=niik_image_free(img);
    exit(0);
  } /* gaussnoise */

  else if(!strncmp(argv[1],"maskthresh",10)) {
    fprintf(stdout,"[niikmath] maskthresh\n");
    if(img==NULL || outname==NULL || niik_check_double_problem(thresh)) {
      fprintf(stdout,"  usage: %s %s -img=<ref.nii> -out=<out.nii> -thresh=<T>\n",argv[0],argv[1]);
      exit(1);
    }
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_voxels_as_double_vector\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      dimg[i]=(dimg[i]<thresh)?0:dimg[i];
    }
    if(!niik_check_double_problem(uthresh)) {
      dimg[i]=(dimg[i]>thresh)?0:dimg[i];
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(1);
    }
    img=niik_image_free(img);
    exit(0);
  } /* maskthresh*/

  else if(!strncmp(argv[1],"grid-lines",10)) {
    fprintf(stdout,"  grid-lines   creates an image with lines with some spacing\n");
    if(refimg==NULL || outname==NULL || ijk[0]==0) {
      fprintf(stdout,"  usage: grid-lines -ref=<ref.nii> -out=<out.nii> -ijk=<spacing_in_voxels>\n");
      exit(1);
    }
    NIIK_EXIT((!niik_image_clear(refimg)),fcname,"niik_image_clear",9);
    for(k=n=0; k<refimg->nz; k++) {
      for(j=0; j<refimg->ny; j++) {
        for(i=0; i<refimg->nx; n++,i++) {
          if(ijk[1]>0) if(!(i%(int)ijk[1])) niik_image_set_voxel(refimg,n,1);
          if(ijk[2]>0) if(!(j%(int)ijk[2])) niik_image_set_voxel(refimg,n,1);
          if(ijk[3]>0) if(!(k%(int)ijk[3])) niik_image_set_voxel(refimg,n,1);
        }
      }
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(refimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,refimg)),fcname,"niik_image_write",9);
    refimg=niik_image_free(refimg);
    exit(0);
  } /* grid-lines */

  else if(!strncmp(argv[1],"heavisideF",10)) {
    fprintf(stdout,"[%s] apply heaviside function\n",fcname);
    if(img==NULL || outname==NULL || vlist==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -out=<out.nii> -vlist=<mean>,<stdv>\n",argv[0],argv[1]);
      exit(1);
    } /* options */
    NIIK_EXIT((numvlist!=2),fcname,"please use '-vlist=<mean>,<stdv>'",9);
    NIIK_EXIT((!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert",9);
    fprintf(stdout,"[%s] heaviside function %.6f %.6f\n",fcname,vlist[0],vlist[1]);
    for(i=0; i<img->nvox; i++) {
      niik_image_set_voxel(img,i,NIIK_Heaviside(niik_image_get_voxel(img,i)-vlist[0],vlist[1]));
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
    exit(0);
  } /* OP = heavisideF */

  else if(!strncmp(argv[1],"niikcortex",10)) {
    if(!strncmp(argv[1],"niikcortex-thick-color",18)) {
      fprintf(stdout,"[niikmath] niikcortex-thick-color\n");
      if(objlist==NULL || numobjlist!=2 || outname==NULL) {
        fprintf(stdout,"  usage: %s %s -objlist=<white.off>,<pial.off> -out=<out_white.off>,<out_pial.off>\n",argv[0],argv[1]);
        fprintf(stdout,"\n  optional usage:\n");
        fprintf(stdout,"  -omin=<min>             minimum thickness for color range [default=auto]\n");
        fprintf(stdout,"  -omax=<max>             maximum thickness for color range [default=auto]\n");
        exit(1);
      }
      if((CSlist=niik_csstring_to_list(outname,&numobjlist))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list\n");
        exit(1);
      }
      if(numobjlist!=2) {
        fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list: not 2\n");
        exit(1);
      }
      vlist = niik_calloc_double_vector(objlist[0]->nvert);
      fprintf(stdout,"  calc thk\n");
      fprintf(stdout,"    thickness = %12.9f\n",niikcortex_calc_area_weighted_thickness(objlist[0],objlist[1],vlist,2,3));
      if(niik_check_double_problem(omax)) {
        omax = niik_get_max_from_double_vector(vlist,objlist[0]->nvert);
      }
      if(niik_check_double_problem(omin)) {
        omin = niik_get_min_from_double_vector(vlist,objlist[0]->nvert);
      }
      fprintf(stdout,"  min/max %8.4f %8.4f\n",omin,omax);
      if(!niikcortex_add_thickness_color(objlist[0],objlist[1],vlist,omin,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niikcortex_add_thickness_color(ics,ocs,thk,omin,omax)\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] write   %s\n",CSlist[0]);
      if(!off_kobj_write_offply(CSlist[0],objlist[0],0)) {
        fprintf(stderr,"[niikmath] ERROR: writing %s\n",CSlist[0]);
        exit(1);
      }
      fprintf(stdout,"[niikmath] write   %s\n",CSlist[1]);
      if(!off_kobj_write_offply(CSlist[1],objlist[1],0)) {
        fprintf(stderr,"[niikmath] ERROR: writing %s\n",CSlist[1]);
        exit(1);
      }
      exit(0);
    } /* niikcortex-thick-color */

    if(!strncmp(argv[1],"niikcortex-color",16)) {
      double *meas=NULL;
      int meas_num;
      fprintf(stdout,"[niikmath] niikcortex-color\n");
      if(obj==NULL || inname==NULL || outname==NULL) {
        fprintf(stdout,"  usage: %s %s -obj=<input.off> -inname=<img.txt> -out=<out.off>\n",argv[0],argv[1]);
        fprintf(stdout,"\n  optional usage:\n");
        fprintf(stdout,"  -omin=<min>             minimum for color range [default=auto]\n");
        fprintf(stdout,"  -omax=<max>             maximum for color range [default=auto]\n");
        exit(1);
      }
      meas=niik_read_double_vector(inname,&meas_num);

      if(meas_num!=obj->nvert) {
        fprintf(stderr,"[%s] Inconsistent number of measurement: %i , expected %i\n",fcname,meas_num,obj->nvert);
        return 1;
      }

      if(niik_check_double_problem(omax)) {
        omax = niik_get_max_from_double_vector(meas,meas_num);
      }
      if(niik_check_double_problem(omin)) {
        omin = niik_get_min_from_double_vector(meas,meas_num);
      }
      fprintf(stdout,"  min/max %8.4f %8.4f\n",omin,omax);

      if(!niikcortex_add_color(obj, meas,omin,omax,NIIK_COLORMAP_SPECTRAL, 50)) {
        fprintf(stderr,"[niikmath] ERROR: niikcortex_add_color(ics,ocs,thk,omin,omax)\n");
        exit(1);
      }
      if(!off_kobj_write_offply(outname,obj,0)) {
        fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
        exit(1);
      }
      obj=off_kobj_free(obj);
      exit(0);
    } /* niikcortex-thick-color */

    else if(!strncmp(argv[1],"niikcortex-edge",15)) {
      fprintf(stdout,"  edge-based cortical atrophy method\n");
      if(img==NULL || maskimg==NULL) {
        fprintf(stdout,"\n  usage: %s %s -imglist=<base.nii>,<FU.nii> -img=<avg.nii> -mask=<edge_mask.nii> -matrixlist=<base.matrix>,<FU.matrix> [-maskimglist=<brain_mask.nii>,<ven_mask.nii>]\n",argv[0],argv[1]);
        fprintf(stdout,"  -FWHM=<FWHM>               : FWHM for blurring average image (calculation of normal direction)\n");
        fprintf(stdout,"  -niikcortex-edge-tissue=<CSF_mean>,<CSF_stdv>,<GM_mean>,<GM_stdv>,<WM_mean>,<WM_stdv>\n");
        fprintf(stdout,"                             : intensity and stdv for CSF, GM and WM\n");
        fprintf(stdout,"  -val=<range>               : range for colormap [default=1]\n");
        exit(1);
      }
      if(niik_check_double_problem(thresh)) thresh=5e3;
      if(niik_check_double_problem(inval)) inval=1;
      if(tstats==NULL) {
        tstats=niikmat_init(3,2);
        if(maskimglist==NULL) {
          fprintf(stderr,"[niikmath] ERROR: need to calculate tissue stats\n");
          fprintf(stderr,"           please use '-maskimglist=<brain_mask.nii>,<ven_mask.nii>' or\n");
          fprintf(stderr,"           '-niikcortex-edge-tissue=<CSF_mean>,<CSF_stdv>,<GM_mean>,<GM_stdv>,<WM_mean>,<WM_stdv>\n");
          exit(1);
        }
        if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,maskimglist[0]);
        if(niik_check_double_problem(imax)) {
          /*imax=niik_image_get_max(img,maskimglist[0]);*/
          imax=niik_image_get_percentile(img,maskimglist[0],0.999999);
        }
        if(histo_num < 0) histo_num=500;
        if(niik_check_double_problem(delta)) delta = (imax-imin)/(histo_num-1);
        fprintf(stdout,"[niikmath] bimodal %9.4f : %9.4f : %9.4f    %i\n",imin,delta,imax,histo_num);
        if(!niik_image_bimodal_fit(img,maskimglist[0],imin,delta,imax,histo_num,bimodal_fit+2,bimodal_fit+4,bimodal_fit,&dd)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_bimodal_fit\n");
          exit(1);
        }
        tstats->m[1][0]=bimodal_fit[2];
        tstats->m[1][1]=bimodal_fit[4];
        tstats->m[2][0]=bimodal_fit[3];
        tstats->m[2][1]=bimodal_fit[5];
        if(maskimglist[1]!=NULL || nummaskimglist!=2) {
          tstats->m[0][0]=niik_image_get_mean(img,maskimglist[1]);
          tstats->m[0][1]=sqrt(niik_image_get_var(img,maskimglist[1]));
          for(n=0; n<numimglist; n++) {
            maskimglist[n]=niik_image_free(maskimglist[n]);
          }
          free(maskimglist);
          nummaskimglist=0;
        } else {
          fprintf(stderr,"[niikmath] ERROR: please use '-maskimglist=<brain_mask>,<ven_mask.nii>'\n");
          exit(1);
        }
      } /* tissue stats */
      for(n=0; n<4; n++) idx[n]=(int)ijk[n];
      if(niik_check_double_problem(FWHM)) FWHM=2.5;
      if((tmpimglist[0]=niik_image_copy(img))==NULL) {
        fprintf(stderr,"ERROR: niik_image_copy\n");
        exit(1);
      }
      if(!niik_image_iscale(tmpimglist[0],0,niik_image_get_percentile(img,NULL,0.98),0,255)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_iscale\n");
        exit(1);
      }
      if(!niik_image_type_convert(tmpimglist[0],NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
      if(!niik_image_convert_to_color_image(tmpimglist[0])) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_convert_to_color_image\n");
        exit(1);
      }
      if(FWHM>0) {
        fprintf(stdout,"[niikmath] gaussian filter %0.2f\n",FWHM);
        if(!niik_image_filter_gaussian_update(img,17,FWHM)) {
          fprintf(stderr,"[niikmath] niik_image_filter_gaussian_update\n");
          exit(1);
        }
      } /* blur */

      atv=niikmat_init(3,img->nvox);
      if(!niik_image_niikcortex_atrophy_edge(img,imglist,matlist,maskimg,tstats,thresh,atv,idx,2)) {
        fprintf(stderr,"ERROR: niik_image_niikcortex_atrophy_edge\n");
        exit(1);
      }

      /*
       * write output image for the atrophy amount
       */
      if(outname!=NULL) {
        if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
          exit(1);
        }
        niik_image_clear(img);
        bimg=(unsigned char*)maskimg->data;
        for(i=k=0; i<img->nvox; i++) {
          if(bimg[i]==0) continue;
          niik_image_set_voxel(img,i,atv->m[0][k++]);
        }
        fprintf(stdout,"[niikmath] writing output    %s\n",outname);
        niik_image_append_history(img,timestamp);
        if(!niik_image_write(outname,img)) {
          fprintf(stderr,"ERROR: niik_image_write\n");
          exit(1);
        }
      } /* write otuput image */

      if(outcheckname!=NULL) {
        cmap = niik_colormap_get(NIIK_COLORMAP_ATROPHY,1000);
        /* niikmat_display(cmap); */
        for(i=k=0; i<img->nvox; i++) {
          if(niik_image_get_voxel(maskimg,i)>0) {
            j = NIIK_IMINMAX(floor((atv->m[0][k++] + inval) * 999.0 / (2.0*inval) + 0.5),0,999);
            niik_image_set_voxel(tmpimglist[0],i,floor(cmap->m[j][0]*255));
            niik_image_set_voxel(tmpimglist[0],i+img->nvox,floor(cmap->m[j][1]*255));
            niik_image_set_voxel(tmpimglist[0],i+img->nvox*2,floor(cmap->m[j][2]*255));
            /*if(i==140+137*img->nx+73*img->nx*img->ny){
              fprintf(stdout,"%i   %8.3f (=%8.3f + %8.3f) -> %i %5.3f %5.3f %5.3f\n",k-1,
                      atv->m[0][k-1],atv->m[1][k-1],atv->m[2][k-1],
                      j,cmap->m[j][0],cmap->m[j][1],cmap->m[j][2]);  }*/
          }
        }
        fprintf(stdout,"[niikmath] writing output    %s\n",outcheckname);
        niik_image_append_history(img,timestamp);
        if(!niik_image_write(outcheckname,tmpimglist[0])) {
          fprintf(stderr,"ERROR: niik_image_write\n");
          exit(1);
        }
        tmpimglist[0]=niik_image_free(tmpimglist[0]);
        n=niik_image_count_mask(maskimg);
        niik_write_double_vector("tmp_atvp.txt",atv->m[1],n);
        niik_write_double_vector("tmp_atvw.txt",atv->m[2],n);
        niik_write_double_vector("tmp_atv.txt",atv->m[0],n);
      } /* writing output check image */
      img=niik_image_free(img);
      maskimg=niik_image_free(maskimg);
      exit(0);
    } /* OP = niikcortex-edge */

    else if(!strncmp(argv[1],"niikcortex-nonctx",18)) {
      fprintf(stdout,"[niikmath] niikcortex-nonctx\n");
      if(maskimg==NULL || objlist==NULL || outname==NULL || numobjlist!=2) {
        fprintf(stdout,"  usage: %s %s -mask=<mask.nii> -objlist=<white_surf.off>,<pial_surf.off> -out=<ctxlabel.txt>\n",argv[0],argv[1]);
        exit(1);
      }
      if(niik_check_double_problem(inval)) inval = 1.0;
      if(niik_check_double_problem(delta)) delta = 0.1;
      if(iter<0) iter=2; /*TODO: potentially make sure niikcortex_ctx_label is NULL*/
      if((niikcortex_ctx_label = niikcortex_non_cortex_label(niikcortex_ctx_label,maskimg,objlist[0],objlist[1],inval,delta,iter))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niikcortex_non_cortex_label\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing output  %s\n",outname);
      if(!niik_write_int_vector(outname,niikcortex_ctx_label,objlist[0]->nvert)) {
        fprintf(stderr,"[niikmath] ERROR: niik_write_int_vector\n");
        exit(1);
      }
      if(outcheckname!=NULL) {
        off_kobj_add_color(objlist[0]);
        for(f=objlist[0]->face,n=0; f!=NULL; f=f->next,n++) {
          if(niikcortex_ctx_label[f->vert[0]->index-1]>0) {
            f->color[0] = f->color[1] = f->color[2] = f->color[3] = 1.0;
          } else if(niikcortex_ctx_label[f->vert[1]->index-1]>0) {
            f->color[0] = f->color[1] = f->color[2] = f->color[3] = 1.0;
          } else if(niikcortex_ctx_label[f->vert[2]->index-1]>0) {
            f->color[0] = f->color[1] = f->color[2] = f->color[3] = 1.0;
          } else {
            f->color[0] = 1.0;
            f->color[1] = f->color[2] = f->color[3] = 0.0;
          }
        }
        fprintf(stdout,"[niikmath] writing output  %s\n",outcheckname);
        if(!off_kobj_write_offply(outcheckname,objlist[0],0)) {
          fprintf(stderr,"[niikmath] ERROR: off_kobj_write_off\n");
          exit(1);
        }
      } /* outcheckname */
    } /* niikcortex-nonctx */

    else if(!strncmp(argv[1],"niikcortex-mimg",15)) {
      fprintf(stdout,"[niikmath] niikcortex-mimg\n");
      if(imglist==NULL || outname==NULL) {
        fprintf(stdout,"  usage: %s %s -in=<img.nii> -imglist=<brain_mask.nii>,<ventricle_mask.nii>,<gwi_mask.nii>,<cerebellum+brainstem.mask.nii>,<dgm_mask.nii> -out=<out.nii>\n",argv[0],argv[1]);
        exit(1);
      }
      if(!niikcortex_estimate_tissue_values(img,imglist[0],imglist[1],imglist[2],
                                            &niikcortex_tissue_val[2],&niikcortex_tissue_val[1],&niikcortex_tissue_val[0],
                                            &niikcortex_tissue_val[3],&niikcortex_tissue_val[4],&niikcortex_tissue_val[5],
                                            &niikcortex_tissue_val[6],&niikcortex_tissue_val[7], 0.5, 0.5)) {
        fprintf(stderr,"[niikmath] ERROR: niikcortex_estimate_tissue_values\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] niikcortex_estimate_tissue_values:\n");
      fprintf(stdout,"           GM WM CSF Brain   %7.3f %7.3f %7.3f %7.3f\n",niikcortex_tissue_val[0],niikcortex_tissue_val[1],niikcortex_tissue_val[2],niikcortex_tissue_val[3]);
      fprintf(stdout,"           WM_Surf           %7.3f (%7.3f)\n",niikcortex_tissue_val[4],niikcortex_tissue_val[6]);
      fprintf(stdout,"           Pial_Surf         %7.3f (%7.3f)\n",niikcortex_tissue_val[5],niikcortex_tissue_val[7]);
      if(niik_check_double_problem(radius )) radius =2.5;
      if(niik_check_double_problem(radius2)) radius2=1.0;
      outimg = niikcortex_modify_image(img,imglist[0],imglist[1],radius,radius,
                                       imglist[3],imglist[4],
                                       radius2,radius2,
                                       maskimg, niikcortex_tissue_val[1]);
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(outimg,timestamp);
      niik_image_write(outname,outimg);
      exit(1);
      exit(0);
    } /* modifying image */

    else if(!strncmp(argv[1],"niikcortex-gwi",14)) {
      fprintf(stdout,"[niikmath] niikcortex-gwi\n");
      if(imglist==NULL || outname==NULL) {
        fprintf(stdout,"  usage: %s %s -in=<img.nii> -imglist=<brain_mask.nii>,<ventricle_mask.nii>,<wm_mask.nii>,<dgm_mask.nii>,<avoid_mask.nii> -warp-loc-map=<warpimg.nii.gz> -out=<out.nii>\n",argv[0],argv[1]);
        fprintf(stdout,"\n");
        fprintf(stdout,"  optional usage:\n");
        fprintf(stdout,"  -radius=<R>         : median filter at the end [default=0.8]\n");
        fprintf(stdout,"  -radius2=<R2>       : raduis for morphologic opening to remove\n");
        fprintf(stdout,"                        bright vessels (etc) [default=1.2]\n");
        fprintf(stdout,"  -mask=<les.nii>     : optional lesion mask\n");
        fprintf(stdout,"  -val=<VDR>          : radius for ventricle dilation [default=1.5]\n");
        exit(1);
      }
      if(niik_check_double_problem(radius))  radius = 0.8;
      if(niik_check_double_problem(radius2)) radius2= 1.2;
      if(niik_check_double_problem(inval))   inval = 1.5;
      if((outimg = niikcortex_gwi_mask(img,
                                       imglist[0],         /* brain mask */
                                       imglist[1], inval,  /* ven mask */
                                       imglist[2],       /* wm mask */
                                       0,                /* dgm dilate */
                                       imglist[4],       /* avoid mask */
                                       maskimg,          /* lesion mask */
                                       warpimg,
                                       radius2,          /* blood vessel opening radius */
                                       radius)           /* median filter */
         ) == NULL ) {
        fprintf(stderr,"[niikmath] ERROR: niikcortex_gwi_mask\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(outimg,timestamp);
      if(!niik_image_write(outname,outimg)) {
        fprintf(stderr,"ERROR: niik_image_write\n");
        exit(1);
      }
      /* write check color-coded boundary image */
      if(outcheckname!=NULL) {
        fprintf(stdout,"    Creating boundary image\n");
        vlist = niik_calloc_double_vector(3);
        vlist[0] = vlist[2] = 0;
        vlist[1] = niik_image_get_percentile(img,NULL,0.98);
        if(!niik_image_boundary(img,outimg,vlist,0,1,0)) {
          fprintf(stderr,"ERROR: niik_image_boundary\n");
          exit(1);
        }
        fprintf(stdout,"[niikmath] writing check     %s\n",outcheckname);
        niik_image_append_history(img,timestamp);
        if(!niik_image_write(outcheckname,img)) {
          fprintf(stderr,"ERROR: niik_image_write\n");
          exit(1);
        }
      }
      exit(1);
    } /* OP = niikcortex-gwi */

    else if(!strncmp(argv[1],"niikcortex-wm",13)) {
      fprintf(stdout,"[niikmath] white matter segmentation for niikcortex\n");
      if(img==NULL || outname==NULL || maskimg==NULL || warpimg==NULL) {
        /* fprintf(stdout,"[niikmath] niikcortex-wm\n"); */
        fprintf(stdout,"  usage: %s %s -in=<img.nii> -mask=<brain_mask.nii> -out=<wm_mask.nii> -warp-loc-map=<warp.nii>\n",
                argv[0],argv[1]);
        fprintf(stdout,"\n");
        fprintf(stdout,"  optional usage:\n");
        fprintf(stdout,"  -percent=<P>        : percentile for threshold [default=0.5]\n");
        fprintf(stdout,"  -inval=<m>          : magnitude of change for regional threshold [default=0.08]\n");
        fprintf(stdout,"  -thresh=<T>         : white matter threshold instead of percent [default from '-percent']\n");
        fprintf(stdout,"  -FWHM=<fwhm>        : blurring kernel FWHM [default=OFF]\n");
        exit(1);
      }
      /* check optional usage */
      if(niik_check_double_problem(percent)) {
        percent=0.50;
      }
      if(niik_check_double_problem(inval)) {
        inval=0.08;
      }
      if(niik_check_double_problem(thresh)) {
        thresh=niik_image_get_percentile(img,maskimg,percent);
      }
      fprintf(stdout,"[niikmath] threshold   %.2f\n",thresh);
      /* blur image */
      if(!niik_check_double_problem(FWHM)) {
        fprintf(stdout,"[niikmath] gaussian filter: FWHM %.2f\n",FWHM);
        if(!niik_image_filter_gaussian_update(img,11,FWHM)) {
          fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
          exit(1);
        }
      }
      /* get the white matter mask */
      if((outimg = niikcortex_wm_mask(img,maskimg,warpimg,inval,thresh))==NULL) {
        fprintf(stderr,"ERROR: niikcortex_wm_mask\n");
        exit(1);
      }
      /* write output white matter mask */
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(outimg,timestamp);
      if(!niik_image_write(outname,outimg)) {
        fprintf(stderr,"ERROR: niik_image_write\n");
        exit(1);
      }
      exit(0);
    } /* OP = niikcortex-wm */

    else {
      fprintf(stderr,"[niikmath] ERROR: unknown operation, %s\n",argv[1]);
      fprintf(stderr,"           niikcortex functions\n");
      fprintf(stdout,"\n\n%s\n\n",*niikmath_op_niikcortex_list);
      exit(1);
    }
  } /* niikcortex functions */

  else if(!strncmp(argv[1],"avgmatrix",9)) {
    if(matlist==NULL&&xfmlist!=NULL) {
      matlist=xfmlist;
      nummatlist=numxfmlist;
    }
    if(matlist==NULL || outname==NULL) {
      fprintf(stdout,"[niikmath] matrix averaging\n");
      fprintf(stdout,"  usage 1: %s %s -matrixlist=<M1>,<M2>... -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if((afmat=niikmat_average_from_affine_param(matlist,nummatlist))==NULL) {
      fprintf(stderr,"ERROR: niikmat_average_from_affine_param\n");
      exit(1);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write %s\n",outname);
      exit(1);
    }
    exit(0);
  } /* OP = avgmatrix */

  else if(!strncmp(argv[1],"applywarp",9)) {
    if(img==NULL || outname==NULL || warpimg==NULL || refimg==NULL) {
      fprintf(stdout,"[%s] applywarp\n",fcname);
      fprintf(stdout,"  usage 1: niikmath applywarp -in=<img.nii> -warp-loc-map=<warp_loc_map.nii.gz> -ref=<ref.nii> -out=<out.nii>\n");
      fprintf(stdout,"  usage 2: niikmath applywarp -in=<img.nii> -warp-disp-map=<warp_disp_map.nii.gz> -ref=<ref.nii> -out=<out.nii>\n");
      fprintf(stdout,"  usage 3: niikmath applywarp -in=<img.nii> -warp-disp-bspline=<warp_disp_bspline_img.nii.gz> -ref=<ref.nii> -out=<out.nii>\n");
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -linear -bspline -nn     choice of interpolation [default=linear]\n");
      exit(1);
    }
    fprintf(stdout,"[%s] apply warp\n",fcname);
    if(!niik_image_apply_3d_warp_update(img,refimg,warpimg,warp_map_type,interp)) {
      fprintf(stderr,"ERROR: niik_image_apply_3d_warp_update\n");
      exit(1);
    }
    if(datatype>=0) {
      fprintf(stdout,"[niikmath] convert type %s\n",nifti_datatype_string(datatype));
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath] writing image:    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write\n");
      exit(1);
    }
    if(outcheckname!=NULL) { /* writing a banded image for checking */
      imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
      imglist[0]=img;
      imglist[1]=refimg;
      if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_combine\n");
        exit(1);
      }
      if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(outimg,datatype)\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing output    %s\n",outcheckname);
      niik_image_append_history(outimg,timestamp);
      if(!niik_image_write(outcheckname,outimg)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outcheckname);
        exit(1);
      }
      outimg=niik_image_free(outimg);
    } /* out-check image */
    img=niik_image_free(img);
    warpimg=niik_image_free(warpimg);
    refimg=niik_image_free(refimg);
    exit(0);
  } /* OP = applywarp */

  else if(!strncmp(argv[1],"make-offe",9)) {
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<in.off> -out=<out.off>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(1);
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* OP = make-offe */

  else if(!strncmp(argv[1],"fsasc2off",9)) {
    if(fsasc_name==NULL || outname==NULL) {
      fprintf(stdout,"  usage: fsasc2off -fsasc=<in.asc> -out=<out.off>\n");
      exit(1);
    }
    NIIK_EXIT((!niik_off_fsasc2off(fsasc_name,outname)),fcname,"niik_off_fsasc2off",1)
    exit(0);
  } /* OP = fsasc2off */

  else if(!strncmp(argv[1],"mergetdim",9) ||
          !strncmp(argv[1],"mergeudim",9) ||
          !strncmp(argv[1],"mergevdim",9)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s [mergetdim | mergeudim | mergevdim] -imglist=<in1.nii>,<in2.nii>... -out=<out.nii>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -merges images in the imglist in 4th (t) or 5th (u) dimension\n");
      fprintf(stdout,"  -writes <out.nii> as output\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -val=<max_value>                : scale approximately to this value [default=none]\n");
      fprintf(stdout,"  -vlist=<s1>,<s2>[,<s3>...]      : scaling factors [default=none]\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(vlist!=NULL) {
      NIIK_EXIT((numvlist!=numimglist),fcname,"#vlist and #imglist are different",9);
      for(n=0; n<numimglist; n++) {
        fprintf(stdout,"[%s] %9.5f    %s\n",fcname,vlist[n],imglist[n]->fname);
        NIIK_EXIT((!niik_image_mul_voxels_ROI(imglist[n],0,0,0,imglist[n]->nx,imglist[n]->ny,imglist[n]->nz,vlist[n])),
                  fcname,"niik_image_mul_voxels_ROI",9);
      }
    }
    if(niik_check_double_problem(inval)) {
      fprintf(stdout,"[%s] image combining (dim %s)  (no intensity normalization)\n",fcname,argv[1]+5);
      inval=0;
    } else {
      fprintf(stdout,"[%s] image combining (dim %s)    target scale %8.4f\n",fcname,argv[1]+5,inval);
    }
    if     (!strncmp(argv[1],"mergetdim",9)) n=4;
    else if(!strncmp(argv[1],"mergeudim",9)) n=5;
    else if(!strncmp(argv[1],"mergevdim",9)) n=6;
    else {
      fprintf(stderr,"[%s] ERROR: unknown dimension, %s\n",fcname,argv[1]);
      exit(1);
    }
    if((outimg=niik_image_combine(imglist,numimglist,n,inval))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_combine\n",fcname);
      exit(1);
    }
    if(inval>0) {
      outimg->cal_min = 0;
      outimg->cal_max = inval;
    }
    /* change the datatype */
    if(datatype<0) {
      datatype=imglist[0]->datatype;
    }
    if(!niik_image_type_convert(outimg,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(outimg,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    outimg=niik_image_free(outimg);
  } /* OP mergeudim */

  else if(!strncmp(argv[1],"invmatrix",9) ||
          !strncmp(argv[1],"invertmatrix",12) ||
          !strncmp(argv[1],"invert-matrix",13) ||
          !strncmp(argv[1],"inversematrix",13) ||
          !strncmp(argv[1],"inverse-matrix",14) ||
          !strncmp(argv[1],"matrixinverse",13) ||
          !strncmp(argv[1],"matrix-inverse",14) ||
          !strncmp(argv[1],"matrixinvert",12) ||
          !strncmp(argv[1],"matrix-invert",13)) {
    if(afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -matrix=<in.matrix> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -inverse <in.matrix> and write <out.matrix>\n");
      fprintf(stdout,"  -for more info on matrix, please '--help-matrix'\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    fprintf(stdout,"input");
    niikmat_display(afmat);
    niikmat_inverse_update(afmat);
    fprintf(stdout,"inverse");
    niikmat_display(afmat);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1);
    }
    afmat=niikmat_free(afmat);
  } /* OP = invmatrix */

  else if(!strncmp(argv[1],"sobelvoxel",9)) {
    if(img==NULL || ijk[0]==0.0) {
      fprintf(stderr,"  usage: %s %s -in=<in.nii> -ijk=<x>,<y>,<z>\n",argv[0],argv[1]);
      fprintf(stderr,"  -sobel filter for a voxel\n");
      fprintf(stderr,"  -input image <in.nii>\n");
      fprintf(stderr,"  -voxel location <x>,<y>,<z>\n");
      fprintf(stderr,"\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   sobel filter (%s)\n",argv[1]);
    n=floor(ijk[1])+floor(ijk[2])*img->nx+floor(ijk[3])*img->nx*img->ny;
    pt = niik_image_sobel_filter_voxel(img,n);
    if(niik_check_double_problem(pt.x)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_sobel_filter\n");
      exit(1);
    }
    fprintf(stdout,"  sobel %15.9f %15.9f %15.9f        %15.9f at %i %i %i    idx=%i\n",pt.x,pt.y,pt.z,pt.w,(int)ijk[1],(int)ijk[2],(int)ijk[3],n);
    img=niik_image_free(img);
  } /* OP = sobelvoxel */

  else if(!strncmp(argv[1],"cuberille",9)) {
    if(maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -mask=<mask.nii> -out=<out.off>\n",argv[0],argv[1]);
      exit(1);
    }
    if((obj=off_cuberille_kobj(maskimg,0))==NULL) {
      fprintf(stderr,"ERROR: off_cuberille_kobj\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_off\n");
      exit(1);
    }
    exit(0);
  } /* cuberille */

  else if(!strncmp(argv[1],"restruct=",9)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<mask.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_EXIT((!niik_image_restructure(img,argv[1]+9)),fcname,"niik_image_restructure",9);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* restructure */

  else if(!strncmp(argv[1],"gabormap",8)) {
    fprintf(stdout,"[%s] gabor filter\n",fcname);
    if(img==NULL || ijk[0]==0.0 || vlist==NULL || numvlist!=9) {
      fprintf(stdout,"[%s] usage: -in=<img.nii> -ijk=<x>,<y>,<z> -vlist=<lambda_x>,<lambda_y>,<lambda_z>,<sigma_x>,<sigma_y>,<sigma_z>,<rotate_x>,<rotate_y>,<rotate_z>\n",fcname);
      exit(1);
    }
    fprintf(stdout,"gabor = %24.20f\n",
            niik_image_voxel_gabor_filter(img,(int)ijk[1],(int)ijk[2],(int)ijk[3],
                                          vlist[0],vlist[1],vlist[2],
                                          vlist[3],vlist[4],vlist[5],
                                          vlist[6],vlist[7],vlist[8]));
    exit(0);
  } /* OP = gabormap */

  else if(!strncmp(argv[1],"cropmask",8)) {
    if(img==NULL || maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -mask=<mask.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -xyz=<add_x>,<add_y>,<add_z>      border size in mm\n");
      fprintf(stdout,"  -ijk=<add_x>,<add_y>,<add_z>      border size in voxels\n");
      exit(1);
    }
    if(!niik_image_get_mask_crop_roi(maskimg,img_dim_min+1,img_dim_max+1,img_dim_min+2,img_dim_max+2,img_dim_min+3,img_dim_max+3)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_mask_crop_roi(maskimg,img_dim_min+1,img_dim_max+1,img_dim_min+2,img_dim_max+2,img_dim_min+3,img_dim_max+3)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i\n",img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    if(xyz[0]>0) {
      img_dim_min[1]=img_dim_min[1]-floor(xyz[1]/maskimg->dx);
      img_dim_min[2]=img_dim_min[2]-floor(xyz[2]/maskimg->dy);
      img_dim_min[3]=img_dim_min[3]-floor(xyz[3]/maskimg->dz);
      img_dim_max[1]=img_dim_max[1]+floor(xyz[1]/maskimg->dx);
      img_dim_max[2]=img_dim_max[2]+floor(xyz[2]/maskimg->dy);
      img_dim_max[3]=img_dim_max[3]+floor(xyz[3]/maskimg->dz);
      fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i   (add edge)\n",
              img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    } else if(ijk[0]>0) {
      img_dim_min[1]=img_dim_min[1]-floor(ijk[1]);
      img_dim_min[2]=img_dim_min[2]-floor(ijk[2]);
      img_dim_min[3]=img_dim_min[3]-floor(ijk[3]);
      img_dim_max[1]=img_dim_max[1]+floor(ijk[1]);
      img_dim_max[2]=img_dim_max[2]+floor(ijk[2]);
      img_dim_max[3]=img_dim_max[3]+floor(ijk[3]);
      fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i   (add edge)\n",
              img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    }
    if(img_dim_min[1]<0) img_dim_min[1]=0;
    if(img_dim_min[2]<0) img_dim_min[2]=0;
    if(img_dim_min[3]<0) img_dim_min[3]=0;
    if(img_dim_max[1]>=img->nx) img_dim_max[1]=img->nx-1;
    if(img_dim_max[2]>=img->ny) img_dim_max[2]=img->ny-1;
    if(img_dim_max[3]>=img->nz) img_dim_max[3]=img->nz-1;
    fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i   (final)\n",img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    if(!niik_image_crop(img,img_dim_min[1],img_dim_max[1],img_dim_min[2],img_dim_max[2],img_dim_min[3],img_dim_max[3])) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_crop(img,img_dim_min[1],img_dim_max[1],img_dim_min[2],img_dim_max[2],img_dim_min[3],img_dim_max[3])\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    maskimg=niik_image_free(maskimg);
    img=niik_image_free(img);
    exit(0);
  } /* OP = cropmask */

  else if(!strncmp(argv[1],"edgeinfo",8)) {
    if(obj==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -val=<edge_index>\n",argv[0],argv[1]);
      exit(1);
    }
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    for(edge=obj->edge; edge!=NULL; edge=edge->next) {
      if(edge->index==(int)inval) {
        fprintf(stdout,"[niikmath] showing edge info for %i\n",edge->index);
        off_display_edge_info(edge,20);
      }
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* edgeinfo */

  else if(!strncmp(argv[1],"faceinfo",8)) {
    if(obj==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -val=<face_index>\n",argv[0],argv[1]);
      fprintf(stdout,"  -object's edge info\n");
      fprintf(stdout,"  -val=<index>         : vertex index\n");
      exit(1);
    }
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    for(face=obj->face; face!=NULL; face=face->next) {
      if(face->index==(int)inval) {
        fprintf(stdout,"[niikmath] showing face info for %i\n",face->index);
        off_display_face_info(face,20);
      }
    }
    obj=off_kobj_free(obj);
    exit(0);
  } /* faceinfo */

  else if(!strncmp(argv[1],"get3dimg",8)) {
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<img_num>\n",argv[0],argv[1]);
      exit(1);
    }
    size = img->nx*img->ny*img->nz;
    n = size * inval;
    for(i=0; i<size; n++,i++) {
      niik_image_set_voxel(img,i,niik_image_get_voxel(img,n));
    }
    img->ndim = img->dim[0] = 3;
    img->nvox = size;
    fprintf(stdout,"[niikmath] writing image:      %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    exit(0);
  } /* get3dimg */

  else if(!strncmp(argv[1],"diffinfo",8) ||
          !strncmp(argv[1],"infodiff",8)) {
    fprintf(stdout,"[%s] diffinfo\n",fcname);
    if(imglist==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii> [-mask=<mask.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"  -statistical info about the difference image = in1 - in2\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(outname!=NULL) {
      fprintf(stderr,"[niikmath] ERROR: there's no output\n");
      exit(1);
    }
    if(numimglist!=2) {
      fprintf(stderr,"[niikmath] ERROR: require 2 images\n");
      exit(1);
    }
    if(niik_image_cmp_dim(imglist[0],imglist[1])!=0) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_cmp_dim\n");
      for(n=0; n<2; n++)
        fprintf(stderr,"       img%i %5i %5i %5i\n",n+1,imglist[n]->nx,imglist[n]->ny,imglist[n]->nz);
      exit(1);
    }
    if(maskimg!=NULL) {
      if(niik_image_cmp_dim(imglist[0],maskimg)!=0) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_cmp_dim\n");
        exit(1);
      }
    }
    for(i=0; i<imglist[0]->nvox; i++) {
      niik_image_add_voxel(imglist[0],i,niik_image_get_voxel(imglist[1],i)*-1);
    }
    nifti_image_free(imglist[1]);
    if(!niik_image_display_stats(imglist[0],maskimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_display_stats\n");
      exit(1);
    }
    nifti_image_free(imglist[0]);
    nifti_image_free(maskimg);
    imglist[0]=maskimg=NULL;
    free(imglist);
    imglist=NULL;
  } /* OP = diffinfo */

  else if(!strncmp(argv[1],"calc_mtr",8)) {
    fprintf(stdout,"[niikmath] calc_mtr\n");
    if(imglist==NULL || outname==NULL || numimglist!=2) {
      fprintf(stdout,"  usage: %s %s -imglist=<MToff.nii>,<MTon.nii> -out=<MTR.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -???\n");
      exit(1);
    }
    if(imglist[0]->nvox != imglist[1]->nvox) {
      fprintf(stderr,"[niikmath] ERROR: nvox is different %i, %i\n",imglist[0]->nvox,imglist[1]->nvox);
      exit(1);
    }
    /*
     * calculate mtr
     */
    if(!niik_image_type_convert(imglist[0],NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(imglist[0],NIFTI_TYPE_FLOAT32)\n");
      exit(1);
    }
    fimg = (float*)imglist[0]->data;
    for(i=0; i<imglist[0]->nvox; i++) {
      dmtoff = fimg[i];
      dmton  = niik_image_get_voxel(imglist[1],i);
      if(fabs(dmtoff)<=1.0) {
        fimg[i]=0;
        continue;
      }
      fimg[i] = (dmtoff - dmton) / dmtoff * 100;
      if(fimg[i]<0) fimg[i]=0;
      if(fimg[i]>100) fimg[i]=100;
    }
    /* write output */
    fprintf(stdout,"[niikmath] writing mtr image:  %s\n",outname);
    niik_image_append_history(imglist[0],timestamp);
    if(!niik_image_write(outname,imglist[0])) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    exit(0);
  } /* OP = calc_mtr */

  else if(!strncmp(argv[1],"seedfill",8)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> [-mask=<mask.nii> -ijk=<i>,<j>,<k> -xyz=X,Y,Z]\n",argv[0],argv[1]);
      fprintf(stdout,"  -ijk for start position in voxel space\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]  mask count = %i\n",niik_image_count_mask(img));
    if(ijk[0]) {
      for(n=1; n<=3; n++) idx[n]=(int)ijk[n];
      fprintf(stdout,"[niikmath]  seed fill  [ijk]= %i %i %i grad=%i\n",idx[1],idx[2],idx[3],flag_grad);
      if(!niik_image_seed_fill_xyz(img, idx, flag_grad)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_seed_fill_xyz\n");
        exit(1);
      }
    } else if(xyz[0]>0) {
      pt = niikpt_val(xyz[1],xyz[2],xyz[3],0);
      if(afmat!=NULL) {
        pt = niikpt_affine_transform(afmat,pt);
      }
      pt = niikpt_affine_transform_m44(img->sto_ijk,pt);
      idx[1]=(int)floor(0.5+pt.x);
      idx[2]=(int)floor(0.5+pt.y);
      idx[3]=(int)floor(0.5+pt.z);
      fprintf(stdout,"[niikmath]  seed fill [XYZ]=%f %f %f -> [ijk]= %i %i %i grad=%i\n",xyz[1],xyz[2],xyz[3], idx[1],idx[2],idx[3],flag_grad);
      if(!niik_image_seed_fill_xyz(img, idx, flag_grad)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_seed_fill_xyz\n");
        exit(1);
      }
    } else if(maskimg!=NULL) {
      if(!niik_image_seed_fill(img,maskimg,flag_grad)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_seed_fill\n");
        exit(1);
      }
      img=niik_image_free(img);
      img=maskimg;
      maskimg=NULL;
    } else {
      if(!niik_image_seed_fill_from_middle(img,flag_grad)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_seed_fill\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath]  mask count = %i\n",niik_image_count_mask(img));
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"ERROR: writing file");
    img=niik_image_free(img);
  } /* OP = seedfill */

  else if(!strncmp(argv[1],"resample",8)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in1.nii> -out=<out.nii> [options]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -xyz=<x>,<y>,<z>     use voxel spacing <x> <y> and <z>\n");
      fprintf(stdout,"  -ijk=<i>,<j>,<k>     use image dimension <i> <j> and <k>\n");
      fprintf(stdout,"  -ref=<ref.nii>       use reference image <ref.nii>\n");
      exit(1);
    }
    /* using ijk as image dimensions */
    if(ijk[0]>0) {
      idx[1]=(int)ijk[1];
      idx[2]=(int)ijk[2];
      idx[3]=(int)ijk[3];
      /* no new pixel spacing -> use old */
      if(xyz[0]<=1e-8) {
        xyz[1]=img->dx;
        xyz[2]=img->dy;
        xyz[3]=img->dz;
      }
      /* fprintf(stdout,"[niikmath] resample\n");
      fprintf(stdout,"[niikmath] voxel size = %7.4f %7.4f %7.4f\n",xyz[1],xyz[2],xyz[3]);
       fprintf(stdout,"[niikmath] image dim  = %7i %7i %7i\n",idx[1],idx[2],idx[3]); */
      datatype=img->datatype;
      if(!niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],idx[1],idx[2],idx[3],interp)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],idx[1],idx[2],idx[3],interp)\n");
        exit(1);
      }
    } /* using image size */
    /* new pixel spacing, update the image size automatically */
    else if(xyz[0]>0) {
      /* fprintf(stdout,"[niikmath] resample\n");
      fprintf(stdout,"[niikmath] voxel size = %7.4f %7.4f %7.4f\n",xyz[1],xyz[2],xyz[3]); */
      datatype=img->datatype;
      if(!niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],-1,-1,-1,interp)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],idx[1],idx[2],idx[3],interp)\n");
        exit(1);
      }
    }
    /* use the reference (target) image */
    else if(refimg!=NULL) {
      fprintf(stdout,"[niikmath] using reference image\n");
      datatype=img->datatype;
      if(!niik_image_affine_transform_3d_update(img,refimg,NULL,interp)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update(img,refimg,NULL,interp)\n");
        exit(1);
      }
    } /* using reference image */
    /* unknown */
    else {
      fprintf(stderr,"[niikmath] ERROR: please use '-ref=<img.nii>' or '-ijk=<xdim>,<ydim>,<zdim>'\n");
      fprintf(stderr,"       '-xyz=<dx>,<dy>,<dz>' to specify the new image size and voxel size\n");
      exit(1);
    }
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    /* write output */
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    img=niik_image_free(img);
  } /* OP = resample */

  else if(!strncmp(argv[1],"permute=",8)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(!niik_image_rotate(img,argv[1]+8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_rotate(img,%s)\n",argv[1]+8);
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    img=niik_image_free(img);
  } /* OP = permute */

  else if(!strncmp(argv[1],"addsform",8)) {
    if(img==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<s.matrix> [-out=<out.nii>]\n",argv[0],argv[1]);
      exit(1);
    }
    if(!niik_image_update_sform(img,afmat)) {
      fprintf(stderr,"[%s] ERROR: niik_image_update_sform\n",fcname);
      exit(1);
    }
    if(outname!=NULL) {
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outname,img)) {
        fprintf(stderr,"ERROR: niik_image_write\n");
        exit(1);
      }
    } else {
      outname=(char *)calloc(4096,sizeof(char));
      sprintf(outname,"%s",img->fname);
      fprintf(stdout,"[niikmath] writing output    %s\n",img->fname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outname,img)) {
        fprintf(stderr,"ERROR: niik_image_write\n");
        exit(1);
      }
      free(outname);
      outname=NULL;
    }
    exit(0);
  } /* addsform */

  else if(!strncmp(argv[1],"tiffallz",8)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out_prefix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -imin=<min>         minimum intensity value [default=min]\n");
      fprintf(stdout,"  -imax=<max>         maximum intensity value [default=max]\n");
      fprintf(stdout,"  -percent=<percent>  percentile for max intensity value\n");
      exit(1);
    }
    if(!niik_check_double_problem(percent)) {
      imax=niik_image_get_percentile(img,NULL,percent);
      fprintf(stdout,"[%s] percentile [%.3f] -> %.4f\n",fcname,percent,imax);
    }
    if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,maskimg);
    if(niik_check_double_problem(imax)) imax=niik_image_get_max(img,maskimg);
    if(xyz[0]>0) { /* if resampling */
      fprintf(stdout,"[%s] resampling        %7.3f %7.3f %7.3f\n",fcname,xyz[1],xyz[2],xyz[3]);
      if(!niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],-1,-1,-1,interp)) {
        fprintf(stderr,"[%s] ERROR: niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],idx[1],idx[2],idx[3],interp)\n",fcname);
        exit(1);
      }
    } /* resampling */
    for(i=0; i<img->nz; i++) {
      sprintf(fname,"%s-%i.tif",outname,i);
      if(!niik_tiff_write_zslice(fname,img,imin,imax,i)) {
        fprintf(stderr,"[niikmath] ERROR: niik_tiff_write_xslice\n");
        exit(1);
      }
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = tiffall */

  else if(!strncmp(argv[1],"thinedge",8)) {
    fprintf(stdout,"[niikmath] Non-maximum Suppression)\n");
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -mask=<mask.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage\n");
      fprintf(stdout,"  -vlist=<mean>,<stdv>       : estimate mean and standard deviation for edge detection [default=off]\n");
      fprintf(stdout,"  -percent=<P>               : percentile for gradient threshold [default=0.5]\n");
      fprintf(stdout,"                             : this parameter determines <T>\n");
      fprintf(stdout,"  -thresh=<T>                : gradient threshold [default=auto using percent]\n");
      fprintf(stdout,"                             : if both <P> and <T> are set, <T> will be used and <P> is discarded.\n");
      exit(1);
    } /* options */
    if(niik_check_double_problem(FWHM))FWHM=3.0;
    if(numvlist==2) {
      dmean=vlist[0];
      dstdv=vlist[1];
    } else {
      dmean=dstdv=-1;
    }
    if(niik_check_double_problem(percent)) {
      percent = 0.5;
    }
    if((outimg=niik_image_non_maximum_suppression(img,maskimg,FWHM,dmean,dstdv,thresh,percent))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_non_maximum_suppression\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    if(outcheckname!=NULL) {
      fprintf(stdout,"[niikmath] preparing output check image %s\n",outcheckname);
      maskimg=niik_image_free(maskimg);
      if((maskimg=niik_image_mask_add_red_color_uint8(img,outimg))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mask_add_red_color_uint8\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing check file %s\n",outcheckname);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outcheckname,maskimg)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outcheckname);
        exit(1);
      }
    } /* outcheckname */
    exit(0);
  } /* OP = thinedge */

  else if(!strncmp(argv[1],"gaussPDF",8)) {
    fprintf(stdout,"[niikmath] calculate gaussian PDF (probability from Gaussian probability density function)\n");
    if(img==NULL || outname==NULL || vlist==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -out=<out.nii> -vlist=<mean>,<stdv>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage\n");
      exit(1);
    } /* options */
    if(numvlist!=2) {
      fprintf(stderr,"[niikmath] ERROR: please use '-vlist=<mean>,<stdv>'\n");
      exit(1);
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      niik_image_set_voxel(img,i,NIIK_GaussPDF(niik_image_get_voxel(img,i)-vlist[0],vlist[1]));
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = gaussPDF */


  /******************************************************
   *
   * voxel-wise comparison
   *
   ******************************************************/

  else if(!strncmp(argv[1],"voxel-ge",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img >= img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) >= niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-ge */

  else if(!strncmp(argv[1],"voxel-gt",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img > img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) > niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-gt */

  else if(!strncmp(argv[1],"voxel-le",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img <= img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) <= niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-le */

  else if(!strncmp(argv[1],"voxel-lt",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img < img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) < niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-lt */

  else if(!strncmp(argv[1],"voxel-eq",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img == img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) == niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-eq */

  else if(!strncmp(argv[1],"voxel-ne",8)) {
    if(img==NULL || img2==NULL) {
      fprintf(stdout,"  voxel-wise comparison: if img != img2, then 1, otherwise 0\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,(niik_image_get_voxel(img,n) != niik_image_get_voxel(img2,n)));
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(0);
  } /* OP = voxel-ne */


  /******************************************************
   * end of voxel-wise comparison
   ******************************************************/



  else if(!strncmp(argv[1],"pseudoT2",8)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stderr,"  usage: %s %s -imglist=<PD.nii>,<T2.nii> -out=<pseudoT2.nii> -TE=<PD_TE_in_msec>,<T2_TE_in_msec> -TR=<PD_TR_in_msec>[,<T2_TR_in_msec>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage\n");
      fprintf(stdout,"  -mask=<ROI.nii>          : region-of-interest mask\n");
      fprintf(stdout,"  -omax=<max_val>          : maximum possible value\n");
      /*fprintf(stdout,"  -val=<epsilon>           : allowable difference in T2\n");*/
      fprintf(stdout,"\n  N.B.\n");
      fprintf(stdout,"  TR must be the same for PD and T2. If two PD_TR and T2_TR are used,\n");
      fprintf(stdout,"  the program checks if they are same within <epsilon> msec.\n");
      exit(0);
    }
    if(num_MRI_TR==2) {
      if(fabs(MRI_TR[0] - MRI_TR[1]) > 5) {
        fprintf(stderr,"[%s] ERROR: TR is very different %12.4f %12.4f\n",fcname,MRI_TR[0],MRI_TR[1]);
        exit(0);
      }
    }
    if((outimg=niik_image_copy_as_type(imglist[0],NIFTI_TYPE_FLOAT32))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
      exit(0);
    }
    if(!niik_image_calculate_pseudoT2(imglist[0],MRI_TE[0],MRI_TR[0],imglist[1],MRI_TE[1],MRI_TR[0],omax,maskimg,outimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_calculate_pseudoT2\n",fcname);
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write\n");
      exit(1);
    }
    exit(0);
  } /* OP = pseudoT2 */



  else if(!strncmp(argv[1],"segcolor",8)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<img.nii>,<brain_mask.nii>,<GM_pve.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage\n");
      fprintf(stdout,"  -mask=<lesion.nii>           lesion mask\n");
      fprintf(stdout,"  -omax=<max>                  linear intensity scaling (new 98%%-value)\n");
      exit(1);
    }
    if(numimglist!=3) {
      fprintf(stderr,"ERROR: wrong # imglist, %i != 3\n",numimglist);
      exit(1);
    }
    if(vlist==NULL) {
      vlist=(double *)calloc(12,sizeof(double));
      vlist[0] = 0.5;
      vlist[1] = 1.0;
      vlist[2] = 2.0;
    } else if(numvlist!=3) {
      fprintf(stderr,"ERROR: wrong # vlist, %i != 3\n",numvlist);
      exit(1);
    }
    imax = niik_image_get_percentile(imglist[0],NULL,0.9);
    if((outimg=niik_image_copy(imglist[0]))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      exit(1);
    }
    outimg->cal_max=imax;
    outimg->cal_min=0;
    outimg->ndim=outimg->dim[0]=6;
    outimg->nt=outimg->dim[4]=outimg->nu=outimg->dim[5]=1;
    outimg->nv=outimg->dim[6]=3;
    outimg->dv=outimg->pixdim[6]=1;
    size=outimg->nx*outimg->ny*outimg->nz;
    outimg->nvox=outimg->nx*outimg->ny*outimg->nz*outimg->nt*outimg->nu*outimg->nv;
    free(outimg->data);
    NIIK_EXIT(((outimg->data=(void *)malloc(outimg->nvox*outimg->nbyper))==NULL),fcname,"malloc(outimg->data)",9);
    dval = niik_image_get_max(imglist[2],NULL);
    if(fabs(dval-255.0)<2.5) {
      fprintf(stdout,"[%s] modifying the gray matter probability range from 0-255 to 0-1\n",fcname);
      NIIK_EXIT((!niik_image_mul_voxels_ROI(imglist[2],0,0,0,imglist[2]->nx,imglist[2]->ny,imglist[2]->nz,1.0/255)),fcname,"niik_image_mul_voxels_ROI",9);
    } else if(fabs(dval-100.0)<0.5) {
      fprintf(stdout,"[%s] modifying the gray matter probability range from 0-100 to 0-1\n",fcname);
      NIIK_EXIT((!niik_image_mul_voxels_ROI(imglist[2],0,0,0,imglist[2]->nx,imglist[2]->ny,imglist[2]->nz,1.0/100.0)),fcname,"niik_image_mul_voxels_ROI",9);
    }
    dimg=(double *)calloc(outimg->nvox,sizeof(double));
    if(!strncmp(argv[1],"segcolor1",9)) { /* 2013-06-11, incomplete */
      for(i=0; i<size; i++) {
        dval=niik_image_get_voxel(imglist[0],i);
        dimg[i]=dimg[i+size]=dimg[i+size*2]=dval;
        if(niik_image_get_voxel(imglist[1],i)==0) {
          dimg[i]=dval*(1+vlist[0]);
        }
        dimg[i       ]+=dval*niik_image_get_voxel(imglist[2],i)*vlist[0];
        dimg[i+size  ]+=dval*niik_image_get_voxel(imglist[2],i)*vlist[1];
        dimg[i+size*2]+=dval*niik_image_get_voxel(imglist[2],i)*vlist[2];
        if(maskimg!=NULL) {
          if(niik_image_get_voxel(maskimg,i)>0)
            dimg[i+size  ]+=dval*vlist[2];
        }
      } /* each voxel */
    }  /* pink? */
    else { /* default = cyan */
      for(i=0; i<size; i++) {
        dval=niik_image_get_voxel(imglist[0],i);
        dimg[i]=dval;
        dimg[i+size]=dval;
        dimg[i+size*2]=dval;
        if(niik_image_get_voxel(imglist[1],i)==0) {
          dimg[i]=dval*(1+vlist[0]);
        }
        dimg[i+size  ]+=dval*niik_image_get_voxel(imglist[2],i)*vlist[1];
        dimg[i+size*2]+=dval*niik_image_get_voxel(imglist[2],i)*vlist[1];
        if(maskimg!=NULL) {
          if(niik_image_get_voxel(maskimg,i)>0)
            dimg[i+size  ]+=dval*vlist[2];
        }
      } /* each voxel */
    }  /* cyan */

    if(!niik_image_set_voxels_from_double_vector(outimg,dimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(0);
    }
    if(!niik_check_double_problem(omax)) {
      imax = niik_image_get_percentile(outimg,NULL,0.98);
      if(!niik_image_iscale(outimg,0,imax,0,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_iscale\n");
        exit(0);
      }
    } /* intensity scaling */
    if(datatype>=0) {
      if(!niik_image_type_convert(outimg,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    } /* type conversion */
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"ERROR: niik_image_write %s\n",outname);
      exit(0);
    }
    exit(0);
  } /* OP = segcolor */

  else if(!strncmp(argv[1],"labelvol",8)) {
    if(img==NULL) {
      fprintf(stdout,"[%s] usage: labelvol -in=<img>\n",fcname);
      exit(1);
    }
    NIIK_EXIT(((limg=niik_image_get_voxels_as_long_vector(img))==NULL),fcname,"niik_image_get_goxels_as_long_vector",9);
    imax=niik_image_get_max(img,NULL);
    dval = niik_image_get_voxel_size(img);
    for(n=0; n<(int)imax+2; n++) {
      for(i=j=0; i<img->nvox; i++) {
        if(limg[i]==(long)n) {
          j++;
        }
      }
      if(j>0) {
        fprintf(stdout,"[%s] label %-6i  %9i  %19.9f  %s\n",fcname,n,j,j*dval,img->fname);
      }
    }
    img=niik_image_free(img);
    free(limg);
    limg=NULL;
    exit(0);
  } /* OP  = labelvol */

  else if(!strncmp(argv[1],"add-rice",8)) {
    fprintf(stdout,"[%s] add-rice: adding Rician noise\n",fcname);
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"[%s] usage: -in=<in.nii> -out=<out.nii> -val=<s>\n",fcname);
      exit(1);
    }
    NIIK_EXIT((niik_check_double_problem(inval)),fcname,"please use -val=<s> option",9);
    NIIK_EXIT(((dimg=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector",9);
    if(!strncmp(argv[1],"add-rice-random",15)) {
      fprintf(stdout,"[%s] real pseudo-random generator     %.6f\n",fcname,inval);
      for(i=0; i<img->nvox; i++) {
        dimg[i]=NIIK_RiceRnd2(dimg[i],inval);
      }
    } else if(!strncmp(argv[1],"add-rice",8)) {
      fprintf(stdout,"[%s] reproducible pseudo-random generator    %.6f\n",fcname,inval);
      for(i=0; i<img->nvox; i++) {
        dimg[i]=NIIK_RiceRnd(dimg[i],inval);
      }
    }
    NIIK_EXIT((!niik_image_set_voxels_from_double_vector(img,dimg)),fcname,"niik_image_set_voxels_from_double_vector",9);
    free(dimg);
    if(datatype>=0) {
      NIIK_EXIT((!niik_image_type_convert(outimg,datatype)),fcname,"niik_image_type_convert",9);
    }
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* add-rice */

  else if(!strncmp(argv[1],"cropimg",7)) {
    if(img==NULL || vlist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: niikmath cropimg -in=<in.nii> -out=<out.nii> -vlist=<xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>\n");
      exit(0);
    }
    NIIK_EXIT((numvlist!=6),fcname,"please use \'-vlist=<xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>\'",9);
    img_dim_min[1]=vlist[0];
    img_dim_min[2]=vlist[1];
    img_dim_min[3]=vlist[2];
    img_dim_max[1]=vlist[3];
    img_dim_max[2]=vlist[4];
    img_dim_max[3]=vlist[5];
    fprintf(stdout,"[%s]   roi %3i %3i %3i  -  %3i %3i %3i   (add edge)\n",fcname,
            img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    if(img_dim_min[1]<0) img_dim_min[1]=0;
    if(img_dim_min[2]<0) img_dim_min[2]=0;
    if(img_dim_min[3]<0) img_dim_min[3]=0;
    if(img_dim_max[1]>=img->nx) img_dim_max[1]=img->nx-1;
    if(img_dim_max[2]>=img->ny) img_dim_max[2]=img->ny-1;
    if(img_dim_max[3]>=img->nz) img_dim_max[3]=img->nz-1;
    fprintf(stdout,"[%s]   roi %3i %3i %3i  -  %3i %3i %3i   (final)\n",fcname,
            img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    NIIK_EXIT((!niik_image_crop(img,img_dim_min[1],img_dim_max[1],img_dim_min[2],img_dim_max[2],img_dim_min[3],img_dim_max[3])),
              fcname,"niik_image_crop",9);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
    exit(0);
  } /* OP = cropimg */

  else if(!strncmp(argv[1],"unmerge",7)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s unmerge -in=<img.nii> -out=<out>\n",fcname);
      fprintf(stdout,"\n  [%s] optional usage\n",fcname);
      fprintf(stdout," -vlist=<v1>[,<v2>...]        : list of indices for writing output\n");
      exit(1);
    }
    fprintf(stdout,"[%s] unmerging volumes\n",fcname);
    if((imglist=niik_image_unmerge3d(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_unmerge3d\n");
      exit(0);
    }
    numimglist=img->nvox/img->nx/img->ny/img->nz;
    if(vlist==NULL) {
      vlist=(double *)calloc(numimglist,sizeof(double));
      for(n=0; n<numimglist; n++) {
        vlist[n]=n;
      }
      numvlist=numimglist;
    }
    for(m=0; m<numvlist; m++) {
      n=(int)vlist[m];
      sprintf(fname,"%s%i.nii.gz",outname,n);
      fprintf(stdout,"[niikmath] writing output  %s\n",fname);
      niik_image_append_history(imglist[n],timestamp);
      NIIK_EXIT((!niik_image_write(fname,imglist[n])),fcname,"niik_image_write",9);
    }
    exit(0);
  } /* OP = unmerge */

  else if(!strncmp(argv[1],"invwarp",7)) {
    if(warpimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -warp-disp-bspline=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  -output is a displacement map\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath]   invert warp\n");
    if(!niik_image_type_convert(warpimg,NIFTI_TYPE_FLOAT64)) {
      fprintf(stderr,"ERROR: niik_image_type_convert\n");
      return 0;
    }
    if(!niik_image_nregister_invert_nonlinear_map(warpimg,warp_map_type)) {
      fprintf(stderr,"ERROR: niik_image_nregister_invert_nonlinear_map\n");
      return 0;
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niik_image_write(outname,warpimg)) {
      fprintf(stderr,"ERROR: niik_image_write %s\n",outname);
      exit(0);
    }
    if(img!=NULL) {
      refimg=niik_image_copy(img);
    }
    exit(0);
  } /* OP = invwarp */

  else if(!strncmp(argv[1],"combine-obj",11)) {
    if(objlist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: combine-obj -objlist=<obj1.off>,<obj2.off>[...] -out=<out.off>\n");
      exit(1);
    }
    obj=objlist[0];
    NIIK_EXIT((obj==NULL),fcname,"first obj is null",1);

    fprintf(stdout,"[%s] #vert = %i\n",fcname,objlist[0]->nvert);
    /* Fast forward to the end of linked list */
    for(vert=objlist[0]->vert; vert->next!=NULL; vert=vert->next) {}
    fprintf(stdout,"[%s] #edge = %i\n",fcname,objlist[0]->nedge);
    for(edge=objlist[0]->edge; edge->next!=NULL; edge=edge->next) {}
    fprintf(stdout,"[%s] #face = %i\n",fcname,objlist[0]->nface);
    for(face=objlist[0]->face; face->next!=NULL; face=face->next) {}

    /* glue additional objects in the end of linked-lists */
    for(n=1; n<numobjlist; n++) {
      int d_vertex= vert->index;
      int d_edge  = edge->index;
      int d_face  = face->index;

      fprintf(stdout,"[%s] object %i\n",fcname,n);
      fprintf(stdout,"[%s] #vert = %i\n",fcname,objlist[n]->nvert);
      fprintf(stdout,"[%s] #edge = %i\n",fcname,objlist[n]->nedge);
      fprintf(stdout,"[%s] #face = %i\n",fcname,objlist[n]->nface);

      vert->next=objlist[n]->vert;
      vert->next->prev=vert;
      edge->next=objlist[n]->edge;
      edge->next->prev=edge;
      face->next=objlist[n]->face;
      face->next->prev=face;
      /*TODO: updated indexes ?*/
      /* Fast forward to the end of linked list */
      for(vert=vert->next; vert->next!=NULL; vert=vert->next) {vert->index+=d_vertex;}
      for(edge=edge->next; edge->next!=NULL; edge=edge->next) {edge->index+=d_edge;}
      for(face=face->next; face->next!=NULL; face=face->next) {face->index+=d_face;}

      /*update the last element*/
      vert->index+=d_vertex;
      edge->index+=d_edge;
      face->index+=d_face;
    }
    /*update counts*/
    obj->nvert=vert->index;
    obj->nface=face->index;
    obj->nedge=edge->index;
    
    fprintf(stdout,"[%s] output object \n",fcname);
    fprintf(stdout,"[%s] #vert = %i\n",fcname,obj->nvert);
    fprintf(stdout,"[%s] #edge = %i\n",fcname,obj->nedge);
    fprintf(stdout,"[%s] #face = %i\n",fcname,obj->nface);
   
    fprintf(stdout,"[%s] writing output %s\n",fcname,outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,objlist[0],0)),fcname,"off_kobj_write_off",1);
    exit(0);
  } /* OP = combine-obj */

  else if(!strncmp(argv[1],"avgobjs",7)) {
    if(objlist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: avgobjs -objlist=<obj1.off>,<obj2.off>[,...] -out=<out.off>\n");
      exit(1);
    }
    if(vlist==NULL) {
      vlist=(double *)calloc(numobjlist,sizeof(double));
      for(n=0; n<numobjlist; n++) vlist[n]=1.0/numobjlist;
    }
    for(n=0; n<numobjlist; n++)  {
      fprintf(stdout,"[%s] averaging %3i   %9.8f\n",fcname,n,vlist[n]);
    }
    vertlist=(kvert **)calloc(numobjlist,sizeof(kvert *));
    for(n=0; n<numobjlist; n++) {
      vertlist[n]=objlist[n]->vert;
    }
    for(vertlist[0]=objlist[0]->vert; vertlist[0]!=NULL; vertlist[0]=vertlist[0]->next) {
      vertlist[0]->v = niikpt_kmul(vertlist[0]->v,vlist[0]);
      for(n=1; n<numobjlist; n++) {
        vertlist[0]->v = niikpt_move_normal(vertlist[0]->v,vertlist[n]->v,vlist[n]);
      }
      for(n=1; n<numobjlist; n++) vertlist[n]=vertlist[n]->next;
    } /* each vertex */
    fprintf(stdout,"[%s] writing output %s\n",fcname,outname);
    NIIK_EXIT((!off_kobj_write_offply(outname,objlist[0],0)),fcname,"off_kobj_write_off",1);
    exit(0);
  } /* OP = avgobjs */

  else if(!strncmp(argv[1],"obj2img",7) ||
          !strncmp(argv[1],"obj2img-qform",13)) {
    if(obj==NULL || img==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -obj=<obj.off> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    fprintf(stdout,"[niikmath] obj2img start\n");
    if(niik_check_double_problem(inval))
      inval=niik_image_get_max(img,NULL)+100;
    if(!strncmp(argv[1],"obj2img-qform",13)) {
      if(!off_obj2img_use_qform(img,obj,inval)) {
        fprintf(stderr,"[niikmath] ERROR: off_obj2img_use_qform(img,obj,inval)\n");
        exit(0);
      }
    } else if(!strncmp(argv[1],"obj2img",7)) {
      if(!off_obj2img_color(img,obj,inval)) {
        fprintf(stderr,"[niikmath] ERROR: off_obj2img(img,obj,inval)\n");
        exit(0);
      }
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    img=niik_image_free(img);
    obj=off_kobj_free(obj);
    exit(0);
  } /* obj2img-qform */

  else if(!strncmp(argv[1],"avgimgs",7)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>,... -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_EXIT((img!=NULL),fcname,"img is not used, use imglist\n",9);
    if(datatype>=0) {
      fprintf(stdout,"[%s] converting type %s\n",fcname,nifti_datatype_string(datatype));
      NIIK_EXIT((!niik_image_type_convert(imglist[0],datatype)),fcname,"niik_image_typeconvert",9);
    }
    NIIK_EXIT(((outimg=niik_image_average_multiple(imglist,numimglist,vlist))==NULL),fcname,"niik_image_average_multiple",9);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_image_write",9);
    outimg=niik_image_free(outimg);
  } /* OP = avgimgs */

  else if(!strncmp(argv[1],"maximgs",7)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>,... -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(img!=NULL) {
      fprintf(stderr,"[niikmath] ERROR: img is not used, use imglist\n");
      exit(1);
    }
    if((outimg=niik_image_maximize_multiple(imglist,numimglist,vlist))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: outimg=niik_image_average_multiple(imglist,imglist)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    outimg=niik_image_free(outimg);
  } /* OP = maximgs */

  else if(!strncmp(argv[1],"addimgs",7)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>,... -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(img!=NULL) {
      fprintf(stderr,"[niikmath] ERROR: img is not used, use imglist\n");
      exit(1);
    }
    if((outimg=niik_image_add_multiple(imglist,numimglist,vlist))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: outimg=niik_image_average_multiple(imglist,imglist)\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(outimg,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(outimg,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    outimg=niik_image_free(outimg);
  } /* OP = addimgs */

  else if(!strncmp(argv[1],"mulimgs",7)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>,... -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(img!=NULL) {
      fprintf(stderr,"[niikmath] ERROR: img is not used, use imglist\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(imglist[0],datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
    }
    if((outimg=niik_image_mul_multiple(imglist,numimglist))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: outimg=niik_image_average_multiple(imglist,imglist)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    outimg=niik_image_free(outimg);
  } /* OP = mulimgs */

  else if(!strncmp(argv[1],"mgz2nii",7)) {
    /*
     * NOT WORKING!!!
     * --I TRIED BUT FAILED MAY 5, 2012
     */
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.mgz> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] reading mgz file\n");
    exit(1);
  } /* OP = mgz2nii */

  else if(!strncmp(argv[1],"tiffimg",7)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"[%s] usage: tiffimg -in=<img.nii> -out=<out.tif> [-tiff-xslice=<x> | -tiff-yslice=<y> | -tiff-zslice=<z>]\n",argv[0]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -xyz=<dx>,<dy>,<dz>            : resample image [default=off]\n");
      fprintf(stdout,"  -imin=<m>                      : minimal intensity [default=auto]\n");
      fprintf(stdout,"  -imax=<M>                      : maximum intensity [default=auto]\n");
      fprintf(stdout,"  -percent=<P>                   : percentile for maximum intensity [default=maximum]\n");
      fprintf(stdout,"  -tiff-xslices=<x1>[,<x2>...]   : write multiple slices [default=none]\n");
      fprintf(stdout,"                                 : output filename should be prefix without extension\n");
      fprintf(stdout,"  -tiff-yslices=<y1>[,<y2>...]   : write multiple slices [default=none]\n");
      fprintf(stdout,"                                 : output filename should be prefix without extension\n");
      fprintf(stdout,"  -tiff-zslices=<z1>[,<z2>...]   : write multiple slices [default=none]\n");
      fprintf(stdout,"                                 : output filename should be prefix without extension\n");
      exit(1);
    }
    if(xyz[0]>0) { /* if resampling */
      fprintf(stdout,"[niikmath] resampling        %7.3f %7.3f %7.3f\n",xyz[1],xyz[2],xyz[3]);
      NIIK_EXIT((!niik_image_resample_3d_update(img,xyz[1],xyz[2],xyz[3],-1,-1,-1,interp)),fcname,"niik_image_resample_3d_update",9);
    } /* resampling */
    if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,NULL);
    if(niik_check_double_problem(imax)) {
      if(niik_check_double_problem(percent)) {
        imax=niik_image_get_max(img,NULL);
      } else {
        imax=niik_image_get_percentile(img,NULL,percent);
      }
    }

    if(tiff_xslice>=0) {
      NIIK_EXIT((!niik_tiff_write_xslice(outname,img,imin,imax,tiff_xslice)),fcname,"niik_tiff_write_xslice",9);
    } else if(tiff_yslice>=0) {
      NIIK_EXIT((!niik_tiff_write_yslice(outname,img,imin,imax,tiff_yslice)),fcname,"niik_tiff_write_yslice",9);
    } else if(tiff_zslice>=0) {
      NIIK_EXIT((!niik_tiff_write_zslice(outname,img,imin,imax,tiff_zslice)),fcname,"niik_tiff_write_zslice",9);
    } else { /* multiple slices */
      if(tiff_xslices!=NULL) {
        for(n=1; n<tiff_xslices[0]; n++) {
          sprintf(fname,"%s_x%i.tif",outname,tiff_xslices[n]);
          fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
          NIIK_EXIT((tiff_xslices[n]<0),fcname,"xslice too small",9);
          NIIK_EXIT((tiff_xslices[n]>=img->nx),fcname,"xslice too large",9);
          NIIK_EXIT((!niik_tiff_write_xslice(fname,img,imin,imax,tiff_xslices[n])),fcname,"niik_tiff_write_xslice",9);
        } /* each specified xslices */
      } /* x-slices */
      if(tiff_yslices!=NULL) {
        for(n=1; n<tiff_yslices[0]; n++) {
          sprintf(fname,"%s_y%i.tif",outname,tiff_yslices[n]);
          fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
          NIIK_EXIT((tiff_yslices[n]<0),fcname,"yslice too small",9);
          NIIK_EXIT((tiff_yslices[n]>=img->ny),fcname,"yslice too large",9);
          NIIK_EXIT((!niik_tiff_write_yslice(fname,img,imin,imax,tiff_yslices[n])),fcname,"niik_tiff_write_yslice",9);
        } /* each specified yslices */
      } /* y-slices */
      if(tiff_zslices!=NULL) {
        for(n=1; n<tiff_zslices[0]; n++) {
          sprintf(fname,"%s_z%i.tif",outname,tiff_zslices[n]);
          fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
          NIIK_EXIT((tiff_zslices[n]<0),fcname,"zslice too small",9);
          NIIK_EXIT((tiff_zslices[n]>=img->nx),fcname,"zslice too large",9);
          NIIK_EXIT((!niik_tiff_write_zslice(fname,img,imin,imax,tiff_zslices[n])),fcname,"niik_tiff_write_zslice",9);
        } /* each specified zslices */
      } /* z-slices */
      /*else {
        fprintf(stderr,"[niikmath] ERROR: please enter the slice by one of these:\n");
        fprintf(stderr,"           -tiff-xslice=<x>     -tiff-yslice=<y>     | -tiff-zslice=<z>\n");
        fprintf(stderr,"           -tiff-xslices=<x>... -tiff-yslices=<y>... | -tiff-zslicex=<z>...\n");
        exit(1); }*/
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = tiffimg */

  else if(!strncmp(argv[1],"bspline",7)) {
    fprintf(stdout,"  calculates bspline coefficients\n");
    fprintf(stdout,"\n\n  *** TO DO ***\n\n");
    exit(0);
  } /* OP = bspline */

  else if(!strncmp(argv[1],"aregimg",7)) {
    if(img==NULL || refimg==NULL || outname==NULL) {
      fprintf(stdout,"\n  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -matrix=<M.matrix>             initial matrix [default=identity]\n");
      fprintf(stdout,"  -affpar=<Rx>,<Ry>,<Rz>,<Tx>,<Ty>,<Tz>[,<Sx>,<Sy>,<Sz>,<Kx>,<Ky>,<Kz>,<Px>,<Py>,<Pz>\n");
      fprintf(stdout,"                                 initial affine parameters for 3 rotations, 3 translations,\n");
      fprintf(stdout,"                                 3 scaling, 3 shears, and 3 pivot [default=none]\n");
      fprintf(stdout,"                                 Scaling, shear, and pivot are optional\n");
      fprintf(stdout,"  -nmi                           use normalized mutual information [default=correlation coefficient]\n");
      fprintf(stdout,"  -nmi-num=<size1>,<size2>       nmi histogram size [default=32,32 for size1 and size2]\n");
      fprintf(stdout,"  -dof=<DOF>                     set dof [default=6, available: 3,6,7,9,12]\n");
      fprintf(stdout,"  -delta=<D>                     set sampling distance [default=2.5]\n");
      fprintf(stdout,"  -fwhm=<FWHM>                   set FWHM for blurring [default=5.5]\n");
      fprintf(stdout,"  -mask=<refmask.nii>            mask image [default=none]\n");
      fprintf(stdout,"  -outcheck=<outcheck.nii>       create a check image [default=no output]\n");
      exit(0);
    }
    NIIK_EXIT((!niik_image_type_convert(   img,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert, img",9);
    NIIK_EXIT((!niik_image_type_convert(refimg,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert, refimg",9);
    if(maskimg!=NULL) {
      NIIK_EXIT((!niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert, maskimg",9);
    }
    switch(dof) {
    case 12:
      daffpar[11]=daffpar[12]=daffpar[13]=0.05;
    case 9:
      daffpar[8]=daffpar[9]=daffpar[10]=0.05;
    case 7:
      daffpar[7]=0.05;
    case 6:
      daffpar[1]=daffpar[2]=daffpar[3]=15;
      daffpar[4]=daffpar[5]=daffpar[6]=15;
      break;
    case 3:
      daffpar[3]=15;
      daffpar[4]=daffpar[5]=15;
      break;
    }
    affpar[7]=affpar[8]=affpar[9]=affpar[10]=1.0;
    if(niik_check_double_problem(delta)) delta=2.5;
    if(niik_check_double_problem(FWHM))  {
      FWHM = 5.5;
      if(areg_cost==NIIK_REGISTER_NMI) {
        FWHM=2.5;
      }
    }
    if(areg_cost==NIIK_REGISTER_NMI) {
      if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,NULL);
      if(niik_check_double_problem(imax)) imax=niik_image_get_max(img,NULL);
      if(niik_check_double_problem(omin)) omin=niik_image_get_min(refimg,NULL);
      if(niik_check_double_problem(omax)) omax=niik_image_get_max(refimg,NULL);
      if(!nmi_num[0] || !nmi_num[1]) {
        nmi_num[0] = nmi_num[1] = 32;
      }
      if(!niik_aregister_nmi_obj_update_var(imin,imax,omin,omax,nmi_num[0],nmi_num[1])) {
        fprintf(stderr,"[%s] ERROR: niik_aregister_nmi_obj_update_var\n",fcname);
        exit(0);
      }
      nmiobj=niik_aregister_g_nmi_obj_get();
      niik_aregister_nmi_obj_display(nmiobj);
    } else if(areg_cost==NIIK_REGISTER_UNKNOWN) {
      areg_cost = NIIK_REGISTER_CC;
    }
    ctr=niikpt_avg(niikpt_image_get_centroid(img,NULL),niikpt_image_get_centroid(refimg,NULL));
    affpar[14]=ctr.x;
    affpar[15]=ctr.y;
    affpar[16]=ctr.z;
    pt=niikpt_sub(niikpt_image_get_centroid(refimg,NULL),niikpt_image_get_centroid(img,NULL));
    if(affpar[4]==0) affpar[4]=pt.x;
    if(affpar[5]==0) affpar[5]=pt.y;
    if(affpar[6]==0) affpar[6]=pt.z;
    niik_aregister_display_affine(affpar);
    /* new affine register2 */
    if(!strncmp(argv[1],"aregimg2",8)) {
      if(!niik_image_aregister2_test1(refimg,maskimg,img,afmat,affpar,dof,areg_cost,nmiobj)) {
        fprintf(stderr,"ERROR: niik_image_aregister2_test1\n");
        exit(0);
      }
    } else {
      if(!niik_image_aregister_imat(refimg,maskimg,img,NULL,afmat,affpar,daffpar,areg_cost,delta,FWHM)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_aregister(refimg,maskimg,img,NULL,affpar,daffpar,areg_cost,delta,FWHM)\n");
        exit(0);
      }
    }
    if(afmat==NULL) {
      if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar(affpar)\n");
        exit(0);
      }
    } else {
      if((afmat2=niik_aregister_matrix_from_affpar(affpar))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar(affpar)\n");
        exit(0);
      }
      niikmat_display(afmat2);
      if((afmat=niikmat_multiply_free12(afmat,afmat2))==NULL) {
        fprintf(stderr,"ERROR: niikmat_multiply_free12\n");
        exit(0);
      }
    }
    fprintf(stdout,"[%s] writing           %s\n",fcname,outname);
    if(!strncmp(outname+(strlen(outname)-4),".xfm",4)) {
      NIIK_RET0((!niikmat_convert_to_xfm(afmat,img,refimg)),fcname,"niik_convert_to_xfm");
      NIIK_RET0((!niikmat_write_as_linear_xfm(outname,afmat)),fcname,"niik_write");
    } else {
      NIIK_RET0((!niikmat_write(outname,afmat)),fcname,"niikmat_write");
    }
    if(outcheckname!=NULL) {
      fprintf(stdout,"    transform image for checking: %s\n",outcheckname);
      if((outimg=niik_image_affine_transform_3d(img,refimg,afmat,NIIK_INTERP_LINEAR))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d\n");
        exit(1);
      }
      imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
      imglist[0]=outimg;
      imglist[1]=refimg;
      if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_combine\n");
        exit(1);
      }
      if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
      outimg->cal_min=0;
      outimg->cal_max=200;
      fprintf(stdout,"[niikmath] writing output  %s\n",outcheckname);
      niik_image_append_history(outimg,timestamp);
      niik_image_write(outcheckname,outimg);
      /* clean up the temporary variables */
      outimg=niik_image_free(outimg);
      /*imglist[0]=niik_image_free(imglist[0]);*/
      /*imglist[1]=niik_image_free(imglist[1]);*/
      free(imglist);
      imglist=NULL;
    }
    afmat=niikmat_free(afmat);
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    refimg=niik_image_free(refimg);
    exit(0);
  } /* OP = aregimg */

  else if(!strncmp(argv[1],"aregmni",7)) {
    /* available operations:
     *   aregmni-nbcsr aregmni-NBCSR
     *   aregmni-sform
     *   aregmni
     */
    if(img==NULL) {
      if(!strncmp(argv[1],"aregmni-nbcsr",13) ||
          !strncmp(argv[1],"aregmni-NBCSR",13)) {
        fprintf(stdout,"\n  usage: %s %s -in=<img.nii> [-out=<out.matrix> -mask=<seg.nii> -outcheck=<check.nii>]\n",
                argv[0],argv[1]);
        fprintf(stdout,"\n  optional usage:\n");
        fprintf(stdout,"  -mask=<seg.nii>          -use brain mask\n");
        fprintf(stdout,"  -out=<out.matrix>        -write an output matrix\n");
        fprintf(stdout,"  -outcheck=<check.nii>    -write an output image for error checking\n");
        fprintf(stdout,"\n  writes to the same input filename\n");
      } else if(!strncmp(argv[1],"aregmni-sform",13)) {
        fprintf(stdout,"\n  usage: %s %s -in=<img.nii> [-mask=<seg.nii> -outcheck=<check.nii>]\n",
                argv[0],argv[1]);
        fprintf(stdout,"\n  optional usage:\n");
        fprintf(stdout,"  -mask=<seg.nii>          -use brain mask\n");
        fprintf(stdout,"  -outcheck=<check.nii>    -write an output image for error checking\n");
        fprintf(stdout,"\n  writes to the same input filename\n");
      } else if(!strncmp(argv[1],"aregmni",7)) {
        fprintf(stdout,"\n  usage: %s %s -in=<img.nii> -out=<out.matrix> [-mask=<seg.nii> -outcheck=<check.nii>]\n",
                argv[0],argv[1]);
        fprintf(stdout,"\n  optional usage:\n");
        fprintf(stdout,"  -mask=<seg.nii>          -use brain mask\n");
        fprintf(stdout,"  -outcheck=<check.nii>    -write an output image for error checking\n");
      }
      exit(1);
    } /* error messages */
    /* remove negative numbers */
    imin=niik_image_get_min(img,NULL);
    for(i=0; i<img->nvox; i++) niik_image_add_voxel(img,i,-imin);
    if(!strncmp(argv[1],"aregmni-sform",13)) {}
    else if(!strncmp(argv[1],"aregmni",7)) {
      if(outname==NULL) {
        fprintf(stderr,"[niikmath] ERROR: please enter '-out=<out.matrix>'\n");
        exit(0);
      }
    }
    /*
     * start affine registration (aregmni)
     */
    datatype=img->datatype;
    fprintf(stdout,"  affine mni registration\n");
    if(!niik_image_type_convert_scl(img,NIFTI_TYPE_FLOAT32,1)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)\n",fcname);
      exit(1);
    }
    if(maskimg!=NULL) {
      if(!niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[%s] ERROR: niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)\n",fcname);
        exit(1);
      }
    }
    /* read mni images */
    NIIK_EXIT((niik_check_fsldir_exists()>0),fcname,"fsldir is not defined",1);
    sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
    /* sprintf(fname,"/lab2/Kunio/kproj/data/mni_icbm152_t1_tal_nlin_sym_09c.mnc"); */
    fprintf(stdout,"[%s] reading mni img   %s\n",fcname,fname);
    if((mni_img=niik_image_read(fname))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
      exit(1);
    }
    sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
    /* sprintf(fname,"/lab2/Kunio/kproj/data/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"); */
    fprintf(stdout,"[%s] reading mni seg   %s\n",fcname,fname);
    if((mni_seg=niik_image_read(fname))==NULL) {
      fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
      exit(1);
    }
    /* use default values */
    if(niik_check_double_problem(FWHM))  {
      FWHM=6.0;
    }
    if(niik_check_double_problem(delta)) {
      delta=3.2;
    }
    if(areg_cost==NIIK_REGISTER_UNKNOWN) {
      areg_cost=NIIK_REGISTER_CC;
    }
    /*
     * if NBCSR, then do that part
     */
    if(!strncmp(argv[1],"aregmni-nbcsr",13) ||
        !strncmp(argv[1],"aregmni-NBCSR",13)) {
      affpar[0]=0;
      fprintf(stdout,"[niikmath] start NBCSR\n");
      fprintf(stdout,"[niikmath]   base FWHM = %7.3f\n",FWHM);
      fprintf(stdout,"[niikmath]   base dist = %7.3f\n",delta);
      fprintf(stdout,"[niikmath]   cost func = %s\n",niik_aregister_method_string(areg_cost));
      if(!niik_aregister_nbcsr(mni_img,mni_seg,img,maskimg,affpar,areg_cost,FWHM,delta,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_nbcsr(mni_img,mni_seg,img,maskimg,affpar,areg_cost,FWHM,delta,0)\n");
        exit(1);
      }
    } /* NBCSR */
    else {
      /* start the align_mni function!! */
      fprintf(stdout,"[%s] start aregister_align_mni\n",fcname);
      fprintf(stdout,"[%s]   base FWHM = %7.3f\n",fcname,FWHM);
      fprintf(stdout,"[%s]   base dist = %7.3f\n",fcname,delta);
      fprintf(stdout,"[%s]   cost func = %s\n",fcname,niik_aregister_method_string(NIIK_REGISTER_NMI));
      if(!niik_aregister_align_mni(mni_img,mni_seg,img,maskimg,affpar,NIIK_REGISTER_NMI,FWHM,delta)) {
        fprintf(stderr,"[%s] ERROR: niik_aregister_align_mni\n",fcname);
        exit(0);
      }
    }
    /* prepare output (image or text) */
    fprintf(stdout,"[%s] writing output\n",fcname);
    if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar\n");
      exit(1);
    }
    niikmat_display(afmat);
    niik_aregister_display_affine(affpar);
    if(!strncmp(argv[1],"aregmni-sform",13)) {
      niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(img->dx,img->dy,img->dz));
      niikmat_multiply_mat2_free1(niikmat_mat44_matrix(mni_img->sto_xyz),afmat);
      img->sform_code=NIFTI_XFORM_MNI_152;
      for(m=0; m<4; m++) /* forward matrix */
        for(n=0; n<4; n++)
          img->sto_xyz.m[m][n]=afmat->m[m][n];
      if(!niikmat_inverse_update(afmat)) {
        fprintf(stderr,"[niikmath] ERROR: niikmat_inverse_update\n");
        exit(1);
      }
      for(m=0; m<4; m++) /* inverse matrix */
        for(n=0; n<4; n++)
          img->sto_ijk.m[m][n]=afmat->m[m][n];
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing image:    %s\n",img->fname);
      nifti_image_write(img);
    } /* writing to image's sform */
    if(outname!=NULL) { /* writing output transformation matrix */
      fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
      if(!strncmp(outname+(strlen(outname)-4),".xfm",4)) {
        afmat2=niikmat_copy(afmat);
        niikmat_convert_to_xfm(afmat2,img,mni_img);
        NIIK_RET0((!niikmat_write(outname,afmat2)),fcname,"niikmat_write");
        afmat2=niikmat_free(afmat2);
      } else {
        NIIK_RET0((!niikmat_write(outname,afmat)),fcname,"niikmat_write");
      }
    }
    /*
     * writing output transformed image
     */
    if(outcheckname!=NULL) {
      if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update\n");
        exit(1);
      }
      fprintf(stdout,"    transform image for checking: %s\n",outcheckname);
      if((outimg=niik_image_affine_transform_3d(img,mni_img,afmat,NIIK_INTERP_LINEAR))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d\n");
        exit(1);
      }
      imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
      imglist[0]=outimg;
      imglist[1]=mni_img;
      if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_combine\n");
        exit(1);
      }
      if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
      outimg->cal_min=0;
      outimg->cal_max=200;
      fprintf(stdout,"[niikmath] writing output  %s\n",outcheckname);
      niik_image_append_history(outimg,timestamp);
      niik_image_write(outcheckname,outimg);
      /* clean up the temporary variables */
      nifti_image_free(outimg);
      nifti_image_free(imglist[0]);
      nifti_image_free(imglist[1]);
      outimg=mni_img=imglist[0]=imglist[1]=NULL;
      free(imglist);
      imglist=NULL;
    }
    img=niik_image_free(img);
    nifti_image_free(mni_img);
    nifti_image_free(mni_seg);
    img=outimg=NULL;
  } /* OP  = aregmni or aregmni-sform */

  else if(!strncmp(argv[1],"nii2img",7)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  nii2img\n");
      fprintf(stdout,"  usage: %s %s -in=<img> -out=<out.img>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if((fp=fopen(outname,"wb"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",outname);
      exit(0);
    }
    i=fwrite(img->data,img->nvox*img->nbyper,1,fp);
    fclose(fp);
    fprintf(stdout,"[%s] MATLAB: fid=fopen('%s','r'); img=fread(fid,%i,'float'); fclose(fid); img=reshape(img,%i,%i,%i);\n",fcname,outname,img->nvox,img->nx,img->ny,img->nz);
    img=niik_image_free(img);
    exit(0);
  } /* OP = nii2img */

  else if(!strncmp(argv[1],"img2nii",7)) {
    if(inname==NULL || outname==NULL || datatype<0) {
      fprintf(stdout,"  img2nii\n");
      fprintf(stdout,"  usage: %s %s -inname=<img> -out=<out.nii> -vlist=<xdim>,<ydim>,<zdim>,... -<datatype>\n",argv[0],argv[1]);
      exit(1);
    }
    switch(numvlist) {
    case 0:
      fprintf(stderr,"ERROR: please  use -vlist=<xdim>,<ydim>,...\n");
      exit(0);
    case 1:
      if((img=niik_image_init(vlist[0],1,1,1,1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 2:
      if((img=niik_image_init(vlist[0],vlist[1],1,1,1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 3:
      fprintf(stdout,"[niikmath] img2nii 3d image %3i %3i %3i\n",(int)vlist[0],(int)vlist[1],(int)vlist[2]);
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],1,1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 4:
      fprintf(stdout,"[%s] 4-dims\n",fcname);
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(1);
      }
      fprintf(stdout,"[%s] 4-dims initialized\n",fcname);
      break;
    case 5:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 6:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],vlist[5],1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 7:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],vlist[5],vlist[6],
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    } /* dimensions */
    fprintf(stdout,"[%s] reading %s, nvox = %i, nbyper = %i\n",fcname,inname,img->nvox,img->nbyper);
    NIIK_EXIT(((fp=fopen(inname,"rb"))==NULL),fcname,"fopen file",1);
    fprintf(stdout,"[%s] fread image, nvox = %i, nbyper = %i\n",fcname,img->nvox,img->nbyper);
    /*NIIK_EXIT(((fread(img->data,img->nvox*img->nbyper,1,fp))!=(img->nvox*img->nbyper)),fcname,"fread",1);*/
    i=fread(img->data,img->nvox*img->nbyper,1,fp);
    fclose(fp);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",1);
    img=niik_image_free(img);
    exit(0);
  } /* OP = img2nii */

  else if(!strncmp(argv[1],"txt2nii",7)) {
    if(inname==NULL || outname==NULL || datatype<0) {
      fprintf(stdout,"  img2nii\n");
      fprintf(stdout,"  usage: %s %s -inname=<img.txt> -out=<out.nii> -vlist=<xdim>,<ydim>,<zdim>,... -<datatype>\n",argv[0],argv[1]);
      exit(1);
    }
    switch(numvlist) {
    case 0:
      fprintf(stderr,"ERROR: please  use -vlist=<xdim>,<ydim>,...\n");
      exit(0);
    case 1:
      NIIK_RET0(((img=niik_image_init(vlist[0],1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,datatype))==NULL),
                fcname,"niik_image_init");
      break;
    case 2:
      if((img=niik_image_init(vlist[0],vlist[1],1,1,1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 3:
      fprintf(stdout,"[niikmath] img2nii 3d image %3i %3i %3i\n",(int)vlist[0],(int)vlist[1],(int)vlist[2]);
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],1,1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 4:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],1,1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 5:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],1,1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 6:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],vlist[5],1,
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    case 7:
      if((img=niik_image_init(vlist[0],vlist[1],vlist[2],vlist[3],vlist[4],vlist[5],vlist[6],
                              1,1,1,1,1,1,1,datatype))==NULL) {
        fprintf(stderr,"ERROR: niik_image_init\n");
        exit(0);
      }
      break;
    } /* dimensions */
    if((fp=fopen(inname,"r"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",inname);
      exit(0);
    }
    for(i=0; i<img->nvox; i++) {
      NIIK_EXIT((!fscanf(fp,"%lf ",&dval)),fcname,"fscanf",9);
      niik_image_set_voxel(img,i,dval);
    }
    fclose(fp);
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"niik_image_write");
    img=niik_image_free(img);
    exit(0);
  } /* OP = txt2nii */

  else if(!strncmp(argv[1],"imginfo",7)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s imginfo -in=<img.nii>\n",argv[0]);
      exit(1);
    }
    fprintf(stdout,"image info %s\n%s\n\n",img->fname,nifti_image_to_ascii(img));
    img=niik_image_free(img);
    img=NULL;
  } /* OP = imginfo */

  else if(!strncmp(argv[1],"calcvol",7)) {
    if(img==NULL && imglist==NULL) {
      fprintf(stdout,"  usage 1: %s %s -in=<img.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  usage 2: %s %s -imglist=<img1.nii>,<img2.nii>...\n",argv[0],argv[1]);
      exit(1);
    }
    if(img!=NULL) {
      fprintf(stdout,"%9.9f ml calcvol    %s\n",niik_image_count_mask(img) * niik_image_get_voxel_size(img)/1000.0,img->fname);
      img=niik_image_free(img);
      img = NULL;
    }
    if(imglist!=NULL) {
      for(n=0; n<numimglist; n++) {
        fprintf(stdout,"%9.9f ml calcvol    %s\n",niik_image_count_mask(imglist[n]) * niik_image_get_voxel_size(imglist[n])/1000.0,
                imglist[n]->fname);
        nifti_image_free(imglist[n]);
        imglist[n] = NULL;
      }
      free(imglist);
      imglist=NULL;
    }
    exit(0);
  } /* OP = interp */

  else if(!strncmp(argv[1],"flipxyz",7)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(img->ndim>3) {
      fprintf(stderr,"[niikmath] ERROR: no support for more than 3D\n");
      exit(1);
    }
    if(!niik_image_rotate(img,"flipxyz")) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_rotate\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(0);
    }
    img=niik_image_free(img);
  } /* OP = flipxyz */

  else if(!strncmp(argv[1],"adirimg",7) ||
          !strncmp(argv[1],"xdirimg",7) ||
          !strncmp(argv[1],"ydirimg",7) ||
          !strncmp(argv[1],"zdirimg",7)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] %c-dir image\n",argv[1][0]);
    fimg = (float*) img->data;
    switch(argv[1][0]) {
    case 'x':
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; n++,i++) {
            fimg[n]=img->dx*i;
          }
        }
      }
      break;
    case 'y':
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; n++,i++) {
            fimg[n]=img->dy*j;
          }
        }
      }
      break;
    case 'z':
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; n++,i++) {
            fimg[n]=img->dz*k;
          }
        }
      }
      break;
    case 'a':
      size=img->nx*img->ny*img->nz;
      fimg=(float *)calloc(size*3,sizeof(float));
      img->nu=3;
      for(k=n=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; n++,i++) {
            fimg[n       ]=img->dx*i;
            fimg[n+size  ]=img->dy*j;
            fimg[n+size*2]=img->dz*k;
          }
        }
      }
      free(img->data);
      img->data=fimg;
      img->nu=3;
      img->du=1;
      img->ndim=5;
      img->nvox=img->nx*img->ny*img->nz*img->nu;
      break;
    default:
      fprintf(stderr,"[niikmath] ERROR: unknown diretion %c\n",argv[1][0]);
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
  } /* OP = xdirimg ydirimg zdirimg */

  else if(!strncmp(argv[1],"nregimg",7)) {
    if(img==NULL || refimg==NULL || maskimg==NULL || maskref==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -mask=<mask.nii> -refmask=<refseg.nii> -ref=<ref.nii> -out=<out_warp.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -matrix=<img2ref.matrix>        : affine registration matrix\n");
      fprintf(stdout,"  -outcheck=<check.nii>           : output check image\n");
      fprintf(stdout,"  -iter=<iter>                    : iteration [default=5]\n");
      fprintf(stdout,"  -FWHM=<fwhm>                    : FWHM [default=2.0]\n");
      exit(0);
    }
    /*
     * prepare output warp image
     */
    fprintf(stdout,"[niikmath]   prepare output warp image\n");
    datatype=refimg->datatype;
    if((outimg=niik_image_copy_as_type(refimg,NIFTI_TYPE_FLOAT32))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
      exit(0);
    }
    outimg->ndim=outimg->dim[0]=5;
    outimg->nt=outimg->dim[4]=1;
    outimg->nu=outimg->dim[5]=3;
    outimg->nvox=outimg->nx*outimg->ny*outimg->nz*outimg->nt*outimg->nu;
    free(outimg->data);
    outimg->data=(float *)calloc(outimg->nvox,outimg->nbyper);
    /*
     * nonlinear registration
     */
    if(iter<0) iter=5;
    afmat2=niikmat_identity(4,4);
    if(niik_check_double_problem(FWHM)) FWHM=2.0;
    if(FWHM_list!=NULL) {
      for(n=0; n<FWHM_list->num; n++) {
        fprintf(stdout,"[%s] making gradient image   %-9.3f\n",fcname,FWHM_list->v[n]);
        if((gradimg=niik_image_gauss_sobel_filters_with_mag_single_output(refimg,FWHM_list->v[n]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag_single_output\n",fcname);
          return 0;
        }
        if(!niik_image_nregister_demons(outimg,refimg,maskref,refimg,afmat2,gradimg,img,afmat,iter,FWHM)) {
          fprintf(stderr,"[%s] ERROR: niik_image_nregister_demons\n",fcname);
          exit(0);
        }
        gradimg=niik_image_free(gradimg);
      }
    } else {
      if(!niik_image_nregister_demons(outimg,refimg,maskref,refimg,afmat2,NULL,img,afmat,iter,FWHM)) {
        fprintf(stderr,"[%s] ERROR: niik_image_nregister_demons\n",fcname);
        exit(0);
      }
    }

    afmat2=niikmat_free(afmat2);
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(0);
    }
    /*
     * image to check
     */
    if(outcheckname!=NULL) {
      imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
      imglist[0]=refimg;
      imglist[1]=niik_image_copy(refimg);
      if(!niik_image_nregister_intrasubject_apply_warp(outimg,imglist[1],img,afmat,NIIK_INTERP_BSPLINE)) {
        fprintf(stderr,"[%s] ERROR: niik_image_nregister_intrasubject_apply_warp\n",fcname);
        exit(0);
      }
      if(!niik_image_type_convert(imglist[0],datatype)) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
        exit(0);
      }
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outcheckname);
      if(!niik_image_combine_and_write(outcheckname,imglist,2,'t',0)) {
        fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write\n",fcname);
        exit(0);
      }
    } /* outcheckname */
    exit(0);
  } /* OP = nregimg */

  else if(!strncmp(argv[1],"off2pts",7) ||
          !strncmp(argv[1],"obj2pts",7)) {
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<in.off> -out=<out.pts>\n",argv[0],argv[1]);
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing   %s\n",outname);
    if((fp=fopen(outname,"w"))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: fopen\n");
      exit(0);
    }
    for(v=obj->vert; v!=NULL; v=v->next) {
      fprintf(fp,"%15.9f %15.9f %15.9f\n",v->v.x,v->v.y,v->v.z);
    }
    fclose(fp);
    obj=off_kobj_free(obj);
    exit(0);
  } /* OP = off2pts */

  else if(!strncmp(argv[1],"offinfo",7)) {
    if(obj==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"  offinfo\n");
    fprintf(stdout,"      Filename :   %s\n",obj->fname);
    fprintf(stdout,"   No Vertices :   %i\n",obj->nvert);
    fprintf(stdout,"  No Triangles :   %i\n",obj->nface);
    fprintf(stdout,"      No Edges :   %i\n",obj->nedge);
    fprintf(stdout,"  Euler's Char :   %i\n",obj->nvert-obj->nedge+obj->nface);
    pt=off_calc_kobj_pmax(obj);
    fprintf(stdout,"     Max Point :   %15.9f %15.9f %15.9f\n",pt.x,pt.y,pt.z);
    pt=off_calc_kobj_pmin(obj);
    fprintf(stdout,"     Min Point :   %15.9f %15.9f %15.9f\n",pt.x,pt.y,pt.z);
    off_display_kobj_edge_stats(obj);
    off_kobj_display_valency_stats(obj);
    obj=off_kobj_free(obj);
    exit(0);
  } /* OP = offinfo */

  else if(!strncmp(argv[1],"maskimg",7)) {
    if(img==NULL || maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img_input.nii> -mask=<mask_input.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -outcheck=<check.nii>        : writes color image\n");
      exit(0);
    }
    if(outcheckname!=NULL) {
      NIIK_EXIT(((outimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy",9);
    }
    NIIK_EXIT((!niik_image_mask(img,maskimg)),fcname,"niik_image_mask",9);
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    if(outcheckname!=NULL) {
      fprintf(stdout,"[niikmath] preparing output check image %s\n",outcheckname);
      img=niik_image_free(img);
      NIIK_EXIT(((img=niik_image_mask_add_red_color_uint8(outimg,maskimg))==NULL),fcname,"niik_image_mask_add_red_color_uint8",9);
      fprintf(stdout,"[niikmath] writing check file %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      NIIK_EXIT((!niik_image_write(outcheckname,img)),fcname,"niik_image_write",9);
      outimg=niik_image_free(outimg);
    } /* outcheckname */
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
  } /* OP = maskimg */

  else if(!strncmp(argv[1],"maskout",7)) {
    if(img==NULL || maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img_input.nii> -mask=<mask_input.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -outcheck=<check.nii>           :  creates a color image to check output\n");
      exit(0);
    }
    if(outcheckname!=NULL) {
      NIIK_RET0(((outimg=niik_image_copy(img))==NULL),fcname,"niik_image_copy");
    }
    NIIK_RET0((!niik_image_maskout(img,maskimg)),fcname,"niik_image_maskout");
    if(datatype>=0) {
      NIIK_RET0((!niik_image_type_convert(img,datatype)),fcname,"niik_image_type_convert");
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"niik_image_write");
    if(outcheckname!=NULL) {
      fprintf(stdout,"[niikmath] preparing output check image %s\n",outcheckname);
      img=niik_image_free(img);
      NIIK_RET0(((img=niik_image_mask_add_red_color_uint8(outimg,maskimg))==NULL),fcname,"niik_image_mask_add_red_color_uint8");
      fprintf(stdout,"[niikmath] writing check file %s\n",outcheckname);
      NIIK_RET0((!niik_image_write(outcheckname,img)),fcname,"niik_image_write");
      outimg=niik_image_free(outimg);
    } /* outcheckname */
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = maskout */

  else if(!strncmp(argv[1],"montage",7)) {
    if(imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<in1.nii>,<in2.nii>... -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -matrixlist=<mat1>,<mat2>...         transformation matrix\n");
      fprintf(stdout,"  -ref=<ref.nii>                       reference image\n");
      fprintf(stdout,"                                       if matrixlist is used\n");
      fprintf(stdout,"  -mask=<mask.nii>                     crop image to include just around the mask\n");
      fprintf(stdout,"  -roi=<xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>\n");
      fprintf(stdout,"                                       crop image\n");
      fprintf(stdout,"  -val=<ncol>                          number of columns\n");
      exit(1);
    }
    if(niik_check_double_problem(inval)) {
      if(numimglist<=4) n = numimglist;
      else {
        n = floor(sqrt(numimglist)+0.5);
      }
    } else {
      n = (int)inval;
    }
    if(maskimg!=NULL) {
      if(!niik_image_get_mask_crop_roi(maskimg,img_dim_min+1,img_dim_max+1,img_dim_min+2,img_dim_max+2,img_dim_min+3,img_dim_max+3)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_get_mask_crop_roi(maskimg,img_dim_min+1,img_dim_max+1,img_dim_min+2,img_dim_max+2,img_dim_min+3,img_dim_max+3)\n");
        exit(1);
      }
      if(ijk[0]>0) {
        for(i=1; i<=3; i++) {
          img_dim_min[i]=NIIK_IMAX(0,img_dim_min[i]-ijk[i]);
        }
        for(i=1; i<=3; i++) {
          img_dim_max[i]=NIIK_IMAX(img_dim_max[i]+ijk[i],maskimg->dim[i]-1);
        }
      }
      fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i\n",img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    } else if(img_dim_min[0]>0) {
      fprintf(stdout,"[niikmath]   roi %3i %3i %3i  -  %3i %3i %3i\n",img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    } else {
      for(i=1; i<=3; i++) {
        img_dim_min[i]=img_dim_max[i]=-1;
      }
    }
    if(nummatlist==numimglist) {
      if(refimg!=NULL) {
        for(i=0; i<numimglist; i++) {
          fprintf(stdout,"[niikmath] image transformation %s  %s\n",niik_interpolate_string(interp),imglist[i]->fname);
          if(!niik_image_affine_transform_3d_update(imglist[i],refimg,matlist[i],interp)) {
            fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
            exit(0);
          }
          /*sprintf(fname,"tmp%i.nii.gz",n);
            niik_image_write(fname,imglist[i]);*/
        } /* nimg */
      } else {
        fprintf(stderr,"ERROR: not implemented yet\n");
        fprintf(stderr,"       please use -ref=<ref.nii>\n");
        exit(0);
      }
    } /* transformation */
    else if(nummatlist==0 && numimglist>0) {
      fprintf(stdout,"[niikmath] no matrix transformation\n");
    } else {
      fprintf(stderr,"ERROR: #matrix and #image did not match\n");
      exit(0);
    }
    if((outimg=niik_image_montage(imglist,numimglist,n,img_dim_min[1],img_dim_min[2],img_dim_min[3],img_dim_max[1],img_dim_max[2],img_dim_max[3]))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_montage\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(outimg,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    } /* convert type */
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,outimg)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    exit(0);
  } /* OP = montage */

  else if(!strncmp(argv[1],"evalreg",7)) {
    fprintf(stdout,"[%s] evaluate registrtaion\n",fcname);
    if(img==NULL || refimg==NULL) {
      fprintf(stdout,"  usage: -in=<img> -ref=<ref>\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage\n");
      fprintf(stdout,"  -mask=<mask>            : volume of interest\n");
      exit(1);
    }

    if(maskimg!=NULL && radius>0.1) {
      NIIK_EXIT((!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_DILATE,radius)),fcname,"niik_image_morph_3d_radius",9);
    }
    NIIK_EXIT((img->nvox!=refimg->nvox),fcname,"different nvox",9);
    xsum=ysum=xssq=yssq=wsum=0;
    for(i=0,xysum=0; i<img->nvox; i++) {
      if(maskimg!=NULL)
        if(niik_image_get_voxel(maskimg,i)<0.5)
          continue;
      d1 = niik_image_get_voxel(img,i);
      d2 = niik_image_get_voxel(refimg,i);
      xsum  += d1;
      ysum  += d2;
      xssq  += d1*d1;
      yssq  += d2*d2;
      xysum += d1*d2;
      wsum  += 1.0;
    } /* each voxel */
    fprintf(stdout,"[%s] CC      %0.12f\n",fcname,(wsum*xysum - xsum*ysum) / sqrt(wsum*xssq-xsum*xsum) / sqrt(wsum*yssq-ysum*ysum));
    img=niik_image_free(img);
    refimg=niik_image_free(refimg);
    exit(0);
  } /* OP = evalreg */

  else if(!strncmp(argv[1],"histo2d",7)) {
    if(img==NULL || refimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: histo2d -in=<img_x> -ref=<img_y> -out=<histo2d.txt>\n");
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -imin=<xmin>           : min for img_x [default=auto]\n");
      fprintf(stdout,"  -imax=<xmax>           : max for img_x [default=auto]\n");
      fprintf(stdout,"  -imin=<ymin>           : min for img_y [default=auto]\n");
      fprintf(stdout,"  -imax=<ymax>           : max for img_y [default=auto]\n");
      fprintf(stdout,"  -histo-num=<N>         : number of bins [default=128]\n");
      fprintf(stdout,"  -mask=<mask.nii>       : mask image\n");
      exit(0);
    }
    if(histo_num<0) histo_num=-histo_num;
    histomat=niikmat_init(histo_num,histo_num);
    if(niik_check_double_problem(imin)) imin = niik_image_get_min(img,maskimg);
    if(niik_check_double_problem(imax)) imax = niik_image_get_max(img,maskimg);
    if(niik_check_double_problem(omin)) omin = niik_image_get_min(refimg,maskimg);
    if(niik_check_double_problem(omax)) omax = niik_image_get_max(refimg,maskimg);
    fprintf(stdout,"[%s] histo2d parameters:\n",fcname);
    fprintf(stdout,"  x_img min/max num   : %12.6f %12.6f %9i\n",imin,imax,histo_num);
    fprintf(stdout,"  y_img min/max num   : %12.6f %12.6f %9i\n",omin,omax,histo_num);
    if(!niik_image_histogram_2d(img,refimg,maskimg,imin,(imax-imin)/(histo_num-1.0),omin,(omax-omin)/(histo_num-1.0),histomat,1)) {
      fprintf(stderr,"[%s] ERROR: niik_image_histogra_2d\n",fcname);
      exit(0);
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niikmat_write(outname,histomat);
    exit(0);
  } /* OP = histo2d */

  else if(!strncmp(argv[1],"ocvobj",6)) {
    if(maskimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -mask=<mask.nii> -out=<out.nii>[,<out.off>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -elen=<edge>           : target edge length [default=3.0]\n");
      fprintf(stdout,"  -matrix=<mni.matrix>   : affine MNI registration matrix\n");
      exit(0);
    }
    if(niik_check_double_problem(elen)) elen=3.0;
    if(obj==NULL && afmat!=NULL) {
      if(((obj=niik_image_bseg_get_brain_surface(afmat,0,0)))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_bseg_get_brain_surface\n",fcname);
        exit(0);
      }
    } /* get object if not found */
    if(!niik_off_outer_contour_object(maskimg,obj,elen)) {
      fprintf(stderr,"[%s] ERROR: niik_off_outer_contour_object\n",fcname);
      exit(0);
    }
    fprintf(stdout,"[%s] off volume: %12.5f\n",fcname,off_get_kobj_volume(obj)/1000.0);
    fprintf(stdout,"[%s] writing output  %s\n",fcname,outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[%s] ERROR: off_kobj_write_off\n",fcname);
      exit(0);
    }
    exit(0);
  } /* OP = ocvobj */

  else if(!strncmp(argv[1],"athresh",6)) {
    if(img==NULL || afmat==NULL) {
      fprintf(stderr,"ERROR: test function athresh\n");
      fprintf(stderr,"  usage -in=<img.nii> -matrix=<mat.matrix>\n");
      exit(0);
    }
    thresh = 1.2;
    if(!niik_image_bseg_basic_thresh(img,afmat,&thresh,0)) {
      fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh\n");
      exit(0);
    }
    fprintf(stdout,"  threshold = %.3f\n",thresh);
    exit(0);
  }

  else if(!strncmp(argv[1],"padimg",6)) {
    fprintf(stdout,"[niikmath] image padding\n");
    if(img==NULL || outname==NULL || ijk[0]==0) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ijk=<x>,<y>,<z> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    NIIK_EXIT((!niik_image_pad_3d(img,(int)ijk[1],(int)ijk[2],(int)ijk[3])),fcname,"niik_image_pad_3d",9);
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = padimg */

  else if(!strncmp(argv[1],"subimg",6)) {
    fprintf(stdout,"[niikmath] image subtraction\n");
    if(img2!=NULL && refimg==NULL) {
      refimg=img2;
      img2=NULL;
    }
    if(img==NULL || refimg==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  out = img - ref\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -vlist=<v1>,<v2>      : out = img * v1 - ref * v2\n");
      fprintf(stdout,"  -val=<v2>             : out = img - ref * v2\n");
      exit(0);
    }
    if(niik_image_cmp_dim(img,refimg)>0) {
      fprintf(stderr,"[niikmath] ERROR: ref image and image have different dimensions\n");
      exit(1);
    }
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_voxels_as_double_vector\n");
      exit(1);
    }
    if(vlist!=NULL) {
      for(i=0; i<img->nvox; i++) {
        dimg[i] = dimg[i] * vlist[0] - vlist[1] * niik_image_get_voxel(refimg,i);
      }
    } else if(!niik_check_double_problem(inval)) {
      for(i=0; i<img->nvox; i++) {
        dimg[i] -= inval * niik_image_get_voxel(refimg,i);
      }
    } else {
      for(i=0; i<img->nvox; i++) {
        dimg[i] -= niik_image_get_voxel(refimg,i);
      }
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    free(dimg);
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = subimg */

  else if(!strncmp(argv[1],"divimg",4)) {
    if(img==NULL || img2==NULL || outname==NULL) {
      fprintf(stdout,"  divide an image <img.nii> by another <img2.nii>\n");
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -img2=<img2.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(0);
    }
    if(datatype<0) datatype=img->datatype;
    if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT64)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,NIFTI_TYPE_FLOAT64)\n");
      exit(1);
    }
    for(n=0; n<img->nvox; n++) {
      niik_image_mul_voxel(img,n,1.0/niik_image_get_voxel(img2,n));
    }
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(1);
    exit(0);
  } /* OP = divimg */

  else if(!strncmp(argv[1],"remesh",6)) {
    fprintf(stdout,"[niikmath] remesh\n");
    if(obj==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -obj=<obj.off> -elen=<L> -out=<out.off>\n",argv[0],argv[1]);
      exit(0);
    }
    if(iter<0) iter=10;
    fprintf(stdout,"[niikmath] elen = %7.4f\n",elen);
    fprintf(stdout,"[niikmath] iter = %i\n",iter);
    if(!off_remesh_kobj(obj,elen,iter,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_remesh_kobj(obj,elen,iter,0)\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    if(!off_kobj_write_offply(outname,obj,0)) {
      fprintf(stderr,"[niikmath] ERROR: off_kobj_write_offply(outname,obj,0)\n");
      exit(0);
    }
    obj = off_kobj_free(obj);
    exit(0);
  } /* OP = remesh */

  else if(!strncmp(argv[1],"interp",6)) {
    if((img==NULL && imglist==NULL)) {
      fprintf(stdout,"  usage: %s interp -in=<img.nii> -ijk=<i>,<j>,<k>\n",argv[0]);
      exit(1);
    }
    if(imglist!=NULL) {
      if(ijk[0]>0) {
        vlist=(double *)calloc(999,sizeof(double));
        for(n=0; n<numimglist; n++) {
          if(interp==NIIK_INTERP_BSPLINE) {
            if(!niik_image_interpolate_convert_3d_bspline_coeff(imglist[n])) {
              fprintf(stderr,"[niikmath] ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n");
              exit(0);
            }
          }
          pt=niikpt_val(ijk[1]*imglist[n]->dx,ijk[2]*imglist[n]->dy,ijk[3]*imglist[n]->dz,0);
          if(!niik_image_interpolate_3d_xyz_update(imglist[n],pt,interp,vlist)) {
            fprintf(stderr,"[niikmath] ERROR: niik_image_interpolate_3d_xyz_update\n");
            exit(0);
          }
          k=imglist[n]->nvox/imglist[n]->nx/imglist[n]->ny/imglist[n]->nz;
          for(m=0; m<k; m++) {
            fprintf(stdout,"  %7.3f %7.3f %7.3f    %15.9f\n",pt.x,pt.y,pt.z,vlist[m]);
          }
        } /* each image */
        free(vlist);
      } /* ijk */
      else if(xyz[0]>0) {
        fprintf(stderr,"[niikmath] ERROR: not implemented yet\n");
        exit(0);
      } else if(afmat!=NULL) {
        fprintf(stderr,"[niikmath] ERROR: not implemented yet\n");
        exit(0);
      } else {
        fprintf(stderr,"[niikmath] ERROR: not implemented yet\n");
        exit(0);
      }
      exit(0);
    } else { /* one image */
      if(ijk[0]>0) {
        pt.x=ijk[1];
        pt.y=ijk[2];
        pt.z=ijk[3];
      } else if(xyz[0]>0) {
        pt.x=xyz[1]/img->dx;
        pt.y=xyz[2]/img->dy;
        pt.z=xyz[3]/img->dz;
      } else if(afmat!=NULL) {
        fprintf(stdout,"[niikmath] using 3d-point list\n");
        if(afmat->col!=3) {
          fprintf(stderr,"[niikmath] ERROR: matrix is used for 3d points\n");
          niikmat_display(afmat);
          exit(0);
        }
        if(interp==NIIK_INTERP_BSPLINE) {
          fprintf(stdout,"[niikmath] calculate b-spline coefficients\n");
          if(!niik_image_interpolate_convert_3d_bspline_coeff(img)) {
            fprintf(stderr,"[niikmath] ERROR: niik_image_interpolate_convert_3d_bspline_coeff(img)\n");
            exit(1);
          }
        }
        m=img->nt*img->nu*img->nv*img->nw;
        if(m>1) {
          fprintf(stderr,"[niikmath] ERROR: not a 3d image\n");
          exit(0);
        }
        for(n=0; n<afmat->row; n++) {
          pt.x=afmat->m[n][0];
          pt.y=afmat->m[n][1];
          pt.z=afmat->m[n][2];
          pt.w=0;
          fprintf(stdout,"  %7.3f %7.3f %7.3f    %15.9f\n",pt.x,pt.y,pt.z,
                  niik_image_interpolate_3d(img,pt,interp));
        }
        exit(0);
      } /* using text file to read points */
      else {
        fprintf(stderr,"[niikmath] ERROR: please use '-ijk=<i>,<j>,<k>' or '-xyz=<i>,<j>,<k>'\n");
        fprintf(stderr,"           ijk for voxel coordinates and xyz for voxel coordinates\n");
        fprintf(stderr,"           with voxel spacing; that is 'ijk*pixdim'\n");
        exit(1);
      }
      if(interp==NIIK_INTERP_BSPLINE) {
        fprintf(stdout,"[niikmath] calculate b-spline coefficients\n");
        if(!niik_image_interpolate_convert_3d_bspline_coeff(img)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_interpolate_convert_3d_bspline_coeff(img)\n");
          exit(1);
        }
      }
      m=img->nt*img->nu*img->nv*img->nw;
      if(m<=0) m=1;
      /* fprintf(stdout,"m=%i\n",m); */
      vlist=(double *)calloc(m,sizeof(double));
      if(!niik_image_interpolate_3d_ijk_update(img,pt,interp,vlist)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_interpolate_3d_ijk_update(img,pt,vlist)\n");
        exit(1);
      }
      for(n=0; n<m; n++) {
        fprintf(stdout,"  %7.3f %7.3f %7.3f    %15.9f\n",ijk[1],ijk[2],ijk[3],vlist[n]);
      }
      free(vlist);
      img=niik_image_free(img);
    }
  } /* OP = interp */

  else if(!strncmp(argv[1],"swapzv",6) ||
          !strncmp(argv[1],"swapvz",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_ISWAP(&img->nz,&img->nv);
    NIIK_ISWAP(&img->dim[3],&img->dim[6]);
    NIIK_FSWAP(&img->dz,&img->dv);
    NIIK_FSWAP(&img->pixdim[3],&img->pixdim[6]);
    if(img->dw>1) img->ndim=img->dim[0]=7;
    else img->ndim=img->dim[0]=6;
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    img=niik_image_free(img);
  } /* OP = swapzv swapvz */

  else if(!strncmp(argv[1],"swaput",6) ||
          !strncmp(argv[1],"swaptu",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_ISWAP(&img->nu,&img->nt);
    NIIK_ISWAP(&img->dim[4],&img->dim[5]);
    NIIK_FSWAP(&img->du,&img->dt);
    NIIK_FSWAP(&img->pixdim[4],&img->pixdim[5]);
    if(img->du>1) img->ndim=img->dim[0]=5;
    else img->ndim=img->dim[0]=4;
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
    img = NULL;
  } /* OP = swaput swaptu */

  else if(!strncmp(argv[1],"swaptv",6) ||
          !strncmp(argv[1],"swapvt",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_ISWAP(&img->nv,&img->nt);
    NIIK_ISWAP(&img->dim[4],&img->dim[6]);
    NIIK_FSWAP(&img->dv,&img->dt);
    NIIK_FSWAP(&img->pixdim[4],&img->pixdim[6]);
    if(img->dv>1) img->ndim=img->dim[0]=6;
    else img->ndim=img->dim[0]=6;
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
    img = NULL;
  } /* OP = swaptv swapvt */

  else if(!strncmp(argv[1],"swapuv",6) ||
          !strncmp(argv[1],"swapvu",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii>\n",argv[0],argv[1]);
      exit(1);
    }
    NIIK_ISWAP(&img->nv,&img->nu);
    NIIK_ISWAP(&img->dim[5],&img->dim[6]);
    NIIK_FSWAP(&img->dv,&img->du);
    NIIK_FSWAP(&img->pixdim[5],&img->pixdim[6]);
    if(img->dv>1) img->ndim=img->dim[0]=6;
    else img->ndim=img->dim[0]=6;
    fprintf(stdout,"[niikmath] writing output  %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
    img = NULL;
  } /* OP = swapuv swapvu */

  else if(!strncmp(argv[1],"thresh",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii> -thresh=<T>\n",argv[0],argv[1]);
      fprintf(stdout,"  -threshold <in.matrix> with threshold = <T>\n");
      fprintf(stdout,"  -writes <out.matrix>\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -uthresh=<U>          : upper threshold [default=no upper threshold]\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(thresh)) {
      fprintf(stderr,"[%s] ERROR: please enter threshold -thresh=<thresh>\n",fcname);
      exit(1);
    }
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_double_vector\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] threshold %.5f\n",fcname,thresh);
    for(i=0; i<img->nvox; i++) {
      dimg[i]=(dimg[i]>=thresh);
    }
    if(!niik_check_double_problem(uthresh)) {
      fprintf(stdout,"[niikmath] upper threshold %.5f\n",uthresh);
      for(i=0; i<img->nvox; i++) {
        if(dimg[i]>0) {
          if(niik_image_get_voxel(img,i)>uthresh)
            dimg[i]=0;
        }
      }
    }
    if(outcheckname!=NULL) {
      if((outimg=niik_image_copy(img))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
        exit(0);
      }
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    free(dimg);
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
      exit(1);
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    if(outcheckname!=NULL) {
      fprintf(stdout,"[niikmath] preparing output check image %s\n",outcheckname);
      if((maskimg=niik_image_mask_add_red_color_uint8(outimg,img))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mask_add_red_color_uint8\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing check file %s\n",outcheckname);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outcheckname,maskimg)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outcheckname);
        exit(0);
      }
      maskimg=niik_image_free(maskimg);
      outimg=niik_image_free(outimg);
    } /* outcheckname */
    img=niik_image_free(img);
    exit(0);
  } /* OP = thresh */

  else if(!strncmp(argv[1],"dilate",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage 1: %s %s -in=<in.nii> -out=<out.nii> -radius=<R>\n",argv[0],argv[1]);
      fprintf(stdout,"  usage 2: %s %s -in=<in.nii> -out=<out.nii> -val=<kernel_size> [-diamond]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -dilate image <in.matrix> with spherical kernel with radius = <R>\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(!niik_check_double_problem(inval)) {
      fprintf(stdout,"[niikmath]   dilation dim %i\n",(int)inval);
      if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_DILATE,morph_kernel_shape,(int)inval)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_morph_3d_mask\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outname,img)) {
        fprintf(stderr,"ERROR: writing file\n");
        exit(0);
      }
      img=niik_image_free(img);
      exit(0);
    }
    if(niik_check_double_problem(radius)) {
      fprintf(stderr,"[niikmath] ERROR: please set radius by -radius=<R>\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   dilation R %8.3f\n",radius);
    NIIK_RET0((!niik_image_morph_3d_radius(img,NIIK_MORPH_DILATE,radius)),fcname,"niik_image_morph_3d_radius");
    if(datatype>=0) {
      NIIK_RET0((!niik_image_type_convert(img,datatype)),fcname,"niik_image_type_convert");
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: writing file\n");
      exit(0);
    }
    img=niik_image_free(img);
  } /* OP = dilate */

  else if(!strncmp(argv[1],"ridler",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii> [-thresh=<T>]\n",argv[0],argv[1]);
      fprintf(stdout,"  -threshold image with Ridler's algorithm\n");
      fprintf(stdout,"  -writes threshold mask <out.matrix>\n");
      fprintf(stdout,"  -<T> is a value used by the algorithm and\n");
      fprintf(stdout,"   ranges approximately [0.5,4.0] [default=2]\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(thresh)) thresh=2;
    if(!niik_image_thresh_ridler(img,maskimg,&thresh)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_thresh_ridler\n");
      exit(1);
    }
    if(maskimg!=NULL) {
      nifti_image_free(maskimg);
      maskimg=NULL;
    }
    fprintf(stdout,"[niikmath]   ridler threshold  %8.3f\n",thresh);
    if(!niik_image_threshold(img,thresh)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_threshold\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    if(nifti_set_filenames(img,outname,0,img->byteorder)) {
      fprintf(stderr,"[niikmath] ERROR: nifti_set_filenames %s\n",outname);
      exit(1);
    }
    nifti_image_write(img);
    img=niik_image_free(img);
    img=NULL;
  } /* OP = ridler */

  else if(!strncmp(argv[1],"sobelm",6) ||
          !strncmp(argv[1],"sobelx",6) ||
          !strncmp(argv[1],"sobely",6) ||
          !strncmp(argv[1],"sobelz",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -directional / magnitude sobel filter\n");
      fprintf(stdout,"  -input image <in.nii>\n");
      fprintf(stdout,"  -output image <out.nii>\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -FWHM=<fwhm>       : blurs image before applying sobel filter\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(!niik_check_double_problem(FWHM)) {
      fprintf(stdout,"[%s]   gaussian filter %7.3f\n",fcname,FWHM);
      if(!niik_image_filter_gaussian_update(img,11,FWHM)) {
        fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian_update\n",fcname);
        exit(0);
      }
    }
    fprintf(stdout,"[%s]   sobel filter (%s)\n",fcname,argv[1]);
    if(!niik_image_sobel_filter_update(img,argv[1][5])) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_sobel_filter\n");
      exit(1);
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(1);
    }
    img=niik_image_free(img);
  } /* OP = sobelx sobely sobelz */

  else if(!strncmp(argv[1],"bounds",6) ||
          !strncmp(argv[1],"bounds-inc",10) ||
          !strncmp(argv[1],"bounds-color",12)) {
    if(img==NULL || imglist==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -imglist=<mask1.nii>,<mask2.nii>... -out=<out.nii> [-vlist=<value1>,<value2>\n",
              argv[0],argv[1]);
      fprintf(stdout,"  -creates boundary image\n");
      fprintf(stdout,"  -<in.nii> is the overlay image\n");
      fprintf(stdout,"  -<mask1.nii>,<mask2.nii>... are the mask(s)\n");
      fprintf(stdout,"  -optionally set the intensities (or color) of the boundary\n");
      fprintf(stdout,"  -writes boundary image <out.matrix>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(maskimg!=NULL) {
      fprintf(stderr,"WARNING: mask is not used in this subroutine\n");
      exit(1);
    }
    if(!strncmp(argv[1],"bounds-color",12)) { /* color */
      if(numimglist!=(numvlist/3)) {
        fprintf(stderr,"[niikmath] ERROR: no of mask files [%i] and no of values [%i->%i] did not match\n",numimglist,numvlist,numvlist/3);
        exit(1);
      }
    } else if(numvlist>0) {
      if(numimglist!=numvlist) {
        fprintf(stderr,"[niikmath] ERROR: no of mask files [%i] and no of values [%i] did not match\n",numimglist,numvlist);
        exit(1);
      }
    }
    dval=niik_image_get_max(img,NULL) + 50;
    for(n=0; n<numimglist; n++) { /* for each mask image */
      if(!strncmp(argv[1],"bounds-color",12)) { /* color */
        fprintf(stdout,"[%s]   boundary %8.2f %8.2f %8.2f  %s\n",fcname,vlist[n*3],vlist[n*3+1],vlist[n*3+2],imglist[n]->fname);
        if(!niik_image_boundary(img,imglist[n],vlist+n*3,1,1,0)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
          exit(1);
        }
      } else if(!strncmp(argv[1],"bounds-inc",10)) { /* inclusive */
        if(numvlist>0) {
          dval=vlist[n];
        }
        if(!niik_image_boundary(img,imglist[n],&dval,1,0,1)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
          exit(1);
        }
      } else { /* exclusive */
        if(numvlist>0) {
          dval=vlist[n];
        }
        if(!niik_image_boundary(img,imglist[n],&dval,1,0,0)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
          exit(1);
        }
      }
      nifti_image_free(imglist[n]);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: writing file\n");
      exit(0);
    }
    img=niik_image_free(img);
    img=NULL;
    free(imglist);
    imglist=NULL;
  } /* OP = bounds bounds-inc bounds-color */

  else if(!strncmp(argv[1],"median",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: median %s -in=<in.nii> -out=<out.nii> -radius=<R>\n",argv[1]);
      fprintf(stdout,"  -filters image <in.nii> with median spherical kernel with radius = <R>in mm\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"\n");
      exit(0);
    }
    if(niik_check_double_problem(radius)) {
      fprintf(stderr,"[%s] ERROR: please set FWHM by -radius=<R>\n",fcname);
      exit(0);
    }
    fprintf(stdout,"[%s]   median filter   %8.3f\n",fcname,radius);
    if(!niik_image_filter_median_radius(img,maskimg,radius)) {
      fprintf(stderr,"[%s] ERROR: niik_image_filter_median_radius\n",fcname);
      exit(0);
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((niik_image_write(outname,img)==0),fcname,"niik_image_write",9);
    img=niik_image_free(img);
    exit(0);
  } /* OP = median */

  else if(!strncmp(argv[1],"setval",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<K> -vlist=<xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>\n",argv[0],argv[1]);
      exit(0);
    }
    if(datatype<0) datatype=img->datatype;
    if(niik_check_double_problem(inval)) {
      fprintf(stderr,"[niikmath] ERROR: please use '-val=<K>' to set a new voxel value\n");
      exit(0);
    }
    if(numvlist!=6) {
      fprintf(stderr,"[niikmath] ERROR: please use '-vlist=<xmin>,<ymin>,<zmin>,<xmax>,<ymax>,<zmax>' to set the ROI\n");
      exit(0);
    }
    for(n=1,m=0; n<=3; n++) {
      img_dim_min[n] = vlist[m++];
    }
    for(n=1; n<=3; n++) {
      img_dim_max[n] = vlist[m++];
    }
    for(n=1; n<=3; n++) {
      if(img_dim_min[n] > img_dim_max[n]) NIIK_ISWAP(&img_dim_min[n],&img_dim_max[n]);
      if(img_dim_min[n]<0) img_dim_min[n]=0;
      if(img_dim_max[n]>=img->dim[n]) img_dim_max[n]=img->dim[n]-1;
    }
    fprintf(stdout,"[niikmath] ROI = %3i %3i %3i   ->   %3i %3i %3i\n",
            img_dim_min[1],img_dim_min[2],img_dim_min[3],
            img_dim_max[1],img_dim_max[2],img_dim_max[3]);
    for(k=img_dim_min[3]; k<=img_dim_max[3]; k++) {
      for(j=img_dim_min[2]; j<=img_dim_max[2]; j++) {
        for(i=img_dim_min[1]; i<=img_dim_max[1]; i++) {
          n = i + j*img->nx + k*img->nx*img->ny;
          niik_image_set_voxel(img,n,inval);
        }
      }
    }
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
    exit(0);
  } /* OP = setval */

  else if(!strncmp(argv[1],"iscale",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<brainmask.nii.gz> [options]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -omin=<omin>         set output min value [default=0]\n");
      fprintf(stdout,"  -omax=<omax>         set output max value [default=140]\n");
      fprintf(stdout,"  -imin=<imin>         set input min value [default=auto]\n");
      fprintf(stdout,"  -imax=<imax>         set input max value [default=auto]\n");
      fprintf(stdout,"  -percent=<p>         set input max from percentile [default=max]\n");
      exit(0);
    }
    if(niik_check_double_problem(imax)) {
      if(!niik_check_double_problem(percent)) {
        imax=niik_image_get_percentile(img,maskimg,percent);
      }
    }
    if(!niik_image_iscale(img,imin,imax,omin,omax)) {
      fprintf(stderr,"[%s] ERROR: niik_image_iscale(img,NIIKMAX,NIIKMAX,omin,omax)\n",fcname);
      exit(0);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
        exit(0);
      }
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outname);
      exit(0);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = iscale */

  else if(!strncmp(argv[1],"binimg",6)) {
    if(img==NULL || outname==NULL || vlist==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<brainmask.nii.gz> -vlist=<v1>[,<v2>...] [options]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      exit(0);
    }
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_voxels_as_double_vector\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      for(n=0; n<numvlist; n++) {
        if(dimg[i]==vlist[n]) break;
      }
      if(n<numvlist) dimg[i]=1;
      else dimg[i]=0;
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    free(dimg);
    if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    niik_image_write(outname,img);
    img=niik_image_free(img);
    exit(0);
  } /* OP = binimg */

  else if(!strncmp(argv[1],"setone",6)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      exit(1);
    }
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_voxels_as_double_vector\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      dimg[i]=1;
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    free(dimg);
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    if(nifti_set_filenames(img,outname,0,img->byteorder)) {
      fprintf(stderr,"[niikmath] ERROR: nifti_set_filenames %s\n",outname);
      exit(1);
    }
    nifti_image_write(img);
    img=NULL;
  }  /* OP = setone */

  else if(!strncmp(argv[1],"inorm2",6)) {
    /* intensity normalization function */
    fprintf(stdout,"[niikmath] intensity normalization using histogram\n");
    if(outname==NULL || img==NULL || refimg==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"    -step=<S>                 number of steps [default=5]\n");
      fprintf(stdout,"    -mask=<img_mask.nii>      mask for <img.nii>\n");
      fprintf(stdout,"    -refmask=<ref_mask.nii>   mask for <ref.nii>\n");
      fprintf(stdout,"\n");
      exit(0);
    }
    if(step<=0) step=5;
    if(!niik_image_histogram_matching_test1(refimg,maskref,img,maskimg,step)) {
      fprintf(stderr,"ERROR: niik_image_histogram_matching_test1\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = inorm2 */

  else if(!strncmp(argv[1],"inorm3",6)) {
    /* intensity normalization function */
    fprintf(stdout,"[%s] intensity normalization using linear regression\n",fcname);
    if(outname==NULL || img==NULL || refimg==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: inorm3 -in=<img.nii> -ref=<ref.nii> -out=<out.nii> -matrix=<img2ref.matrix>\n");
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"    -mask=<img_mask.nii>        : input mask for <img.nii>\n");
      fprintf(stdout,"\n");
      exit(0);
    }
    if(xfmflag>0 && afmat!=NULL) {
      NIIK_EXIT((!niikmat_convert_from_xfm(afmat,img,refimg)),
                fcname,"niikmat_convert_from_xfm",9);
    }
    NIIK_EXIT((!niik_image_linear_normalization(refimg,img,maskimg,afmat)),fcname,"niik_image_linear_normalization",9);
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
    exit(0);
  } /* OP = inorm3 */

  else if(!strncmp(argv[1],"qcimg",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"[%s] usage: qcimg -img=<img.nii> -out=<outname>\n",fcname);
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -xfm=<stx.xfm>            : stx registration\n");
      fprintf(stdout,"  -mat=<stx.mat>            : stx registration\n");
      fprintf(stdout,"  -ref=<ref.nii>            : reference image in stx space\n");
      fprintf(stdout,"  -bspline -nn -linear      : choose one of interpolation types\n");
      exit(1);
    }

    if(xfmflag>0 && refimg!=NULL) {
      NIIK_EXIT((!niikmat_convert_from_xfm(afmat,img,refimg)),fcname,"niikmat_convert_from_xfm",9);
    }
    NIIK_EXIT((!niik_image_affine_transform_3d_update(img,refimg,afmat,interp)),fcname,"niik_image_affine_transform_3d_update",9);
    tiff_xslices=(int *)calloc(33,sizeof(int));
    tiff_xslices[0]=33;
    for(n=20,i=1; i<tiff_xslices[0]; i++,n+=5) tiff_xslices[i]=n;
    tiff_yslices=(int *)calloc(37,sizeof(int));
    tiff_yslices[0]=37;
    for(n=20,i=1; i<tiff_yslices[0]; i++,n+=5) tiff_yslices[i]=n;
    tiff_zslices=(int *)calloc(33,sizeof(int));
    tiff_zslices[0]=33;
    for(n=20,i=1; i<tiff_zslices[0]; i++,n+=5) tiff_zslices[i]=n;
    for(n=1; n<tiff_xslices[0]; n++) {
      sprintf(fname,"%s_x%i.tif",outname,tiff_xslices[n]);
      fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
      NIIK_EXIT((tiff_xslices[n]<0),fcname,"xslice too small",9);
      NIIK_EXIT((tiff_xslices[n]>=img->nx),fcname,"xslice too large",9);
      NIIK_EXIT((!niik_tiff_write_xslice(fname,img,imin,imax,tiff_xslices[n])),fcname,"niik_tiff_write_xslice",9);
    } /* each specified xslices */
    for(n=1; n<tiff_yslices[0]; n++) {
      sprintf(fname,"%s_y%i.tif",outname,tiff_yslices[n]);
      fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
      NIIK_EXIT((tiff_yslices[n]<0),fcname,"yslice too small",9);
      NIIK_EXIT((tiff_yslices[n]>=img->ny),fcname,"yslice too large",9);
      NIIK_EXIT((!niik_tiff_write_yslice(fname,img,imin,imax,tiff_yslices[n])),fcname,"niik_tiff_write_yslice",9);
    } /* each specified yslices */
    for(n=1; n<tiff_zslices[0]; n++) {
      sprintf(fname,"%s_z%i.tif",outname,tiff_zslices[n]);
      fprintf(stdout,"[%s] %i  writing %s\n",fcname,n,fname);
      NIIK_EXIT((tiff_zslices[n]<0),fcname,"zslice too small",9);
      NIIK_EXIT((tiff_zslices[n]>=img->nx),fcname,"zslice too large",9);
      NIIK_EXIT((!niik_tiff_write_zslice(fname,img,imin,imax,tiff_zslices[n])),fcname,"niik_tiff_write_zslice",9);
    } /* each specified zslices */


  } /* OP = qcimg */

  else if(!strncmp(argv[1],"NBCSR",5) ||
          !strncmp(argv[1],"nbcsr",5)) {
    if(imglist==NULL || outname==NULL || numimglist<3 || numimglist>4) {
      fprintf(stdout,"  usage: %s %s -imglist=<refimg.nii>,<refseg.nii>,<movimg.nii>[,<movseg.nii>] -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -delta=<dist>          sampling distance in mm [default=2.4]\n");
      fprintf(stdout,"  -radius=<R>            radius for brain filling [default=3.2]\n");
      fprintf(stdout,"  -FWHM=<f>              blurring kernel FWHM [default=6.0]\n");
      fprintf(stdout,"  -cc                    cost function is correlation coefficient [default]\n");
      fprintf(stdout,"  -nmi                   cost function is normalized mutual info\n");
      fprintf(stdout,"  -nmi-num=<s1>,<s2>     nmi histogram size [default=32,32 for s1 and s2]\n");
      fprintf(stdout,"  -outcheck=<check.nii>  writes check image\n");
      exit(1);
    }
    /* set default values: FWHM delta areg_cost */
    if(!(numimglist==3 || numimglist==4)) {
      fprintf(stderr,"[niikmath] ERROR: 3 or 4 images needed, %i\n",numimglist);
      exit(1);
    }
    if(niik_image_cmp_dim(imglist[0],imglist[1])>0) {
      fprintf(stderr,"[niikmath] ERROR: ref image and mask have different dimensions\n");
      exit(1);
    }
    if(!niik_image_type_convert(imglist[1],NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_convert_type\n");
      exit(0);
    }
    if(numimglist==4) {
      if(!niik_image_type_convert(imglist[3],NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_convert_type\n");
        exit(0);
      }
    }
    if(niik_image_count_mask(imglist[1])==0) {
      fprintf(stderr,"[niikmath] ERROR: no mask in %s\n",imglist[1]->fname);
      exit(0);
    }
    if(numimglist==4) {
      if(niik_image_cmp_dim(imglist[2],imglist[3])>0) {
        fprintf(stderr,"[niikmath] ERROR: moving image and mask have different dimensions\n");
        exit(1);
      }
      if(niik_image_count_mask(imglist[3])==0) {
        fprintf(stderr,"[niikmath] ERROR: no mask in %s\n",imglist[3]->fname);
        exit(0);
      }
    }
    if(niik_check_double_problem(FWHM)) {
      FWHM=6.0;
    }
    if(niik_check_double_problem(delta)) {
      delta=2.4;
    }
    if(niik_check_double_problem(radius)) {
      radius=3.2;
    }
    if(areg_cost==NIIK_REGISTER_UNKNOWN) {
      areg_cost=NIIK_REGISTER_CC;
    }
    if     (!strncmp(argv[1],"NBCSR-rigid",11)) dof=6;
    else if(!strncmp(argv[1],"nbcsr-rigid",11)) dof=6;
    else if(!strncmp(argv[1],"NBCSR",       5)) dof=12;
    else if(!strncmp(argv[1],"nbcsr",       5)) dof=12;
    fprintf(stdout,"  Non-Brain Constrained Symmetric Registration (NBCSR)\n");
    fprintf(stdout,"  basic FWHM     %8.4f\n",FWHM);
    fprintf(stdout,"  basic delta    %8.4f\n",delta);
    fprintf(stdout,"  radius         %8.4f\n",radius);
    fprintf(stdout,"  DOF            %4i\n",dof);
    fprintf(stdout,"  cost           %s\n",niik_aregister_method_string(areg_cost));
    if(radius>0) {
      fprintf(stdout,"[niikmath]   close brain %5.2f\n",radius);
      if(!niik_image_morph_close_brain(imglist[1],radius,radius)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_morph_close_brain(imglist[1],radius,radius)\n");
        exit(1);
      }
      if(debug) {
        fprintf(stdout,"[niikmath] writing tmp_nbcsr_seg1.nii.gz for debugging\n");
        if(!niik_image_write("tmp_nbcsr_seg1.nii.gz",imglist[1])) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_write(\"tmp_nbcsr_seg1.nii.gz\",imglist[1])\n");
          exit(1);
        }
      }
    }
    fprintf(stdout,"  vol ref        %9.4f\n",niik_image_get_mask_vol(imglist[1]));
    if(numimglist==4) {
      if(radius>0) {
        fprintf(stdout,"[niikmath]   close brain %5.2f\n",radius);
        if(!niik_image_morph_close_brain(imglist[3],radius,radius)) {
          fprintf(stderr,"[niikmath] ERROR: niik_image_morph_close_brain(imglist[3],radius,radius)\n");
          exit(1);
        }
        if(debug) {
          fprintf(stdout,"[niikmath] writing tmp_nbcsr_seg2.nii.gz for debugging\n");
          if(!niik_image_write("tmp_nbcsr_seg2.nii.gz",imglist[3])) {
            fprintf(stderr,"[niikmath] ERROR: niik_image_write(\"tmp_nbcsr_seg2.nii.gz\",imglist[3])\n");
            exit(1);
          }
        }
      }
      if(numimglist==4)
        fprintf(stdout,"  vol mov        %9.4f\n",niik_image_get_mask_vol(imglist[3]));
    }
    /*
     * start NBCSR
     */
    if(debug) nifti1_kunio_nbcsr_turn_on_debug();
    fprintf(stdout,"[niikmath]   start NBCSR\n");
    if(numimglist==4) {
      if(!niik_aregister_nbcsr(imglist[0],imglist[1],imglist[2],imglist[3],affpar,areg_cost,FWHM,delta,(dof==6))) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_nbcsr(imglist[0],imglist[1],imglist[2],imglist[3],affpar,areg_cost,FWHM,delta,(dof==6))\n");
        exit(1);
      }
    } else if(numimglist==3) {
      if(!niik_aregister_nbcsr(imglist[0],imglist[1],imglist[2],NULL,affpar,areg_cost,FWHM,delta,(dof==6))) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_nbcsr\n");
        exit(1);
      }
    } else if(numimglist>4) {
      fprintf(stderr,"[niikmath] ERROR: imglist too large\n");
      exit(1);
    } else if(numimglist<3) {
      fprintf(stderr,"[niikmath] ERROR: imglist too small\n");
      exit(1);
    }
    /*
     * create the affine matrix from parameters
     */
    if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar\n");
      exit(1);
    }
    niikmat_display(afmat);
    /*
     * write output matrix
     */
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
      exit(1);
    }
    /*
     * write output transformed image
     */
    if(outcheckname!=NULL) {
      if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update\n");
        exit(1);
      }
      fprintf(stdout,"    transform image for checking: %s\n",outcheckname);
      if(!niik_image_affine_transform_3d_update(imglist[2],imglist[0],afmat,NIIK_INTERP_BSPLINE)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update\n");
        exit(1);
      }
      /* fprintf(stdout,"      modify imglist\n"); */
      imglist[1]=niik_image_free(imglist[1]);
      if(numimglist==4) {
        if(imglist[3]!=NULL) {
          imglist[3]=niik_image_free(imglist[3]);
        }
      }
      imglist[1]=imglist[2];
      imglist[2]=NULL;
      /* fprintf(stdout,"      create output banded image\n"); */
      if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_combine\n");
        exit(1);
      }
      /* fprintf(stdout,"      convert\n"); */
      if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
      /* fprintf(stdout,"      cal_min cal_max\n"); */
      outimg->cal_min=0;
      outimg->cal_max=200;
      fprintf(stdout,"[niikmath] writing output    %s\n",outcheckname);
      niik_image_append_history(outimg,timestamp);
      niik_image_write(outcheckname,outimg);
      /* clean up the temporary variables */
      outimg=niik_image_free(outimg);
      imglist[0]=niik_image_free(imglist[0]);
      imglist[1]=niik_image_free(imglist[1]);
      free(imglist);
      imglist=NULL;
    }
    exit(0);
  } /* OP = NBCSR */

  else if(!strncmp(argv[1],"sobel",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -sobel filter for all x,y,z-directions\n");
      fprintf(stdout,"  -<out.nii> each direction is updated in 5th dimension\n");
      fprintf(stdout,"  -for sobel filter in each directon or magnitude,\n");
      fprintf(stdout,"   use 'sobelx', 'sobely', 'sobelz', or 'sobelm' \n");
      fprintf(stdout,"\n");
      exit(1);
    }
    fprintf(stdout,"[%s]   sobel filter (4-dimensional with x,y,z)\n",fcname);
    NIIK_EXIT(((imglist=niik_image_sobel_filters(img))==NULL),fcname,"niik_image_sobel_filters",9);
    img = niik_image_free(img);
    NIIK_EXIT(((outimg=niik_image_combine(imglist,3,5,0))==NULL),fcname,"niik_image_combine",9);
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(outimg,timestamp);
    NIIK_EXIT((!niik_image_write(outname,outimg)),fcname,"niik_image_write",9);
    outimg=niik_image_free(outimg);
    for(n=0; n<3; n++) {
      imglist[n] = niik_image_free(imglist[n]);
    }
    free(imglist);
    exit(0);
  } /* OP = sobel */

  else if(!strncmp(argv[1],"clear",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s clear -in=<in.nii> -out=<out.nii>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -clears image <in.matrix> (i.e., all voxels have intensity zero\n");
      fprintf(stdout,"  -writes <out.matrix>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   clear image\n");
    if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_get_voxels_as_double_vector\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      dimg[i]=0;
    }
    if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    free(dimg);
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
  }  /* OP = clear */

  else if(!strncmp(argv[1],"erode",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage 1: %s %s -in=<in.nii> -out=<out.nii> -radius=<R>\n",argv[0],argv[1]);
      fprintf(stdout,"  usage 2: %s %s -in=<in.nii> -out=<out.nii> -val=<kernel_size>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -erodes image <in.matrix> with spherical kernel with radius = <R>\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(!niik_check_double_problem(inval)) {
      fprintf(stdout,"[niikmath]   erosion dim %i\n",(int)inval);
      if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_ERODE,morph_kernel_shape,(int)inval)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_morph_3d_mask\n");
        exit(1);
      }
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outname,img)) {
        fprintf(stderr,"ERROR: writing file\n");
        exit(0);
      }
      img=niik_image_free(img);
      exit(0);
    } /* kernel based erosion */
    if(niik_check_double_problem(radius)) {
      fprintf(stderr,"[niikmath] ERROR: please set radius by -radius=<R>\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   erosion R %8.3f\n",radius);
    if(!niik_image_morph_3d_radius(img,NIIK_MORPH_ERODE,radius)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    if(datatype>=0) {
      NIIK_RET0((!niik_image_type_convert(img,datatype)),fcname,"niik_image_type_convert");
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: writing file\n");
      exit(0);
    }
    img=niik_image_free(img);
    exit(0);
  } /* OP = erode */

  else if(!strncmp(argv[1],"close",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s erode -in=<in.nii> -out=<out.nii> -radius=<R> [-mask=<mask.nii>]\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -closes (dilates and erodes) image <in.nii> with spherical\n");
      fprintf(stdout,"   kernel with radius = <R>\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"  -optionally use '-mask=<mask.nii>'\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    NIIK_RET0((niik_check_double_problem(radius)),fcname,"please set radius by -radius=<R>");
    fprintf(stdout,"[%s]   closing R %8.3f\n",fcname,radius);
    NIIK_RET0((!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_CLOSE,radius)),fcname,"niik_image_morph_3d_radius_mask");
    if(datatype>=0) {
      NIIK_RET0((!niik_image_type_convert(img,datatype)),fcname,"niik_image_type_convert");
    }
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"ERROR: writing file");
    img=niik_image_free(img);
  } /* OP = close */

  else if(!strncmp(argv[1],"gauss",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii> -FWHM=<fwhm>\n",argv[0],argv[1]);
      fprintf(stdout,"  -filters image <in.nii> with gaussian square kernel with FWHM = <fwhm>\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(FWHM)) {
      fprintf(stderr,"[niikmath] ERROR: please set FWHM by -FWHM=<FWHM>\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath]   gaussian filter\n");
    if(!niik_image_filter_gaussian_update(img,NIIK_DMAX(1,floor(FWHM*2.5+0.5)),FWHM)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_filter_gaussian_update\n");
      exit(1);
    }
    if(datatype>=0) {
      if(!niik_image_type_convert(img,datatype)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
        exit(1);
      }
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"niik_image_write");
    img=niik_image_free(img);
  } /* OP = gauss */

  else if(!strncmp(argv[1],"histo",5)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> [-out=<out.txt>]\n",argv[0],argv[1]);
      fprintf(stdout,"  -write to screen or to output file with '-out=<out.txt>'\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -histo-num=<N>      : histogram bin number\n");
      fprintf(stdout,"  -delta=<D>          : histogram bin size\n");
      fprintf(stdout,"  -imin=<min>         : histogram min value\n");
      fprintf(stdout,"  -imax=<max>         : histogram max value\n");
      fprintf(stdout,"  -histo-avg=<a>      : running average size [default=0]\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  see also:\n");
      fprintf(stdout,"  niik_plot_histo     : script to plot the histogram using octave\n");
      exit(0);
    }
    if(niik_check_double_problem(imin)) imin=niik_image_get_min(img,maskimg);
    if(niik_check_double_problem(imax)) imax=niik_image_get_max(img,maskimg);
    if(niik_check_double_problem(delta) && histo_num<0) {
      /* fprintf(stdout,"  estimating both delta and num\n"); */
      dd = niik_image_histogram_optim_bin_size(img,maskimg,2);
      histo_num = (imax-imin)/dd;
    } else if(!niik_check_double_problem(delta)) {
      /* fprintf(stdout,"  estimating num\n"); */
      dd = delta;
      histo_num = (imax-imin)/delta+1;
    } else if(histo_num>0) {
      /* fprintf(stdout,"  estimating delta\n"); */
      dd = delta = (imax-imin)/(histo_num-1);
    }
    histomat=niikmat_init(2,histo_num);
    dd=(imax-imin)/(histo_num-1.0);
    for(n=0; n<histo_num; n++) {
      histomat->m[0][n]=imin + dd*n;
    }
    fprintf(stdout,"[niikmath]   histogram binning: %5.5f %5.5f %5.5f   | %i\n",imin,dd,imax,histo_num);
    if(!niik_image_histogram(img,maskimg,histomat->m[0],histomat->m[1],histo_num)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_histogram(img,maskimg,histo_hx,histo_hy,histo_num)\n");
      exit(0);
    }
    if(histo_anum>0) {
      fprintf(stdout,"[niikmath]   running average %i\n",histo_anum);
      if(!niik_runavg_double_vector(histomat->m[1],histo_num,histo_anum)) {
        fprintf(stderr,"[niikmath] ERROR: niik_runavg_double_vector\n");
        exit(0);
      }
    }
    if(outname!=NULL) {
      if(!niikmat_write(outname,histomat)) {
        fprintf(stderr,"[niikmath] ERROR: niikmat_write\n");
        exit(0);
      }
    } else {
      for(n=0; n<histo_num; n++) {
        fprintf(stdout,"%-15.9f %-15.9f\n",histomat->m[0][n],histomat->m[1][n]);
      }
    }
    histomat=niikmat_free(histomat);
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    exit(0);
  } /* OP = histo */

  else if(!strncmp(argv[1],"inorm",5)) {
    /* intensity normalization function */
    if(outname==NULL || img==NULL || afmat==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<areg-mni.matrix> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -omin=<omin>         set output min value [default=0]\n");
      fprintf(stdout,"  -omax=<omax>         set output max value [default=2500]\n");
      fprintf(stdout,"  -imin=<imin>         set input min value [default=auto]\n");
      fprintf(stdout,"  -imax=<imax>         set input max value [default=auto]\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  -linearly maps imax to omax and imin to omin\n");
      fprintf(stdout,"  -values can be smaller than omin and larger than omax\n");
      fprintf(stdout,"  -to remove negative numbers, run niikmath 'maskthresh'\n");
      exit(0);
    }
    if(niik_check_double_problem(omax)) omax=2500;
    if(niik_check_double_problem(omin)) omin=0;
    if(niik_check_double_problem(imax)) {
      if(niik_check_double_problem(inval)) inval=0.95;
      imax = niik_image_get_percentile(img,NULL,inval);
    }
    if(niik_check_double_problem(imin)) {
      /* read mni images */
      sprintf(fname,"%s/data/CLADA/ .nii.gz",NIIKDIR);
      fprintf(stdout,"[niikmath] reading mni bg     %s\n",fname);
      if((mni_seg=niik_image_read(fname))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
        exit(1);
      }
      if(!niik_image_affine_transform_3d_update(mni_seg,img,afmat,NIIK_INTERP_NN)) {
        fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
        exit(0);
      }
      imin = niik_image_get_median(img,mni_seg);
    }
    if(!niik_image_iscale(img,imin,imax,omin,omax)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_iscale(img,NIIKMAX,NIIKMAX,omin,omax)\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    img=niik_image_free(img);
    exit(0);
  } /* inorm */

  else if(!strncmp(argv[1],"bseg0",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<regmni.matrix> -out=<brainmask.nii.gz>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -radius=<R>        set radius for erosion/dilation [default=3.6]\n");
      fprintf(stdout,"  -thresh=<T>        set threshold [default=auto]\n");
      exit(0);
    }
    if(niik_check_double_problem(radius)) radius=3.6;
    if(iter<0) iter=3;
    if((maskimg=niik_image_bseg_simplest(img,afmat,thresh,iter,radius))==NULL) {
      fprintf(stderr,"ERROR: niik_image_bseg_simplest\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath]   vol %8.5f ml\n",niik_image_get_mask_vol(maskimg));
    fprintf(stdout,"[niikmath] write %s\n",outname);
    niik_image_append_history(maskimg,timestamp);
    if(!niik_image_write(outname,maskimg)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    if(outcheckname!=NULL) {
      omax=140 / niik_image_get_percentile(img,NULL,0.95);
      if(!niik_image_mul_voxels_ROI(img,0,0,0,img->nx-1,img->ny-1,img->nz-1,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mul_voxels_ROI\n");
        exit(0);
      }
      vlist=(double *)calloc(3,sizeof(double));
      /*omax=niik_image_get_max(img,NULL);*/
      vlist[0]=vlist[1]=255;
      vlist[2]=0;
      if(!niik_image_boundary(img,maskimg,vlist,1,1,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
        exit(0);
      }
      if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"ERROR: niik_image_type_convert\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing check   %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outcheckname,img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
        exit(1);
      }
    }
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    afmat=niikmat_free(afmat);
    exit(0);
  } /* OP = bseg0 */

  else if(!strncmp(argv[1],"bseg1",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<regmni.matrix> -out=<brain_area.nii.gz>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -radius=<R>           : radius for opening [default=2.4]\n");
      fprintf(stdout,"  -ithresh=<ithresh>    : initial threshold [default=auto]\n");
      fprintf(stdout,"  -thresh=<thresh>      : threshold [default=auto]\n");
      fprintf(stdout,"  -uthresh=<uthresh>    : upper threshold [default=auto]\n");
      fprintf(stdout,"  -iter=<iter>          : number of iterations [default=10]\n");
      exit(0);
    }
    if(niik_check_double_problem(radius))  radius=2.4;
    if(niik_check_double_problem(ithresh)) ithresh=-1.5;
    if(niik_check_double_problem(thresh))  thresh=-1.8;
    if(niik_check_double_problem(uthresh)) uthresh=-0.99;
    if(iter<0) iter=15;
    if(afmat==NULL) {
      if(!niik_aregister_align_mni_predefined_imgs(img,maskimg,affpar,NIIK_REGISTER_NMI,6.0,3.2)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_align_mni\n");
        exit(1);
      }
      if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar\n");
        exit(1);
      }
    }
    if(obj==NULL) {
      if((obj=niik_image_bseg_get_brain_surface(afmat,0,0))==NULL) {
        fprintf(stderr,"ERROR: niik_image_bseg_get_brain_surface\n");
        exit(0);
      }
    }
    if((maskimg=niik_image_bseg_test1(img,afmat,obj,ithresh,thresh,uthresh,radius,1.5,12,iter,1,6.5,3.0,5.0,1.0))==NULL) {
      fprintf(stderr,"ERROR: niik_image_bseg_test1\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath]   vol %8.5f ml\n",niik_image_get_mask_vol(maskimg));
    fprintf(stdout,"[niikmath] write %s\n",outname);
    niik_image_append_history(maskimg,timestamp);
    if(!niik_image_write(outname,maskimg)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    if(outcheckname!=NULL) {
      omax=140 / niik_image_get_percentile(img,NULL,0.95);
      if(!niik_image_mul_voxels_ROI(img,0,0,0,img->nx-1,img->ny-1,img->nz-1,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mul_voxels_ROI\n");
        exit(0);
      }
      vlist=(double *)calloc(3,sizeof(double));
      vlist[0]=vlist[1]=255;
      vlist[2]=0;
      if(!niik_image_boundary(img,maskimg,vlist,1,1,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
        exit(0);
      }
      if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"ERROR: niik_image_type_convert\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing check   %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outcheckname,img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
        exit(1);
      }
    }
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    afmat=niikmat_free(afmat);
    exit(0);
  } /* OP = bseg1 */

  else if(!strncmp(argv[1],"bseg2",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -matrix=<regmni.matrix> -out=<brain_area.nii.gz>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -radius=<R>           : radius for opening [default=2.4]\n");
      fprintf(stdout,"  -ithresh=<ithresh>    : initial threshold [default=auto]\n");
      fprintf(stdout,"  -thresh=<thresh>      : threshold [default=auto]\n");
      fprintf(stdout,"  -uthresh=<uthresh>    : upper threshold [default=auto]\n");
      fprintf(stdout,"  -iter=<iter>          : number of iterations [default=10]\n");
      exit(0);
    }
    if(niik_check_double_problem(uthresh)) uthresh=-0.99;
    if(niik_check_double_problem(thresh))  thresh=-1.8;
    if(niik_check_double_problem(radius))  radius=5.0;
    if((maskimg=niik_image_bseg_test2(img,afmat,thresh,uthresh,radius,10.0))==NULL) {
      fprintf(stderr,"ERROR: niik_image_bseg_test1\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath]   vol %8.5f ml\n",niik_image_get_mask_vol(maskimg));
    fprintf(stdout,"[niikmath] write %s\n",outname);
    niik_image_append_history(maskimg,timestamp);
    if(!niik_image_write(outname,maskimg)) {
      fprintf(stderr,"ERROR: niik_image_write\n");
      exit(0);
    }
    if(outcheckname!=NULL) {
      omax=140 / niik_image_get_percentile(img,NULL,0.95);
      if(!niik_image_mul_voxels_ROI(img,0,0,0,img->nx-1,img->ny-1,img->nz-1,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mul_voxels_ROI\n");
        exit(0);
      }
      vlist=(double *)calloc(3,sizeof(double));
      vlist[0]=vlist[1]=255;
      vlist[2]=0;
      if(!niik_image_boundary(img,maskimg,vlist,1,1,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
        exit(0);
      }
      if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"ERROR: niik_image_type_convert\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing check   %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outcheckname,img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
        exit(1);
      }
    }
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    afmat=niikmat_free(afmat);
    exit(0);
  } /* OP = bseg2 */

  else if(!strncmp(argv[1],"PABIC",5)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -mask=<mask.nii> -out=<cor.nii>[,<bias.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -out=<cor.nii>,<bias.nii>      : in addition to bias-corrected image <cor.nii>,\n");
      fprintf(stdout,"                                   bias field is also written\n");
      fprintf(stdout,"                                 : use comma-separated entry\n");
      fprintf(stdout,"  -ndegree=<nd>                  : number of degrees (order) [default=2]\n");
      fprintf(stdout,"  -delta=<sample>                : sampling distance [default=2.27]\n");
      exit(0);
    }
    if(ndegree<0) ndegree=2;
    if(niik_check_double_problem(delta)) delta=2.27;
    if((outimg=niik_image_pabic_bias(img,maskimg,ndegree,delta))==NULL) {
      fprintf(stderr,"EROR: niik_image_pabic_bias\n");
      exit(0);
    }
    if(!niik_image_multiply_2_images(img,outimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_multiply_2_images\n");
      exit(0);
    }
    if((CSlist=niik_csstring_to_list(outname,&n))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",outname);
      exit(0);
    }
    switch(n) {
    case 2:
      fprintf(stdout,"[niikmath] write %s\n",CSlist[1]);
      niik_image_append_history(outimg,timestamp);
      if(!niik_image_write(CSlist[1],outimg)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",CSlist[1]);
        exit(0);
      }
    case 1:
      fprintf(stdout,"[niikmath] write %s\n",CSlist[0]);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(CSlist[0],img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",CSlist[0]);
        exit(0);
      }
      break;
    default:
      fprintf(stderr,"[niikmath] ERROR: unknown outname %s\n",outname);
      exit(0);
    }
    for(m=0; m<n; m++) {
      free(CSlist[m]);
    }
    free(CSlist);
    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
    outimg=niik_image_free(outimg);
    exit(0);
  } /* OP = bseg2 */

  else if(!strncmp(argv[1],"NBCR",4)) {
    if( img==NULL || refimg==NULL || outname==NULL || maskref==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<refimg.nii> -refmask=<refmask.nii> -out=<out.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -delta=<dist>          : sampling distance in mm [default=2.4]\n");
      fprintf(stdout,"  -radius=<R>            : radius for brain filling [default=3.2]\n");
      fprintf(stdout,"  -FWHM=<f>              : blurring kernel FWHM [default=6.0]\n");
      fprintf(stdout,"  -cc                    : cost function is correlation coefficient [default]\n");
      fprintf(stdout,"  -nmi                   : cost function is normalized mutual info\n");
      fprintf(stdout,"  -nmi-num=<s1>,<s2>     : nmi histogram size [default=32,32 for s1 and s2]\n");
      fprintf(stdout,"  -outcheck=<check.nii>  : writes check image\n");
      fprintf(stdout,"  -matrix=<MNI.matrix>   : MNI affine registration matrix for <refimg.nii>\n");
      exit(0);
    }
    if(niik_check_double_problem(FWHM))  {
      FWHM=6.0;
    }
    if(niik_check_double_problem(delta)) {
      delta=2.4;
    }
    if(niik_check_double_problem(radius)) {
      radius=3.2;
    }
    if(areg_cost==NIIK_REGISTER_UNKNOWN) {
      areg_cost=NIIK_REGISTER_CC;
    }
    dof=12;
    fprintf(stdout,"  Non-Brain Constrained Registration (NBCR)\n");
    fprintf(stdout,"  basic FWHM     %8.4f\n",FWHM);
    fprintf(stdout,"  basic delta    %8.4f\n",delta);
    fprintf(stdout,"  radius         %8.4f\n",radius);
    fprintf(stdout,"  DOF            %4i\n",dof);
    fprintf(stdout,"  cost           %s\n",niik_aregister_method_string(areg_cost));
    if(radius>0) {
      fprintf(stdout,"[%s]   close brain %5.2f %s\n",fcname,radius,maskref->fname);
      if(!niik_image_morph_close_brain(maskref,radius,radius)) {
        fprintf(stderr,"[%s] ERROR: niik_image_morph_close_brain\n",fcname);
        exit(0);
      }
      if(maskimg!=NULL) {
        fprintf(stdout,"[%s]   close brain %5.2f %s\n",fcname,radius,maskimg->fname);
        if(!niik_image_morph_close_brain(maskimg,radius,radius)) {
          fprintf(stderr,"[%s] ERROR: niik_image_morph_close_brain\n",fcname);
          exit(0);
        }
        fprintf(stdout,"  vol mov        %9.4f\n",niik_image_get_mask_vol(maskimg));
      }
    } /* close brain */
    fprintf(stdout,"  vol ref        %9.4f\n",niik_image_get_mask_vol(maskref));
    if(afmat!=NULL) {
      fprintf(stdout,"[%s] updating sform\n",fcname);
      if(!niik_image_update_sform(refimg,afmat)) {
        fprintf(stderr,"[%s] ERROR: niik_image_update_sform\n",fcname);
        exit(0);
      }
    } /* affine matrix */
    /*
     * start NBCR
     */
    if(debug>=1) nifti1_kunio_nbcsr_turn_on_debug();
    ctr=niikpt_image_get_centroid(refimg,maskref);
    affpar[14]=ctr.x;
    affpar[15]=ctr.y;
    affpar[16]=ctr.z;
    if(!niik_aregister_nbcr(refimg,maskref,img,maskimg,affpar,areg_cost,FWHM,delta,(dof==6))) {
      fprintf(stderr,"[%s] ERROR: niik_aregister_nbcr\n",fcname);
      exit(0);
    }
    /* create the affine matrix from parameters */
    if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_aregister_matrix_from_affpar\n",fcname);
      exit(1);
    }
    niikmat_display(afmat);
    /* write output matrix */
    fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
    if(!niikmat_write(outname,afmat)) {
      fprintf(stderr,"[%s] ERROR: niikmat_write\n",fcname);
      exit(0);
    }
    /* write output transformed image */
    if(outcheckname!=NULL) {
      if(!niik_aregister_matrix_from_affpar_update(afmat,affpar)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar_update\n");
        exit(0);
      }
      fprintf(stdout,"    transform image for checking: %s\n",outcheckname);
      if(!niik_image_affine_transform_3d_update(img,refimg,afmat,NIIK_INTERP_BSPLINE)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update\n");
        exit(1);
      }
      imglist=(nifti_image **)calloc(2,sizeof(nifti_image *));
      imglist[0]=refimg;
      imglist[1]=img;
      /* fprintf(stdout,"      create output banded image\n"); */
      if((outimg=niik_image_combine(imglist,2,4,140))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_combine\n");
        exit(1);
      }
      /* fprintf(stdout,"      convert\n"); */
      if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert\n");
        exit(1);
      }
      /* fprintf(stdout,"      cal_min cal_max\n"); */
      outimg->cal_min=0;
      outimg->cal_max=200;
      fprintf(stdout,"[niikmath] writing output    %s\n",outcheckname);
      niik_image_append_history(outimg,timestamp);
      if(!niik_image_write(outcheckname,outimg)) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outcheckname);
        exit(0);
      }
      /* clean up the temporary variables */
      outimg=niik_image_free(outimg);
      imglist[0]=niik_image_free(imglist[0]);
      imglist[1]=niik_image_free(imglist[1]);
      free(imglist);
      imglist=NULL;
    }
    exit(0);
  } /* OP = NBCR */

  else if(!strncmp(argv[1],"otsu",4)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"  -threshold image with otsu's algorithm\n");
      fprintf(stdout,"  -writes threshold mask <out.matrix>\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(!niik_image_thresh_otsu(img,maskimg,&thresh)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_thresh_otsu\n");
      exit(1);
    }
    if(maskimg!=NULL) {
      nifti_image_free(maskimg);
      maskimg=NULL;
    }
    fprintf(stdout,"[niikmath]   otsu threshold  %8.3f\n",thresh);
    if(!niik_image_threshold(img,thresh)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_threshold\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    img=niik_image_free(img);
  } /* OP = otsu */

  else if(!strncmp(argv[1],"img2txt",7)) {
    fprintf(stdout,"  image to text file\n");
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -img=<img.nii> -out=<out.txt>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"  -mask=<mask.nii>      : writes only within mask\n");
      exit(0);
    } /* usage */
    fprintf(stdout,"[%s] writing output   %s\n",fcname,outname);
    if((fp=fopen(outname,"w"))==NULL) {
      fprintf(stderr,"ERROR: fopen %s\n",outname);
      exit(0);
    }
    if(maskimg!=NULL) {
      bimg=niik_image_get_voxels_as_uint8_vector(maskimg);
      for(i=0; i<img->nvox; i++) {
        if(bimg[i]>0) {
          fprintf(fp,"%-15.9f\n",niik_image_get_voxel(img,i));
        }
      }
      free(bimg);
    } else {
      /*fprintf(stderr,"[niikmath] ERROR: under construction\n");*/
      for(i=0; i<img->nvox; i++) {
        fprintf(fp,"%-15.9f\n",niik_image_get_voxel(img,i));
      }
    }
    fclose(fp);
    exit(0);
  } /* OP = img2txt */

  else if(!strncmp(argv[1],"stats",5)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s stats -in=<in.nii> [-mask=<mask.nii>]\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -statistical info about the image\n");
      fprintf(stdout,"\n");
      exit(0);
    }
    NIIK_RET0((outname!=NULL),fcname,"There's no output");
    if(!niik_check_double_problem(FWHM)) {
      NIIK_EXIT((!niik_image_filter_gaussian_update(img,FWHM*3.5,FWHM)),fcname,"niik_image_filter_gausssian_update",1);
    }
    if(!niik_check_double_problem(thresh)) {
      fprintf(stdout,"[%s] threshold = %.5f\n",fcname,thresh);
      if(maskimg!=NULL) {
        NIIK_EXIT((!niik_image_threshold(maskimg,thresh)),fcname,"niik_image_threshold",1);
      } else {
        NIIK_EXIT(((maskimg=niik_image_threshold_new(img,thresh))==NULL),fcname,"niik_image_threshold",1);
      }
    }
    NIIK_RET0((!niik_image_display_stats(img,maskimg)),fcname,"niik_image_display_stats");
    if(!niik_check_double_problem(imax) && maskimg!=NULL) {
      NIIK_EXIT(((vec=niik_image_get_voxels_as_double_vector_within_mask(img,maskimg))==NULL),fcname,"vec=niik_image_get_voxels_as_double_vector_within_mask",1);
      NIIK_EXIT((!niik_get_trimmed_average_from_double_vector(vec->v,vec->num,imax,&dval)),fcname,"niik_get_trimmed_average_from_double_vector",1);
      fprintf(stdout," trimmed average      %15.9f   %5.3f%%\n",dval,imax*100.0);
      vec=niikvec_free(vec);
    } /* trimmed average */

    /* principal axes */
    vec=niikvec_init(3);
    mat=niikmat_init(3,3);
    niik_image_principal_axes(img,maskimg,mat,vec);
    fprintf(stdout,"[%s]: eigen values:   %12.6f %12.6f %12.6f\n",fcname,vec->v[0],vec->v[1],vec->v[2]);
    for(i=0; i<3; i++) {
      fprintf(stdout,"[%s]: eigen vector %i: %12.6f %12.6f %12.6f\n",fcname,i,mat->m[i][0],mat->m[i][1],mat->m[i][2]);
    }

    img=niik_image_free(img);
    maskimg=niik_image_free(maskimg);
  } /* OP = stats */

  else if(!strncmp(argv[1],"siena",5)) {
    fprintf(stdout,"niik implmenetaito of SIENA (not tested/validated)\n");
    if(imglist==NULL || maskimg==NULL || img==NULL) {
      fprintf(stdout,"  usage: %s %s -imglist=<img1.nii>,<img2.nii> -matrixlist=<reg1.matrix>,<reg2.matrix> -vlist=<ulim>,<uran> -out=<atrophy.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  optional usage:\n");
      fprintf(stdout,"\n");
      fprintf(stdout,"  example 1: ... \n");
      exit(0);
    }
    if(numimglist!=nummatlist) {
      fprintf(stderr,"[niikmath] ERROR: # image list and # matrix list did not match\n");
      exit(0);
    }
    n=niik_image_count_mask(maskimg);
    fprintf(stdout,"  mask count %i\n",n);
    atv=niikmat_init(numimglist,n);
    if(niik_check_double_problem(dmean) || niik_check_double_problem(dstdv)) {
      if(brainmask==NULL) {
        fprintf(stderr,"[niikmath] ERROR: missing brain mask for calculation of white matter mean and stdv\n");
        exit(0);
      }
      if(niik_check_double_problem(imax)) imax=niik_image_get_max(img,brainmask);
      if(histo_num<0) histo_num=imax/5;
      if(!niik_image_bimodal_fit(img,brainmask,0,10,imax,histo_num,bimodal_fit+2,bimodal_fit+4,bimodal_fit,&dd)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_bimodal_fit\n");
        exit(0);
      }
      dmean=bimodal_fit[3];
      dstdv=bimodal_fit[5];
      fprintf(stdout,"[niikmath]: mean/stdv for WM = %9.3f %9.3f\n",dmean,dstdv);
    }
    if(!niik_image_siena(img,maskimg,imglist,matlist,numimglist,dmean,dstdv,1,atv)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_sinea\n");
      exit(0);
    }
    for(n=1; n<numimglist; n++) {
      niik_display_stats_for_double_vector(atv->m[n],atv->col);
    }
    exit(0);
  } /* OP siena */

  else if(!strncmp(argv[1],"open",4)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s open -in=<in.nii> -out=<out.nii> -radius=<R>\n",argv[0]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -opens (erodes and dilates) image <in.nii> with spherical\n");
      fprintf(stdout,"   kernel with radius = <R>\n");
      fprintf(stdout,"  -writes <out.nii>\n");
      fprintf(stdout,"  -optionally use '-mask=<mask.nii>'\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(niik_check_double_problem(radius)) {
      fprintf(stderr,"[niikmath] ERROR: please set radius by -radius=<R>\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   opening R %8.3f\n",radius);
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_OPEN,radius)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_set_voxels_from_double_vector\n");
      exit(1);
    }
    if(datatype>=0) {
      NIIK_RET0((!niik_image_type_convert(img,datatype)),fcname,"niik_image_type_convert");
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_RET0((!niik_image_write(outname,img)),fcname,"ERROR: writing file");
    img=niik_image_free(img);
    img=NULL;
  } /* OP = open */

  else if(!strncmp(argv[1],"info",4)) {
    if(img==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<in.nii> [-mask=<mask.nii>]\n",argv[0],argv[1]);
      fprintf(stdout,"\n");
      fprintf(stdout,"  -statistical info about the image\n");
      fprintf(stdout,"\n");
      exit(1);
    }
    if(outname!=NULL) {
      fprintf(stderr,"[niikmath] ERROR: there's no output\n");
      exit(1);
    }
    /*fprintf(stdout,"  bin size = %9.5f\n",niik_image_histogram_optim_bin_size(img,maskimg,2));
      exit(0);*/
    if(!niik_image_display_stats(img,maskimg)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_display_stats\n");
      exit(1);
    }
    img=niik_image_free(img);
    nifti_image_free(maskimg);
    img=maskimg=NULL;
  } /* OP = info */

  else if(!strncmp(argv[1],"kmul",4)) {
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<K>\n",argv[0],argv[1]);
      exit(1);
    }
    fprintf(stdout,"[niikmath]   val = %12.9f\n",inval);
    for(n=0; n<img->nvox; n++) {
      niik_image_mul_voxel(img,n,inval);
    }
    if(datatype<0) datatype=img->datatype;
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(1);
    exit(0);
  } /* OP = kmul */

  else if(!strncmp(argv[1],"ksub",4)) {
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<K>\n",argv[0],argv[1]);
      exit(1);
    }
    if(datatype<0) datatype=img->datatype;
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   val = %12.9f\n",inval);
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,inval-niik_image_get_voxel(img,n));
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(1);
    exit(0);
  } /* OP = ksub */

  else if(!strncmp(argv[1],"kadd",4)) {
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<K>\n",argv[0],argv[1]);
      exit(1);
    }
    if(datatype<0) datatype=img->datatype;
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   val = %12.9f\n",inval);
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,inval+niik_image_get_voxel(img,n));
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(1);
    exit(0);
  } /* OP = kadd */

  else if(!strncmp(argv[1],"kdiv",4)) {
    if(img==NULL || outname==NULL || niik_check_double_problem(inval)) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<out.nii> -val=<K>\n",argv[0],argv[1]);
      exit(1);
    }
    if(datatype<0) datatype=img->datatype;
    if(!niik_image_type_convert(img,datatype)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_type_convert(img,datatype)\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath]   val = %12.9f\n",inval);
    for(n=0; n<img->nvox; n++) {
      niik_image_set_voxel(img,n,inval/niik_image_get_voxel(img,n));
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
      exit(1);
    }
    exit(1);
    exit(0);
  } /* OP = kdiv */

  else if(!strncmp(argv[1],"vseg",4)) {
    if(img==NULL || outname==NULL || warpimg==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -mask=<brain_mask.nii> -warp-loc-map=<warp.nii> -out=<ventricle_mask.nii.gz> [options]\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -iter=<i>            change number of dilation iteration [default=auto]\n");
      fprintf(stdout,"  -thresh=<T>          manually set threshold [default=auto]\n");
      exit(1);
    }
    /* checking inputs and options */
    if(warp_map_type!=NIIK_WARP_MAP_LOC)  {
      fprintf(stderr,"[niikmath] ERROR: please use -warp-loc-map=<warp.nii>\n");
      exit(1);
    }
    /* resample the brain mask to low resolution */
    if((maskroi=niik_image_resample_3d(maskimg,3,3,3,-1,-1,-1,NIIK_INTERP_NN))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_resample_3d(maskimg,3,3,3-1,-1,-1,NIIK_INTERP_NN)\n");
      exit(1);
    }
    /*   -dilate the brain mask few times */
    for(i=0; i<3; i++) {
      if(!niik_image_morph_3d_radius(maskroi,NIIK_MORPH_DILATE,3.6)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_morph_3d_radius(maskroi,NIIK_MORPH_DILATE,3.2)\n");
        exit(1);
      }
    }
    /*   -close holes for the dilated brain mask */
    if(!niik_image_close_holes(maskroi)) {
      fprintf(stderr,"[niikmath] ERROR: niik_close_holes\n");
      return 0;
    }
    /*   -erode the closed brain mask */
    for(i=0; i<2; i++) {
      if(!niik_image_morph_3d_radius(maskroi,NIIK_MORPH_ERODE,3.2)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_morph_3d_radius(maskroi,NIIK_MORPH_DILATE,3.2)\n");
        exit(1);
      }
    }
    /* -resample back to the original resolution */
    if(!niik_image_affine_transform_3d_update(maskroi,maskimg,NULL,NIIK_INTERP_NN)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_affine_transform_3d_update(maskroi,maskimg,NULL,NIIK_INTERP_NN)\n");
      exit(1);
    }
    if(debug) { /* write intermediate image */
      fprintf(stdout,"[niikmath] writing tmp_niikmath_vseg_segroi.nii.gz");
      niik_image_write("tmp_niikmath_vseg_segroi.nii.gz",maskroi);
    }
    /* read and warp mni ventricle masks */
    for(n=1; n<=3; n++) {
      if((mni_ven_list[n]=niik_image_segment_ventricle_prep_from_mni_warp(img,warpimg,n))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_segment_ventricle_prep_from_mni_warp(img,warpimg,%i)\n",n);
        exit(1);
      }
      if(!niik_image_mask(mni_ven_list[n],maskroi)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mask\n");
        exit(1);
      }
    }
    if(debug) {
      niik_image_write("tmp_niikmath_vseg_mask.nii.gz",mni_ven_list[1]);
      niik_image_write("tmp_niikmath_vseg_lothresh_roi.nii.gz",mni_ven_list[2]);
      niik_image_write("tmp_niikmath_vseg_dilate_roi.nii.gz",mni_ven_list[3]);
    }
    /* default iteration is 30 */
    if(iter<0) {
      dvol=niik_image_get_mask_vol(mni_ven_list[1]);
      fprintf(stdout,"[niikmath]     vol %8.5f ml\n",dvol);
      iter=30+dvol/img->dx/2;
    }
    if(niik_check_double_problem(thresh)) {
      NIIK_EXIT((niik_check_fsldir_exists()>0),fcname,"fsldir is not defined",1);
      sprintf(fname,"%s/data/standard/MNI152_T1_2mm_brain_mask.nii.gz",FSLDIR);
      fprintf(stdout,"[niikmath] reading mni seg   %s\n",fname);
      if((mni_seg=niik_image_read(fname))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
        exit(1);
      }
      if(!niik_image_apply_3d_warp_update(mni_seg,img,warpimg,NIIK_WARP_MAP_LOC,NIIK_INTERP_NN)) {
        fprintf(stderr,"ERROR: niik_image_apply_3d_warp_update\n");
        exit(0);
      }
      thresh=2.2;
      /* use bseg's threshold calculation to estimate brain intensity threshold */
      if(!niik_image_bseg_basic_thresh_with_mask(img,mni_seg,&thresh,0)) {
        fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh_with_mask\n");
        exit(0);
      }
      mni_seg=niik_image_free(mni_seg);
    }
    if(!niik_image_segment_ventricle(img,maskimg,mni_ven_list[1],mni_ven_list[2],mni_ven_list[3],thresh,iter)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_segment_ventricle(img,maskimg,mni_ven_list[1],mni_ven_list[2],mni_ven_list[3],thresh,iter)\n");
      exit(1);
    }
    /* 2012-08-22 Kunio
     * -added to automatically remove with brain roi */
    fprintf(stdout,"[niikmath] masking out with brain roi\n");
    if(!niik_image_mask(mni_ven_list[1],maskroi)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_mask\n");
      exit(0);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(mni_ven_list[1],timestamp);
    if(!niik_image_write(outname,mni_ven_list[1])) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,mni_ven_list[1])\n");
      exit(1);
    }
    for(n=1; n<=3; n++) {
      mni_ven_list[n]=niik_image_free(mni_ven_list[n]);
    }
    exit(1);
    exit(0);
  } /* OP = vseg */

  else if(!strncmp(argv[1],"none",4)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: none -in=<img> -out=<out>\n");
      exit(1);
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    exit(0);
  } /* OP = none */

  else if(!strncmp(argv[1],"bseg",4)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<brain_mask.nii>,<brain_mask_roi.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -radius=<R>        set radius for erosion/dilation [default=3.6]\n");
      fprintf(stdout,"  -thresh=<T>        set threshold [default=auto]\n");
      fprintf(stdout,"  -uthresh=<UT>      set upper threshold (global) [default=auto]\n");
      fprintf(stdout,"  -uthresh2=<UT2>    set upper threshold (regional-near brain edge) [default=auto]\n");
      fprintf(stdout,"  -ithresh=<IT>      set initial threshold [default=auto]\n");
      fprintf(stdout,"  -iter=<N>          set number of dilation for MNI mask [default=2]\n");
      exit(1);
    }
    if(niik_check_double_problem(radius)) {  /* default radius */
      dvol=niik_image_get_voxel_size(img);
      if(dvol>2) {
        radius=4.5;
      } else {
        radius = 3.6;
      }
    } /* default radius */
    if(iter<0) {
      iter=2;
    }
    if(!niik_image_type_convert_scl(img,NIFTI_TYPE_FLOAT32,1)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert(img,NIFTI_TYPE_FLOAT32)\n",fcname);
      exit(1);
    }
    if((maskimg=niik_image_bseg_basic(img,iter,thresh,ithresh,uthresh,uthresh2,radius,afmat,0))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_bseg_basic\n");
      exit(1);
    }
    if(niik_image_count_mask(maskimg)==0) {
      fprintf(stderr,"ERROR: no mask, zero-volume\n");
      exit(0);
    }
    if((outnamelist=niik_csstring_to_list(outname,&n))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",argv[nc]);
      exit(0);
    }
    switch(n) {
    case 1:
      fprintf(stdout,"[niikmath] writing output    %s\n",outname);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outname,maskimg)) {
        fprintf(stderr,"ERROR: niik_image_write %s\n",outname);
        exit(0);
      }
      break;
    case 2:
      fprintf(stdout,"[niikmath] writing output    %s\n",outnamelist[0]);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outnamelist[0],maskimg)) {
        fprintf(stderr,"ERROR: niik_image_write %s\n",outnamelist[0]);
        exit(0);
      }
      free(outnamelist[0]);
      fprintf(stdout,"[niikmath] close brain-5.5mm\n");
      if(!niik_image_morph_close_brain(maskimg,5.5,5.5)) {
        fprintf(stderr,"ERROR: niik_image_morph_close_brain\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing output    %s\n",outnamelist[1]);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outnamelist[1],maskimg)) {
        fprintf(stderr,"ERROR: niik_image_write %s\n",outnamelist[1]);
        exit(0);
      }
      free(outnamelist[1]);
      break;
    }
    if(outcheckname!=NULL) {
      omax=140 / niik_image_get_percentile(img,NULL,0.95);
      if(!niik_image_mul_voxels_ROI(img,0,0,0,img->nx-1,img->ny-1,img->nz-1,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mul_voxels_ROI\n");
        exit(0);
      }
      vlist=(double *)calloc(3,sizeof(double));
      /*omax=niik_image_get_max(img,NULL);*/
      vlist[0]=vlist[1]=255;
      vlist[2]=0;
      if(!niik_image_boundary(img,maskimg,vlist,1,1,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
        exit(0);
      }
      if(!niik_image_type_convert(img,NIFTI_TYPE_UINT8)) {
        fprintf(stderr,"ERROR: niik_image_type_convert\n");
        exit(0);
      }
      fprintf(stdout,"[niikmath] writing check     %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outcheckname,img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
        exit(1);
      }
    }
    maskimg=niik_image_free(maskimg);
    img=niik_image_free(img);
  } /* OP = bseg */

  else if(!strncmp(argv[1],"BPF",3)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -out=<brainmask.nii.gz> -matrix=<mni.matrix>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -matrix=<mni.matrix>  : affine MNI registration matrix\n");
      fprintf(stdout,"  -radius=<R>           : radius for erosion/dilation [default=3.6]\n");
      fprintf(stdout,"  -thresh=<T>           : threshold [default=auto]\n");
      fprintf(stdout,"  -uthresh=<UT>         : upper threshold (global) [default=auto]\n");
      fprintf(stdout,"  -uthresh2=<UT2>       : upper threshold (regional-near brain edge) [default=auto]\n");
      fprintf(stdout,"  -ithresh=<IT>         : initial threshold [default=auto]\n");
      fprintf(stdout,"  -iter=<N>             : number of dilation for MNI mask [default=2]\n");
      fprintf(stdout,"  -elen=<E>             : edge length for outer contour surface [default=3.5]\n");
      exit(1);
    }
    /*********************************
     * affine MNI matrix
     **********************************/
    if(afmat==NULL) {
      /* read mni images */
      NIIK_EXIT((niik_check_fsldir_exists()>0),fcname,"fsldir is not defined",1);
      sprintf(fname,"%s/data/standard/MNI152_T1_1mm.nii.gz",FSLDIR);
      fprintf(stdout,"[niikmath] reading mni img   %s\n",fname);
      if((mni_img=niik_image_read(fname))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
        exit(1);
      }
      sprintf(fname,"%s/data/standard/MNI152_T1_1mm_brain_mask.nii.gz",FSLDIR);
      fprintf(stdout,"[niikmath] reading mni seg   %s\n",fname);
      if((mni_seg=niik_image_read(fname))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: nifti_image_read %s\n",fname);
        exit(1);
      }
      /* use default values */
      if(niik_check_double_problem(FWHM)) {
        FWHM=6.0;
      }
      if(niik_check_double_problem(delta)) {
        delta=3.2;
      }
      if(!niik_aregister_align_mni(mni_img,mni_seg,img,maskimg,affpar,NIIK_REGISTER_NMI,FWHM,delta)) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_align_mni\n");
        exit(0);
      }
      if((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL) {
        fprintf(stderr,"[niikmath] ERROR: niik_aregister_matrix_from_affpar\n");
        exit(1);
      }
    } /* affine matrix */
    if(niik_check_double_problem(elen)) {  /* default edgelength */
      elen=3.5;
    }
    if(niik_check_double_problem(radius)) {  /* default radius */
      dvol=niik_image_get_voxel_size(img);
      if(dvol>2) {
        radius=4.5;
      } else {
        radius=3.6;
      }
    } /* default radius */
    if(iter<0) {
      iter=2;
    }


    /* BRAIN SEGMENTATION */
    if((brainmask=niik_image_bseg_basic(img,iter,thresh,ithresh,uthresh,uthresh2,radius,afmat,0))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_bseg_basic\n",fcname);
      exit(1);
    }

    if(niik_check_double_problem(imin)) imin=0;
    if(niik_check_double_problem(imax)) imax=-1;
    if(histo_num<0) {
      if(niik_image_get_voxel_size(img)<2.0)
        histo_num=500;
      else
        histo_num=150;
    }
    histomat=niikmat_init(5,5);
    if(!niik_image_fit_gaussian_tissues(img,brainmask,imin,imax,histo_num,0,histomat->m[0],histomat->m[1],histomat->m[2],&dval)) {
      fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian_tissues\n",fcname);
      exit(0);
    }

    exit(0);

    /* BRAIN SEGMENTATION */
    if((brainmask=niik_image_bseg_basic(img,iter,thresh,ithresh,uthresh,uthresh2,radius,afmat,1))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_bseg_basic\n",fcname);
      exit(1);
    }
    if(niik_image_count_mask(brainmask)==0) {
      fprintf(stderr,"[%s] ERROR: no mask, zero-volume\n",fcname);
      exit(0);
    }
    /* BRAIN ROI FOR OUTER CONTOUR */
    fprintf(stdout,"[%s] close brain-5.5mm\n",fcname);
    if((maskimg=niik_image_threshold_new(brainmask,127))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_morph_close_brain\n",fcname);
      exit(0);
    }
    if(!niik_image_morph_close_brain(maskimg,5.5,5.5)) {
      fprintf(stderr,"ERROR: niik_image_morph_close_brain\n");
      exit(0);
    }
    /* surface */
    if(obj==NULL && afmat!=NULL) {
      if(((obj=niik_image_bseg_get_brain_surface(afmat,0,0)))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_bseg_get_brain_surface\n",fcname);
        exit(0);
      }
    } /* get object if not found */
    if(!niik_off_outer_contour_object(maskimg,obj,elen)) {
      fprintf(stderr,"[%s] ERROR: niik_off_outer_contour_object\n",fcname);
      exit(0);
    }
    dvol = off_get_kobj_volume(obj)/1000.0;
    fprintf(stdout,"[%s] outer contour surface volume: %12.5f\n",fcname,dvol);
    fprintf(stdout,"[%s] BPF: %12.5f\n",fcname,niik_image_get_voxel_size(brainmask)*niik_image_get_sum(brainmask,brainmask)*0.001/255.0 / dvol);
    /*
     * WRITING OUTPUT
     */
    if((outnamelist=niik_csstring_to_list(outname,&n))==NULL) {
      fprintf(stderr,"[niikmath] ERROR: niik_csstring_to_list %s\n",argv[nc]);
      exit(0);
    }
    switch(n) {
    default:
    case 3:
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outnamelist[2]);
      if(!off_kobj_write_offply(outnamelist[2],obj,0)) {
        fprintf(stderr,"[%s] ERROR: off_kobj_write_off\n",fcname);
        exit(1);
      }
      free(outnamelist[2]);
    case 2:
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outnamelist[1]);
      niik_image_append_history(maskimg,timestamp);
      if(!niik_image_write(outnamelist[1],maskimg)) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outnamelist[1]);
        exit(1);
      }
      free(outnamelist[1]);
    case 1:
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outnamelist[0]);
      niik_image_append_history(brainmask,timestamp);
      if(!niik_image_write(outnamelist[0],brainmask)) {
        fprintf(stderr,"[%s] ERROR: niik_image_write %s\n",fcname,outnamelist[0]);
        exit(1);
      }
      free(outnamelist[0]);
      break;
    }
    if(outcheckname!=NULL) {
      omax=140 / niik_image_get_percentile(img,NULL,0.95);
      if(!niik_image_mul_voxels_ROI(img,0,0,0,img->nx-1,img->ny-1,img->nz-1,omax)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_mul_voxels_ROI\n");
        exit(1);
      }
      vlist=(double *)calloc(3,sizeof(double));
      /*omax=niik_image_get_max(img,NULL);*/
      vlist[0]=vlist[1]=255;
      vlist[2]=0;
      if(!niik_image_boundary(img,maskimg,vlist,1,1,0)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_boundary\n");
        exit(1);
      }
      NIIK_EXIT((!niik_image_type_convert(img,NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert",9);
      fprintf(stdout,"[niikmath] writing check     %s\n",outcheckname);
      niik_image_append_history(img,timestamp);
      if(!niik_image_write(outcheckname,img)) {
        fprintf(stderr,"[niikmath] ERROR: niik_image_write(outname,img)\n");
        exit(1);
      }
    }
    maskimg=niik_image_free(maskimg);
    img=niik_image_free(img);
  } /* OP = bseg */

  else if(!strncmp(argv[1],"log",3)) {
    fprintf(stdout,"[niikmath natural log\n");
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -val=<v>         : add this value to prevent log(zero)\n");
      exit(1);
    }
    for(i=0; i<img->nvox; i++) {
      niik_image_set_voxel(img,i,log(niik_image_get_voxel(img,i)));
    }
    fprintf(stdout,"[niikmath] writing output    %s\n",outname);
    niik_image_append_history(img,timestamp);
    if(!niik_image_write(outname,img)) {
      fprintf(stderr,"[niikmath] ERROR: niik_image_write %s\n",outname);
      exit(1);
    }
    exit(0);
  } /* OP = log */

  else if(!strncmp(argv[1],"dbc",3) ||
          !strncmp(argv[1],"DBC",3)) {
    if(img==NULL || outname==NULL) {
      fprintf(stdout,"  usage: %s %s -in=<img.nii> -ref=<ref.nii> -out=<out.nii>\n",argv[0],argv[1]);
      fprintf(stdout,"\n  optional usage:\n");
      fprintf(stdout,"  -matrix=<img2ref.matrix>     affine registration matrix [default=none]\n");
      fprintf(stdout,"  -radius=<R>                  radius [default=7.5]\n");
      exit(1);
    }
    fprintf(stdout,"[%s] differential bias correction\n",fcname);
    if(niik_check_double_problem(inval))  inval=3;
    if(niik_check_double_problem(radius)) radius=7.5;
    if(maskimg!=NULL) {
      if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
        fprintf(stdout,"[%s] converting mask to uint8\n",fcname);
        NIIK_EXIT((!niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert",9);
      }
    }
    NIIK_EXIT(((outimg=niik_image_dbc_with_scaling(refimg,img,maskimg,afmat,radius,inval,0))==NULL),fcname,
              "niik_image_dbc_with_scaling",9);
    /*
     *  output may be
     * -out=<out.nii> or -out=<out.nii>,<bias.nii>
     */
    NIIK_EXIT(((outnamelist=niik_csstring_to_list(outname,&n))==NULL),fcname,"niik_csstring_to_list",9);
    if(n==1) {
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outname);
      niik_image_append_history(img,timestamp);
      NIIK_EXIT((!niik_image_write(outname,img)),fcname,"niik_image_write",9);
    } else if(n==2) {
      fprintf(stdout,"[%s] writing output    %s\n",fcname,outnamelist[0]);
      niik_image_append_history(img,timestamp);
      NIIK_EXIT((!niik_image_write(outnamelist[0],img)),fcname,"niik_image_write",9);
      fprintf(stdout,"[%s] writing bias      %s\n",fcname,outnamelist[1]);
      niik_image_append_history(outimg,timestamp);
      NIIK_EXIT((!niik_image_write(outnamelist[1],outimg)),fcname,"niik_image_write",9);
    }
    img=niik_image_free(img);
    outimg=niik_image_free(outimg);
    refimg=niik_image_free(refimg);
    exit(0);
  } /* OP = dbc */

  else {
    fprintf(stderr,"[%s] ERROR: unknown operation, %s\n",fcname,argv[1]);
    exit(1);
  }
  free(timestamp);
  exit(0);
} /* main */
/*
 kate: space-indent on; indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
 */
