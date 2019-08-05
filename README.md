# FALCON cortical surface extraction tool
![Example Thickness map](https://github.com/vfonov/falcon/blob/gh-pages/examples/OASIS-TRT-20-18_thickness-s66_lt.png?raw=true)

## Preliminary results
* [Google drive slides from the meeting on 2019-01-28](https://docs.google.com/presentation/d/1z2YPbEM3VdAtisJ61f-7iaqEFusUhtN4l2r1_5jiWT0/edit?usp=sharing)

## Installation 
### Using conda:
```conda install -c vfonov falcon``` - will work on linux
### Compiling from the source:
Prerequisites: minc-toolkit-v2 (http://bic-mni.github.io/), libtiff, povray, imagemagick
```
cmake <source_dir> \
 -DCMAKE_INSTALL_PREFIX:PATH=<install prefix> \
 -DMINC_TOOLKIT_DIR:PATH=/opt/minc/1.9.17 \
 -DCMAKE_BUILD_TYPE:STRING=Release \
 -DHAVE_POVRAY:BOOL=ON \
 -DUSE_OPENMP:BOOL=ON

make && make install
```

## Included higlevel scripts
* falcon_run.sh - execute FALCON surface extraction pipeline 
```
falcon_run.sh <input.mnc> <output_base> -brain <brain_mask.mnc>

  parameters:
  -help                          :  show this usage

  ---- Required parameters ---
  -brain <brain mask.mnc>        :  brain mask

  --- Recomended parameters ---
  -nl <nl.xfm>                   :  nonlinear registration to icbm [default = None, will run ANTs]
  -omp <omp>                     :  change number of processors [default=1]
  -use_icbm                      :  use mesh from ICBM for initialization , default - use shrink-wrap
  -anlm                          :  apply anlm filter to input t1w scan (if not done before)
  -postprocess                   :  apply post-procesiing: resample mesh to atlas, calculate thickness

  --- Optional parameters  ---
  -vent <vent.mnc>               :  ventricle mask [default = None]
  -cerebellum <cerebellum.mnc>   :  cerebellum mask [default = None]
  -brainstem <brainstem.mnc>     :  brainstem mask [default = None]
  -cls <cls.mnc>                 :  tissue classification map
  -priors <WM> <GM> <CSF>        :  tissue priors, ( default: none )
  -nopriors                      :  don't use tissue priors, for surface fitting

  --- Optional parameters (don't touch, if you don't know what you are doing)   ---
  -gwimask <mask.mnc>            :  gray matter - white matter interface mask [default = None]
  -cerebellummask <mask.mnc>     :  cerebellum mask [default = None]
  -wmmask <mask.mnc>             :  white matter mask [default = None]
  -csfmask <csf.mnc>             :  CSF mask [default = None]
  -left <left.mnc>               :  left hemisphere mask [default = None]
  -right <right.mnc>             :  right hemisphere mask [default = None]
  -sides <left.mnc> <right.mnc>  :  left and right hemisphere masks [default = None]
  -variant  <var>                :  for debugging
  -trace                         :  produce traces of intermediate surfaces
  -hr                            :  produce high-resolution thickness measurement in ICBM space
  -smooth <fwhm>                 :  smooth thickness maps (in addition)
  -noremesh                      :  don't remesh initial mesh (for DEBUGGING only)
  -debug                         :  run in debug mode (keep temp, echo commands)
  -verbose                       :  echo commands
```
* falcon_off_qc_2.sh - generate surface QC image for two hemispheres
```  
falcon_off_qc_2.sh <input_lt.ply> <input_measurement_lt.csv> <input_rt.ply> <input_measurement_rt.csv> <output.png>   [foreground] [background] 
  --- Optional parameters  ---
  -min <m>  - minimal value for the measurement (default: minimum)
  -max <m>  - maximal value for the measurement (default: maximum)
  -title <title> - plot title
  -spectral - use spectral colour map (default)
  -atrophy  - use atrophy colour map
  -summer   - use summer colour map
  -jacobian - use jacobian colour map
  -gray     - use gray colour map
  -sphere   - output on the spherical map, using spherical coordinates instead of x,y,z
  -column  <n> - specify column name from csv file to use (default will use the first column)
```
* falcon_off_qc.sh   - generrate surface QC image for one hemisphere 
```
falcon_off_qc.sh <input.off> <input_measurement.txt> <output.png> [foreground] [background]
  --- Optional parameters  ---
  -min <m>  - minimal value for the measurement (default: minimum)
  -max <m>  - maximal value for the measurement (default: maximum)
  -title <title> - plot title
  -spectral - use spectral colour map (default)
  -atrophy  - use atrophy colour map
  -summer   - use summer colour map
  -jacobian - use jacobian colour map
  -gray     - use gray colour map
  -sphere   - output on the spherical map, using spherical coordinates instead of x,y,z
```
* falcon_slice_qc.sh - generate QC image showing intersection of surface and volume
```
  falcon_slice_qc.sh <input.mnc> <ics.ply> <ocs.ply> <output.jpg/png/tiff>
  --- Optional parameters  ---
  --min <f> - set image min (default: minimum)
  --max <f> - set image max (default: maximum)
  --pct <f> - set image max using percentile
```


## Included low level tools
* falcon_math - mathematical operations
    * primitive operations
        * typeconvert  -converts the image data type
        * info         -displays image stats
        * infodiff     -displays difference image stats
        * infoabsdiff  -displays absolute difference image stats
        * clear        -clears an image
        * setone       -sets one for an image
        * setval       -sets voxel value
        * calcvol      -calculates volumes
        * labelvol     -calculates volumes of label file
    * voxel-wise comparisons
        * voxel-gt     -greater than comparison
        * voxel-ge     -greater than or equal to comparison
        * voxel-lt     -less than comparison
        * voxel-le     -less than or equal to comparison
        * voxel-eq     -equal to comparison
        * voxel-ne     -not equal to comparison
    * arithmetic operations
        * maskimg       -mask image with a mask
        * maskout       -mask-out image with a mask
        * resample      -resamples image with new voxel spacing or image size
        * avgimgs       -average multiple images (voxel-by-voxel)
        * mulimgs       -multiply multiple images (voxel-by-voxel)
        * maximgs       -maximize multiple images (voxel-by-voxel)
        * addimgs       -add multiple images (voxel-by-voxel)
        * subimg        -subtract two images (voxel-by-voxel)
        * absdiffimg    -absolute difference of two images (voxel-by-voxel)
        * divimg        -subtract two images (voxel-by-voxel)
        * ksub          -subtract image voxel value from a constant
        * heavisideF    -heaviside function
        * map-spectral  -put spectral colors to an image
        * map-jacobian  -put jacobian colors to an image
        * log           -natural log
        * avgimgs-fov   -average multiple images with field-of-view correction
        * kmul          -multiply a constant (-val=<K>)
        * iscale        -intensity scaling
        * gaussPDF      -gaussian normal PDF
        * histo         -creates histogram
        * binimg        -binarize image
        * gaussnoise    -image of Gaussian random
        * fill-lesion   -fill lesions as NAWM mean/stdv
        * tiffimg       -creates tiff of a slice
        * bimodalfit    -fits the intensity histogram to bimodal Gaussian distributions
        * brainhistofit -fits intensity histogram to brain tissues (GM,WM,CSF)
        * pseudoT2     -calculates T2 map from dual echo images
    * changing image size
        * montage      -image montage
        * cropimg      -crop image
        * cropmask     -crop image according to mask
        * padimg       -pads image with zeros
    * threshold operations
        * thresh       -threshold image
        * otsu         -threshold image using Otsu's algorithm
        * ridler       -threshold image using Ridler's algorithm
        * maskthresh   -threshold image and mask image (output is an image)
    * segmentation operations
        * bseg         -brain segmentation
        * vseg         -ventricle segmentation
        * segcolor     -write colored segmentation results
    * intensity normalization (and related functions)
        * iscale       -intensity scaling
        * inorm        -intensity normalization
        * inorm2       -intensity normalization using histograms
                    -registration is not required
        * inorm3       -intensity normalization using linear regression
                    -registration is required
        * dbc          -differential bias correction
                    [Lewis and Fox 2004 NeuroImage]
    *  boundary images
        * bounds       -creates a boundary image
        * bounds-inc   -creates a boundary image (boundary within mask)
        * bounds-color -creates a boundary image with colors
    *  create maps
        *   laplacemap   -creates laplace map
        *   distancemap  -creates signed distance map
        *   xdirimg      -x-directional map (voxel intensity for voxel position)
        *   ydirimg      -y-directional map (voxel intensity for voxel position)
        *   zdirimg      -z-directional map (voxel intensity for voxel position)
        *   adirimg      -directional maps (voxel intensity for voxel position, 3D)
    *  image permutation
        * permute-norm         -converts to more or less RAS directions
        * permute-left         -turns left  [x'=-y; y'=x;  z'=z]
        * permute-right        -turns right [x'=y;  y'=-x; z'=z]
        * permute-cor2sag      -rotates     [x'=z;  y'=y;  z'=x]
        * permute-cor2axi      -rotates     [x'=x;  y'=z;  z'=y]
        * permute-sag2cor      -rotates     [x'=z;  y'=z;  z'=x]
        * permute-sag2axi      -rotates     [x'=y;  y'=z;  z'=x]
        * permute-cor2axi      -rotates     [x'=x;  y'=z;  z'=y]
        * permute-axi2sag      -rotates     [x'=y;  y'=z;  z'=x]
        * permute-axi2cor      -rotates     [x'=x;  y'=z;  z'=y]
        * permute-flipxyz      -flips in all x,y,z directions [x'=-x; y'=-y; z'=-z]
    *  morphologic operations
                -uses -radius=<R> to set the kernel radius
        * open         -morphologic opening
        * close        -morphologic closing
        * dilate       -morphologic dilation
        * erode        -morphologic erosion
        * closeholes   -morphologic closing by seed-filling background from
                        image edge
        * closebrain   -series of morphologic operations to find brain ROI:
                    (1) dilation, (2) close holes, and (3) erosion
    *  simple filters
        * gauss        -gaussian filter
        * median       -median filter
        * sobelx       -x-directional sobel filter
        * sobely       -y-directional sobel filter
        * sobelz       -z-directional sobel filter
        * sobelm       -sobel filter magnitude sqrt(x^2+y^2+z^2)
        * sobel        -sobel filter for all directions
        * sobelvoxel   -sobel filter for a voxel
    * more than 3 dimensions
        -v (6th)-dimension is used as color in this set of programs
        * mergetdim         -merges images in t(4th)-dim
        * mergeudim         -merges images in u(5th)-dim
    * swap routines assume the largest dimension is unused
        * swaput            -swaps u and t dimensions
        * swapuv            -swaps u and v dimensions
        * swapvt            -swaps v and t dimensions
    * changing header info
        * qformradio        -puts arbitrary qform with center in middle
        * removesqform      -removes sform and qform
        * copysqform        -copies sform and qform to output image
        * addsform          -adds sform to an image
        * modify-header-scl_slope -changes scl_slope
        * modify-header-dz  -changes dz
    * matrix operations
        * invmatrix         -inverts matrix
        * multiplymatrix    -multiplies (multiple) matrix(matrices)
        * affine-decompose decompose-affine -decomposes matrix into affine parameters
        * affinematrix      -creates affine matrix
        * halfway-matrix    -calculates halfway matrix
        * avgmatrix         -calculates average of matrices (approximate)
        * determinant       -calculates determinant of a matrix
    * affine registration / transformation
        * aregmni           -affine registration to MNI152 in FSL
        * aregmni-sform     -affine registration to MNI152 in FSL
                            -output matrix is put into image header sform
        * applysform        -applies affine transformation using sform
        * writesform        -writes sform affine transformation
        * NBCSR             -non-brain constrained symmetric registration
        * NBCR              -non-brain constrained registration (non-symmetric)
        * NBCSR-rigid       -non-brain constrained symmetric registration
                            rigid (6-dof) version
        * applyaffine       -applies affine transformation using input matrix
        * aregimg           -affine registration
        * aregimg2          -multi-seed/multi-resolution affine registration
    * nonlinear registration / transformation
        *   applywarp         -nonlinearly transforms images
        *   nregimg           -nonlinearly registers images (under construction)
        *   intratemplate     -nonlinearly registers images from the same subject and create a template
    *  OFF object processing
        *   cuberille         -creates a surface model from a binary image
        *   applyaffine-obj   -applies affine transformation using input
        *                      matrix on an off object
        *   fsasc2off         -converts FreeSurfer's ASCII format to off object
        *   shrinkwrap        -shrink-wraps a mask image (not well validated)
        *   sphere-off        -creates a sphere
        *   obj2img           -creates boundary at the object surface
        *   avgobjs           -averages objects, same object config required
        *   obj2pts off2pts   -creates point file with a list of 3d coordinates
                            for each vertex
        *   make-offe         -write complete off file with offe file
        *   combine-obj       -combine multiple off files
        *   add-off-color
        *   add-obj-color
                            -adds one color to an object
        *   add-off-color-curvature
        *   add-obj-color-curvature
                            -adds colors to an object based on curvatures
        *   add-off-color-sphere-theta
                            -creates a color-map of spherical coordinates (theta)
        *   add-off-color-sphere-psi
                            -creates a color-map of spherical coordinates (psi)
        *   offinfo           -displays off information
        *   vertexinfo        -displays vertex information
        *   edgeinfo          -displays edge information
        *   faceinfo          -displays face information
    * Cortex object processing (niikcortex)
        *   niikcortex-initocs - creates initial pial surface

        *   niikcortex-thick-color 
                              -assigns colour according to the cortical thickness
                              -objlist=<white.off>,<pial.off> -out=<out_white.off>,<out_pial.off>
                              -omin=<min>             minimum thickness for color range [default=auto]
                              -omax=<max>             maximum thickness for color range [default=auto]
                        
    * inter-packet (interleave) processing
        * interpacket-correction     -corrects interpacket patient motion by registering
        * interpacket-split          -splits interpacket (interleave) images by specified number
        * interpacket-merge          -merges interpacket (interleave) images
    * Other (less common)operations:
        * grid-lines            -creates image with grid-lines
                            e.g. to show nonlinear deformation fields
        * offinfo
    * general optional usage
        * -u -help          -help and general usage
        * --version         -version info
        * -list             -list of operations
    *  common options
        * -in=<img.nii>     -set input image
        * -out=<out>        -set output filename
        * -obj=<obj.off>    -set input object filename
        * -imglist=<in1.nii>,<in2.nii>,...
                        -list of input images
    *  matrix options
        The matrix format is same as FSL
        The matrix is in ascii format and should not have empty rows
        * -matrix=<mat>     -set input matrix
        * -matrixlist=<mat1>,<mat2>,...
        *                 -list of of input matrices
        *                 -separate by comma
        * -invmatrix=<mat>  -set and use inverse of input matrix
    *  XFM options
       The transformation format is linear MINC XFM with single transform
        * -xfm=<trans.xfm>     -set input matrix
        * -invxfm=<trans.xfm>  -set and use inverse of input matrix
    *  interpolation options
        * -nn                -nearest neighbor
        * -linear            -trilinear
        * -bspline           -b-spline
    *  data types
        * -uint8             -unsigned 8-bit integer
        * -uint16            -unsigned 16-bit integer
        * -uint32            -unsigned 32-bit integer
        * -uint64            -unsigned 64-bit integer
        * -int8              -8-bit signed integer
        * -int16             -16-bit signed integer
        * -int32             -32-bit signed integer
        * -int64             -64-bit signed integer
        * -float32           -32-bit floating point
        * -float64           -64-bit floating point
* falcon_cortex_shrinkwrap - create initial surface, starting from a sphere using fast marching shrink-wrap algorithm
* falcon_cortex_initics - create initial inner surface , from an initial one
* falcon_cortex_initocs - create initial outer surface by expanding inner surface
* falcon_cortex_check_off - check two surfaces self-intersections and correct it
* falcon_cortex_refine - refine inner and outer mesh simultaneously
* falcon_cortex_calc_thickness - calculate distance between two meshes
* falcon_cortex_smooth - smooth two surfaces simultaneously
* falcon_cortex_measure - calculate surface parameters (mostly for debugging)
* falcon_surface_check - check surface for self-intersections
* falcon_surface_refine - modify surface 
* falcon_surface_colour - colorize surface using measurements from another file
* falcon_surface_signals - extrace values along the surface from a volume
* falcon_surface_split - split surface into disconnected components
* falcon_transform_surface - apply transformation (.xfm file) - to a surface
* falcon_field_smooth - apply blurring along the surface
* falcon_obj2off - convert MNI .obj to .off or .ply file
* falcon_off2asc - convert .off or .ply file to .asc file
* falcon_off2ply - convert .off to .ply file
* falcon_ply2off - convert .ply to .off
* falcon_slice_extract - extract slices from a 3D volume with intersection lines from surfaces
* falcon_midsag - split white matter mask into hemishpheres

## Build
 * Using cmake , need to specify location of minc-toolkit

## Running
  * To extract cortical surfaces from T1w scan in stereotaxic space:

## Environment variables
* Verbosity (and generation of debug files) is controlled by environment variable NIIKVERBOSE, set to integer to specify verbosity level, default is 0
* Location of data directory is specified by NIIKDIR , set to location of <PREFIX>/data 
* If OpenMP is used, OMP_NUM_THREADS specifies maximum number of concurrent threads, default is number of CPU cores
