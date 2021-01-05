## Simplified MINC2 API for C,Python and LUA
### More interfaces to come

## Goal
To provide access to most functionality of MINC2 file format with minimal effort, with consistent interface

## Installing

### Common Requirements
 * MINC2 library with headers, either by itself or as part of minc-toolkit or minc-toolkit-v2, see http://bic-mni.github.io/
 

### C
 * Requirements: cmake
 * Installation:
    ```
    mkdir build
    cd build
    cmake .. 
    make 
    ```
 * If location of minc2 library is not found:
    ```
    cmake -DLIBMINC_DIR:PATH=<location of libminc>
    ```
 
### LUA
 * Requirements: torch ( http://torch.ch/ )
 * Installation:
    ```
    cd lua
    luarocks make
    ```
  * If location of minc2 library is not found:
    ```
    cd build.luarocks/
    cmake -DLIBMINC_DIR:PATH=<location of libminc>
    cd ..
    luarocks make
    ```
    
### Python
 * Requirements: cffi, numpy, six
 * Optional requirements: scipy
 * Installation:
    ```
    python python/setup.py build
    python python/setup.py install 
    ```
 * If libminc is not found: set environment `MINC_TOOLKIT`  variable to point to the base on installation
