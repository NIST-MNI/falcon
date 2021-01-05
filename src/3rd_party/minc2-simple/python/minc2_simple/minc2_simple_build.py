import os
import sys
from cffi import FFI
from sys import platform
  
ffi = FFI()

minc_prefix=os.environ.get('MINC_TOOLKIT',"/opt/minc/1.9.17")
source_path=os.path.join(os.path.dirname(__file__), "../../src")

_extra_link_args_debug=['-g']
_extra_compile_args_debug=['-fsanitize=address','-fno-omit-frame-pointer','-g']

_extra_link_args = []
_extra_compile_args = []

minc2_simple_src=""
minc2_simple_defs=""

with open(os.path.join(source_path,"minc2-simple.c"),'r') as f:
    minc2_simple_src=f.read()
with open(os.path.join(source_path,"minc2-matrix-ops.c"),'r') as f:
    minc2_simple_src+=f.read()
with open(os.path.join(source_path,"minc2-simple-int.h"),'r') as f:
  minc2_simple_defs+=f.read()

# add free system call
minc2_simple_defs+="""
void free(void *ptr);
"""

if platform == "linux" or platform == "linux2":
  _extra_link_args=['-Wl,-rpath={}'.format(os.path.join(minc_prefix,"lib"))]
elif platform == "darwin":
  _extra_link_args=['-Xlinker','-rpath','-Xlinker',os.path.join(minc_prefix,"lib")]

ffi.set_source("minc2_simple._simple",
    minc2_simple_src,
    # The important thing is to include libc in the list of libraries we're
    # linking against:
    libraries=["minc2","c"],
    include_dirs=[os.path.join(minc_prefix,"include"),source_path],
    library_dirs=[os.path.join(minc_prefix,"lib")],
    extra_compile_args=_extra_compile_args,
    extra_link_args=_extra_link_args
)

ffi.cdef(minc2_simple_defs )


if __name__ == "__main__":
    ffi.compile(verbose=True)
