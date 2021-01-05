package = "minc2_simple"
version = "scm-1"

source = {
   url = "git://github.com/vfonov/minc2-simple.git",
   dir = "minc2-simple/lua"
   
}
description = {
   homepage = "https://github.com/vfonov/minc2-simple",
   license = "BSD",
   summary = "Interface for reading/writing MINC2 files"
}
dependencies = {}

build = {
   type = "cmake",
   modules = {},
   variables = {
    CMAKE_BUILD_TYPE="Release",
    CMAKE_PREFIX_PATH="$(LUA_BINDIR)/..",
    CMAKE_INSTALL_PREFIX="$(PREFIX)"
  }
}
