#ifndef IGL_STB_IMAGE_H
#define IGL_STB_IMAGE_H

namespace igl {
  unsigned char * stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
  void stbi_image_free(void *retval_from_stbi_load);
  int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
} // namespace igl

#endif
