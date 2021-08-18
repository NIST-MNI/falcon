/*Implementation file for stb */
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

//#define STBI_WRITE_NO_STDIO
#ifdef HAVE_ZLIB
/* ZLIB */
#include <zlib.h>
#include <stdlib.h>

#pragma warning "Using ZLIB"

unsigned char * _zlib_compress(unsigned char *data, int data_len, int *out_len, int quality)
{
    /*quick and dirty*/
    unsigned long out_len_=compressBound(data_len);
    unsigned char *data_out=malloc(out_len_);

    if(data_out!=NULL)
    {
        if(compress2(data_out,&out_len_, data,data_len,quality)==Z_OK)
        {
            *out_len=(int)out_len_;
            return data_out;
        }
        else {
            /*TODO: handle compression error?*/
            free(data_out);
            return NULL;
        }
    } else {
        /*Out of memory?*/
        return NULL;
    }
}
#define STBIW_ZLIB_COMPRESS _zlib_compress
#endif




#include "stb_image.h"
#include "stb_image_write.h"