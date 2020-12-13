/*Based on sconvert from rply */

#include <stdio.h> 
#include <locale.h>
#include "rply.h"


static const char *const ply_type_list[] = {
    "int8", "uint8", "int16", "uint16",
    "int32", "uint32", "float32", "float64",
    "char", "uchar", "short", "ushort",
    "int", "uint", "float", "double",
    "list", NULL
};     /* order matches e_ply_type enum */


static void error_callback(p_ply ply, const char *message) {
    fprintf(stderr,"PLY:%s\n",message);
}


static int callback(p_ply_argument argument) {
    void *pdata;
    /* just pass the value from the input file to the output file */
    /*ply_get_argument_user_data(argument, &pdata, NULL);
    ply_write((p_ply) pdata, ply_get_argument_value(argument));*/
    return 1;
}

static int setup_callbacks(p_ply iply) {
    p_ply_element element = NULL;
    /* iterate over all elements in input file */
    fprintf(stdout,"Elements:\n");
    while ((element = ply_get_next_element(iply, element))) {
        p_ply_property property = NULL;
        long ninstances = 0;
        const char *element_name;
        ply_get_element_info(element, &element_name, &ninstances);
        fprintf(stdout,"\t%s:%ld\n", element_name, ninstances);

        /* iterate over all properties of current element */
        while ((property = ply_get_next_property(element, property))) {
            const char *property_name;
            e_ply_type type, length_type, value_type;
            ply_get_property_info(property, &property_name, &type, 
                    &length_type, &value_type);
            
            fprintf(stdout,"\t\t%s type:%s length type:%s value_type:%s\n",
                    property_name,
                    type>=0&&type<=PLY_LIST? ply_type_list[type]       :"NA",
                    length_type>=0&&length_type<=PLY_LIST ? ply_type_list[length_type]:"NA",
                    value_type>=0&&value_type<=PLY_LIST  ? ply_type_list[value_type] :"NA" );
        }
    }
    return 1;
}

int main(int argc, char *argv[]) {
    const char *value;
    const char *old_locale;
    p_ply iply; 

    if(argc<2)
    {
        fprintf(stdout,"Usage:\n"
        " falcon_ply_info <input.ply>\n"
        );
        return 1;
    }
    old_locale = setlocale(LC_NUMERIC, NULL);
    setlocale(LC_NUMERIC, "C");

    iply = ply_open(argv[1], error_callback, 0, NULL);
    if (!iply) {
        fprintf(stderr,"Error opeining %s\n",argv[1]);
        return 1; 
    }
    if (!ply_read_header(iply)){
        fprintf(stderr,"Error parsing header %s\n",argv[1]);
        return 1; 
    }

    if (!setup_callbacks(iply)) return 1; 
    /* pass comments and obj_infos from input to output */

    value = NULL;
    fprintf(stdout,"Comments:\n");
    while ((value = ply_get_next_comment(iply, value)))
    {
        fprintf(stdout,"\t%s:\n",value);
    }
    value = NULL;
    fprintf(stdout,"Object:\n");
    while ((value = ply_get_next_obj_info(iply, value)))
    {
        fprintf(stdout,"\t%s:\n",value);
    }
    if (!ply_read(iply)) return 1; 
    /* close up, we are done */
    if (!ply_close(iply)) return 1; 
    return 0;
}
