#! /bin/bash

# <extract surfaces on icbm template> 

# create mid surface
falcon_igl_mesh_avg \
                    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_ocs-0.9.9.ply \
                    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_ics-0.9.9.ply \
                    -o mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9.ply

# split left-right
falcon_igl_mesh_split \
    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9.ply \
    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-%d.ply 


for s in 0 1;do 
# sample
../build/src/igl/falcon_igl_field_sampler \
    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-${s}.ply \
    icbm152_model_09c/mni_icbm152_CerebrA_tal_nlin_sym_09c.mnc  \
    cerebra_atlas_${s}_pre.csv.gz \
    --header cerebra --labels --clobber --distance 15.0
done

# remove non cortex labels

R --vanilla <<END
library("tidyverse")

#select cortical labels

ctx<-read_csv("icbm152_model_09c/cerebra_surface_labels.csv")|>
   filter(mindboggle_id>=2000) |> select(cerebra_label)

for(s in c(0,1)) {
i=paste0("cerebra_atlas_",s,"_pre.csv.gz")
o=paste0("cerebra_atlas_",s,".csv.gz")

df<-read_csv(i) |> 
      mutate(cerebra=if_else(cerebra %in% ctx\$cerebra_label,cerebra,0) )

write_csv(df,o)
}
END

for s in 0 1 ;do
falcon_igl_mesh_render \
      mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-${s}.ply --csv cerebra_atlas_${s}.csv.gz  \
      --field cerebra \
      --zoom 2.0 \
      --six --relevel --double \
      --output mni_icbm152_mid_cerebra_${s}.png
done


# resample to atlas space
falcon_igl_field_resample \
        -i cerebra_atlas_0.csv.gz \
        mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-0.ply \
        icbm152_model_09c/mni_icbm152_ics_sm_lt.ply \
        -o icbm152_model_09c/mni_icbm152_ics_sm_lt_atlas_cerebra.csv.gz  \
        --majority_invexp --SO3 --knn 3 --clobber

falcon_igl_field_resample \
        -i cerebra_atlas_1.csv.gz \
        mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-1.ply \
        icbm152_model_09c/mni_icbm152_ics_sm_rt.ply \
        -o icbm152_model_09c/mni_icbm152_ics_sm_rt_atlas_cerebra.csv.gz  \
        --majority_invexp --SO3 --knn 3 --clobber



# show result rendering
falcon_igl_mesh_render \
      icbm152_model_09c/mni_icbm152_ics_sm_lt.ply --csv icbm152_model_09c/mni_icbm152_ics_sm_lt_atlas_cerebra.csv.gz  \
      --field cerebra \
      --zoom 2.0  \
      --six --relevel \
      --output mni_icbm152_cerebra_sm_lt.png

falcon_igl_mesh_render \
      icbm152_model_09c/mni_icbm152_ics_sm_rt.ply --csv icbm152_model_09c/mni_icbm152_ics_sm_rt_atlas_cerebra.csv.gz  \
      --field cerebra \
      --zoom 2.0  \
      --six --relevel \
      --output mni_icbm152_cerebra_sm_rt.png

