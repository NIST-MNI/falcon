require 'minc2_simple'
require 'xlua'

--qqq=minc2_file.new('/home/vfonov/data/viola03/models/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc')
--qqq=minc2_file.new('/home/vfonov/mni/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc')
qqq=minc2_file.new('test_in.mnc')


xfm=minc2_xfm.new()

print(xfm:get_n_concat())

xfm:open('test.xfm')
print(string.format("Loaded minc %dD file",qqq:ndim()))
-- going to read and write in standard order (XYZ)
qqq:setup_standard_order()
-- let's check coordinate transfer
my_dims=qqq:representation_dims()


ooo=minc2_file.new()
ooo:define(my_dims, minc2_file.MINC2_BYTE, minc2_file.MINC2_FLOAT)
ooo:create('test_res.mnc')

ooo:setup_standard_order()

-- load into a c buffer
data=qqq:load_complete_volume(minc2_file.MINC2_FLOAT)
local sz=data:size()
out_data=torch.Tensor(sz)
out_data:fill(0.0)


print("Running naive resampler...")
for i=0,sz[1]-1 do
    xlua.progress(i,sz[1]-1)
    for j=0,sz[2]-1 do
        for k=0,sz[3]-1 do
            local out_xyz=ooo:voxel_to_world({i,j,k})
            local in_xyz=xfm:inverse_transform_point(out_xyz)
            local in_ijk=qqq:world_to_voxel(in_xyz)
            -- using nearest neigbour interpolation
            if in_ijk[1]>=0 and in_ijk[1]<=sz[1] and in_ijk[2]>=0 and in_ijk[2]<=sz[2] and in_ijk[3]>=0 and in_ijk[3]<=sz[3] then
                out_data[{i+1,j+1,k+1}]=data[{math.floor(in_ijk[1])+1,math.floor(in_ijk[2])+1,math.floor(in_ijk[3])+1}] -- indexes in minc are 0-based
            end
        end
    end
end


-- save from buffer to volume
ooo:save_complete_volume(out_data)

-- not strictly needed , but will flush the data to disk immedeately
ooo:close()
-- not strictly needed 
qqq:close()
