require 'minc2_simple'

--qqq=minc2_file.new('/home/vfonov/data/viola03/models/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc')
--qqq=minc2_file.new('/home/vfonov/mni/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc')
qqq=minc2_file.new('test_in.mnc')
print(string.format("Loaded minc %dD file",qqq:ndim()))
dims=qqq:store_dims()

for i=0,(qqq:ndim()-1) do -- contrary to common LUA convention, it is 0-based
    print(string.format('Dimension %d length:%d id:%s start:%f step:%f',
        i,dims[i].length,dims[i].id,dims[i].start,dims[i].step))
end

print("Representation volume size:")
print(qqq:volume_size())


ooo=minc2_file.new()
-- will create file with same dimensions

my_dims={
     {id=minc2_file.MINC2_DIM_X,  length=73,start=96.0,  step=-1.0},
     {id=minc2_file.MINC2_DIM_Y,  length=87,start=-132.0,step=1.0},
     {id=minc2_file.MINC2_DIM_Z,  length=91,start=-78.0, step=1.0}
}

print("History:",qqq:read_attribute("","history"))
print("Test number attribute:",qqq:read_attribute("test","att1"))
print("Reading metadata")
m=qqq:metadata()
print("Metadata:")
print(m)

ooo:define(my_dims, minc2_file.MINC2_BYTE, minc2_file.MINC2_FLOAT)
ooo:create('test_out.mnc')
ooo:setup_standard_order()

print("Writing metadata")
ooo:write_metadata(m)

-- going to read and write in standard order (XYZ)
qqq:setup_standard_order()
print("Representation volume size again (standard mode):")
print(qqq:volume_size())


-- let's check coordinate transfer
my_dims=qqq:representation_dims()

ijk=torch.Tensor({0,0,0})
xyz=qqq:voxel_to_world(ijk)
print(string.format("Voxel %s corresponds to %s",ijk,xyz))

ijk=torch.Tensor({10,10,1.5})
xyz=qqq:voxel_to_world(ijk)
print(string.format("Voxel %s corresponds to %s",ijk,xyz))

function str(t) 
    local s="{"..t[1]
    local i,j
    for i=2,#t do
        s=s..","..t[i]
    end
    s=s.."}"
    return s
end

ijk={2,3,1.5}
xyz=qqq:voxel_to_world(ijk)
print(string.format("Voxel %s corresponds to %s",ijk,xyz))


ooo:setup_standard_order()
print("Loading from file...")

-- load into a c buffer
data=qqq:load_complete_volume(minc2_file.MINC2_FLOAT)

print(string.format("Complete tensor: %s of size :%dx%dx%d",torch.type(data),data:size(1),data:size(2),data:size(3)))

data=qqq:load_hyperslab(minc2_file.MINC2_FLOAT, {{1,91},{2,88},{3,75}})

print(string.format("Hyperslab tensor %s of size :%dx%dx%d",torch.type(data),data:size(1),data:size(2),data:size(3)))


-- save from buffer to volume
ooo:save_complete_volume(data)

rep=qqq:representation_dims()
ttt=minc2_file.new()
ttt:define(rep, minc2_file.MINC2_FLOAT, minc2_file.MINC2_FLOAT)
ttt:create('test_hyperslab.mnc')
ttt:setup_standard_order()

local z=torch.Tensor(rep[2].length,rep[1].length,rep[0].length):zero()

--print(z:size(),data:size())
--ttt:save_complete_volume(z)
print("-----")
ttt:set_volume_range(0,100)
ttt:save_hyperslab(data,{1,1,1})

-- not strictly needed , but will flush the data to disk immedeately
ooo:close()
qqq:close()
ttt:close()
