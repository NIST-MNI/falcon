#! /usr/bin/env th

require 'minc2_simple'
-- require 'torch'
require 'math'

local xfm_file1=arg[1]
local xfm_file2=arg[2]


local xfm=minc2_xfm.new(xfm_file1)


if xfm_file2 then
    local xfm2=minc2_xfm.new(xfm_file2)
    
    --concatenate inverted xfm2
    xfm:invert()
    xfm:concat_xfm(xfm2)
--    xfm2:invert()
    
--    xfm:concat_xfm(xfm2)
end


local edges={ {-60,-94, -52},
              { 60, 50, 78}
            }



local x,y,z
local max_dist=0.0
for x=1,2 do
    for y=1,2 do
        for z=1,2 do
            local p_in={ edges[x][1],edges[y][2],edges[z][3] }
            local p_out=xfm:transform_point(p_in)
            
            local dist=math.sqrt( (p_in[1]-p_out[1])^2+
                                  (p_in[2]-p_out[2])^2+
                                  (p_in[3]-p_out[3])^2)
            if dist>max_dist then
                max_dist=dist
            end
            
        end
    end
end

print(max_dist)


