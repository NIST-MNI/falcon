require 'minc2_simple'
require 'xlua'

function prin_transform(tr) 
    print(string.format("center: %f %f %f",tr.center[1],tr.center[2],tr.center[3]))
    print(string.format("translations: %f %f %f",tr.translations[1],tr.translations[2],tr.translations[3]))
    print(string.format("scales: %f %f %f",tr.scales[1],tr.scales[2],tr.scales[3]))
    print(string.format("rotations: %f %f %f",tr.rotations[1],tr.rotations[2],tr.rotations[3]))
    print(string.format("shears: %f %f %f",tr.shears[1],tr.shears[2],tr.shears[3]))
end

xfm=minc2_xfm.new()
xfm:open('test.xfm')
transform=xfm:get_linear_transform_param()

prin_transform(transform)