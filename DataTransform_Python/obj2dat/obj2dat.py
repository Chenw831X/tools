#
# obj2dat.py
# obj2dat
#
# Created by Wei Chen on 5/22/21
#

# transform 'inputPath+start+.obj' ~ 'inputPath+end+.obj' to
# 'outputPath+start+.dat' ~ 'outputPath+end+.dat'

import meshio

inputPath = 'E:\\项目\\清锋\\result_obj\\lattice_type501__10_2_10\\' # need modify
outputPath = 'E:\\项目\\清锋\\result_dat\\lattice_type501__10_2_10\\' # need modify
start = 0 # need modify
end = 100 # need modify

for i in range(start, end+1):
    if i%10 == 0:
        print(i)
    mesh = meshio.read(inputPath+str(i)+'.obj')
    mesh.write(outputPath+str(i)+'.dat')

print('finished!')