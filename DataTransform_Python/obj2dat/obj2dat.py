#
# obj2dat.py
# obj2dat
#
# Created by Wei Chen on 5/22/21
#

# transform 'inputPath+start+.obj' ~ 'inputPath+end+.obj' to
# 'outputPath+start+.dat' ~ 'outputPath+end+.dat'

import meshio

inputPath = '/home/cw/project/2021-12-9----shoe_simulation/T90_sphere/obj/' # need modify
outputPath = '/home/cw/project/2021-12-9----shoe_simulation/T90_sphere/dat/' # need modify
start = 0 # need modify
end = 200 # need modify

for i in range(start, end+1):
    if i%10 == 0:
        print(i)
    mesh = meshio.read(inputPath+str(i)+'.obj')
    mesh.write(outputPath+str(i)+'.dat')

print('finished!')