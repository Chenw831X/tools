#
# msh_bbox.py
# msh_bbox
#
# Created by Wei Chen on 5/22/21
#

# read vertices in inputPath(.msh) to ndarray
# print bbox of the vertices

import numpy as np

inputPath = '/home/cw/MyCode/IPC/input/tetMeshes/T90_upper.msh' # need modifiy

with open(inputPath, 'r') as f:
    line = f.readline()
    while line:
        if line[:6] == '$Nodes':
            break
        line = f.readline()
    line = f.readline()
    line_split = line.split()
    vAmt = int(line_split[1])
    V = np.empty((vAmt, 3), dtype=np.float64)
    line = f.readline()

    for i in range(vAmt):
        line = f.readline()

    for i in range(vAmt):
        line = f.readline()
        line_split = line.split()
        V[i, 0] = float(line_split[0])
        V[i, 1] = float(line_split[1])
        V[i, 2] = float(line_split[2])

print('min:', np.min(V, axis=0))
print('max:', np.max(V, axis=0))
