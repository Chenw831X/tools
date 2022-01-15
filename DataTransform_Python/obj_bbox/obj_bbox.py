#
# obj_bbox.py
# obj_bbox
#
# Created by Wei Chen on 5/21/21
#

# read vertices in inputPath(.obj) to np.array
# print bbox of the vertices

import numpy as np

inputPath = 'T90.obj' # need modify
V = np.empty((0, 3), dtype=np.float64)

with open(inputPath, 'r') as f:
    line = f.readline()
    while line:
        line_split = line.split()
        if line_split[0] != 'v':
            break
        V = np.concatenate([V, [[float(line_split[1]), float(line_split[2]), float(line_split[3])]]], axis=0)
        line = f.readline()

# print(V)
print('max:', np.max(V, axis=0))
print('min:', np.min(V, axis=0))