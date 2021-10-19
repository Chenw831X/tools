#
# ReadMat_scipy.py
# ReadMat
#
# Created by Wei Chen on 5/31/21
#

# load data in .mat, transform it to np.array, and the obtained array is transpose of data in matlab
# notice:
# if version of .mat is v7.3, use h5py to load .mat
# else, use scipy

import scipy.io as scio
import numpy as np

inputPath = 'input/ia.mat'
outputPath = 'output/Ab.txt'

f = scio.loadmat(inputPath)

print('keys: ', f.keys())
data = f['ia']
print('shape:', data.shape)
print('type:', data.dtype)

n = 41472
nnz = data.shape[0]
print('n: %d, nnz: %d' %(n, nnz))
#print('0:', data[0][0])
#print('%d: %d' %(nnz-1, data[nnz-1][0]))

csrR = np.empty((n+1), dtype=np.int32)
csrR[0] = 0
row = 1
cnt = 0

for i in range(nnz):
    #if i%10000000 == 0:
    #    print('i: ', i)
    print(i)
    
    now = data[i][0]
    if now == row:
        cnt = cnt + 1
    else:
        csrR[row] = csrR[row-1] + cnt
        for j in range(row+1, now):
            csrR[j] = csrR[j-1]
        row = now
        cnt = 1

csrR[row] = csrR[row-1] + cnt
for j in range(row+1, n+1):
    csrR[j] = csrR[j-1]

output = open(outputPath, 'a')
print(n+1, file=output)

for i in range(n+1):
    print(csrR[i], file=output)

print('write iK_ma finished!')