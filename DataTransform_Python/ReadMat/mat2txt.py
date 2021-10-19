
import h5py
import numpy as np

inputPath = '3d/prepare_ma.mat'
outputPath = '3d/3d.txt'

f = h5py.File(inputPath, 'r')

print('keys: ', f.keys())
data = np.array(f['iK_ma'], dtype=np.int32).transpose()
print('shape:', data.shape)
n = 2855883
nnz = data.shape[0]
print('n: %d, nnz: %d' %(n, nnz))
# print('0:', data[0][0])
# print('%d: %d' %(nnz-1, data[nnz-1][0]))

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
