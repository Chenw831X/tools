#
# ReadMat.py
# ReadMat
#
# Created by Wei Chen on 5/20/21
#

# load data in .mat, transform it to np.array, and the obtained array is transpose of data in matlab
# notice:
# if version of .mat is v7.3, use h5py to load .mat
# else, use scipy

import h5py
import numpy as np

inputPath = input("Please enter path of .mat: ")

f = h5py.File(inputPath, 'r')

print('keys: ', f.keys())
data = np.array(f['jK_ma'], dtype=np.int32) # key in f[key] differs in different .mat, dtype similar
data = data.transpose() # transpose to obtain right data
print('shape: ', data.shape)
print('data')
print(data)