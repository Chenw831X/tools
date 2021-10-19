#
# msh2ipcmsh.py
# msh2ipcmsh
#
# Created by Wei Chen on 5/22/21
#

# transform inputPath(.msh) to outputPath(.msh), which is the format used in ipc

import numpy as np

inputPath = '/home/cw/project/2021-09-01/msh/' # need modify
outputPath = '/home/cw/project/2021-09-01/ipcmsh/' # need modify
output = open(outputPath, 'w')

with open(inputPath, 'r') as f:
    line = f.readline().strip('\n')
    while line:
        print(line, file=output)
        if line[:6] == '$Nodes':
            break
        line = f.readline().strip('\n')
    line = f.readline().strip('\n')
    line_split = line.split()
    print(line_split[0], line_split[1], file=output)
    vAmt = int(line_split[1])
    line = f.readline().strip('\n')
    print(line, file=output)

    for i in range(vAmt):
        line = f.readline().strip('\n')
    
    for i in range(vAmt):
        line = f.readline().strip('\n')
        print(i+1, line, file=output)
    
    line = f.readline().strip('\n')
    print(line, file=output)
    line = f.readline().strip('\n')
    print(line, file=output)

    line = f.readline().strip('\n')
    line_split = line.split()
    print(line_split[0], line_split[1], file=output)
    fAmt = int(line_split[1])
    line = f.readline().strip('\n')
    print(line, file=output)

    for i in range(fAmt):
        line = f.readline().strip('\n')
        print(line, file=output)
    
    line = f.readline().strip('\n')
    print(line, file=output)

    print('finished!')

output.close()