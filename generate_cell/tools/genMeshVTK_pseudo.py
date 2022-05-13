import numpy as np
import meshio
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()
    inputPath = args.input
    outputPath = args.output

    V = np.zeros((0, 3), dtype=np.float64)
    F = np.zeros((0, 1), dtype=np.int32)
    m = 1 # 1, 2, 4, 8
    with open(inputPath, 'r') as f:
        line = f.readline()
        line_split = line.split()
        numV = int(line_split[0])

        line = f.readline()
        line_split = line.split()
        numF = int(line_split[0])

        m = int(numV / numF)
        V.resize((numV, 3))
        F.resize((numF, m))

        for vI in range(numV): # read V
            line = f.readline()
            line_split = line.split()
            V[vI, 0] = float(line_split[0])
            V[vI, 1] = float(line_split[1])
            V[vI, 2] = float(line_split[2])

        for fI in range(numF): # read F
            line = f.readline()
            line_split = line.split()
            for j in range(m):
                F[fI, j] = int(line_split[j])

    points = V
    cells = [] # m == 1: point
    if m == 2:
        cells = [("line", F)]
    elif m == 4:
        cells = [("quad", F)]
    elif m ==8:
        cells = [("hexahedron", F)]
            
    meshio.write_points_cells(outputPath, points, cells)
