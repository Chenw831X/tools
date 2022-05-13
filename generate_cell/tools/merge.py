import numpy as np
import argparse

if __name__ == '__main__':
    inputPath = []
    inputPath.append("../output/Initial.txt")
    inputPath.append("../output/solidL1.txt")
    inputPath.append("../output/solidL2.txt")
    inputPath.append("../output/solidL3.txt")
    outputPath = "../output/merge.txt"

    numV = 0
    numF = 0
    V = np.zeros((0, 3), dtype=np.float64)
    F = np.zeros((0, 8), dtype=np.int32)

    for i in range(4):
        with open(inputPath[i], 'r') as f:
            line = f.readline()
            line_split = line.split()
            nV = int(line_split[0])
            line = f.readline()
            line_split = line.split()
            nF = int(line_split[0])

            V = np.resize(V, (numV+nV, 3))
            F = np.resize(F, (numF+nF, 8))

            for vI in range(nV):
                line = f.readline()
                line_split = line.split()
                V[numV+vI, 0] = float(line_split[0])
                V[numV+vI, 1] = float(line_split[1])
                V[numV+vI, 2] = float(line_split[2])

            for fI in range(nF):
                line = f.readline()
                line_split = line.split()
                F[numF+fI, 0] = numV + int(line_split[0])
                F[numF+fI, 1] = numV + int(line_split[1])
                F[numF+fI, 2] = numV + int(line_split[2])
                F[numF+fI, 3] = numV + int(line_split[3])
                F[numF+fI, 4] = numV + int(line_split[4])
                F[numF+fI, 5] = numV + int(line_split[5])
                F[numF+fI, 6] = numV + int(line_split[6])
                F[numF+fI, 7] = numV + int(line_split[7])

            numV += nV
            numF += nF

    with open(outputPath, 'w') as f:
        f.write("{}\n".format(numV))
        f.write("{}\n".format(numF))
        
        for vI in range(numV):
            f.write("{} {} {}\n".format(V[vI, 0], V[vI, 1], V[vI, 2]))

        for fI in range(numF):
            f.write("{} {} {} {} {} {} {} {}\n".format(F[fI, 0], F[fI, 1], F[fI, 2], F[fI, 3],
                F[fI, 4], F[fI, 5], F[fI, 6], F[fI, 7]))

    # with open(inputPath, 'r') as f:
    #     for vI in range(numV): # read V
    #         line = f.readline()
    #         line_split = line.split()
    #         V[vI, 0] = float(line_split[0])
    #         V[vI, 1] = float(line_split[1])
    #         V[vI, 2] = float(line_split[2])

    #     for fI in range(numF): # read F
    #         line = f.readline()
    #         line_split = line.split()
    #         for j in range(m):
    #             F[fI, j] = int(line_split[j])

    # meshio.write_points_cells(outputPath, points, cells)