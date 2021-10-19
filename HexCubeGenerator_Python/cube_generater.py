import meshio
import numpy as np


target_vtk = "bbox.vtk"
x_min, x_max, x_n = 0, 100, 10
y_min, y_max, y_n = 0, 20, 2
z_min, z_max, z_n = 0, 100, 10

'''
target_vtk = "HexMesh4Test.vtk"
x_min, x_max, x_n = 21, 29, 2
y_min, y_max, y_n = 12, 15, 1
z_min, z_max, z_n = 12, 15, 1
'''
if __name__ == '__main__':
    xs = np.linspace(x_min, x_max, x_n + 1)
    ys = np.linspace(y_min, y_max, y_n + 1)
    zs = np.linspace(z_min, z_max, z_n + 1)
    
    points = []
    cells = []
    ids = {}
    
    
    def gen_box(last_x, x, last_y, y, last_z, z):
        corners = [
            (x, y, z),
            (last_x, y, z),
            (last_x, last_y, z),
            (x, last_y, z),
            (x, y, last_z),
            (last_x, y, last_z),
            (last_x, last_y, last_z),
            (x, last_y, last_z)
        ]
        cube = []
        for p in corners:
            if p in ids:
                id = ids[p]
            else:
                id = len(points)
                points.append(p)
                ids[p] = id
            cube.append(id)
        cells.append(cube)
    
    
    for zi in range(1, len(zs)):
        z = zs[zi]
        last_z = zs[zi - 1]
        for yi in range(1, len(ys)):
            y = ys[yi]
            last_y = ys[yi - 1]
            for xi in range(1, len(xs)):
                x = xs[xi]
                last_x = xs[xi - 1]
                gen_box(last_x, x, last_y, y, last_z, z)
    
    points = np.array(points)
    cells = [("hexahedron", np.array(cells))]
    meshio.write_points_cells(
        target_vtk,
        points,
        cells
    )
