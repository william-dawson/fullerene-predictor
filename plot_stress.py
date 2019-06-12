''' @package a script to visualize a fullerene
'''
from sys import argv
from fullerene_curvature.fullerene import Fullerene
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def onpick(event):
    ind = event.ind[0]
    x, y, z = event.artist._offsets3d
    print(ind, x[ind], y[ind], z[ind])


if __name__ == "__main__":
    file_name = argv[1]
    atom_number = int(argv[2]) - 1

    input_fullerene = Fullerene(file_name)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection=Axes3D.name)

    for ring in input_fullerene.ring_list:
        for i in range(1, len(ring)):
            x_values = [input_fullerene.atoms_array[ring[i - 1]]
                        [0], input_fullerene.atoms_array[ring[i]][0]]
            y_values = [input_fullerene.atoms_array[ring[i - 1]]
                        [1], input_fullerene.atoms_array[ring[i]][1]]
            z_values = [input_fullerene.atoms_array[ring[i - 1]]
                        [2], input_fullerene.atoms_array[ring[i]][2]]
            ax1.plot3D(x_values, y_values, z_values, 'k')
        x_values = [input_fullerene.atoms_array[ring[-1]]
                    [0], input_fullerene.atoms_array[ring[0]][0]]
        y_values = [input_fullerene.atoms_array[ring[-1]]
                    [1], input_fullerene.atoms_array[ring[0]][1]]
        z_values = [input_fullerene.atoms_array[ring[-1]]
                    [2], input_fullerene.atoms_array[ring[0]][2]]
        ax1.plot3D(x_values, y_values, z_values, 'k')

    x_values = []
    y_values = []
    z_values = []
    for point in input_fullerene.atoms_array:
        x_values.append(point[0])
        y_values.append(point[1])
        z_values.append(point[2])
    ax1.scatter(x_values, y_values, z_values, s=20, c='b', picker=5)

    ax1.set_xlim3d(-4, 4)
    ax1.set_ylim3d(-4, 4)
    ax1.set_zlim3d(-4, 4)

    plt.show()
