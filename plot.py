''' @package a script to visualize a fullerene
'''
from sys import argv
from fullerene_curvature.plot_curvature import plot_k, plot_g
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
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax2 = fig2.add_subplot(111, projection='3d')

    for ring in input_fullerene.ring_list:
        x_values = []
        y_values = []
        z_values = []
        for i in ring:
            x_values.append(input_fullerene.atoms_array[i][0])
            y_values.append(input_fullerene.atoms_array[i][1])
            z_values.append(input_fullerene.atoms_array[i][2])
        x_values.append(input_fullerene.atoms_array[ring[0]][0])
        y_values.append(input_fullerene.atoms_array[ring[0]][1])
        z_values.append(input_fullerene.atoms_array[ring[0]][2])
        ax1.plot_wireframe(x_values, y_values, z_values)
        ax2.plot_wireframe(x_values, y_values, z_values)

    x_values = []
    y_values = []
    z_values = []
    for point in input_fullerene.atoms_array:
        x_values.append(point[0])
        y_values.append(point[1])
        z_values.append(point[2])
    ax1.scatter(x_values, y_values, z_values, s=20, c='b', picker=5)
    ax2.scatter(x_values, y_values, z_values, s=20, c='b', picker=5)

    plot_k(input_fullerene, ax1, atom_number)
    plot_g(input_fullerene, ax2, atom_number)

    fig1.canvas.mpl_connect('pick_event', onpick)
    fig2.canvas.mpl_connect('pick_event', onpick)

    ax1.set_xlim3d(-4, 4)
    ax1.set_ylim3d(-4, 4)
    ax1.set_zlim3d(-4, 4)
    ax2.set_xlim3d(-4, 4)
    ax2.set_ylim3d(-4, 4)
    ax2.set_zlim3d(-4, 4)

    plt.show()
