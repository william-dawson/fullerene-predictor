''' @package a script to visualize a fullerene
'''
from sys import argv
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_energy, compute_energy_novel, compute_euler_characteristic
from fullerene_curvature.fullerene import Fullerene
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":
    file_name = argv[1]

    input_fullerene = Fullerene(file_name)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

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
        ax.plot_wireframe(x_values,y_values,z_values)

    plt.show()
