''' @package plot_curvature.py
Visualize the calculation of the curvature
'''
import numpy
import random
from fullerene_curvature.polygon import sort_points
from fullerene_curvature.sphere import compute_sphere
from fullerene_curvature.triangle import compute_angle


def plot_k(fullerene, ax, atom_number):
    '''!
    plot_k: visualization of the calculation of K

    fullerene: to compute
    ax: to plot on
    atom_number: which atom to focus on
    '''
    print("Plot_K_Legend:")
    print('''\tGreen: "Point on the Fullerene"''')
    print('''\tRed: "That point's nearest neighbors"''')
    print('''\tYellow: "Center of the circle"''')
    i = atom_number
    neighbors = fullerene.connectivity[i]
    a = numpy.array(fullerene.atoms_array[i])
    b = numpy.array(fullerene.atoms_array[neighbors[0]])
    c = numpy.array(fullerene.atoms_array[neighbors[1]])
    d = numpy.array(fullerene.atoms_array[neighbors[2]])

    R, center = compute_sphere(a, b, c, d)

    x_values = []
    y_values = []
    z_values = []

    ax.scatter(a[0], a[1], a[2], s=80, c='g')
    for point in ([b, c, d]):
        x_values.append(point[0])
        y_values.append(point[1])
        z_values.append(point[2])

    ax.scatter(x_values, y_values, z_values, s=80, c='r')

    ax.scatter(center[0], center[1], center[2], s=80, c='y', picker=5)


def plot_g(fullerene, ax, atom_number):
    '''!
    plot_k: visualization of the calculation of G

    fullerene: to compute
    ax: to plot on
    atom_number: which atom to focus on
    '''
    print("Plot_G_Legend:")
    print('''\tGreen: "Point on the Fullerene"''')
    print('''\tColored_Dots: "Center of rings it participates in"''')
    print('''\tColored_Lines: "Center of neighbor's of those rings"''')
    i = atom_number
    # Who are my three neighboring rings
    v0 = fullerene.ring_lookup[i][0]
    v1 = fullerene.ring_lookup[i][1]
    v2 = fullerene.ring_lookup[i][2]

    # And what is the center of those neighbor rings
    a = numpy.array(fullerene.atoms_array[i])
    b = numpy.array(fullerene.ring_center[v0])
    c = numpy.array(fullerene.ring_center[v1])
    d = numpy.array(fullerene.ring_center[v2])

    x_values = []
    y_values = []
    z_values = []

    ax.scatter(a[0], a[1], a[2], s=80, c='g')
    for point in ([b, c, d]):
        x_values.append(point[0])
        y_values.append(point[1])
        z_values.append(point[2])

    # ax.scatter(x_values, y_values, z_values, s=80, c='r')

    # What rings are connected to each of my neighbor rings
    for v, c in zip([v0, v1, v2], ['black', 'yellow', 'pink']):
        ring_center = fullerene.ring_center[v]
        ax.scatter(ring_center[0], ring_center[1], ring_center[2], s=80, c=c)
        neighbor_vertex = numpy.array(fullerene.rings_connectivity[v])
        neighbor_vertex_positions = numpy.array([
            fullerene.ring_center[x] for x in neighbor_vertex])

        x_values = []
        y_values = []
        z_values = []

        for point in (neighbor_vertex_positions):
            x_values.append(point[0])
            y_values.append(point[1])
            z_values.append(point[2])

        # ax.scatter(x_values, y_values, z_values, s=80, c='y')

        shortest_rings = sort_points(
            neighbor_vertex, neighbor_vertex_positions)
        x_values = []
        y_values = []
        z_values = []
        for j in range(0, len(fullerene.ring_list[v])):
            p = fullerene.ring_center[shortest_rings[j]]
            x_values.append(p[0])
            y_values.append(p[1])
            z_values.append(p[2])
        p = fullerene.ring_center[shortest_rings[0]]
        x_values.append(p[0])
        y_values.append(p[1])
        z_values.append(p[2])
        ax.plot_wireframe(x_values, y_values, z_values, color=c)
