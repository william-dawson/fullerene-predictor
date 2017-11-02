''' @package curvature.py
Main routines to compute the curvature energy of a fullerene.
'''
import numpy
from fullerene_curvature.polygon import sort_points
from fullerene_curvature.sphere import compute_sphere
from fullerene_curvature.triangle import compute_angle


def compute_energy(k_array, g_array):
    '''!
    compute_energy compute the curvature energy using the K and G values.
    Equation 5.

    @param k_array: K values (equation 6)
    @param g_array: G values (equation 8)

    @return: curvature energy
    '''
    A = 2.62
    D = 1.41
    alpha = 0.165
    sum_value = 0
    for i in range(0, len(k_array)):
        sum_value = sum_value + 2 * k_array[i]**2 - ((1 - alpha) * g_array[i])
    sum_value = D * A * sum_value
    return sum_value


def compute_energy_novel(fullerene, k_array, g_array):
    '''!
    compute_energy_novel experimental ways to compute the curvature energy
    that take into account the different ring sizes.

    @param fullerene: the fullerne to calculate.
    @param k_array: K values (equation 6)
    @param g_array: G values (equation 8)

    @return: curvature energy
    '''
    A5 = 1.7292
    A6 = 2.62
    D5 = 0.423
    D6 = 1.41
    alpha5 = 0.1155
    alpha6 = 0.165

    sum_value = 0
    for i in range(0, len(fullerene.atoms_array)):
        # Determine the ring mixing coefficients
        atom_i_rings = fullerene.ring_lookup[i]
        length1 = len(fullerene.ring_list[atom_i_rings[0]])
        length2 = len(fullerene.ring_list[atom_i_rings[1]])
        length3 = len(fullerene.ring_list[atom_i_rings[2]])

        # Compute mixing coefficients
        A, D, alpha = 0, 0, 0
        for length in [length1, length2, length3]:
            if length == 5:
                A = A + A5 / 3.0
                D = D + D5 / 3.0
                alpha = alpha + alpha5 / 3.0
            else:
                A = A + A6 / 3.0
                D = D + D6 / 3.0
                alpha = alpha + alpha6 / 3.0

        # Add to the sum
        sum_value = sum_value + D * A * (2 *
                                         k_array[i]**2 - ((1 - alpha) * g_array[i]))

    return sum_value


def compute_k_values(fullerene):
    '''!
    compute_k_values (equation 5)

    @param fullerene: Fullerene to process.

    return: k value for each atom.
    '''
    k_values = []
    for i in range(0, len(fullerene.atoms_array)):
        neighbors = fullerene.connectivity[i]
        point_a = numpy.array(fullerene.atoms_array[i])
        point_b = numpy.array(fullerene.atoms_array[neighbors[0]])
        point_c = numpy.array(fullerene.atoms_array[neighbors[1]])
        point_d = numpy.array(fullerene.atoms_array[neighbors[2]])

        R = compute_sphere(point_a, point_b, point_c, point_d)
        k_values.append(1.0 / R)

    return k_values


def compute_g_values(fullerene):
    '''!
    compute_g_values (equation 8)

    @param fullerene: Fullerene to process.

    return: g value for each atom.
    '''
    A = 2.62
    g_values = []

    for i in range(0, len(fullerene.atoms_array)):
        # Who are my three neighboring rings
        v0 = fullerene.ring_lookup[i][0]
        v1 = fullerene.ring_lookup[i][1]
        v2 = fullerene.ring_lookup[i][2]

        # What rings are connected to each of my neighbor rings
        neighbor_vertex_v0 = numpy.array(fullerene.rings_connectivity[v0])
        neighbor_vertex_v1 = numpy.array(fullerene.rings_connectivity[v1])
        neighbor_vertex_v2 = numpy.array(fullerene.rings_connectivity[v2])

        # And what is the center of those neighbor's neighbor rings
        neighbor_vertex_positions_v0 = numpy.array([
            fullerene.ring_center[x] for x in neighbor_vertex_v0])
        neighbor_vertex_positions_v1 = numpy.array([
            fullerene.ring_center[x] for x in neighbor_vertex_v1])
        neighbor_vertex_positions_v2 = numpy.array([
            fullerene.ring_center[x] for x in neighbor_vertex_v2])

        # How many triangles surround each vertex neighbor
        n0 = len(fullerene.ring_list[v0])
        n1 = len(fullerene.ring_list[v1])
        n2 = len(fullerene.ring_list[v2])

        # For each vertex neighbor, what is its delta value
        Del_V0 = 2 * numpy.pi
        Del_V1 = 2 * numpy.pi
        Del_V2 = 2 * numpy.pi

        shortest_rings = sort_points(
            neighbor_vertex_v0, neighbor_vertex_positions_v0)
        for j in range(0, len(fullerene.ring_list[v0])):
            ring1 = v0
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(fullerene.ring_list[v0])]
            ring1_point = fullerene.ring_center[ring1]
            ring2_point = fullerene.ring_center[ring2]
            ring3_point = fullerene.ring_center[ring3]

            Del_V0 = Del_V0 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        shortest_rings = sort_points(
            neighbor_vertex_v1, neighbor_vertex_positions_v1)
        for j in range(0, len(fullerene.ring_list[v1])):
            ring1 = v1
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(fullerene.ring_list[v1])]
            ring1_point = fullerene.ring_center[ring1]
            ring2_point = fullerene.ring_center[ring2]
            ring3_point = fullerene.ring_center[ring3]

            Del_V1 = Del_V1 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        shortest_rings = sort_points(
            neighbor_vertex_v2, neighbor_vertex_positions_v2)
        for j in range(0, len(fullerene.ring_list[v2])):
            ring1 = v2
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(fullerene.ring_list[v2])]
            ring1_point = fullerene.ring_center[ring1]
            ring2_point = fullerene.ring_center[ring2]
            ring3_point = fullerene.ring_center[ring3]

            Del_V2 = Del_V2 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        Del_P = Del_V0 / n0 + Del_V1 / n1 + Del_V2 / n2
        g_values.append(Del_P / A)

    return g_values
