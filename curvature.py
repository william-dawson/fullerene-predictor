import numpy
import math
import numpy.linalg
from polygon import sort_points
from sphere import compute_sphere
from triangle import compute_angle, compute_area

##########################################################################


def compute_ring_connectivity(rings, rings_lookup):
    ring_connectivity = []
    for i in range(0, len(rings)):
        temp_list = []
        for j in range(0, len(rings[i])):
            neighbor_ring = rings_lookup[rings[i][j]]
            temp_list.extend(neighbor_ring)
        temp_set = set(temp_list)
        temp_set.remove(i)
        ring_connectivity.append(list(temp_set))

    return ring_connectivity

##########################################################################


def compute_energy(k_array, g_array):
    A = 2.62
    D = 1.41
    alpha = 0.165
    sum_value = 0
    for i in range(0, len(k_array)):
        sum_value = sum_value + 2 * k_array[i]**2 - ((1 - alpha) * g_array[i])
    sum_value = D * A * sum_value
    return sum_value

##########################################################################


def compute_k_values(atoms, connect):
    k_values = []
    for i in range(0, len(atoms)):
        neighbors = connect[i]
        point_a = numpy.array(atoms[i])
        point_b = numpy.array(atoms[neighbors[0]])
        point_c = numpy.array(atoms[neighbors[1]])
        point_d = numpy.array(atoms[neighbors[2]])

        R = compute_sphere(point_a, point_b, point_c, point_d)
        k_values.append(1.0 / R)

    return k_values

##########################################################################


def compute_g_values(atoms, connect, rings, rings_center, rings_lookup):
    A = 2.62
    g_values = []

    rings_connectivity = compute_ring_connectivity(rings, rings_lookup)

    for i in range(0, len(atoms)):
        # Who are my three neighboring rings
        v0 = rings_lookup[i][0]
        v1 = rings_lookup[i][1]
        v2 = rings_lookup[i][2]

        # What rings are connected to each of my neighbor rings
        neighbor_vertex_v0 = numpy.array(rings_connectivity[v0])
        neighbor_vertex_v1 = numpy.array(rings_connectivity[v1])
        neighbor_vertex_v2 = numpy.array(rings_connectivity[v2])

        # And what is the center of those neighbor's neighbor rings
        neighbor_vertex_positions_v0 = numpy.array([
            rings_center[x] for x in neighbor_vertex_v0])
        neighbor_vertex_positions_v1 = numpy.array([
            rings_center[x] for x in neighbor_vertex_v1])
        neighbor_vertex_positions_v2 = numpy.array([
            rings_center[x] for x in neighbor_vertex_v2])

        # How many triangles surround each vertex neighbor
        n0 = len(rings[v0])
        n1 = len(rings[v1])
        n2 = len(rings[v2])

        # For each vertex neighbor, what is its delta value
        Del_V0 = 2 * numpy.pi
        Del_V1 = 2 * numpy.pi
        Del_V2 = 2 * numpy.pi

        shortest_rings = sort_points(
            neighbor_vertex_v0, neighbor_vertex_positions_v0)
        for j in range(0, len(rings[v0])):
            ring1 = v0
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(rings[v0])]
            ring1_point = rings_center[ring1]
            ring2_point = rings_center[ring2]
            ring3_point = rings_center[ring3]

            Del_V0 = Del_V0 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        shortest_rings = sort_points(
            neighbor_vertex_v1, neighbor_vertex_positions_v1)
        for j in range(0, len(rings[v1])):
            ring1 = v1
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(rings[v1])]
            ring1_point = rings_center[ring1]
            ring2_point = rings_center[ring2]
            ring3_point = rings_center[ring3]

            Del_V1 = Del_V1 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        shortest_rings = sort_points(
            neighbor_vertex_v2, neighbor_vertex_positions_v2)
        for j in range(0, len(rings[v2])):
            ring1 = v2
            ring2 = shortest_rings[j]
            ring3 = shortest_rings[(j + 1) % len(rings[v2])]
            ring1_point = rings_center[ring1]
            ring2_point = rings_center[ring2]
            ring3_point = rings_center[ring3]

            Del_V2 = Del_V2 - (compute_angle(ring1_point,
                                             ring2_point, ring3_point))

        Del_P = Del_V0 / n0 + Del_V1 / n1 + Del_V2 / n2
        g_values.append(Del_P / A)

    return g_values
