import numpy
import math
import numpy.linalg
from polygon import sort_points
from sphere import compute_sphere
from triangle import compute_angle, compute_area

##########################################################################
def compute_energy(k_array, g_array):
    A = 2.62
    D = 1.41
    alpha = 0.165
    sum_value = 0
    for i in range(0, len(k_array)):
        sum_value = sum_value + 2 * k_array[i]**2 - (1 - alpha) * g_array[i]
        # print(k_array[i]**2, g_array[i])
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
        k_values.append(1.0/R)

    return k_values

##########################################################################
def compute_g_values(atoms, connect, rings, rings_center, rings_lookup):
    A = 2.62
    g_values = []

    for i in range(0, len(atoms)):
        # Who are my vertex neighbors
        v0 = rings_lookup[i][0]
        v1 = rings_lookup[i][1]
        v2 = rings_lookup[i][2]
        ring_positions_v0 = numpy.array([atoms[x] for x in rings[v0]])
        ring_positions_v1 = numpy.array([atoms[x] for x in rings[v1]])
        ring_positions_v2 = numpy.array([atoms[x] for x in rings[v2]])

        # How many triangles surround each vertex neighbor
        n0 = len(rings[v0])
        n1 = len(rings[v1])
        n2 = len(rings[v2])

        # For each vertex neighbor, what is its delta value
        Del_V0 = 2*numpy.pi
        Del_V1 = 2*numpy.pi
        Del_V2 = 2*numpy.pi

        shortest_rings = sort_points(rings[v0], ring_positions_v0)
        for i in range(0,len(rings[v0])):
            point1 = v0
            point2 = shortest_rings[i]
            point3 = shortest_rings[(i+1)%len(rings[v0])]
            point1 = rings_center[v0]
            point2 = atoms[point2]
            point3 = atoms[point3]

            Del_V0 = Del_V0 - (compute_angle(point1,point2,point3))

        shortest_rings = sort_points(rings[v1], ring_positions_v1)
        for i in range(0,len(rings[v1])):
            point1 = v1
            point2 = shortest_rings[i]
            point3 = shortest_rings[(i+1)%len(rings[v1])]
            point1 = rings_center[v1]
            point2 = atoms[point2]
            point3 = atoms[point3]

            Del_V1 = Del_V1 - (compute_angle(point1,point2,point3))

        shortest_rings = sort_points(rings[v2], ring_positions_v2)
        for i in range(0,len(rings[v2])):
            point1 = v2
            point2 = shortest_rings[i]
            point3 = shortest_rings[(i+1)%len(rings[v2])]
            point1 = rings_center[v2]
            point2 = atoms[point2]
            point3 = atoms[point3]

            Del_V2 = Del_V2 - (compute_angle(point1,point2,point3))

        Del_P = Del_V0/n0 + Del_V1/n1 + Del_V2/n2
        g_values.append(0.25*Del_P/A)
    print(g_values)
    return g_values
