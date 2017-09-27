import numpy
import numpy.linalg
import scipy
import scipy.sparse.csgraph
from tsp_solver.greedy import solve_tsp

##########################################################################


def sort_points(points, points_positions):
    shortest_rings = points.copy()
    ring_array = numpy.zeros((len(shortest_rings), len(shortest_rings)))
    for j in range(0, len(shortest_rings)):
        for i in range(0, len(shortest_rings)):
            ring_array[j, i] = numpy.linalg.norm(
                points_positions[j] - points_positions[i])
    
    local_path = solve_tsp(ring_array)
    path = [points[i] for i in local_path]
    return path
