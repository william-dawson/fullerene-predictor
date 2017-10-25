import numpy
import numpy.linalg
import scipy
import scipy.sparse.csgraph
import itertools

##########################################################################
def path_length(point_position_pair):
    length = 0
    for i in range(1,len(point_position_pair)):
        position1 = point_position_pair[i-1][1]
        position2 = point_position_pair[i][1]
        length = length + numpy.linalg.norm(position1-position2)
    position1 = point_position_pair[-1][1]
    position2 = point_position_pair[0][1]
    length = length + numpy.linalg.norm(position1-position2)
    return length

##########################################################################
def sort_points(points, points_positions):
    point_position_pair = []
    for i in range(0,len(points)):
        point_position_pair.append([points[i], points_positions[i]])

    perm_list = itertools.permutations(point_position_pair)
    sorted_list = points.copy()
    shortest_length = 1e16
    for perm in perm_list:
        length = path_length(perm)
        if length < shortest_length:
            shortest_length = length
            sorted_list = []
            for pp in perm:
                sorted_list.append(pp[0])

    return sorted_list
