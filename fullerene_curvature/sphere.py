'''@package sphere helper routines related to spheres.
'''
import numpy
import numpy.linalg

def compute_sphere(point_a, point_b, point_c, point_d):
    '''
    compute_sphere Given four points, this computes the radius on a sphere
    containing them.

    point_a:
    point_b:
    point_c:
    point_d:

    return: radius of the sphere.
    '''
    a0 = point_a[0]
    a1 = point_a[1]
    a2 = point_a[2]
    b0 = point_b[0]
    b1 = point_b[1]
    b2 = point_b[2]
    c0 = point_c[0]
    c1 = point_c[1]
    c2 = point_c[2]
    d0 = point_d[0]
    d1 = point_d[1]
    d2 = point_d[2]

    print(point_a, point_b, point_c, point_d)

    left_matrix = [[a0 - b0, a1 - b1, a2 - b2],
                   [b0 - c0, b1 - c1, b2 - c2],
                   [c0 - d0, c1 - d1, c2 - d2]]

    right_0 = (a0**2 + a1**2 + a2**2 - b0**2 - b1**2 - b2**2)/2.0
    right_1 = (b0**2 + b1**2 + b2**2 - c0**2 - c1**2 - c2**2)/2.0
    right_2 = (c0**2 + c1**2 + c2**2 - d0**2 - d1**2 - d2**2)/2.0
    right_array = [[right_0],[right_1],[right_2]]

    left_array = numpy.dot(numpy.linalg.inv(left_matrix), right_array)
    left_array = left_array.T[0]

    R = numpy.linalg.norm(point_a - left_array)

    return R
