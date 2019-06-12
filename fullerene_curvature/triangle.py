'''@package triangle helper routines related to triangles.
'''
import numpy
import numpy.linalg


def compute_angle(p1, p2, p3):
    '''
    compute_angle: compute angle 23 of triangle 123.

    p1:
    p2:
    p3:

    return: angle value
    '''
    from numpy import array, arccos
    from numpy.linalg import norm

    d12 = norm(array(p1) - array(p2))
    d13 = norm(array(p1) - array(p3))
    d23 = norm(array(p2) - array(p3))

    return_value = arccos((d12**2 + d13**2 - d23**2) / (2 * d12 * d13))

    return return_value


def compute_area(p1, p2, p3):
    '''
    compute_area: compute the area of triangle 123.

    p1:
    p2:
    p3:

    return: area
    '''
    from numpy import array
    from numpy.linalg import norm

    d12 = norm(array(p1) - array(p2))
    d13 = norm(array(p1) - array(p3))
    d23 = norm(array(p2) - array(p3))

    s = 0.5 * (d12 + d13 + d23)
    return (s * (s - d12) * (s - d13) * (s - d23)) ** 0.5
