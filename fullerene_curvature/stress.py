''' @package curvature.py
Main routines to compute the bond stress.
'''


def apply_kernel(matrix, kernel, beta):
    '''!
    Apply a kernel to the matrix representing the fullerene.


    @param matrix: the matrix representing the fullerene.
    @param kernel: name fo th ekernel to use ("I", "EXP", "INV", "Har")
    @param beta: the inverse temperature

    @return: a matrix representing the fullerne after the kernel has been
      applied.
    '''
    from copy import deepcopy
    from scipy.linalg import funm, pinv
    from numpy import exp

    if kernel == "EXP":
        kmat = funm(matrix, lambda x: exp(beta*x))
    elif kernel == "INV":
        kmat = funm(matrix, lambda x: -1.0/(1.0 - beta*x))
    elif kernel == "HAR":
        lapl = deepcopy(matrix)
        for i in range(0, matrix.shape[0]):
            lapl[i, i] = 0 - sum(lapl[i, :])
        kmat = beta * pinv(-1.0 * lapl)
    elif kernel == "I":
        kmat = deepcopy(matrix)

    return kmat


def curvature_stress(fullerene, k_array, g_array, kernel, beta):
    '''!
    compute_bond_stress estimate the stress between all pairs of bonds.

    @param k_array: K values (equation 6)
    @param g_array: G values (equation 8)
    @param kernel: which kernel to use ("I", "EXP", "INV")
    @param beta: the inverse temperature.

    @return: a dictionary mapping tuples of atoms to a stress value.
    '''
    from numpy import zeros
    site_array = []
    A = 2.62
    D = 1.41
    alpha = 0.165
    for i in range(0, len(k_array)):
        site_value = 2 * k_array[i]**2 - ((1 - alpha) * g_array[i])
        site_array.append(D * A * site_value)

    strain_matrix = zeros((len(k_array), len(k_array)))

    for i in range(0, len(k_array)):
        for neigh in fullerene.connectivity[i]:
            strain_matrix[i, neigh] = 0.5*(site_array[neigh] + site_array[i])

    strain_matrix = apply_kernel(strain_matrix, kernel, beta)

    strain_dict = {}
    for i in range(0, len(k_array)):
        for neigh in fullerene.connectivity[i]:
            if (neigh, i) not in strain_dict:
                strain_dict[(i, neigh)] = strain_matrix[i, neigh]
                if kernel == "INV":
                    strain_dict[(i, neigh)] = 1.0/strain_dict[(i, neigh)]

    return strain_dict


def distance_stress(fullerene, kernel, beta):
    '''!
    compute_bond_stress estimate the stress between all pairs of bonds.

    @param kernel: which kernel to use ("I", "EXP", "INV")
    @param beta: the inverse temperature.

    @return: a dictionary mapping tuples of atoms to a stress value.
    '''
    from numpy import zeros

    strain_matrix = fullerene.get_adjacency_matrix(distance=True)
    strain_matrix = apply_kernel(strain_matrix, kernel, beta)

    strain_dict = {}
    for i in range(0, strain_matrix.shape[0]):
        for neigh in fullerene.connectivity[i]:
            if (neigh, i) not in strain_dict:
                strain_dict[(i, neigh)] = strain_matrix[i, neigh]
                if kernel == "INV":
                    strain_dict[(i, neigh)] = 1.0/strain_dict[(i, neigh)]

    return strain_dict
