''' @package curvature.py
Main routines to compute the curvature energy of a fullerene.
'''


def compute_bond_stress(fullerene, k_array, g_array, kernel, beta):
    '''!
    compute_bond_stress estimate the stress between all pairs of bonds.

    \f[
    \Delta E_C = DA\sum_i [ 2k_i^2 - (1 - \alpha)G_i ] .
    \f]

    @param k_array: K values (equation 6)
    @param g_array: G values (equation 8)
    @param kernel: which kernel to use ("I", "EXP", "INV")
    @param beta: the inverse temperature.

    @return: a dictionary mapping tuples of atoms to a stress value.
    '''
    from numpy import zeros, exp
    from scipy.linalg import funm, pinv
    from copy import deepcopy
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

    if kernel == "EXP":
        strain_matrix = funm(strain_matrix, lambda x: exp(beta*x))
    elif kernel == "INV":
        strain_matrix = funm(strain_matrix, lambda x: -1.0/(1.0 - beta*x))
    elif kernel == "HAR":
        lapl = deepcopy(strain_matrix)
        for i in range(0, len(k_array)):
            lapl[i, i] = 0 - sum(lapl[i, :])
        strain_matrix = beta * pinv(-1.0 * lapl)

    strain_dict = {}
    for i in range(0, len(k_array)):
        for neigh in fullerene.connectivity[i]:
            if (neigh, i) not in strain_dict:
                strain_dict[(i, neigh)] = strain_matrix[i, neigh]
                if kernel == "INV":
                    strain_dict[(i, neigh)] = 1.0/strain_dict[(i, neigh)]

    return strain_dict
    