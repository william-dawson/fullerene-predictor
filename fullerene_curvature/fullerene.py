''' @package fullerene.py
A class that stores information about a Fullerene.
'''
import copy
from fullerene_curvature.process_input import process_atoms_array, \
    process_connectivity, process_5_rings, process_6_rings
import numpy


class Fullerene:
    '''!
    Stores a description of a given Fullerene system.
    '''

    def __init__(self, file_name):
        '''!
        __init__ initialize a Fullerene object from file.

        self: Fullerne to initialize.
        file_name: Filename to initialize from.
        '''
        ## List of Atoms.
        self.atoms_array = process_atoms_array(file_name)

        ## Which atoms are connected to which other atoms.
        self.connectivity = process_connectivity(file_name,
                                                 len(self.atoms_array))
        fiverings, fiverings_center = process_5_rings(file_name)
        sixrings, sixrings_center = process_6_rings(file_name)

        ## List of rings in the fullerne
        self.ring_list = copy.copy(fiverings)
        self.ring_list.extend(sixrings)

        ## Where the centers of each ring are.
        self.ring_center = copy.copy(fiverings_center)
        self.ring_center.extend(sixrings_center)

        ## Lookup the rings associated with an atom
        self.ring_lookup = []
        for i in range(0, len(self.atoms_array)):
            self.ring_lookup.append([])
        for i in range(0, len(self.ring_list)):
            for atom in self.ring_list[i]:
                self.ring_lookup[atom].append(i)
        ## Which rings are connected to which other rings.
        self.rings_connectivity = compute_ring_connectivity(self.ring_list,
                                                            self.ring_lookup)

    def get_adjacency_matrix(self):
        '''
        get_adjacency_matrix get an adjacency matrix representation of the
        fullerene.

        @param self: Fullerene to process

        return: an adjacency matrix
        '''
        matrix_dimension = len(self.atoms_array)
        matrix = numpy.zeros((matrix_dimension, matrix_dimension))
        for j in range(0, matrix_dimension):
            for i in range(0, len(self.connectivity[j])):
                index = self.connectivity[j][i]
                matrix[j, index] = 1
        return matrix


def compute_ring_connectivity(rings, rings_lookup):
    '''
    compute_ring_connectivity determine which rings are connected to each other.

    @param rings list of rings in the
    '''
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
