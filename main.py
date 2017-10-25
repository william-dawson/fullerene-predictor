##########################################################################
from sys import argv
from curvature import compute_k_values, compute_g_values, compute_energy
from process_input import process
from sphere import compute_sphere
import copy
import numpy

##########################################################################
if __name__ == "__main__":
    file_name = argv[1]

    atoms_array, connectivity, fiverings, sixrings, fiverings_center, \
        sixrings_center = process(file_name)

    ring_list = copy.copy(fiverings)
    ring_list.extend(sixrings)

    ring_center_list = copy.copy(fiverings_center)
    ring_center_list.extend(sixrings_center)

    ring_lookup = []
    for i in range(0, len(atoms_array)):
        ring_lookup.append([])
    for i in range(0, len(ring_list)):
        for atom in ring_list[i]:
            ring_lookup[atom].append(i)

    k_values = compute_k_values(atoms_array, connectivity)
    g_values = compute_g_values(
        atoms_array, connectivity, ring_list, ring_center_list, ring_lookup)
    energy_value = compute_energy(k_values, g_values)

    print("Energy:", energy_value)
