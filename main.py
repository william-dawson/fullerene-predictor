''' @package main driver for the fullerene energy computer.
'''
from sys import argv
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_energy, compute_energy_novel
from fullerene_curvature.fullerene import Fullerene

if __name__ == "__main__":
    file_name = argv[1]

    input_fullerene = Fullerene(file_name)

    k_values = compute_k_values(input_fullerene)
    g_values = compute_g_values(input_fullerene)

    energy_value = compute_energy(k_values, g_values)
    energy_value_novel = compute_energy_novel(k_values, g_values)

    print("Energy:", energy_value)
    print("Energy_Novel:", energy_value_novel)
