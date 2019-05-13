''' @package main driver for the fullerene energy computer.
'''
from sys import argv
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_energy, compute_euler_characteristic, compute_bond_stress
from fullerene_curvature.fullerene import Fullerene

if __name__ == "__main__":
    file_name = argv[1]
    try:
        site_output = argv[2]
    except:
        site_output = None

    try:
        input_fullerene = Fullerene(file_name)

        k_values = compute_k_values(input_fullerene)
        g_values = compute_g_values(input_fullerene)

        euler_characteristic = compute_euler_characteristic(g_values)
        energy_value = compute_energy(k_values, g_values)

        print("Energy:", energy_value)
        print("Euler_Characteristic:", euler_characteristic)

    except:
        print("Failed:", file_name)

    if site_output:
        site_values = compute_bond_stress(input_fullerene, k_values, g_values)
        if site_output == "Print":
            print(site_values)
        elif site_output == "Plot":
            pass
