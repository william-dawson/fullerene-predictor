''' @package main driver for the fullerene energy computer.
'''
from sys import argv
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_energy, compute_euler_characteristic
from fullerene_curvature.fullerene import Fullerene

if __name__ == "__main__":
    file_name = argv[1]

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
