''' @package driver which computes the stress on a bond
'''
from sys import argv
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_energy, compute_euler_characteristic, compute_bond_stress
from fullerene_curvature.fullerene import Fullerene

if __name__ == "__main__":
    file_name = argv[1]
    num_bonds = int(argv[2])
    try:
        plot = int(argv[3])
        if plot == 0:
            plot = False
    except:
        plot = False

    input_fullerene = Fullerene(file_name)

    k_values = compute_k_values(input_fullerene)
    g_values = compute_g_values(input_fullerene)

    stressdict = compute_bond_stress(input_fullerene, k_values, g_values)

    # Find the largest values
    largest = sorted(stressdict, key=stressdict.get, reverse=True)[:num_bonds]
    print("Largest values:")
    for i in range(0, num_bonds):
        print(largest[i], ":", stressdict[largest[i]])

    if plot:
        from matplotlib import pyplot as plt
        sval = sorted(stressdict.values(), reverse=True)
        plt.plot(sval, 'bx--')
        plt.plot(sval[:num_bonds], 'ro--', label=str(num_bonds) + " largest")
        plt.legend()
        plt.show()
