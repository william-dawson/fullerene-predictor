''' @package driver which computes the stress on a bond
'''
from __future__ import print_function
from argparse import ArgumentParser
from fullerene_curvature.curvature import compute_k_values, compute_g_values
from fullerene_curvature.stress import compute_bond_stress
from fullerene_curvature.fullerene import Fullerene


def visualize(mol, stressdict, largest, check):
    import matplotlib.pyplot as plt
    from matplotlib import cm

    cmap = cm.get_cmap('cividis')
    left = min(stressdict.values())
    lr = max(stressdict.values()) - left
    scaledict = {x: (y - left) / lr for x, y in stressdict.items()}

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')

    if check:
        checkval = [int(check[0]) - 1, int(check[1]) - 1]
        x_values = [mol.atoms_array[checkval[0]]
                    [0], mol.atoms_array[checkval[1]][0]]
        y_values = [mol.atoms_array[checkval[0]]
                    [1], mol.atoms_array[checkval[1]][1]]
        z_values = [mol.atoms_array[checkval[0]]
                    [2], mol.atoms_array[checkval[1]][2]]
        ax1.plot3D(x_values, y_values, z_values,
                   'k', label="check", linewidth=6)

    for link in largest:
        x_values = [mol.atoms_array[link[0]][0], mol.atoms_array[link[1]][0]]
        y_values = [mol.atoms_array[link[0]][1], mol.atoms_array[link[1]][1]]
        z_values = [mol.atoms_array[link[0]][2], mol.atoms_array[link[1]][2]]
        ax1.plot3D(x_values, y_values, z_values, 'r', linewidth=6,
                   label=str(len(largest)) + " largest")

    for ring in mol.ring_list:
        for i in range(1, len(ring)):
            vals = (mol.atoms_array[ring[i - 1]], mol.atoms_array[ring[i]])
            x_values = [vals[0][0], vals[1][0]]
            y_values = [vals[0][1], vals[1][1]]
            z_values = [vals[0][2], vals[1][2]]
            idx = (ring[i - 1], ring[i])
            if idx in scaledict:
                color = cmap(scaledict[idx])
            else:
                color = cmap(scaledict[(idx[1], idx[0])])
            ax1.plot3D(x_values, y_values, z_values,
                       color=color, marker="o", linestyle="--")
        vals = (mol.atoms_array[ring[-1]], mol.atoms_array[ring[0]])
        x_values = [vals[0][0], vals[1][0]]
        y_values = [vals[0][1], vals[1][1]]
        z_values = [vals[0][2], vals[1][2]]
        idx = (ring[-1], ring[0])
        if idx in scaledict:
            color = cmap(scaledict[idx])
        else:
            color = cmap(scaledict[(idx[1], idx[0])])
        ax1.plot3D(x_values, y_values, z_values,
                   marker="o", color=color, linestyle="--")

    ax1.legend()
    plt.show()


if __name__ == "__main__":
    # Setup the argument parser
    parser = ArgumentParser(description="Process the input parameters")
    parser.add_argument("filename", help="the name of the input file")
    parser.add_argument("num_bonds", help="the number of bonds to print",
                        type=int)
    parser.add_argument("--plot", help="whether to plot or not",
                        action="store_true")
    parser.add_argument("--plot3d", help="whether to plot or not",
                        action="store_true")
    parser.add_argument("--kernel", help="which kernel to use", default="I",
                        choices=["INV", "EXP", "HAR", "I"])
    parser.add_argument("--beta", help="the scaling factor", type=float,
                        default=1.0)
    parser.add_argument("--check", help="the actual bond to check against",
                        action='append', default=None)
    args = parser.parse_args()

    if args.check and len(args.check) < 2:
        print("You need to specify two check values")
        quit()

    input_fullerene = Fullerene(args.filename)

    k_values = compute_k_values(input_fullerene)
    g_values = compute_g_values(input_fullerene)

    stressdict = compute_bond_stress(
        input_fullerene, k_values, g_values, args.kernel, args.beta)

    # Find the largest values
    largest = sorted(stressdict, key=stressdict.get,
                     reverse=True)[:args.num_bonds]
    print("Largest values:")
    for i in range(0, args.num_bonds):
        print(str(largest[i][0] + 1) + "," +
              str(largest[i][1] + 1), ":", stressdict[largest[i]])

    if args.check:
        check1 = int(args.check[0]) - 1
        check2 = int(args.check[1]) - 1
        checkval = stressdict[(check1, check2)]
        print("Check value:", checkval)

    if args.plot:
        from matplotlib import pyplot as plt
        sval = sorted(stressdict.values(), reverse=True)
        plt.plot(sval, 'bx--')
        plt.plot(sval[:args.num_bonds], 'ro--',
                 label=str(args.num_bonds) + " largest")
        if args.check:
            plt.hlines(checkval, 0, len(sval), label="check",
                       linestyle="dotted", colors="y")
        plt.legend()
        plt.show()

    if args.plot3d:
        visualize(input_fullerene, stressdict, largest, args.check)
