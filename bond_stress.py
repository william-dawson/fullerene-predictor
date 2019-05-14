''' @package driver which computes the stress on a bond
'''
from argparse import ArgumentParser
from fullerene_curvature.curvature import compute_k_values, compute_g_values, \
    compute_bond_stress
from fullerene_curvature.fullerene import Fullerene

if __name__ == "__main__":
    # Setup the argument parser
    parser = ArgumentParser(description="Process the input parameters")
    parser.add_argument("filename", help="the name of the input file")
    parser.add_argument("num_bonds", help="the number of bonds to print",
                        type=int)
    parser.add_argument("--plot", help="whether to plot or not",
                        action="store_true")
    parser.add_argument("--kernel", help="which kernel to use", default="I")
    parser.add_argument("--beta", help="the scaling factor", type=float,
                        default=1.0)
    args = parser.parse_args()

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
        print(largest[i], ":", stressdict[largest[i]])

    if args.plot:
        from matplotlib import pyplot as plt
        sval = sorted(stressdict.values(), reverse=True)
        plt.plot(sval, 'bx--')
        plt.plot(sval[:args.num_bonds], 'ro--',
                 label=str(args.num_bonds) + " largest")
        plt.legend()
        plt.show()
