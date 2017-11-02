''' @package process_input helpers to build a fullerene from file.
'''
# distance_unit = 0.529177
distance_unit = 1.0

def process_atoms_array(file_name):
    '''
    process_atoms_array build a list of atoms.

    file_name: name of the file to build from.

    returns: list of atoms.
    '''
    atoms_array = []
    with open(file_name, 'r') as ifile:
        found_size = False
        while (not found_size):
            temp = next(ifile)
            if "&General" in temp:
                found_size = True
        split = temp.split()
        number_of_atoms = int(split[1].split("=")[1])
        found_coords = False
        while(not found_coords):
            temp = next(ifile)
            if "New Coordinates:" in temp:
                found_coords = True
        temp = next(ifile)
        for i in range(0, number_of_atoms):
            temp = next(ifile)
            split = temp.split()
            value_list = split[3:6]
            value_list[0] = value_list[0].replace('D', 'E')
            value_list[1] = value_list[1].replace('D', 'E')
            value_list[2] = value_list[2].replace('D', 'E')
            atoms_array.append([distance_unit * float(i) for i in value_list])
    return atoms_array


def process_connectivity(file_name, number_of_atoms):
    '''!
    process_connectivity: determine what atoms are connected to each other.

    file_name: name of the file to build from.
    number_of_atoms: number of atoms in the Fullerene.

    return: for each atom, which other atom is it connected to.
    '''
    connectivity = []
    with open(file_name, 'r') as ifile:
        found_conn = False
        while (not found_conn):
            temp = next(ifile)
            if "Calculate all" in temp:
                found_conn = True
        temp = next(ifile)
        for i in range(0, number_of_atoms):
            temp = next(ifile)
            split = temp.split()
            value_list = split[2:5]
            value_list[2] = value_list[2][:-1]
            connectivity.append([int(i) - 1 for i in value_list])

    return connectivity


def process_5_rings(file_name):
    '''!
    process_5_rings: read in information about the rings with 5 atoms.

    file_name: name of the file to build from.

    return: a list of rings and the contained atoms, centers of the rings.
    '''
    fiverings = []
    fiverings_center = []
    with open(file_name, 'r') as ifile:
        found_ring = False
        while (not found_ring):
            try:
                temp = next(ifile)
            except:
                number_of_five_rings = 0
                exit
            if "five-membered-rings" in temp:
                found_ring = True
        split = temp.split()
        number_of_five_rings = int(split[0])

    if number_of_five_rings > 0:
        with open(file_name, 'r') as ifile:
            found_ring = False
            while (not found_ring):
                temp = next(ifile)
                if "Center for 5-rings" in temp:
                    found_ring = True
            temp = next(ifile)
            temp = next(ifile)
            for i in range(0, number_of_five_rings):
                temp = next(ifile)
                split = temp.split()
                value_list = split[1:6]
                fiverings.append([int(i) - 1 for i in value_list])
                value_list = split[6:9]
                value_list[0] = value_list[0].replace('D', 'E')
                value_list[1] = value_list[1].replace('D', 'E')
                value_list[2] = value_list[2].replace('D', 'E')
                fiverings_center.append(
                    [distance_unit * float(i) for i in value_list])
    return fiverings, fiverings_center


def process_6_rings(file_name):
    '''!
    process_6_rings: read in information about the rings with 6 atoms.

    file_name: name of the file to build from.

    return: a list of rings and the contained atoms, centers of the rings.
    '''
    sixrings = []
    sixrings_center = []
    with open(file_name, 'r') as ifile:
        found_ring = False
        while (not found_ring):
            try:
                temp = next(ifile)
            except:
                number_of_six_rings = 0
                exit
            if "six-membered-rings" in temp:
                found_ring = True
        split = temp.split()
        number_of_six_rings = int(split[0])
    if number_of_six_rings > 0:
        with open(file_name, 'r') as ifile:
            found_ring = False
            while (not found_ring):
                temp = next(ifile)
                if "Center for 6-rings" in temp:
                    found_ring = True
            temp = next(ifile)
            temp = next(ifile)
            for i in range(0, number_of_six_rings):
                temp = next(ifile)
                split = temp.split()
                value_list = split[1:7]
                sixrings.append([int(i) - 1 for i in value_list])
                value_list = split[7:10]
                value_list[0] = value_list[0].replace('D', 'E')
                value_list[1] = value_list[1].replace('D', 'E')
                value_list[2] = value_list[2].replace('D', 'E')
                sixrings_center.append(
                    [distance_unit * float(i) for i in value_list])
    return sixrings, sixrings_center
