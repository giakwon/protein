import subprocess
import matplotlib.pyplot as plt

# directions
DIRS = {
    "U": (0, 1),
    "D": (0, -1),
    "L": (-1, 0),
    "R": (1, 0),
}


def orthogonal_transforms(x, y):
    '''
    Applies all transformation on a point (reflections).

    params:
        x (int): x value on coordinate
        y (int): y value on coorfinate

    return (set): set of all transformations
    '''

    transforms = set()
    points = [
        (x, y),
        (-y, x),
        (-x, -y),
        (y, -x),
    ]
    for px, py in points:
        transforms.add((px, py))
        transforms.add((-px, py))  # reflect over y-axis
        transforms.add((px, -py))  # reflect over x-axis
        transforms.add((-px, -py)) # reflect over origin
    return transforms


def all_structures(n):
    '''
    Finds all possible conformations/structures of a sequence of length n
    based on directions U, R, D, L.

    params:
        n (int): length of sequence

    return (lst): list of all possible conformations

    note: from CodingMaster (Jayden?)), but has one extra "step" 
            and first step always R
    '''

    if n == 0:
        return []
    structures = []
    def backtrack(x=1, y=0, structure=["R"], occupied=set()):
        if len(structure) == n:
            structures.append("".join(structure))
            return
        points = orthogonal_transforms(x, y)
        occupied.update(points)
        for direction, (dx, dy) in DIRS.items():
            if (x + dx, y + dy) in occupied:
                continue
            structure.append(direction)
            backtrack(x + dx, y + dy)
            structure.pop()
        occupied.difference_update(points)
    backtrack()
    return structures


def adjacent(index1, index2, struc):
    '''
    Checks if two residues based on the folded structure. 

    params:
        index1 (int): index of the first residue on the sequence
        index2 (int): index of the second residue on the sequence
        struc (str): the structure of the folded protein

    return (bool): returns True if residues are adjacent, False otherwise
    '''
    coord1 = (0,0)
    for i in range(index1):
        coord1 = tuple(map(lambda x, y: x + y, coord1, DIRS[struc[i]]))

    coord2 = (0,0)
    for i in range(index2):
        coord2 = tuple(map(lambda x, y: x + y, coord2, DIRS[struc[i]]))

    if (abs(coord1[0] - coord2[0]) + abs(coord1[1] - coord2[1]) == 1) \
        and (abs(index1-index2)>1):
        return True 
    
    return False


def num_adj(seq, struc):
    ''''
    Finds the number of total HH contacts in the structure from the sequence.

    params: 
        seq (str): the HP sequence of the protein
        struc (str): the structure of the folded protein

    return (int): returns the number of HH contacts in the structure
    '''
    count = 0
    for i in range(len(seq)): 
        for j in range(len(seq)): 
            if adjacent(i, j, struc) and \
                (seq[i] == 'H') and (seq[j] == 'H'):
                count += 1

    return count//2


def fitness(seq):
    '''
    Fitness function that finds all possible structures of the sequence,
    makes a dictionary with the number 

    param:
        seq (str): the HP sequence of the protein

    return (int): total number of distinct number of HH contacts for the sequence
    '''
    structures = all_structures(len(seq))

    level_diff = {}
    for struc in structures:
        adj = num_adj(seq, struc)
        if level_diff.get(adj) == None: 
            level_diff[adj] = []   
        level_diff[adj].append(struc)

    return len(level_diff)


def enum_dict(seq):
    """
    Parses the output from C file ('../lattice-enum/lattice-enum') and turns it 
    into a dictionary. 

    param:
        seq (str): the HP sequence of the protein

    return (dict): key: total number of HH contacts, value: list of corresponding structures
    """
    process = subprocess.run(['../lattice-enum/lattice-enum'] + [seq], 
                             capture_output=True, text=True, check=True)
    c_output = process.stdout.strip().split('\n')

    conformations = {}
    for line in c_output: 
        kv_pair = line.strip().split(" ")
        key = kv_pair[0]
        value = kv_pair[1]
        if conformations.get(key) == None: 
            conformations[key] = [value]   
        conformations[key].append(value)

    return conformations


def make_histogram(enum_dict):
    counts = {k: len(v) for k, v in enum.items()}

    x_axis = list(counts.keys())
    y_axis = list(counts.values())

    plt.bar(x_axis, y_axis)
    plt.yscale('log')
    plt.ylim(1, 10e6)
    plt.ylabel('Count', fontsize = 13)
    plt.xlabel('Number of HH Contacts in Conformation', fontsize = 13)
    plt.title('High-Fitness Level for an 16-Residue Chain​', fontsize = 16)
    plt.show()

if __name__ == "__main__": 

    # good examples 
    # seq = "HPHPPH" # len=6
    # seq = "HPPHHPHP" # len=8
    seq = "HHHHHHHHPPPPPPPP"
    # seq = "HPHHHHHHPHHPHHPH" # len=16 (in paper)
    # seq = "HPHPPHHPHHPHHPPHPH" # len=18
    # seq = "HHPHHHPHPPHHPHPPHP" # len=18

    # bad example
    # seq = "HHHHHHHHHPPPPPPPPP"

    enum = enum_dict(seq)
    make_histogram(enum)