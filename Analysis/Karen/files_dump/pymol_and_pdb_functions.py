def start_pymol():
    import sys
    import pymol
    pymol.pymol_argv = ['pymol','-qc'] #+ sys.argv[1😏
    stdout = sys.stdout
    stderr = sys.stderr
    pymol.finish_launching()
    sys.stdout = stdout
    sys.stderr = stderr
    cmd = pymol.cmd


def calc_dist_matrix(chain1, chain2, stripping=True, mode='min'):
    """Returns a matrix of C-alpha distances between two chains"""
    chain1_residues = rstrip_non_aa_residues(list(chain1.get_residues()))
    chain2_residues = rstrip_non_aa_residues(list(chain2.get_residues()))
    answer = np.zeros((len(chain1_residues), len(chain2_residues)), np.float)
    for row, residue_one in enumerate(chain1_residues):
        for col, residue_two in enumerate(chain2_residues) :
            answer[row, col] = get_distance_between_residues(residue_one, residue_two, mode=mode)
    return answer


def get_minimal_distances_in_a_complex(chains):
    distances = []
    chains = [chain for chain in chains if chain.id in ascii_letters]
    for chain1 in chains:
        for chain2 in chains:
            distances.append(calc_dist_matrix(chain1, chain2))
    return np.amin(np.array(distances), axis=0)

def strip_non_aa_residues(residues):
    for index, residue in enumerate(residues):
        if residue.resname not in aa3:
            pass
        else:
            break
    residues = residues[index:]
    return residues


def rstrip_non_aa_residues(residues):
    for index, residue in enumerate(residues):
        if residue.resname not in aa3:
            break
    residues = residues[:index]
    return residues