from biotite.structure import AtomArray
import numpy as np
from scipy.spatial import cKDTree


def clash_score(structure: AtomArray, overlap_cutoff=0.4) -> float:
    """
    Calculate the MolProbity clash score of a structure.
    This function checks for steric clashes between atoms in the structure
    and returns a score based on the number of clashes found.

    A clash is defined as two atoms being closer than a certain threshold
    to each other.

    Reference:
        MolProbity: all-atom contacts and structure validation for proteins and nucleic acids
        https://doi.org/10.1093/nar/gkm216

    Args:
        structure : AtomArray
            The structure to calculate the clash score for.
        overlap_cutoff : float, optional
            The cutoff distance for considering two atoms to be in a clash.
            The default is 0.4.

    Returns:
        float: The clash score of the structure.
    """
    vdw_distances = {
        "C": 1.7,
        "N": 1.55,
        "O": 1.52,
        "H": 1.2,
        "S": 1.8,
    }

    # Find all pairs closer than 2 times the maximum vdw distance
    neighbor_tree = cKDTree(structure.coord)
    pairs = neighbor_tree.query_pairs(2*max(vdw_distances.values()))

    # Calculate the number of clashes.
    clashes = 0
    for i, j in pairs:
        d = np.linalg.norm(structure.coord[i] - structure.coord[j])
        if (structure.element[i] not in vdw_distances or structure.element[j] not in vdw_distances):
            # Ignore elements with no default vdw distance
            continue
        if d < vdw_distances[structure.element[i]] + vdw_distances[structure.element[j]] - overlap_cutoff:
            clashes += 1

    # MolProbity clash score: Number of clashes per 1000 atoms
    return clashes / 1000 / len(structure.coord)

