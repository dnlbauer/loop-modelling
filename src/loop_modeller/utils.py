import numpy as np
from biotite.structure import AtomArray

def find_gaps_in_structure(structure: AtomArray) -> np.ndarray:
    """
    Find gaps in the structure by checking for missing residues.
    A gap is defined by a jump in the residue id numbering within a chain.
    A change in the chain id is not considered a gap.
    Missing residues at the beginning or end of a chain are not considered gaps.

    Parameters
    ----------
    structure : AtomArray
        The structures to check for gaps.

    Returns
    -------
    np.ndarray
        The indices of the residue before the gap in the structure.
    """
    cas = structure[structure.atom_name == "CA"]
    
    # find jumps in res id
    jumps = np.diff(cas.res_id)

    # mask jumps that are between different chains
    mask = cas.chain_id[1:] != cas.chain_id[0:-1]
    jumps[mask] = 0

    # find the start of each gap
    gap_starts = np.where(jumps > 1)[0]
    return gap_starts