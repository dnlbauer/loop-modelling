import warnings
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import pdbx
from biotite.database import rcsb

def load_structure(pdb_id: str) -> AtomArrayStack:
    """
    Load a protein structure from the RCSB PDB database.

    Parameters
    ----------
    pdb_id : str
        The PDB ID of the structure to load.

    Returns
    -------
    AtomArrayStack
        The loaded structure as an AtomArrayStack object.
    """
    # Fetch the PDB file from the RCSB PDB database and load it
    pdb_file = rcsb.fetch(pdb_id, format="mmcif", target_path="structures")
    file = pdbx.CIFFile.read(pdb_file)
    structure = pdbx.get_structure(file)

    # If the structure is a single AtomArray, convert it to an AtomArrayStack
    if isinstance(structure, AtomArray):
        structure = AtomArrayStack(structure)

    return structure
