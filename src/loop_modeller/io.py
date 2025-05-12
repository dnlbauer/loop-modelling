from pathlib import Path
from biotite.sequence import ProteinSequence
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io import pdbx
from biotite.database import rcsb

def load_structure(pdb_id_or_path: str|Path) -> AtomArrayStack:
    """
    Load a protein structure from the RCSB PDB database.

    Args:
        pdb_id_or_path (str): The PDB ID or path of the structure
            to load.

    Returns:
        AtomArrayStack: The loaded structure.
    """
    if Path(pdb_id_or_path).is_file():
        # If a file path is provided, load the structure from the file
        file = pdbx.CIFFile.read(pdb_id_or_path)
        structure = pdbx.get_structure(file)
    else:
        # Fetch the PDB file from the RCSB PDB database and load it
        pdb_file = rcsb.fetch(pdb_id_or_path, format="mmcif", target_path="structures")
        file = pdbx.CIFFile.read(pdb_file)
        structure = pdbx.get_structure(file)

    # If the structure is a single AtomArray, convert it to an AtomArrayStack
    if isinstance(structure, AtomArray):
        structure = AtomArrayStack(structure)

    return structure

def load_sequence(pdb_id_or_path: str|Path) -> dict[str, ProteinSequence]:
    """
    Load the sequence of a protein from the RCSB PDB database.

    Args:
    pdb_id_or_path (str): The PDB ID or path of the sequence to load.

    Returns
    -------
    str
        The sequence of the protein structure.
    """

    if Path(pdb_id_or_path).is_file():
        # If a file path is provided, load the sequence from the file
        file = pdbx.CIFFile.read(pdb_id_or_path)
        return pdbx.get_sequence(file)
        file = pdbx.CIFFile.read(pdb_id_or_path)
    else:
        # Fetch the PDB file from the RCSB PDB database and load it
        pdb_file = rcsb.fetch(pdb_id_or_path, format="mmcif", target_path="structures")
        file = pdbx.CIFFile.read(pdb_file)
    return pdbx.get_sequence(file)
