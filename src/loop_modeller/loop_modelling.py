import numpy as np
import biotite.interface.openmm as openmm_interface
import openmm
from openmm.app.pdbfile import PDBFile
from biotite.structure import AtomArray
import biotite.structure.io.pdb as pdb
from pdbfixer import PDBFixer
from tempfile import NamedTemporaryFile


def fill_gaps_in_structure(structure: AtomArray, gaps: dict[int, tuple[int, list[str]]], random_seed: int|None = None) -> AtomArray:
    """
    Fills gaps in a protein structure by adding missing residues.
    This function uses PDBFixer to add missing residues to a protein structure
    based on the provided gaps.

    Args:
        structure (AtomArray): The protein structure to fill gaps in.
        gaps (dict[int, tuple[int, list[str]]]): A dictionary mapping chain IDs
            to tuples, where each tuple contains the index of a gap and a list
            of amino acid residues to be added at that position.
        random_seed (int, optional): A random seed for the structure generation.

    Returns:
        AtomArray: The protein structure with gaps filled.
    """

    # PDBFixer requires a file to read the structure from,
    # so we write it to a temporary file and read it back
    with NamedTemporaryFile(delete=True, suffix=".pdb") as temp_file:
        file = pdb.PDBFile()
        file.set_structure(structure)
        file.write(temp_file.name)
        temp_file.flush()
        fixer = PDBFixer(temp_file.name)


    # PDBFixer addMissingResidues() adds all missing residues (including termini)
    # Therefore, we add our own calculated gaps
    missing_residues = {}
    
    # chain ids might not be sorted in the structure
    unique_chain_indeces = np.unique(structure.chain_id, return_index=True)[1]
    chain_ids = structure.chain_id[np.sort(unique_chain_indeces)]
    for chain_idx, chain_id in enumerate(chain_ids):
        if chain_id in gaps:
            for gap_idx, gap_residues in gaps[chain_id]:
                missing_residues[(chain_idx, gap_idx)] = gap_residues

    # second_fixer = PDBFixer("structures/3IDP.mmcif")
    # second_fixer.findMissingResidues()
    # second_fixer.findMissingAtoms()

    fixer.missingResidues = missing_residues
    # we are only interested in filling the missing residues
    fixer.missingAtoms = {}
    fixer.missingTerminals = {}

    # Run fixer to fill missing gaps
    fixer.addMissingAtoms(seed=random_seed)

    # Parse output back to AtomArray
    structure = openmm_interface.from_topology(fixer.topology)
    structure.coord = np.array(fixer.positions.value_in_unit(openmm.unit.angstrom))
    return structure
