from biotite.sequence import ProteinSequence
import numpy as np
from biotite.structure import AtomArray

def find_missing_residues(structure: AtomArray, chain_sequences: dict[str, list[str]]) -> dict[str, list[tuple[int, ProteinSequence]]]:
    """
    Identifies position and names of missing residues in a protein structure
    by comparing the structure's sequence with the reference sequences for each chain.
    
    Args:
        structure (AtomArray): The atomic structure of the protein
        chain_sequences (dict[str, ProteinSequence]): A dictionary mapping chain 
            IDs to their corresponding full amino acid sequences.
    Returns:
        dict[str, list[tuple[int, ProteinSequence]]]: A dictionary where each key 
        is a chain ID mappint to a list of tuples.
        Each tuple contains an integer and a `ProteinSequence` object
        representing the start position of the gap in the sequence
        and the missing residues.
    """

    # create a gapped sequences for each chain from the structure
    gapped_sequences = {}
    for chain in np.unique(structure.chain_id):
        chain_structure = structure[structure.chain_id == chain]
        chain_cas = chain_structure[chain_structure.atom_name == "CA"]

        chain_sequence = []
        prev_idx = -1
        for i in range(len(chain_cas)):
            if prev_idx > 0 and chain_cas.res_id[i] > prev_idx + 1:
                gap_size = chain_cas.res_id[i] - prev_idx - 1
                if gap_size > 0:
                    chain_sequence += [None] * gap_size
            chain_sequence.append(str(chain_cas.res_name[i]))
            prev_idx = chain_cas.res_id[i]

        gapped_sequences[chain] = chain_sequence


    missing_residues = {}
    for chain in np.unique(structure.chain_id):
        # There might be an offset between the sequence and the structure
        # This finds the offset and cuts the sequence by removing the first n residues
        ungapped_sequence = None
        chain_sequence = [ProteinSequence.convert_letter_1to3(letter) for letter in chain_sequences[chain]]
        chain_structure = gapped_sequences[chain]
        for offset in range(len(chain_sequence)):
            found = True
            for seq_res, struct_res in zip(chain_sequence[offset:], chain_structure):
                if struct_res != None and struct_res != seq_res:
                    found = False
                    break
            if found:
                ungapped_sequence = chain_sequence[offset:]
                break
        if ungapped_sequence == None:
            raise ValueError(f"Failed to align reference sequence with structure for chain {chain}.")
            
        # for each gap in the structure,
        # find the corresponding residues in the reference sequence
        gapped_sequence = gapped_sequences[chain]
        missing_residues[str(chain)] = []

        # iterate over the gapped sequence and find the position and length of the gaps
        # for each found gap, add the corresponding residues from the reference sequence
        offset = 0  # offset in the sequence due to gaps
        gap_start = None # start of the gap in the sequence
        gap_size = 0
        for idx, residue in enumerate(gapped_sequence):
            if residue == None:
                if gap_start == None:
                    # new gap
                    gap_start = idx
                    gap_size = 1
                else:
                    # continue gap
                    gap_size += 1
            else:
                if gap_start != None:
                    # gap ended, add the missing residues
                    residues = ungapped_sequence[gap_start:gap_start + gap_size]
                    missing_residues[str(chain)].append((gap_start-offset, residues))
                    offset += gap_size
                    gap_start = None

    return missing_residues