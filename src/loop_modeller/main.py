from importlib import metadata
import json
import sys

from biotite.setup_ccd import CIFFile
from biotite.structure import AtomArrayStack
from biotite.structure.io import pdbx

from loop_modeller.io import load_structure, load_sequence
from loop_modeller.loop_modelling import fill_gaps_in_structure
from loop_modeller.scoreing import clash_score
from loop_modeller.utils import find_missing_residues

def main() -> None:
    # Pass command line arguments
    if len(sys.argv) < 2:
        print("Usage: lm <file or pdb id> <number of models to produce>")
        sys.exit(1)
    
    file_or_pdb_id = sys.argv[1]
    if len(sys.argv) > 2:
        number_of_models = int(sys.argv[2])
    else:
        number_of_models = 1

    # Load the structure and sequence
    structure = load_structure(file_or_pdb_id)
    sequences = load_sequence(file_or_pdb_id)
    if isinstance(structure, AtomArrayStack):
        structure = structure[0]

    # Find gaps in the structure
    gaps = find_missing_residues(structure, sequences)
    if len(gaps) == 0:
        print("No gaps found in the structure.")
        sys.exit(0)

    # Model missing residues to structures and scrore them 
    fixed_structures = []
    scores = []
    for i in range(number_of_models):
        # Fill gaps in the structure
        fixed_structure = fill_gaps_in_structure(structure, gaps, random_seed=i)
        fixed_structures.append(fixed_structure)
        score = clash_score(fixed_structure)
        scores.append(score)
        print(f"Model {i+1} score: ", score)

    # Write structures to disk
    for i, structure in enumerate(fixed_structures):
        output_path = f"model_{i+1}.mmcif"
        file = CIFFile()
        pdbx.set_structure(file, structure)
        with open(output_path, "w") as out:
            file.write(out)
        print(f"Model {i+1} written to {output_path}")

    # Output a file with the scores and modeled residues
    model_metadata = []
    for i, score in enumerate(scores):
        model_metadata.append({
            "model_number": i + 1,
            "file_name": f"model_{i+1}.mmcif",
            "score": score,
        })
    metadata = {
        "input": file_or_pdb_id,
        "output": model_metadata,
        "residues": gaps
    }
    json.dump(metadata, open("metadata.json", "w"), indent=4)
