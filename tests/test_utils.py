from loop_modeller.io import load_sequence, load_structure
from loop_modeller.utils import find_missing_residues

def test_find_missing_residues():
    structure = load_structure("3IDP")[0]
    sequences = load_sequence("3IDP")

    gaps = find_missing_residues(structure, sequences)
    assert len(gaps["A"]) == 2
    assert gaps["A"][0][0] == 151
    assert gaps["A"][0][1][0] == "GLU"
    assert gaps["A"][0][1][-1] == "SER"

    assert gaps["A"][1][0] == 165
    assert gaps["A"][1][1] == ["ASP", "LYS"]

    assert len(gaps["B"]) == 1
    assert gaps["B"][0][0] == 149
    assert gaps["B"][0][1][0] == "ALA"
    assert gaps["B"][0][1][-1] == "LEU"
