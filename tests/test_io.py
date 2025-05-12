import pytest
from loop_modeller.io import load_sequence, load_structure
import numpy as np

@pytest.mark.parametrize(
    "pdb_id, num_structures, num_chains",
    [
        ("3IDP", 1, 2),
        ("8RX0", 1, 9),
        ("6X18", 1, 6),
    ],
)
def test_load_structure(pdb_id: str, num_structures: int, num_chains: int):
    structure = load_structure(pdb_id)

    assert structure is not None
    assert structure.shape[0] == num_structures
    assert np.unique(structure.chain_id).size == num_chains

@pytest.mark.parametrize(
    "pdb_id, expected_sequences_startswith",
    [
        ("3IDP", {
            "A": "HHHHHHDRNRM",
            "B": "HHHHHHDRNRM"
            }
        ),
        ("6X18", {
            "A": "MGCL",
            "B": "QSEL",
            "G": "NTASI",
            "N": "QVQL",
            "P": "HAEG",
            "R": "MKTII"
            }
        ),
    ],
)
def test_load_sequence(pdb_id: str, expected_sequences_startswith: dict[str, str]):
    seq = load_sequence(pdb_id)
    assert len(seq) == len(expected_sequences_startswith)
    for chain_id, expected_start in expected_sequences_startswith.items():
        assert str(seq[chain_id]).startswith(expected_start)
