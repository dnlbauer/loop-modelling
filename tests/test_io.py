import pytest
from loop_modeller.io import load_structure
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