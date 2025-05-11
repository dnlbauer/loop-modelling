from re import I
import numpy as np
from numpy.testing import assert_array_equal
import pytest
from biotite.structure import AtomArrayStack
from loop_modeller.io import load_sequence, load_structure
from loop_modeller.loop_modelling import fill_gaps_in_structure
from loop_modeller.utils import find_missing_residues


@pytest.mark.slow
@pytest.mark.parametrize(
    "pdb_id",
    [
        "3IDP",
    ]
)
def test_model_loop(pdb_id: str):
    # load input structure
    structure = load_structure(pdb_id)
    if isinstance(structure, AtomArrayStack):
        structure = structure[0]

    sequences = load_sequence(pdb_id)
    gaps = find_missing_residues(structure, sequences)
    filled_structure = fill_gaps_in_structure(structure, gaps)

    # Make sure that all gaps are filled in the resulting structure
    gaps_after_filling = find_missing_residues(filled_structure, sequences)
    assert len(gaps_after_filling) == 0

    # Make sure that the chain ids did not change
    assert_array_equal(
        np.unique(structure.chain_id),
        np.unique(filled_structure.chain_id),
        err_msg="Chain ids should not change after filling gaps"
    )

    # Make sure that the residue ids did not change
    # Only compare the first and last residue ids, because the filled structure
    # is different from the original structure in the middle
    assert structure.res_id[0] == filled_structure.res_id[0]
    assert structure.res_name[0] == filled_structure.res_name[0]
    assert structure.res_id[-1] == filled_structure.res_id[-1]
    assert structure.res_name[-1] == filled_structure.res_name[-1]