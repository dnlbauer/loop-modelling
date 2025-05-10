from biotite.structure import array, Atom

from loop_modeller.utils import find_gaps_in_structure

def test_find_gaps():
    test_structure = array([
        Atom([0,0,0], res_id=3, chain_id="A", atom_name="CA"),
        Atom([0,0,0], res_id=4, chain_id="A", atom_name="CA"),
        Atom([0,0,0], res_id=5, chain_id="A", atom_name="CA"),
        # gap
        Atom([0,0,0], res_id=7, chain_id="A", atom_name="CA"),
        # gap
        Atom([0,0,0], res_id=12, chain_id="A", atom_name="CA"),
        Atom([0,0,0], res_id=13, chain_id="A", atom_name="CA"),
        # no an actual gap because of chain change
        Atom([0,0,0], res_id=15, chain_id="B", atom_name="CA"),
        # gap
        Atom([0,0,0], res_id=23, chain_id="B", atom_name="CA"),
        Atom([0,0,0], res_id=24, chain_id="B", atom_name="CA"),
    ])
    print(test_structure)
    gaps = find_gaps_in_structure(test_structure)
    assert gaps.shape == (3,)
    assert list(gaps) == [2, 3, 6]