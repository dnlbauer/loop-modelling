from biotite.structure import AtomArrayStack
from loop_modeller.io import load_structure
from loop_modeller.scoreing import clash_score


def test_energy_score():
    structure = load_structure("3IDP")
    if isinstance(structure, AtomArrayStack):
        structure = structure[0]

    score = clash_score(structure)

    # Create an obvious clash:
    structure.coord[0] = structure.coord[1] + 0.1
    score_with_clash = clash_score(structure)

    # since we have one additional clash, the score should increase
    # by a multiple of 1/(1000*number_of_atoms)
    # depending on the number of produced clashes
    expected_increase_per_clash = 1 / (1000 * len(structure.coord))

    match = False
    for i in range(1, 5):
        expected_score = score + i * expected_increase_per_clash
        # make sure the expected score is not too far from the actual score
        # by 10 digits.
        if abs(score_with_clash - expected_score) < 1e-10:
            match = True
            break 
    assert match, f"Expected score: {expected_score}, Actual score: {score_with_clash}"

