import numpy as np

from .divider import Divider


def check_divider(L: int, O: int):
    """
    Asserts that the fragment generation function:
        * Generates fragments of equal size
        * Covers the entirety of the full sequence
        * Does not incur into unintended overlaps
    """
    divider = Divider(L=L, O=O, nmer=2)

    # generate a test sequence
    seqAA = "".join(list("ACDEFGHIKL"))
    sequence = seqAA * 4 + seqAA[::-1]

    print(f"Testing fragment generation with L={L} and O={O}")
    divider._sanitize_parameters(sequence)

    # visually check that size of all fragments are identical
    list_segment_idx = divider._generate_fragment_idx(sequence)
    print(list_segment_idx)

    # check that the fragments cover 100% of the sequence
    vec_test = np.zeros(len(sequence))

    for i, j in list_segment_idx:
        segment = sequence[i:j]
        print(" " * i + segment)

        # check that each fragment is of the correct size
        assert len(segment) == L
        vec_test[i:j] = 1

    # check that there is no *unintended* overlap
    assert vec_test.sum() == len(sequence)
    print(f"No errors found with L={L} and O={O}")


def run_divider_tests():
    """
    Tests the sequence divider with several combinations of L (length) and O (overlap) values
    Some of these combinations should result in a raised error
    """
    list_parameters = [
        ["normal operation", "ok", [20, 10]],
        ["normal operation", "ok", [15, 10]],
        ["normal operation", "ok", [10, 5]],
        ["normal operation", "ok", [10, 0]],
        ["big L limit (no error expected)", "ok", [49, 10]],
        ["big L exception", "error", [50, 0]],
        ["small L exception", "error", [7, 5]],
        ["big O exception", "error", [23, 23]],
    ]

    list_results = []

    for concept, expected_result, [L, O] in list_parameters:
        print("-" * 60)
        print(f"Running test: {concept}\n")
        result = None

        try:
            check_divider(L, O)
            result = "ok"
        except Exception as e:
            print("Error was found")
            print(e)
            result = "error"

        test_status = expected_result == result
        list_results.append(test_status)

    print("#########################################\n")
    print(f"{sum(list_results)} out of {len(list_results)} tests were succesful")


if __name__ == "__main__":
    run_divider_tests()
