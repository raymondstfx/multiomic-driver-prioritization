import numpy as np

from multiomic_driver.effects.background_corrected import difference_in_differences


def test_difference_in_differences() -> None:
    result = difference_in_differences(
        treated_perturbed=10,
        treated_control=4,
        untreated_perturbed=5,
        untreated_control=3,
    )
    assert result == 4


def test_vector_difference_in_differences() -> None:
    result = difference_in_differences(
        np.array([10, 2]), np.array([4, 1]), np.array([5, 4]), np.array([3, 2])
    )
    np.testing.assert_array_equal(result, np.array([4, -1]))

