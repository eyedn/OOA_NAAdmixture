###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_normalize_multik_admixture_rows.py
###############################################################################

# overview: test neutral validation and formatting of multi-K ADMIXTURE rows.


##### set up ##################################################################
from pathlib import Path
import math
import sys
import unittest


PYTHON_UTILS_DIR = Path(__file__).parents[1] / "job_scripts" / "python_utils"
sys.path.insert(0, str(PYTHON_UTILS_DIR))

from onekg_utils.normalize_multik_admixture_rows import (
    normalize_multik_admixture_rows,
)


##### tests ###################################################################
class NormalizeMultikAdmixtureRowsTests(unittest.TestCase):
    def test_preserves_component_order_and_pads_through_five(self):
        values_by_k = {
            2: [0.15, 0.85],
            3: [0.15, 0.25, 0.60],
            4: [0.10, 0.20, 0.30, 0.40],
            5: [0.05, 0.10, 0.15, 0.30, 0.40],
        }

        for k, q_values in values_by_k.items():
            with self.subTest(k=k):
                rows = [
                    {
                        "sample": "sample_1",
                        "pop": "ASW",
                        "q_values": q_values,
                    }
                ]

                normalized = normalize_multik_admixture_rows(rows, k)

                self.assertEqual(len(normalized), 1)
                output_row = normalized[0]
                component_columns = [
                    f"component_{index}_q" for index in range(1, 6)
                ]
                self.assertTrue(
                    all(column in output_row for column in component_columns)
                )
                self.assertEqual(
                    [
                        output_row[f"component_{index}_q"]
                        for index in range(1, k + 1)
                    ],
                    q_values,
                )
                self.assertEqual(
                    [
                        output_row[f"component_{index}_q"]
                        for index in range(k + 1, 6)
                    ],
                    ["NA"] * (5 - k),
                )
                self.assertAlmostEqual(
                    sum(
                        output_row[f"component_{index}_q"]
                        for index in range(1, k + 1)
                    ),
                    1.0,
                )

    def test_rejects_k_outside_supported_range(self):
        for k in (1, 6):
            with self.subTest(k=k):
                with self.assertRaisesRegex(
                    ValueError,
                    "between 2 and 5",
                ):
                    normalize_multik_admixture_rows([], k)

    def test_rejects_incorrect_component_count(self):
        rows = [{"sample": "sample_1", "q_values": [0.2, 0.3, 0.5]}]

        with self.assertRaisesRegex(ValueError, "exactly 2 components"):
            normalize_multik_admixture_rows(rows, 2)

    def test_rejects_nonfinite_components(self):
        for invalid_value in (math.nan, math.inf, -math.inf):
            with self.subTest(invalid_value=invalid_value):
                rows = [
                    {
                        "sample": "sample_1",
                        "q_values": [invalid_value, 0.5],
                    }
                ]

                with self.assertRaisesRegex(
                    ValueError,
                    "valid proportions",
                ):
                    normalize_multik_admixture_rows(rows, 2)

    def test_rejects_components_outside_unit_interval(self):
        for q_values in ([-0.1, 1.1], [1.1, -0.1]):
            with self.subTest(q_values=q_values):
                rows = [{"sample": "sample_1", "q_values": q_values}]

                with self.assertRaisesRegex(
                    ValueError,
                    "valid proportions",
                ):
                    normalize_multik_admixture_rows(rows, 2)

    def test_rejects_components_that_do_not_sum_to_one(self):
        rows = [{"sample": "sample_1", "q_values": [0.2, 0.7]}]

        with self.assertRaisesRegex(ValueError, "valid proportions"):
            normalize_multik_admixture_rows(rows, 2)


if __name__ == "__main__":
    unittest.main()
