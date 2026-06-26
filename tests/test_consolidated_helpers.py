###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_consolidated_helpers.py
###############################################################################

import tempfile
import unittest
from pathlib import Path

from job_scripts.python_utils.sim_utils.run_simulation import _build_metadata
from job_scripts.python_utils.stats_utils.calc_stats import _build_1d_sfs_rows
from job_scripts.python_utils.stats_utils.calc_stats import _build_2d_sfs_rows
from job_scripts.python_utils.stats_utils.calc_stats import _build_pi_theta_rows
from job_scripts.python_utils.stats_utils.read_fam_order import read_q_rows


# Small matrix object that mimics the 2D SFS indexing contract.
class FakeSpectrum2D:
    shape = (2, 2)

    # Store deterministic spectrum values for tuple-index lookups.
    def __init__(self):
        self.values = [[0.0, 1.0], [2.0, 3.0]]

    # Return the requested 2D spectrum cell.
    def __getitem__(self, key):
        row_idx, col_idx = key
        return self.values[row_idx][col_idx]


# Minimal tree-sequence-like object for row builder contract tests.
class FakeTreeSequence:
    sequence_length = 100.0

    # Return two haploid nodes for a fake diploid individual.
    def individual(self, ind_id):
        node_start = ind_id * 2
        return type("Individual", (), {"nodes": [node_start, node_start + 1]})

    # Return 1D or 2D spectra based on the requested sample-set count.
    def allele_frequency_spectrum(self, sample_sets, **kwargs):
        if len(sample_sets) == 1:
            return [0.0, 2.0, 1.0]
        return FakeSpectrum2D()

    # Return a deterministic segregating-site count.
    def segregating_sites(self, sample_sets, **kwargs):
        return [3.0]

    # Return a deterministic pi estimate.
    def diversity(self, sample_sets, **kwargs):
        return [0.25]


# Check consolidated helpers preserve supervised labels and row keys.
class ConsolidatedHelpersTest(unittest.TestCase):
    # Verify ADX receives the unsupervised ADMIXTURE label.
    def test_build_metadata_marks_adx_as_unsupervised(self):
        rows = _build_metadata(["AFR", "EUR", "ADX"], 2)
        labels = {row["iid"]: row["supervised_label"] for row in rows}

        self.assertEqual("AFR", labels["AFR_1"])
        self.assertEqual("EUR", labels["EUR_1"])
        self.assertEqual("-", labels["ADX_1"])
        self.assertEqual(6, rows[-1]["original_order"])

    # Verify Q rows follow the filtered FAM order.
    def test_read_q_rows_uses_filtered_fam_order(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fam_path = temp_path / "filtered.fam"
            q_path = temp_path / "filtered.2.Q"
            fam_path.write_text(
                "0 ADX_2 0 0 0 -9\n0 AFR_1 0 0 0 -9\n",
                encoding="utf-8",
            )
            q_path.write_text("0.7 0.3\n0.9 0.1\n", encoding="utf-8")

            rows = read_q_rows(4, q_path, fam_path)

        self.assertEqual("ADX_2", rows[0]["vcf_sample_id"])
        self.assertEqual("ADX", rows[0]["pop"])
        self.assertEqual(0.7, rows[0]["afr_q"])
        self.assertEqual("AFR_1", rows[1]["vcf_sample_id"])

    # Verify stats row builders preserve replicate and chromosome fields.
    def test_stats_builders_keep_rep_and_chrom(self):
        ts = FakeTreeSequence()
        pops = ["AFR", "EUR"]

        one_d = _build_1d_sfs_rows(ts, 7, pops, 1, "22")
        two_d = _build_2d_sfs_rows(ts, 7, pops, 1, "22")
        pi_theta = _build_pi_theta_rows(ts, 7, pops, 1, 2.5e-8, "22")

        self.assertTrue(all(row["rep"] == 7 for row in one_d))
        self.assertTrue(all(row["chrom"] == "22" for row in one_d))
        self.assertTrue(all(row["rep"] == 7 for row in two_d))
        self.assertTrue(all(row["chrom"] == "22" for row in two_d))
        self.assertTrue(all(row["rep"] == 7 for row in pi_theta))
        self.assertTrue(all(row["chrom"] == "22" for row in pi_theta))


if __name__ == "__main__":
    unittest.main()
