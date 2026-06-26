###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_helper_consolidation.py
###############################################################################

import unittest
from pathlib import Path

import job_scripts.python_utils.sim_utils as sim_utils
import job_scripts.python_utils.stats_utils as stats_utils


PROJECT_DIR = Path(__file__).resolve().parents[1]


# Verify private helper files are folded into their workflow modules.
class HelperConsolidationTest(unittest.TestCase):
    # Verify simulation-only helper files were removed.
    def test_simulation_private_helper_modules_removed(self):
        removed_paths = [
            PROJECT_DIR / "job_scripts/python_utils/sim_utils/build_demography.py",
            PROJECT_DIR / "job_scripts/python_utils/sim_utils/build_metadata.py",
        ]
        self.assertEqual([], [path for path in removed_paths if path.exists()])

    # Verify chromosome-stat helper files were removed.
    def test_chromosome_stats_private_helper_modules_removed(self):
        removed_paths = [
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "build_1d_sfs_rows.py",
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "build_2d_sfs_rows.py",
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "build_ancestry_table.py",
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "build_ld_decay_rows.py",
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "build_pi_theta_rows.py",
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/"
            "sample_nodes_for_pop.py",
        ]
        self.assertEqual([], [path for path in removed_paths if path.exists()])

    # Verify package exports expose only workflow and shared helpers.
    def test_package_exports_hide_private_helpers(self):
        sim_init = (
            PROJECT_DIR / "job_scripts/python_utils/sim_utils/__init__.py"
        ).read_text(encoding="utf-8")
        stats_init = (
            PROJECT_DIR / "job_scripts/python_utils/stats_utils/__init__.py"
        ).read_text(encoding="utf-8")

        self.assertIn('"run_simulation"', sim_init)
        self.assertNotIn("build_demography", sim_init)
        self.assertNotIn("build_metadata", sim_init)

        private_names = [
            "build_1d_sfs_rows",
            "build_2d_sfs_rows",
            "build_ancestry_table",
            "build_ld_decay_rows",
            "build_pi_theta_rows",
            "sample_nodes_for_pop",
        ]
        for private_name in private_names:
            self.assertNotIn(private_name, stats_init)
            self.assertFalse(hasattr(stats_utils, private_name))
        self.assertFalse(hasattr(sim_utils, "build_demography"))
        self.assertFalse(hasattr(sim_utils, "build_metadata"))


if __name__ == "__main__":
    unittest.main()
