###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_pipeline_shell_contracts.py
###############################################################################

# overview: test explicit environment and inference-tool shell interfaces.


##### set up ##################################################################
from pathlib import Path
import unittest


PROJECT_DIR = Path(__file__).parents[1]
SUBMITTERS = [
    PROJECT_DIR / f"{stage}_submit_{name}.sh"
    for stage, name in (
        (1, "sim_model"),
        (2, "calc_sim_chrom_stats"),
        (3, "calc_sim_genome_stats"),
        (4, "calc_onekg_chrom_stats"),
        (5, "calc_onekg_genome_stats"),
    )
]
WORKERS = sorted((PROJECT_DIR / "job_scripts").glob("*.sh"))
ADMIXTURE_WORKERS = [
    "calc_sim_chr_rep_stats.sh",
    "calc_sim_genome_rep_stats.sh",
    "combine_onekg_chr_stats.sh",
    "combine_onekg_genome_stats.sh",
]
EMPIRICAL_COMBINE_WORKERS = [
    "combine_onekg_chr_stats.sh",
    "combine_onekg_genome_stats.sh",
]


##### tests ###################################################################
class PipelineShellContractTests(unittest.TestCase):
    def test_submitters_pass_configured_environment_to_every_worker(self):
        for submitter in SUBMITTERS:
            text = submitter.read_text(encoding="utf-8")
            worker_count = text.count('"job_scripts/')
            with self.subTest(submitter=submitter.name):
                self.assertGreater(worker_count, 0)
                self.assertEqual(
                    text.count(
                        '"job_scripts/'
                    ),
                    text.count(
                        '"${OOA_NAADMIXTURE_CONDA}"'
                    ),
                )

    def test_workers_activate_the_passed_environment(self):
        for worker in WORKERS:
            text = worker.read_text(encoding="utf-8")
            with self.subTest(worker=worker.name):
                self.assertNotIn(
                    "conda activate OOA_NAAdmixture",
                    text,
                )
                self.assertIn('conda activate "${conda_env}"', text)

    def test_admixture_executable_is_explicitly_routed(self):
        for worker_name in ADMIXTURE_WORKERS:
            text = (
                PROJECT_DIR / "job_scripts" / worker_name
            ).read_text(encoding="utf-8")
            with self.subTest(worker=worker_name):
                self.assertIn('admixture_exec="$2"', text)
                self.assertNotIn("admixture_linux-1.4.0", text)

    def test_empirical_combiners_receive_faststructure_contract(self):
        for worker_name in EMPIRICAL_COMBINE_WORKERS:
            text = (
                PROJECT_DIR / "job_scripts" / worker_name
            ).read_text(encoding="utf-8")
            with self.subTest(worker=worker_name):
                for token in (
                    "faststructure_conda_env",
                    "faststructure_structure_py",
                    "faststructure_choose_k_py",
                    "faststructure_prior",
                    "faststructure_cv",
                    "faststructure_dir",
                    "unsupervised_ks",
                ):
                    self.assertIn(token, text)
                self.assertIn(
                    'conda activate "${faststructure_conda_env}"',
                    text,
                )
                self.assertIn("--prior=", text)
                self.assertIn("--cv=", text)
                self.assertNotIn("--full", text)
                self.assertNotIn("--tol=", text)
                self.assertEqual(
                    text.count(
                        'python "${faststructure_choose_k_py}"'
                    ),
                    1,
                )
                for suffix in ("meanQ", "meanP", "log"):
                    self.assertIn(suffix, text)

        for stage in (4, 5):
            submitter = next(
                path for path in SUBMITTERS
                if path.name.startswith(f"{stage}_")
            )
            text = submitter.read_text(encoding="utf-8")
            with self.subTest(submitter=submitter.name):
                for constant in (
                    "FASTSTRUCTURE_CONDA_ENV",
                    "FASTSTRUCTURE_STRUCTURE_PY",
                    "FASTSTRUCTURE_CHOOSE_K_PY",
                    "FASTSTRUCTURE_PRIOR",
                    "FASTSTRUCTURE_CV",
                    "ONEKG_OUT_FASTSTRUCTURE_DIR",
                    "ONEKG_UNSUPERVISED_KS",
                ):
                    self.assertIn(f'"${{{constant}', text)

    def test_empirical_seeds_are_chromosome_and_zero(self):
        chrom_text = (
            PROJECT_DIR / "job_scripts" / "combine_onekg_chr_stats.sh"
        ).read_text(encoding="utf-8")
        genome_text = (
            PROJECT_DIR / "job_scripts" / "combine_onekg_genome_stats.sh"
        ).read_text(encoding="utf-8")
        self.assertIn('faststructure_seed="${chr}"', chrom_text)
        self.assertIn('faststructure_seed="0"', genome_text)

    def test_environment_files_are_split_and_cython_is_pinned(self):
        primary_path = PROJECT_DIR / "environment_OOA_NAAdmixture.yml"
        faststructure_path = PROJECT_DIR / "environment_fastStructure.yml"

        self.assertFalse((PROJECT_DIR / "environment.yml").exists())
        self.assertIn(
            "name: OOA_NAAdmixture",
            primary_path.read_text(encoding="utf-8"),
        )
        faststructure_text = faststructure_path.read_text(encoding="utf-8")
        for dependency in (
            "name: fastStructure",
            "python=2.7",
            "numpy",
            "scipy",
            "cython=0.24",
            "gsl",
        ):
            self.assertIn(dependency, faststructure_text)


if __name__ == "__main__":
    unittest.main()
