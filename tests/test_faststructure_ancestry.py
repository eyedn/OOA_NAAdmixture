###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_faststructure_ancestry.py
###############################################################################

# overview: test empirical ADMIXTURE and fastStructure output contracts.


##### set up ##################################################################
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import os
import subprocess
import sys
import unittest

import polars as pl


PROJECT_DIR = Path(__file__).parents[1]
PYTHON_UTILS_DIR = PROJECT_DIR / "job_scripts" / "python_utils"
sys.path.insert(0, str(PYTHON_UTILS_DIR))

from onekg_utils.parse_faststructure_choose_k import (
    parse_faststructure_choose_k,
)
from onekg_utils.write_stats_table import EMPTY_TABLE_COLUMNS


BUILD_SCRIPT = PYTHON_UTILS_DIR / "build_onekg_ancestry.py"
SUPER_COLUMNS = [
    "rep",
    "chrom",
    "pop",
    "sample_id",
    "vcf_sample_id",
    "afr_tspop",
    "eur_tspop",
    "afr_q",
    "eur_q",
    "span",
]
MULTIK_COLUMNS = [
    "rep",
    "chrom",
    "pop",
    "sample_id",
    "vcf_sample_id",
    "k",
    "component_1_q",
    "component_2_q",
    "component_3_q",
    "component_4_q",
    "component_5_q",
    "span",
]
CHOOSE_K_COLUMNS = [
    "rep",
    "chrom",
    "prior",
    "seed",
    "k_min",
    "k_max",
    "max_marginal_likelihood_k",
    "model_components_k",
]


##### helpers #################################################################
def write_text(path, text):
    path.write_text(text, encoding="utf-8")


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))


##### tests ###################################################################
class FastStructureAncestryTests(unittest.TestCase):
    def setUp(self):
        self.temp_context = TemporaryDirectory()
        self.temp_dir = Path(self.temp_context.name)
        self.input_dir = self.temp_dir / "input"
        self.input_dir.mkdir()
        self.fake_modules_dir = self.temp_dir / "fake_modules"
        self.fake_modules_dir.mkdir()
        self._write_inputs()
        self._write_fake_pandas()

    def tearDown(self):
        self.temp_context.cleanup()

    def _write_inputs(self):
        self.unrels_path = self.input_dir / "unrelated.tsv"
        write_text(
            self.unrels_path,
            "IID\tpopulation\n"
            "sample_1\tYRI\n"
            "sample_2\tCEU\n"
            "sample_3\tASW\n",
        )
        self.source_fam_path = self.input_dir / "source.fam"
        self.inference_fam_path = self.input_dir / "inference.fam"
        fam_text = (
            "ASW sample_3 0 0 0 -9\n"
            "YRI sample_1 0 0 0 -9\n"
            "CEU sample_2 0 0 0 -9\n"
        )
        write_text(self.source_fam_path, fam_text)
        write_text(self.inference_fam_path, fam_text)
        self.supervised_q_path = self.input_dir / "supervised.2.Q"
        write_text(
            self.supervised_q_path,
            "0.40 0.60\n"
            "0.90 0.10\n"
            "0.10 0.90\n",
        )
        values_by_k = {
            2: "0.35 0.65\n0.80 0.20\n0.15 0.85\n",
            3: "0.10 0.30 0.60\n0.70 0.20 0.10\n0.20 0.70 0.10\n",
            4: (
                "0.10 0.20 0.30 0.40\n"
                "0.40 0.30 0.20 0.10\n"
                "0.25 0.25 0.25 0.25\n"
            ),
            5: (
                "0.05 0.10 0.15 0.30 0.40\n"
                "0.40 0.30 0.15 0.10 0.05\n"
                "0.20 0.20 0.20 0.20 0.20\n"
            ),
        }
        self.admixture_paths = {}
        self.faststructure_paths = {}
        for k, values in values_by_k.items():
            admixture_path = self.input_dir / f"admixture.{k}.Q"
            faststructure_path = self.input_dir / f"faststructure.{k}.meanQ"
            write_text(admixture_path, values)
            write_text(faststructure_path, values)
            self.admixture_paths[k] = admixture_path
            self.faststructure_paths[k] = faststructure_path
        self.choose_k_path = self.input_dir / "chooseK.txt"
        write_text(
            self.choose_k_path,
            "Model complexity that maximizes marginal likelihood = 3\n"
            "Model components used to explain structure in data = 2\n",
        )

    def _write_fake_pandas(self):
        write_text(
            self.fake_modules_dir / "pandas.py",
            "import polars as pl\n"
            "\n"
            "class DataFrame:\n"
            "    def __init__(self, table):\n"
            "        self.table = table\n"
            "\n"
            "    def to_parquet(self, path, index=False):\n"
            "        if index:\n"
            "            raise ValueError('index output is unsupported')\n"
            "        self.table.write_parquet(path)\n"
            "\n"
            "def read_csv(path, sep=','):\n"
            "    table = pl.read_csv(\n"
            "        path,\n"
            "        separator=sep,\n"
            "        null_values=['NA'],\n"
            "    )\n"
            "    return DataFrame(table)\n",
        )

    def _run_builder(self, stats_dir, chrom=None):
        stats_dir.mkdir()
        suffix = f".chr{chrom}" if chrom is not None else ""
        for old_family in ("ancestry", "ancestry_multik"):
            for extension in ("tsv", "parquet"):
                write_text(
                    stats_dir
                    / f"{old_family}.rep_0{suffix}.{extension}",
                    "superseded\n",
                )
        command = [
            sys.executable,
            str(BUILD_SCRIPT),
            "--unrels-path",
            str(self.unrels_path),
            "--source-fam-path",
            str(self.source_fam_path),
            "--admixture-fam-path",
            str(self.inference_fam_path),
            "--supervised-q-path",
            str(self.supervised_q_path),
            "--faststructure-choose-k-path",
            str(self.choose_k_path),
            "--faststructure-prior",
            "simple",
            "--faststructure-seed",
            "22" if chrom is not None else "0",
            "--stats-dir",
            str(stats_dir),
            "--afr-pop",
            "YRI",
            "--eur-pop",
            "CEU",
            "--admixed-pop",
            "ASW",
        ]
        if chrom is not None:
            command.extend(["--chrom", str(chrom)])
        for k in range(2, 6):
            command.extend(
                [
                    "--admixture-q-path",
                    f"{k}={self.admixture_paths[k]}",
                    "--faststructure-q-path",
                    f"{k}={self.faststructure_paths[k]}",
                ]
            )
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env={
                **os.environ,
                "PYTHONPATH": str(self.fake_modules_dir),
            },
        )

    def test_writes_renamed_chromosome_and_genome_contracts(self):
        for chrom in ("22", None):
            with self.subTest(chrom=chrom):
                label = f"chrom_{chrom}" if chrom is not None else "genome"
                stats_dir = self.temp_dir / label

                result = self._run_builder(stats_dir, chrom=chrom)

                self.assertEqual(result.returncode, 0, msg=result.stderr)
                suffix = f".chr{chrom}" if chrom is not None else ""
                expected_columns = {
                    "ancestry_ADMIXTURE_super": SUPER_COLUMNS,
                    "ancestry_ADMIXTURE_multik": MULTIK_COLUMNS,
                    "ancestry_fastStructure_multik": MULTIK_COLUMNS,
                    "fastStructure_chooseK": CHOOSE_K_COLUMNS,
                }
                for family, columns in expected_columns.items():
                    prefix = stats_dir / f"{family}.rep_0{suffix}"
                    rows = read_tsv(Path(f"{prefix}.tsv"))
                    expected = list(columns)
                    if chrom is None:
                        expected.remove("chrom")
                    self.assertEqual(list(rows[0]), expected)
                    self.assertTrue(Path(f"{prefix}.parquet").is_file())

                for old_family in ("ancestry", "ancestry_multik"):
                    for extension in ("tsv", "parquet"):
                        self.assertFalse(
                            (
                                stats_dir
                                / (
                                    f"{old_family}.rep_0{suffix}."
                                    f"{extension}"
                                )
                            ).exists()
                        )

                admixture_rows = read_tsv(
                    stats_dir
                    / f"ancestry_ADMIXTURE_multik.rep_0{suffix}.tsv"
                )
                faststructure_rows = read_tsv(
                    stats_dir
                    / f"ancestry_fastStructure_multik.rep_0{suffix}.tsv"
                )
                expected_order = [
                    (str(k), sample)
                    for k in range(2, 6)
                    for sample in ("sample_3", "sample_1", "sample_2")
                ]
                for rows in (admixture_rows, faststructure_rows):
                    self.assertEqual(
                        [(row["k"], row["vcf_sample_id"]) for row in rows],
                        expected_order,
                    )
                    self.assertEqual(
                        [
                            rows[3][f"component_{index}_q"]
                            for index in range(1, 6)
                        ],
                        ["0.1", "0.3", "0.6", "NA", "NA"],
                    )

                parquet_path = (
                    stats_dir
                    / f"ancestry_fastStructure_multik.rep_0{suffix}.parquet"
                )
                parquet_row = pl.read_parquet(parquet_path).row(3, named=True)
                self.assertIsNone(parquet_row["component_4_q"])
                self.assertIsNone(parquet_row["component_5_q"])
                choose_rows = read_tsv(
                    stats_dir / f"fastStructure_chooseK.rep_0{suffix}.tsv"
                )
                self.assertEqual(
                    choose_rows[0]["max_marginal_likelihood_k"],
                    "3",
                )
                self.assertEqual(choose_rows[0]["model_components_k"], "2")

    def test_registered_empty_schemas_use_new_family_names(self):
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["ancestry_ADMIXTURE_super"],
            SUPER_COLUMNS,
        )
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["ancestry_ADMIXTURE_multik"],
            MULTIK_COLUMNS,
        )
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["ancestry_fastStructure_multik"],
            MULTIK_COLUMNS,
        )
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["fastStructure_chooseK"],
            CHOOSE_K_COLUMNS,
        )
        self.assertNotIn("ancestry", EMPTY_TABLE_COLUMNS)
        self.assertNotIn("ancestry_multik", EMPTY_TABLE_COLUMNS)


class ParseFastStructureChooseKTests(unittest.TestCase):
    def test_parses_valid_report(self):
        row = parse_faststructure_choose_k(
            (
                "Model complexity that maximizes marginal likelihood = 4\n"
                "Model components used to explain structure in data = 3\n"
            ),
            "simple",
            22,
            [2, 3, 4, 5],
            "22",
        )

        self.assertEqual(
            row,
            {
                "rep": 0,
                "chrom": "22",
                "prior": "simple",
                "seed": 22,
                "k_min": 2,
                "k_max": 5,
                "max_marginal_likelihood_k": 4,
                "model_components_k": 3,
            },
        )

    def test_rejects_invalid_reports(self):
        valid_first = (
            "Model complexity that maximizes marginal likelihood = 4\n"
        )
        valid_second = (
            "Model components used to explain structure in data = 3\n"
        )
        cases = {
            "missing": valid_first,
            "duplicate": valid_first + valid_first + valid_second,
            "noninteger": valid_first + valid_second.replace("3", "three"),
            "out_of_range": valid_first.replace("4", "7") + valid_second,
        }
        for label, report in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(ValueError):
                    parse_faststructure_choose_k(
                        report,
                        "simple",
                        0,
                        [2, 3, 4, 5],
                    )


if __name__ == "__main__":
    unittest.main()
