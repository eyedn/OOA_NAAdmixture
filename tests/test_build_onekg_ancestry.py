###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_build_onekg_ancestry.py
###############################################################################

# overview: test empirical supervised and neutral multi-K ancestry outputs.


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

from onekg_utils.write_stats_table import EMPTY_TABLE_COLUMNS


BUILD_SCRIPT = (
    PYTHON_UTILS_DIR / "build_onekg_ancestry.py"
)
ANCESTRY_COLUMNS = [
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


##### helpers #################################################################
def write_text(path, text):
    path.write_text(text, encoding="utf-8")


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))


##### tests ###################################################################
class BuildOnekgAncestryTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.input_dir = self.temp_dir / "input"
        self.input_dir.mkdir()
        self.fake_modules_dir = self.temp_dir / "fake_modules"
        self.fake_modules_dir.mkdir()
        self._write_inputs()
        self._write_fake_pandas()

    def tearDown(self):
        self.temp_dir_context.cleanup()

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
        self.admixture_fam_path = self.input_dir / "admixture.fam"
        fam_text = (
            "ASW sample_3 0 0 0 -9\n"
            "YRI sample_1 0 0 0 -9\n"
            "CEU sample_2 0 0 0 -9\n"
        )
        write_text(self.source_fam_path, fam_text)
        write_text(self.admixture_fam_path, fam_text)
        self.supervised_q_path = self.input_dir / "supervised.2.Q"
        write_text(
            self.supervised_q_path,
            "0.40 0.60\n"
            "0.90 0.10\n"
            "0.10 0.90\n",
        )
        q_values_by_k = {
            2: (
                "0.35 0.65\n"
                "0.80 0.20\n"
                "0.15 0.85\n"
            ),
            3: (
                "0.10 0.30 0.60\n"
                "0.70 0.20 0.10\n"
                "0.20 0.70 0.10\n"
            ),
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
        self.q_paths = {}
        for k, text in q_values_by_k.items():
            q_path = self.input_dir / f"unsupervised.{k}.Q"
            write_text(q_path, text)
            self.q_paths[k] = q_path

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

    def _base_command(self, stats_dir):
        return [
            sys.executable,
            str(BUILD_SCRIPT),
            "--unrels-path",
            str(self.unrels_path),
            "--source-fam-path",
            str(self.source_fam_path),
            "--admixture-fam-path",
            str(self.admixture_fam_path),
            "--supervised-q-path",
            str(self.supervised_q_path),
            "--stats-dir",
            str(stats_dir),
            "--afr-pop",
            "YRI",
            "--eur-pop",
            "CEU",
            "--admixed-pop",
            "ASW",
        ]

    def _run_parser(self, stats_dir, chrom=None, q_specs=None):
        stats_dir.mkdir()
        command = self._base_command(stats_dir)
        if chrom is not None:
            command.extend(["--chrom", str(chrom)])
        if q_specs is None:
            q_specs = [
                f"{k}={self.q_paths[k]}" for k in sorted(self.q_paths)
            ]
        for q_spec in q_specs:
            command.extend(["--unsupervised-q-path", q_spec])
        environment = {
            **os.environ,
            "PYTHONPATH": str(self.fake_modules_dir),
        }
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )

    def test_writes_chromosome_and_genome_contracts(self):
        for chrom in ("22", None):
            with self.subTest(chrom=chrom):
                label = f"chrom_{chrom}" if chrom is not None else "genome"
                stats_dir = self.temp_dir / label

                result = self._run_parser(stats_dir, chrom=chrom)

                self.assertEqual(result.returncode, 0, msg=result.stderr)
                suffix = f".chr{chrom}" if chrom is not None else ""
                ancestry_prefix = stats_dir / f"ancestry.rep_0{suffix}"
                multik_prefix = (
                    stats_dir / f"ancestry_multik.rep_0{suffix}"
                )
                ancestry_rows = read_tsv(
                    Path(f"{ancestry_prefix}.tsv")
                )
                multik_rows = read_tsv(
                    Path(f"{multik_prefix}.tsv")
                )
                expected_ancestry_columns = list(ANCESTRY_COLUMNS)
                expected_multik_columns = list(MULTIK_COLUMNS)
                if chrom is None:
                    expected_ancestry_columns.remove("chrom")
                    expected_multik_columns.remove("chrom")
                self.assertEqual(
                    list(ancestry_rows[0]),
                    expected_ancestry_columns,
                )
                self.assertEqual(
                    list(multik_rows[0]),
                    expected_multik_columns,
                )
                self.assertNotIn(
                    "afr_unsupervised_q",
                    ancestry_rows[0],
                )
                self.assertNotIn(
                    "eur_unsupervised_q",
                    ancestry_rows[0],
                )
                self.assertEqual(
                    [row["vcf_sample_id"] for row in ancestry_rows],
                    ["sample_3", "sample_1", "sample_2"],
                )
                self.assertEqual(
                    [row["afr_q"] for row in ancestry_rows],
                    ["0.4", "0.9", "0.1"],
                )
                self.assertEqual(
                    [
                        (row["k"], row["vcf_sample_id"])
                        for row in multik_rows
                    ],
                    [
                        (str(k), sample)
                        for k in range(2, 6)
                        for sample in (
                            "sample_3",
                            "sample_1",
                            "sample_2",
                        )
                    ],
                )
                first_k3_row = multik_rows[3]
                self.assertEqual(
                    [
                        first_k3_row[f"component_{index}_q"]
                        for index in range(1, 6)
                    ],
                    ["0.1", "0.3", "0.6", "NA", "NA"],
                )
                parquet_table = pl.read_parquet(
                    Path(f"{multik_prefix}.parquet")
                )
                self.assertEqual(
                    parquet_table.columns,
                    expected_multik_columns,
                )
                first_k3_parquet_row = parquet_table.row(
                    3,
                    named=True,
                )
                self.assertEqual(
                    [
                        first_k3_parquet_row[
                            f"component_{index}_q"
                        ]
                        for index in range(1, 6)
                    ],
                    [0.1, 0.3, 0.6, None, None],
                )

    def test_registered_empty_table_schemas_match_output_contracts(self):
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["ancestry"],
            ANCESTRY_COLUMNS,
        )
        self.assertEqual(
            EMPTY_TABLE_COLUMNS["ancestry_multik"],
            MULTIK_COLUMNS,
        )

    def test_rejects_duplicate_and_invalid_k_specifications(self):
        invalid_cases = {
            "duplicate": [
                f"2={self.q_paths[2]}",
                f"2={self.q_paths[2]}",
            ],
            "malformed": [str(self.q_paths[2])],
            "below_limit": [
                f"1={self.q_paths[2]}",
                f"2={self.q_paths[2]}",
            ],
            "above_limit": [
                f"2={self.q_paths[2]}",
                f"6={self.q_paths[2]}",
            ],
        }

        for label, q_specs in invalid_cases.items():
            with self.subTest(label=label):
                stats_dir = self.temp_dir / f"invalid_{label}"

                result = self._run_parser(stats_dir, q_specs=q_specs)

                self.assertNotEqual(result.returncode, 0)
                if label == "duplicate":
                    self.assertIn(
                        "Duplicate unsupervised ADMIXTURE K=2",
                        result.stderr,
                    )
                elif label == "malformed":
                    self.assertIn(
                        "Unsupervised Q paths must use the format K=PATH",
                        result.stderr,
                    )
                else:
                    self.assertIn(
                        "ADMIXTURE K must be between 2 and 5",
                        result.stderr,
                    )


if __name__ == "__main__":
    unittest.main()
