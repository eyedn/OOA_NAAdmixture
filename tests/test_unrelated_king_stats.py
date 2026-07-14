###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_unrelated_king_stats.py
###############################################################################


from pathlib import Path
import csv
import os
import subprocess
import tempfile
import unittest


PROJECT_DIR = Path(__file__).resolve().parents[1]
WRITER_PATH = (
    PROJECT_DIR / "job_scripts" / "python_utils" /
    "write_unrelated_kinship.py"
)
FAKE_PANDAS_SOURCE = '''
from pathlib import Path


class DataFrame:
    def __init__(self, data=None, columns=None):
        data = [] if data is None else data
        if isinstance(data, dict):
            lengths = [
                len(value)
                for value in data.values()
                if isinstance(value, (list, tuple))
            ]
            row_count = max(lengths, default=1)
            self._rows = []
            for row_idx in range(row_count):
                row = {}
                for name, value in data.items():
                    if isinstance(value, (list, tuple)):
                        row[name] = value[row_idx]
                    else:
                        row[name] = value
                self._rows.append(row)
        else:
            self._rows = [dict(row) for row in data]
        self.columns = list(columns or (
            self._rows[0].keys() if self._rows else []
        ))

    def __getitem__(self, name):
        return [row[name] for row in self._rows]

    def to_dict(self, orient="dict"):
        if orient != "records":
            raise ValueError("only records are supported")
        return self._rows

    def to_parquet(self, path, index=False):
        Path(path).write_text(repr(self._rows), encoding="utf-8")


def read_csv(path, sep=None):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0].split()
    rows = []
    for line in lines[1:]:
        values = line.split()
        row = dict(zip(header, values))
        row["KINSHIP"] = float(row["KINSHIP"])
        rows.append(row)
    return DataFrame(rows, columns=header)
'''


class UnrelatedKinshipWriterTests(unittest.TestCase):
    def _write_kin0(self, path, rows):
        path.write_text(
            "#IID1 IID2 KINSHIP\n" +
            "".join(
                f"{sample1} {sample2} {kinship}\n"
                for sample1, sample2, kinship in rows
            ),
            encoding="utf-8",
        )

    def _run_writer(self, king_dir, stats_dir, *extra_args):
        with tempfile.TemporaryDirectory() as module_dir:
            module_path = Path(module_dir)
            (module_path / "pandas.py").write_text(
                FAKE_PANDAS_SOURCE,
                encoding="utf-8",
            )
            (module_path / "allel.py").write_text("", encoding="utf-8")
            (module_path / "tszip.py").write_text("", encoding="utf-8")
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["PYTHONPATH"] = str(module_path)
            return subprocess.run(
                [
                    "python",
                    str(WRITER_PATH),
                    "--rep",
                    "2",
                    "--king-dir",
                    str(king_dir),
                    "--stats-dir",
                    str(stats_dir),
                    "--genetic-map",
                    "Map",
                    "--pops",
                    "AFR",
                    "EUR",
                    *extra_args,
                ],
                cwd=PROJECT_DIR,
                env=env,
                capture_output=True,
                text=True,
            )

    def test_writer_combines_chromosome_populations(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            king_dir = tmp_path / "king"
            stats_dir = tmp_path / "stats"
            king_dir.mkdir()
            self._write_kin0(
                king_dir / "Map_2_chr7_AFR_unrelated.kin0",
                [("AFR_1", "AFR_2", "0.125")],
            )
            self._write_kin0(
                king_dir / "Map_2_chr7_EUR_unrelated.kin0",
                [("EUR_1", "EUR_2", "-0.03125")],
            )

            result = self._run_writer(
                king_dir,
                stats_dir,
                "--chr",
                "7",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            tsv_path = stats_dir / "kinship_unrelated.rep_2.chr7.tsv"
            parquet_path = (
                stats_dir / "kinship_unrelated.rep_2.chr7.parquet"
            )
            with open(tsv_path, encoding="utf-8", newline="") as in_file:
                rows = list(csv.DictReader(in_file, delimiter="\t"))
            self.assertEqual(
                rows,
                [
                    {
                        "rep": "2",
                        "chrom": "7",
                        "pop": "AFR",
                        "sample1": "AFR_1",
                        "sample2": "AFR_2",
                        "kinship": "0.125",
                    },
                    {
                        "rep": "2",
                        "chrom": "7",
                        "pop": "EUR",
                        "sample1": "EUR_1",
                        "sample2": "EUR_2",
                        "kinship": "-0.03125",
                    },
                ],
            )
            self.assertGreater(parquet_path.stat().st_size, 0)

    def test_writer_omits_chromosome_for_genome_output(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            king_dir = tmp_path / "king"
            stats_dir = tmp_path / "stats"
            king_dir.mkdir()
            for pop in ("AFR", "EUR"):
                self._write_kin0(
                    king_dir / f"Map_2_genome_{pop}_unrelated.kin0",
                    [(f"{pop}_1", f"{pop}_2", "0.0625")],
                )

            result = self._run_writer(king_dir, stats_dir)

            self.assertEqual(result.returncode, 0, result.stderr)
            tsv_path = stats_dir / "kinship_unrelated.rep_2.tsv"
            with open(tsv_path, encoding="utf-8", newline="") as in_file:
                reader = csv.DictReader(in_file, delimiter="\t")
                rows = list(reader)
            self.assertNotIn("chrom", reader.fieldnames)
            self.assertEqual([row["pop"] for row in rows], ["AFR", "EUR"])
            self.assertTrue(
                (stats_dir / "kinship_unrelated.rep_2.parquet").is_file()
            )

    def test_writer_reports_exact_missing_kin0_path(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            king_dir = tmp_path / "king"
            stats_dir = tmp_path / "stats"
            king_dir.mkdir()
            missing_path = king_dir / "Map_2_chr7_AFR_unrelated.kin0"

            result = self._run_writer(
                king_dir,
                stats_dir,
                "--chr",
                "7",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(str(missing_path), result.stderr)


class BackfillWorkerTests(unittest.TestCase):
    def _write_fake_plink(self, bin_dir):
        plink_path = bin_dir / "plink2"
        plink_path.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "printf '%s\\n' \"$*\" >> \"$PLINK_LOG\"\n"
            "out_prefix=''\n"
            "while (( $# )); do\n"
            "    if [[ \"$1\" == '--out' ]]; then\n"
            "        out_prefix=\"$2\"\n"
            "        shift 2\n"
            "    else\n"
            "        shift\n"
            "    fi\n"
            "done\n"
            "printf '#IID1 IID2 KINSHIP\\nA_1 A_2 0.1\\n' > "
            "\"${out_prefix}.kin0\"\n"
            ": > \"${out_prefix}.king\"\n"
            ": > \"${out_prefix}.king.id\"\n",
            encoding="utf-8",
        )
        plink_path.chmod(0o755)

    def _run_worker(self, script_name, args, task_id):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            bin_dir = tmp_path / "bin"
            bin_dir.mkdir()
            module_dir = tmp_path / "modules"
            module_dir.mkdir()
            (module_dir / "pandas.py").write_text(
                FAKE_PANDAS_SOURCE,
                encoding="utf-8",
            )
            (module_dir / "allel.py").write_text("", encoding="utf-8")
            (module_dir / "tszip.py").write_text("", encoding="utf-8")
            self._write_fake_plink(bin_dir)
            plink_log = tmp_path / "plink.log"
            bash_env = tmp_path / "bash_env.sh"
            bash_env.write_text(
                "module() { :; }\n"
                "ml() { :; }\n"
                "conda() { :; }\n"
                "source() {\n"
                "    if [[ \"$1\" == /apps/* ]]; then return 0; fi\n"
                "    builtin source \"$@\"\n"
                "}\n",
                encoding="utf-8",
            )
            env = os.environ.copy()
            env.update(
                {
                    "BASH_ENV": str(bash_env),
                    "PATH": f"{bin_dir}:{env['PATH']}",
                    "PLINK_LOG": str(plink_log),
                    "PYTHONPATH": str(module_dir),
                    "PYTHONDONTWRITEBYTECODE": "1",
                    "SLURM_ARRAY_TASK_ID": str(task_id),
                    "SLURM_CPUS_PER_TASK": "3",
                }
            )
            result = subprocess.run(
                ["bash", str(PROJECT_DIR / "job_scripts" / script_name), *args],
                cwd=PROJECT_DIR,
                env=env,
                capture_output=True,
                text=True,
            )
            log_text = plink_log.read_text(encoding="utf-8")
            return result, log_text

    def test_chromosome_worker_reuses_retained_ids_and_original_bed(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            bed_dir = tmp_path / "beds"
            king_dir = tmp_path / "king"
            stats_dir = tmp_path / "stats"
            bed_dir.mkdir()
            king_dir.mkdir()
            for extension in ("bed", "bim", "fam"):
                (bed_dir / f"Map_1_chr7_all.{extension}").write_text("x")
            for pop in ("AFR", "EUR"):
                retained = (
                    king_dir / f"Map_1_chr7_{pop}.king.cutoff.in.id"
                )
                retained.write_text(f"0\t{pop}_1\n", encoding="utf-8")

            result, log_text = self._run_worker(
                "temp_calc_chrom_kinship_unrelated.sh",
                [
                    str(bed_dir),
                    str(king_dir),
                    str(stats_dir),
                    "1",
                    "Map",
                    "--",
                    "7",
                    "--",
                    "AFR",
                    "EUR",
                ],
                1,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("--king-cutoff", log_text)
            self.assertNotIn("--king-table-filter", log_text)
            self.assertEqual(log_text.count("--make-king "), 2)
            self.assertEqual(log_text.count("--make-king-table"), 2)
            bed_argument = f"--bfile {bed_dir}/Map_1_chr7_all"
            self.assertEqual(log_text.count(bed_argument), 2)
            for pop in ("AFR", "EUR"):
                retained = king_dir / (
                    f"Map_1_chr7_{pop}.king.cutoff.in.id"
                )
                self.assertIn(f"--keep {retained}", log_text)
                self.assertIn(
                    f"--out {king_dir}/Map_1_chr7_{pop}_unrelated",
                    log_text,
                )
                for extension in ("kin0", "king", "king.id"):
                    raw_path = king_dir / (
                        f"Map_1_chr7_{pop}_unrelated.{extension}"
                    )
                    self.assertTrue(raw_path.is_file())
            self.assertTrue(
                (stats_dir / "kinship_unrelated.rep_1.chr7.tsv").is_file()
            )

    def test_genome_worker_reuses_retained_ids_and_merged_bed(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            bed_dir = tmp_path / "beds"
            king_dir = tmp_path / "king"
            stats_dir = tmp_path / "stats"
            bed_dir.mkdir()
            king_dir.mkdir()
            for extension in ("bed", "bim", "fam"):
                (bed_dir / f"Map_1_genome_all.{extension}").write_text("x")
            for pop in ("AFR", "EUR"):
                retained = (
                    king_dir / f"Map_1_genome_{pop}.king.cutoff.in.id"
                )
                retained.write_text(f"0\t{pop}_1\n", encoding="utf-8")

            result, log_text = self._run_worker(
                "temp_calc_genome_kinship_unrelated.sh",
                [
                    str(bed_dir),
                    str(king_dir),
                    str(stats_dir),
                    "1",
                    "Map",
                    "--",
                    "AFR",
                    "EUR",
                ],
                1,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("--king-cutoff", log_text)
            self.assertNotIn("--king-table-filter", log_text)
            self.assertEqual(log_text.count("--make-king "), 2)
            self.assertEqual(log_text.count("--make-king-table"), 2)
            self.assertEqual(
                log_text.count(f"--bfile {bed_dir}/Map_1_genome_all"),
                2,
            )
            for pop in ("AFR", "EUR"):
                for extension in ("kin0", "king", "king.id"):
                    raw_path = king_dir / (
                        f"Map_1_genome_{pop}_unrelated.{extension}"
                    )
                    self.assertTrue(raw_path.is_file())
            self.assertTrue(
                (stats_dir / "kinship_unrelated.rep_1.tsv").is_file()
            )


class BackfillSubmissionTests(unittest.TestCase):
    def _assert_common_slurm_flags(self, call):
        for flag in (
            "--parsable",
            "--chdir=",
            "--job-name=",
            "--cpus-per-task=",
            "--mem=",
            "--time=1-00:00:00",
            "--partition=qcb",
            "--account=jazlynmo_738",
            "--nodes=1",
            "--output=",
            "--error=",
            "--mail-type=ALL",
            "--mail-user=karatas@usc.edu",
        ):
            self.assertIn(flag, call)

    def _run_submitter(self, script_name):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            bin_dir = tmp_path / "bin"
            bin_dir.mkdir()
            sbatch_log = tmp_path / "sbatch.log"
            sbatch_state = tmp_path / "sbatch.state"
            mkdir_path = bin_dir / "mkdir"
            mkdir_path.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$*\" == *'/home1/karatas/logs/'* ]]; then\n"
                "    exit 0\n"
                "fi\n"
                "exec /bin/mkdir \"$@\"\n",
                encoding="utf-8",
            )
            mkdir_path.chmod(0o755)
            sbatch_path = bin_dir / "sbatch"
            sbatch_path.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "printf '%s\\n' \"$*\" >> \"$SBATCH_LOG\"\n"
                "jid=100\n"
                "if [[ -f \"$SBATCH_STATE\" ]]; then\n"
                "    jid=$(( $(cat \"$SBATCH_STATE\") + 1 ))\n"
                "fi\n"
                "printf '%s\\n' \"$jid\" > \"$SBATCH_STATE\"\n"
                "printf '%s\\n' \"$jid\"\n",
                encoding="utf-8",
            )
            sbatch_path.chmod(0o755)
            env = os.environ.copy()
            env.update(
                {
                    "HOME": str(tmp_path / "home"),
                    "PATH": f"{bin_dir}:{env['PATH']}",
                    "SBATCH_LOG": str(sbatch_log),
                    "SBATCH_STATE": str(sbatch_state),
                    "SCRATCH": str(tmp_path / "scratch"),
                }
            )
            result = subprocess.run(
                ["bash", str(PROJECT_DIR / script_name)],
                cwd=PROJECT_DIR,
                env=env,
                capture_output=True,
                text=True,
            )
            return result, sbatch_log.read_text(encoding="utf-8").splitlines()

    def test_chromosome_submitter_gates_filtered_combine_array(self):
        result, calls = self._run_submitter(
            "temp_submit_calc_chrom_kinship_unrelated.sh"
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(calls), 2)
        self._assert_common_slurm_flags(calls[0])
        self._assert_common_slurm_flags(calls[1])
        self.assertIn("--array=1-1100%100", calls[0])
        self.assertIn("--dependency=afterok:100", calls[1])
        self.assertIn("--array=1-22%100", calls[1])
        self.assertIn("/%A.%x.out", calls[1])
        self.assertIn("/%A.%x.err", calls[1])
        self.assertIn("--tables kinship_unrelated", calls[1])

    def test_genome_submitter_gates_filtered_combine_job(self):
        result, calls = self._run_submitter(
            "temp_submit_calc_genome_kinship_unrelated.sh"
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(calls), 2)
        self._assert_common_slurm_flags(calls[0])
        self._assert_common_slurm_flags(calls[1])
        self.assertIn("--array=1-50%100", calls[0])
        self.assertIn("--dependency=afterok:100", calls[1])
        self.assertNotIn("--array", calls[1])
        self.assertIn("--tables kinship_unrelated", calls[1])


class PermanentWorkflowTests(unittest.TestCase):
    def test_filtered_combine_does_not_rewrite_existing_tables(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            stats_dir = Path(tmp_dir)
            header = "rep\tpop\tsample1\tsample2\tkinship\n"
            for rep in (1, 2):
                path = stats_dir / f"kinship_unrelated.rep_{rep}.tsv"
                path.write_text(
                    header + f"{rep}\tAFR\tA_1\tA_2\t0.1\n",
                    encoding="utf-8",
                )
            ancestry_path = stats_dir / "ancestry.tsv"
            ancestry_path.write_text("do not rewrite\n", encoding="utf-8")
            module_dir = stats_dir / "modules"
            module_dir.mkdir()
            (module_dir / "pandas.py").write_text(
                FAKE_PANDAS_SOURCE,
                encoding="utf-8",
            )
            (module_dir / "allel.py").write_text("", encoding="utf-8")
            (module_dir / "tszip.py").write_text("", encoding="utf-8")
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["PYTHONPATH"] = str(module_dir)

            result = subprocess.run(
                [
                    "python",
                    str(
                        PROJECT_DIR / "job_scripts" / "python_utils" /
                        "combine_sim_stats.py"
                    ),
                    "--stats-dir",
                    str(stats_dir),
                    "--num-reps",
                    "2",
                    "--tables",
                    "kinship_unrelated",
                ],
                cwd=PROJECT_DIR,
                env=env,
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            combined_path = stats_dir / "kinship_unrelated.tsv"
            with open(combined_path, encoding="utf-8", newline="") as in_file:
                rows = list(csv.DictReader(in_file, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertTrue(
                (stats_dir / "kinship_unrelated.parquet").is_file()
            )
            self.assertEqual(
                ancestry_path.read_text(encoding="utf-8"),
                "do not rewrite\n",
            )

    def test_permanent_workers_run_second_pass_after_retained_validation(self):
        for relative_path in (
            "job_scripts/calc_chrom_stats.sh",
            "job_scripts/calc_genome_stats.sh",
        ):
            text = (PROJECT_DIR / relative_path).read_text(encoding="utf-8")
            validation_pos = text.index("if [[ ! -s \"${retained_path}\" ]]")
            unrelated_pos = text.index("unrelated_out_prefix=")
            writer_pos = text.index("write_unrelated_kinship.py")
            self.assertLess(validation_pos, unrelated_pos)
            self.assertLess(unrelated_pos, writer_pos)
            unrelated_block = text[unrelated_pos:writer_pos]
            self.assertIn('--keep "${retained_path}"', unrelated_block)
            self.assertIn("--make-king-table", unrelated_block)
            self.assertNotIn("--king-cutoff", unrelated_block)
            self.assertNotIn("--king-table-filter", unrelated_block)

    def test_combiner_defaults_include_unrelated_and_accepts_tables(self):
        text = (
            PROJECT_DIR / "job_scripts" / "python_utils" /
            "combine_sim_stats.py"
        ).read_text(encoding="utf-8")
        self.assertIn('"kinship_unrelated"', text)
        self.assertIn('"--tables"', text)
        self.assertIn("args.tables", text)


if __name__ == "__main__":
    unittest.main()
