###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           test_python_import_contracts.py
###############################################################################

import ast
import subprocess
import sys
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
ENTRYPOINTS = [
    PROJECT_DIR / "job_scripts/python_utils/sim_model.py",
    PROJECT_DIR / "job_scripts/python_utils/calc_sim_stats.py",
    PROJECT_DIR / "job_scripts/python_utils/aggregate_genome_stats.py",
    PROJECT_DIR / "job_scripts/python_utils/combine_sim_stats.py",
    PROJECT_DIR / "job_scripts/python_utils/write_admixture_pop.py",
    PROJECT_DIR / "job_scripts/python_utils/write_pop_subset.py",
]


# Find Python files invoked by shell scripts and require them to exist.
class PythonImportContractsTest(unittest.TestCase):
    # Verify shell-invoked Python paths resolve to files in this checkout.
    def test_shell_python_entrypoints_exist(self):
        shell_paths = sorted(PROJECT_DIR.glob("**/*.sh"))
        missing_paths = []
        for shell_path in shell_paths:
            text = shell_path.read_text(encoding="utf-8")
            for line in text.splitlines():
                if "job_scripts/" not in line or ".py" not in line:
                    continue
                fragment = line.split("job_scripts/", 1)[1].split(".py", 1)[0]
                python_path = PROJECT_DIR / f"job_scripts/{fragment}.py"
                if not python_path.exists():
                    missing_paths.append(
                        f"{shell_path.relative_to(PROJECT_DIR)} -> "
                        f"{python_path.relative_to(PROJECT_DIR)}"
                    )
        self.assertEqual([], missing_paths)

    # Verify entrypoint parsers can load from the project root.
    def test_entrypoint_help_uses_project_root_imports(self):
        for entrypoint in ENTRYPOINTS:
            with self.subTest(entrypoint=entrypoint.name):
                result = subprocess.run(
                    [sys.executable, str(entrypoint), "--help"],
                    cwd=PROJECT_DIR,
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(
                    0,
                    result.returncode,
                    msg=result.stderr,
                )

    # Verify imports use the package root expected under Slurm.
    def test_no_stale_python_utils_absolute_imports(self):
        stale_imports = []
        for python_path in (PROJECT_DIR / "job_scripts").rglob("*.py"):
            tree = ast.parse(python_path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                module = None
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name == "python_utils":
                            module = alias.name
                elif isinstance(node, ast.ImportFrom):
                    module = node.module
                if module and module.startswith("python_utils"):
                    stale_imports.append(
                        f"{python_path.relative_to(PROJECT_DIR)}:{node.lineno}"
                    )
        self.assertEqual([], stale_imports)

    # Verify imports are not hidden inside functions or main guards.
    def test_imports_stay_at_module_top_level(self):
        nested_imports = []
        for python_path in (PROJECT_DIR / "job_scripts").rglob("*.py"):
            tree = ast.parse(python_path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    for child in ast.walk(node):
                        if isinstance(child, (ast.Import, ast.ImportFrom)):
                            nested_imports.append(
                                f"{python_path.relative_to(PROJECT_DIR)}:"
                                f"{child.lineno}"
                            )
                if isinstance(node, ast.If):
                    is_main_guard = (
                        isinstance(node.test, ast.Compare)
                        and isinstance(node.test.left, ast.Name)
                        and node.test.left.id == "__name__"
                    )
                    if not is_main_guard:
                        continue
                    for child in ast.walk(node):
                        if isinstance(child, (ast.Import, ast.ImportFrom)):
                            nested_imports.append(
                                f"{python_path.relative_to(PROJECT_DIR)}:"
                                f"{child.lineno}"
                            )
        self.assertEqual([], nested_imports)


if __name__ == "__main__":
    unittest.main()
