#!/usr/bin/env python3
"""Orchestrate standalone pytest suites and Snakemake workflow regression tests."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = PACKAGE_ROOT / "workflow" / "Snakefile"
CONFIGFILE = PACKAGE_ROOT / "tests" / "test-config.yaml"
DATABASES_FULL = PACKAGE_ROOT / "tests" / "test-databases.yaml"
DATABASES_PYHMMER = PACKAGE_ROOT / "tests" / "test-pyhmmer.yaml"
OUTPUT_DIR = PACKAGE_ROOT / "tests" / "output"
PUBLISHED_TRUTHS = PACKAGE_ROOT / "tests" / "test_published_biological_truths.py"
PROKANOTA_TRUTHS = PACKAGE_ROOT / "tests" / "test_prokanota_truths.py"
WORKFLOW_OUTPUT_TESTS = PACKAGE_ROOT / "tests" / "test_workflow_outputs.py"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Prokanota tests")
    parser.add_argument(
        "--mode",
        choices=("minimal", "pyhmmer", "full"),
        default="pyhmmer",
        help="Select which test mode to run",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show full tool output in the console",
    )
    parser.add_argument(
        "--very-verbose",
        action="store_true",
        help="Show tool output and verbose Snakemake execution details",
    )
    return parser.parse_args(argv)


def build_snakemake_config(
    config_path: Path, databases_path: Path, mode_label: str
) -> Path:
    """
    Create the effective test configuration for the selected database mode.

    The committed base configuration is left untouched. The generated copy is
    written under the test logs directory so the exact workflow configuration
    remains available alongside the run output.
    """
    config_data = yaml.safe_load(config_path.read_text())
    global_config = config_data.get("global")
    if not isinstance(global_config, dict):
        raise ValueError(f"Invalid config file {config_path}: missing global section")

    global_config["databases"] = str(databases_path)
    logs_dir = OUTPUT_DIR / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    target_path = logs_dir / f"test-config-{mode_label}.yaml"
    target_path.write_text(yaml.safe_dump(config_data, sort_keys=False))
    return target_path


def run_pytest(test_files: list[Path], env: dict[str, str] | None = None) -> int:
    """Run one or more pytest suites with visible, biologically named tests."""
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-v",
        "--no-header",
        *(str(path) for path in test_files),
    ]
    print("Running:", " ".join(cmd))
    return subprocess.run(cmd, cwd=PACKAGE_ROOT, env=env).returncode


def run_snakemake(configfile: Path, verbose: bool, very_verbose: bool) -> int:
    """
    Generate fresh workflow outputs in isolated rule-specific Conda environments.

    ``--forceall`` is intentional for golden regression testing: every rule is
    rerun so the output tests cannot pass by reusing files from an older run.
    """
    if very_verbose:
        verbose = True

    cmd = [
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(configfile),
        "--notemp",
        "--cores",
        "1",
        "--forceall",
        "--sdm",
        "conda",
    ]
    if verbose:
        cmd.extend(["--config", "prokanota_verbose=true"])
    if very_verbose:
        cmd.extend(["--printshellcmds", "--show-failed-logs"])

    print("Running:", " ".join(cmd))
    return subprocess.run(cmd, cwd=PACKAGE_ROOT).returncode


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    # Phase 1: published biological truths and stable Prokanota contracts require
    # only pytest plus the dependencies available in the calling environment.
    standalone_result = run_pytest([PUBLISHED_TRUTHS, PROKANOTA_TRUTHS])
    if standalone_result != 0 or args.mode == "minimal":
        return standalone_result

    # Phase 2: pyhmmer and full modes produce a completely fresh Snakemake run.
    databases_path = DATABASES_FULL if args.mode == "full" else DATABASES_PYHMMER
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)

    config_path = build_snakemake_config(CONFIGFILE, databases_path, args.mode)
    snakemake_result = run_snakemake(config_path, args.verbose, args.very_verbose)
    if snakemake_result != 0:
        return snakemake_result

    # Phase 3: validate the newly generated files. The environment variable both
    # selects the expected database columns and proves that outputs were prepared
    # by this runner rather than discovered accidentally by a plain pytest call.
    workflow_env = os.environ.copy()
    workflow_env["PROKANOTA_TEST_MODE"] = args.mode
    return run_pytest([WORKFLOW_OUTPUT_TESTS], env=workflow_env)


if __name__ == "__main__":
    raise SystemExit(main())
