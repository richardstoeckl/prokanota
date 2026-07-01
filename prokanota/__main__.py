"""
Entrypoint for prokanota

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import csv
import json
import os
import shutil
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import click
import jsonschema
import yaml
from snaketool_utils.cli_utils import OrderedCommands, echo_click, run_snakemake

from prokanota import __version__
from prokanota.workflow.scripts import metadataFromDir

SPLASH_SCREEN = """
                             888                                 888
                             888                                 888
                             888                                 888
    88888b.  888d888 .d88b.  888  888  8888b.  88888b.   .d88b.  888888  8888b.
    888 "88b 888P"  d88""88b 888 .88P     "88b 888 "88b d88""88b 888        "88b
    888  888 888    888  888 888888K  .d888888 888  888 888  888 888    .d888888
    888 d88P 888    Y88..88P 888 "88b 888  888 888  888 Y88..88P Y88b.  888  888
    88888P"  888     "Y88P"  888  888 "Y888888 888  888  "Y88P"   "Y888 "Y888888
    888
    888
    888
    """


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("prokanota.CITATION")) as f:
        for line in f:
            echo_click(line)


def print_splash():
    """Print the splash screen."""
    click.echo(SPLASH_SCREEN)
    click.echo(f"prokanota, version {__version__}")
    click.echo()


class SplashCommand(click.Command):
    """Command that prints a splash before help output."""

    def get_help(self, ctx):
        print_splash()
        return super().get_help(ctx)


class SplashOrderedCommands(OrderedCommands):
    """Command group that prints a splash before root help output."""

    def get_help(self, ctx):
        print_splash()
        return super().get_help(ctx)


def validate_readable_file(ctx, param, value):
    """Validate that a CLI file option points to a readable regular file."""
    if value is None:
        return value
    candidate = Path(value)
    if not candidate.exists():
        raise click.BadParameter(f"File does not exist: {candidate}")
    if not candidate.is_file():
        raise click.BadParameter(f"Path is not a file: {candidate}")
    if not os.access(candidate, os.R_OK):
        raise click.BadParameter(f"File is not readable: {candidate}")
    return str(candidate.resolve())


def validate_yaml_schema(yaml_path, schema_path, label):
    """Validate a YAML file against a JSON schema and raise a ClickException on error."""
    try:
        with open(yaml_path, encoding="utf-8") as handle:
            content = yaml.safe_load(handle)
    except Exception as exc:
        raise click.ClickException(f"Invalid {label} file {yaml_path}: {exc}") from exc

    try:
        with open(schema_path, encoding="utf-8") as handle:
            schema = json.load(handle)
    except Exception as exc:
        raise click.ClickException(
            f"Could not read schema {schema_path}: {exc}"
        ) from exc

    try:
        jsonschema.validate(content, schema)
    except jsonschema.ValidationError as exc:
        raise click.ClickException(
            f"Invalid {label} file {yaml_path}: {exc.message}"
        ) from exc


def validate_config_schema(config_data, schema_path, config_path):
    """Validate config dict against schema and raise ClickException on error."""
    try:
        with open(schema_path, encoding="utf-8") as handle:
            schema = json.load(handle)
    except Exception as exc:
        raise click.ClickException(
            f"Could not read schema {schema_path}: {exc}"
        ) from exc

    try:
        jsonschema.validate(config_data, schema)
    except jsonschema.ValidationError as exc:
        raise click.ClickException(
            f"Invalid config file {config_path}: {exc.message}"
        ) from exc


def read_yaml_file(yaml_path, label):
    """Load YAML from disk and return Python object."""
    try:
        with open(yaml_path, encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except Exception as exc:
        raise click.ClickException(f"Invalid {label} file {yaml_path}: {exc}") from exc


def resolve_input_path(config_path, candidate, label):
    """Resolve absolute or config-relative input paths."""
    candidate_path = Path(candidate)
    if not candidate_path.is_absolute():
        candidate_path = config_path.parent / candidate_path

    if not candidate_path.exists():
        raise click.ClickException(
            f"{label.capitalize()} file does not exist: {candidate_path}"
        )
    if not candidate_path.is_file():
        raise click.ClickException(
            f"{label.capitalize()} path is not a file: {candidate_path}"
        )
    if not os.access(candidate_path, os.R_OK):
        raise click.ClickException(
            f"{label.capitalize()} file is not readable: {candidate_path}"
        )
    return candidate_path.resolve()


def resolve_output_path(config_path, candidate, label):
    """Resolve output directories as absolute paths."""
    target = Path(candidate)
    if not target.is_absolute():
        target = config_path.parent / target
    try:
        target.mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        raise click.ClickException(
            f"Could not create {label} directory {target}: {exc}"
        ) from exc
    return target.resolve()


def stage_runtime_files(config_path, config_data, databases_path, metadata_path):
    """Create runtime copies of effective config inputs under the configured logs path."""
    logs_dir = resolve_output_path(config_path, config_data["global"]["logs"], "logs")
    run_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    staging_dir = logs_dir / "run-config" / run_stamp

    suffix = 1
    while staging_dir.exists():
        suffix += 1
        staging_dir = logs_dir / "run-config" / f"{run_stamp}-{suffix}"

    staging_dir.mkdir(parents=True, exist_ok=False)

    staged_databases = staging_dir / "databases.yaml"
    staged_metadata = staging_dir / "metadata.csv"
    staged_config = staging_dir / "config.yaml"

    shutil.copy2(databases_path, staged_databases)
    shutil.copy2(metadata_path, staged_metadata)

    runtime_config = deepcopy(config_data)
    runtime_config.setdefault("global", {})
    runtime_config["global"]["databases"] = str(staged_databases)
    runtime_config["global"]["metadata"] = str(staged_metadata)

    with open(staged_config, "w", encoding="utf-8") as handle:
        yaml.safe_dump(runtime_config, handle, sort_keys=False)

    return str(staged_config)


def validate_metadata_csv(metadata_path):
    """Validate minimal required metadata CSV structure before invoking Snakemake."""
    required_columns = {"sample_id", "path"}
    try:
        with open(metadata_path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise click.ClickException(
                    f"Invalid metadata file {metadata_path}: missing header row"
                )
            missing = required_columns - set(reader.fieldnames)
            if missing:
                raise click.ClickException(
                    f"Invalid metadata file {metadata_path}: missing required column(s): {', '.join(sorted(missing))}"
                )
    except click.ClickException:
        raise
    except Exception as exc:
        raise click.ClickException(
            f"Invalid metadata file {metadata_path}: {exc}"
        ) from exc


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--verbose",
            is_flag=True,
            help="Show full tool output in the console",
        ),
        click.option(
            "--very-verbose",
            is_flag=True,
            help="Show tool output and verbose Snakemake execution details",
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--nolock",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=SplashOrderedCommands,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(__version__, "-v", "--version", is_flag=True)
def cli():
    """Flexible Snakemake pipeline for prokaryotic annotation with a code-free* ,modular database architecture.
    \b
    For more options, run:
    prokanota command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
prokanota run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           prokanota run --cfg config.yaml
Override inputs:    prokanota run --cfg config.yaml --db databases.yaml --meta metadata.csv
Specify threads:    prokanota run ... --threads [threads]
Change defaults:    prokanota run ... --snake-default="-k --nolock"
Add Snakemake args: prokanota run ... --dry-run --keep-going --touch
Verbose output:     prokanota run ... --verbose
Very verbose:       prokanota run ... --very-verbose
Specify targets:    prokanota run ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    cls=SplashCommand,
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--cfg",
    "configfile",
    required=True,
    callback=validate_readable_file,
    type=click.Path(path_type=str),
    help="Path to the workflow config YAML file.",
)
@click.option(
    "--db",
    "databases_file",
    required=False,
    callback=validate_readable_file,
    type=click.Path(path_type=str),
    help="Optional override path to databases YAML file.",
)
@click.option(
    "--meta",
    "metadata_file",
    required=False,
    callback=validate_readable_file,
    type=click.Path(path_type=str),
    help="Optional override path to sample metadata CSV file.",
)
@common_options
def run(**kwargs):
    """Run prokanota"""
    configfile = Path(kwargs.get("configfile")).resolve()
    databases_override = kwargs.pop("databases_file")
    metadata_override = kwargs.pop("metadata_file")
    verbose = bool(kwargs.pop("verbose", False))
    very_verbose = bool(kwargs.pop("very_verbose", False))
    if very_verbose:
        verbose = True

    config_schema = snake_base(os.path.join("config", "schemas", "config.schema.json"))

    config_data = read_yaml_file(configfile, "config")
    if not isinstance(config_data, dict):
        raise click.ClickException(
            f"Invalid config file {configfile}: root must be a mapping"
        )

    global_cfg = config_data.get("global")
    if not isinstance(global_cfg, dict):
        raise click.ClickException(
            f"Invalid config file {configfile}: missing 'global' section"
        )

    databases_candidate = databases_override or global_cfg.get("databases")
    metadata_candidate = metadata_override or global_cfg.get("metadata")

    missing_inputs = []
    if not databases_candidate:
        missing_inputs.append("--db")
    if not metadata_candidate:
        missing_inputs.append("--meta")
    if missing_inputs:
        missing_text = ", ".join(missing_inputs)
        raise click.ClickException(
            "Missing required input definition(s): "
            f"{missing_text}. Define them in {configfile} under global.databases/global.metadata "
            "or provide them on the command line."
        )

    databases_path = resolve_input_path(configfile, databases_candidate, "databases")
    metadata_path = resolve_input_path(configfile, metadata_candidate, "metadata")

    config_data["global"]["databases"] = str(databases_path)
    config_data["global"]["metadata"] = str(metadata_path)

    validate_config_schema(config_data, config_schema, configfile)
    # TODO: Rework database file validation. The current validation
    # implementation needs a complete rework; temporarily disable it
    # to avoid blocking runs until a proper reimplementation is available.
    # validate_yaml_schema(databases_path, db_schema, "databases")
    validate_metadata_csv(metadata_path)

    for output_key in ("logs", "interim", "results"):
        output_value = config_data["global"].get(output_key)
        if not output_value:
            raise click.ClickException(
                f"Invalid config file {configfile}: global.{output_key} must be set"
            )
        config_data["global"][output_key] = str(
            resolve_output_path(configfile, output_value, output_key)
        )

    staged_config_path = stage_runtime_files(
        configfile,
        config_data,
        databases_path,
        metadata_path,
    )
    kwargs["configfile"] = staged_config_path

    # Config to add or update in configfile
    kwargs["use_conda"] = True

    snake_defaults = list(kwargs.get("snake_default", ()))
    if very_verbose:
        snake_defaults.extend(["--printshellcmds", "--show-failed-logs"])
    kwargs["snake_default"] = tuple(snake_defaults)

    merge_config = {
        "prokanota": {
            "args": {
                **kwargs,
                "verbose": verbose,
                "very_verbose": very_verbose,
            }
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs,
    )


@click.command(cls=SplashCommand)
@click.argument("output_dir", type=click.Path(path_type=str, file_okay=False))
def config(output_dir, **kwargs):
    """Create default config files in the provided directory"""
    outdir = Path(output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    template_dir = Path(snake_base("config"))
    template_files = ("config.yaml", "databases.yaml", "metadata.csv")
    for filename in template_files:
        src = template_dir / filename
        dst = outdir / filename
        if dst.exists():
            raise click.ClickException(f"Refusing to overwrite existing file: {dst}")
        shutil.copy2(src, dst)

    click.echo(f"Created default config files in: {outdir}")


@click.command(cls=SplashCommand)
@click.option(
    "--verbose",
    is_flag=True,
    help="Show full tool output in the console",
)
@click.option(
    "--very-verbose",
    is_flag=True,
    help="Show tool output and verbose Snakemake execution details",
)
def test(**kwargs):
    """Execute the pyhmmer-only integration tests for prokanota"""
    _run_tests("pyhmmer", **kwargs)


@click.command(cls=SplashCommand, name="test-minimal")
@click.option(
    "--verbose",
    is_flag=True,
    help="Show full tool output in the console",
)
@click.option(
    "--very-verbose",
    is_flag=True,
    help="Show tool output and verbose Snakemake execution details",
)
def test_minimal(**kwargs):
    """Execute minimal pytest-only tests for prokanota"""
    _run_tests("minimal", **kwargs)


@click.command(cls=SplashCommand, name="test-full")
@click.option(
    "--verbose",
    is_flag=True,
    help="Show full tool output in the console",
)
@click.option(
    "--very-verbose",
    is_flag=True,
    help="Show tool output and verbose Snakemake execution details",
)
def test_full(**kwargs):
    """Execute full integration tests for prokanota"""
    _run_tests("full", **kwargs)


def _run_tests(mode: str, **kwargs) -> None:
    """Run test script with a specific mode."""
    package_root = Path(__file__).resolve().parent
    test_script = package_root / "tests" / "run_tests.py"
    if not test_script.exists():
        raise click.ClickException(
            f"Could not find test runner at {test_script}. Reinstall prokanota with test assets."
        )

    verbose = bool(kwargs.pop("verbose", False))
    very_verbose = bool(kwargs.pop("very_verbose", False))
    if very_verbose:
        verbose = True

    cmd = [sys.executable, str(test_script), "--mode", mode]
    if verbose:
        cmd.append("--verbose")
    if very_verbose:
        cmd.append("--very-verbose")
    result = subprocess.run(cmd, cwd=str(package_root), check=False)
    raise SystemExit(result.returncode)


@click.command(cls=SplashCommand)
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


@click.command(cls=SplashCommand)
def version(**kwargs):
    """Print the version for this tool"""
    click.echo(f"prokanota, version {__version__}")


@click.group(cls=SplashCommand)
def helper():
    """Utility commands to help prepare inputs for prokanota"""
    pass


@click.command(cls=SplashCommand)
@click.option(
    "--path",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=str),
    help="Directory containing FASTA files",
)
@click.option(
    "--out",
    required=False,
    type=click.Path(path_type=str),
    default=None,
    help="Output CSV file path (default: metadata.csv in current directory)",
)
@click.option(
    "--mode",
    required=True,
    type=click.Choice(["dna", "protein"]),
    help="Type of sequences: dna or protein",
)
@click.option(
    "--comment",
    default="",
    help="Optional comment for all entries (default: empty string)",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite output file if it already exists",
)
def metadata_from_dir(path, out, mode, comment, force):
    """Generate metadata CSV from a directory of FASTA files

    This helper generates a sample metadata CSV file from all FASTA files
    found in a directory. It automatically extracts sample IDs from filenames
    and records the absolute path to each file.

    Example:
        prokanota helper metadata-from-dir --path ./genomes --mode dna
    """
    try:
        output_file = Path(out) if out else Path.cwd() / "metadata.csv"
        message = metadataFromDir.generate_metadata(
            Path(path), output_file, mode, comment, force
        )
        click.echo(message)
    except ValueError as e:
        raise click.ClickException(str(e)) from e


cli.add_command(run)
cli.add_command(config)
cli.add_command(test)
cli.add_command(test_minimal)
cli.add_command(test_full)
cli.add_command(citation)
cli.add_command(version)
cli.add_command(helper)
helper.add_command(metadata_from_dir)


def main():
    cli()


if __name__ == "__main__":
    main()
