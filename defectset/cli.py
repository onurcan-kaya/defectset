from pathlib import Path

import click

from .config import Config
from .pipeline import run_pipeline


@click.group()
@click.version_option()
def cli() -> None:
    """defectset: filter defectool structures and DFT-label them for MLIP training."""


@cli.command("run")
@click.argument("config_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--work-dir",
    "-w",
    type=click.Path(file_okay=False),
    default="defectset_work",
    show_default=True,
    help="Working directory for logs, QE runs, cache, and output datasets.",
)
@click.option(
    "--skip-defectool",
    is_flag=True,
    help="Do not invoke `defectool`; just read the existing defectool DB.",
)
def run_cmd(config_path: str, work_dir: str, skip_defectool: bool) -> None:
    """Run the full pipeline from a defectset YAML config."""
    cfg = Config.from_yaml(config_path)
    if skip_defectool:
        cfg.skip_defectool = True
    run_pipeline(cfg, work_dir=Path(work_dir))


if __name__ == "__main__":
    cli()
