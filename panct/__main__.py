#!/usr/bin/env python

from enum import Enum
from pathlib import Path
from typing_extensions import Annotated

import typer

from . import __version__
from .complexity import AVAILABLE_METRICS

app = typer.Typer()


class Verbosity(str, Enum):
    critical = "CRITICAL"
    error = "ERROR"
    waring = "WARNING"
    info = "INFO"
    debug = "DEBUG"
    notset = "NOTSET"


verbose = Annotated[
    Verbosity,
    typer.Option(
        "-v",
        "--verbosity",
        case_sensitive=False,
        help="The level of verbosity desired",
    ),
]


def version_callback(value: bool = False):
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option(
            "-v",
            "--version",
            callback=version_callback,
            is_eager=True,
            help="Show the application's version and exit.",
        ),
    ] = False,
):
    """
    panct: A collection of tools for working with pangenomes
    """
    pass


@app.command()
def complexity(
    graph: Annotated[
        Path,
        typer.Argument(
            exists=True,
            readable=True,
            help="Path to the .gfa or .gbz file of a pangenome graph",
        ),
    ],
    region: Annotated[
        str, typer.Option("--region", help="Region to compute complexity over")
    ] = "",
    region_file: Annotated[
        Path,
        typer.Option(
            "--region-file",
            exists=True,
            readable=True,
            help="Bed file of regions to compute complexity over",
        ),
    ] = None,
    metrics: Annotated[
        str,
        typer.Option(
            "--metrics",
            help="Comma-separated list of which "
            "complexity metrics to compute. "
            "Options: " + ",".join(AVAILABLE_METRICS),
        ),
    ] = "sequniq-normwalk",
    reference: Annotated[
        str,
        typer.Option(
            "--reference", help="The ID of the reference sequence in the GFA file"
        ),
    ] = "GRCh38",
    output_file: Annotated[
        Path, typer.Option("--out", help="Name of output file")
    ] = Path("-"),
    verbosity: verbose = Verbosity.info,
):
    """
    Compute complexity scores
    """
    from .complexity import main as complexity_main
    from .logging import getLogger

    log = getLogger(name="complexity", level=verbosity.value)
    retcode = complexity_main(
        graph, str(output_file), region, region_file, metrics, reference, log
    )
    if retcode != 0:
        typer.Exit(code=retcode)


# Adding dummy command for now
# Removing this breaks Typer commands?
@app.command()
def walks(
    graph: Annotated[
        Path,
        typer.Argument(
            help="Path to the .gfa file of a pangenome graph",
            exists=True,
            readable=True,
        ),
    ],
    output_file: Annotated[
        Path, typer.Option("-o", "--out", help="Name of output file")
    ] = None,
    verbosity: verbose = Verbosity.info,
):
    """
    Extract walks to a file
    """
    from .walks import extract_walks
    from .logging import getLogger

    log = getLogger(name="walks", level=verbosity.value)

    extract_walks(graph, output_file, log)


typer_click_object = typer.main.get_command(app)


if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panct
    app(prog_name="panct")
