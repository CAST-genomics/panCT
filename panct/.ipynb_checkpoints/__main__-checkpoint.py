#!/usr/bin/env python

from enum import Enum
from pathlib import Path
from typing_extensions import Annotated

import typer

from . import __version__
from .complexity import AVAILABLE_METRICS as complexity_metrics
from .population_uniqueness import AVAILABLE_METRICS as pop_uniq_metrics

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
    panct: A collection of tools for working with pangenomes graphs.
    Designed for use with the human pangenome.
    """
    pass


@app.command()
def complexity(
    graph: Annotated[
        Path,
        typer.Option("-g",
            "--graph_file",
            metavar="PATH",
            show_default=False,
            help="Path to the .gfa or .gbz file of a pangenome graph",
        ),
    ] = ..., 
    walk_file: Annotated[
        Path,
        typer.Option("-w",
            "--walk_file",
            help="Path to the walk file of a pangenome graph",
        ),
    ] = None,
    region: Annotated[
        str,
        typer.Option("-reg",
            "--region",
            help="A region in which to compute complexity, or a BED file of regions",
        ),
    ] = "",
    exclude_samples: Annotated[
        str,
        typer.Option("-e",
            "--exclude",
            help="Comma separated list of samples to exclude from analysis. ",
        ),
    ] = "GRCh38,CHM13",
    metrics: Annotated[
        str,
        typer.Option("-m",
            "--metrics",
            help="Comma-separated list of which "
            "complexity metrics to compute. "
            "Options: " + ", ".join(complexity_metrics),
        ),
    ] = "sequniq-normwalk",
    reference: Annotated[
        str,
        typer.Option(
            "-ref", "--reference", help="The ID of the reference sequence in the GFA file"
        ),
    ] = "GRCh38",
    output_file: Annotated[
        Path, typer.Option("-o", "--out", help="Name of output file")
    ] = Path("/dev/stdout"),
    verbosity: verbose = Verbosity.info,
):
    """
    Compute complexity scores
    """
    from .complexity import main as complexity_main
    from .logging import getLogger

    log = getLogger(name="complexity", level=verbosity.value)
    region_str = region
    if region == "":
        region_str = None
    elif Path(region).exists():
        region_str = Path(region)
    retcode = complexity_main(graph, output_file, region_str, metrics, reference, exclude_samples, walk_file,log)
    if retcode != 0:
        typer.Exit(code=retcode)


        
@app.command(name="population_uniqueness")
def population_uniqueness(
    graph: Annotated[
        Path,
        typer.Option("-g",
            "--graph_file",
            metavar="PATH",
            show_default=False,
            help="Path to the .gfa or .gbz file of a pangenome graph",
        ),
    ] = ...,    
    assemblies_file: Annotated[
        Path,
        typer.Option("-a",
            "--assemblies_file",
            metavar="PATH",
            show_default=False,
            help="Path to the .tsv file that maps samples to population. "
            "Must contain columns 'Sample ID', 'Haplotype' and 'Population Abbreviation'. "
            "Can be downloaded from pangenome consortium data explorer.",
        ),
    ] = ...,
    walk_file: Annotated[
        Path,
        typer.Option("-w",
            "--walk_file",
            help="Path to the walk file of a pangenome graph. "
            "Required if using gbz file.",
        ),
    ] = None,
    region: Annotated[
        str,
        typer.Option("-reg",
            "--region",
            help="A region in which to compute population "
            "specific uniqueness, or a BED file of regions",
        ),
    ] = "",
    exclude_samples: Annotated[
        str,
        typer.Option("-e",
            "--exclude",
            help="Comma separated list of samples to exclude from analysis."
            "Use in case where assembly file lists samples not in graph file."
            "GRCh38,CHM13,HG00272,HG03492 recommended for pangenome v2.0",
        ),
    ] = "GRCh38,CHM13",
    metrics: Annotated[
        str,
        typer.Option("-m",
            "--metrics",
            help="Comma-separated list of which "
            "population specific uniqueness metrics to compute. "
            "Options: " + ", ".join(pop_uniq_metrics),
        ),
    ] = "popuniq-normwalk",
    reference: Annotated[
        str,
        typer.Option(
            "-ref", "--reference", help="The ID of the reference sequence in the GFA file"
        ),
    ] = "GRCh38",
    output_file: Annotated[
        Path, typer.Option("-o", "--out", help="Name of output file")
    ] = Path("/dev/stdout"),
    verbosity: verbose = Verbosity.info,
):
    """
    Compute population specific sequence uniqueness scores
    """
    from .population_uniqueness import main as population_uniqueness_main
    from .logging import getLogger

    log = getLogger(name="population_uniqueness", level=verbosity.value)
    region_str = region
    if region == "":
        region_str = None
    elif Path(region).exists():
        region_str = Path(region)
    retcode = population_uniqueness_main(graph, output_file, region_str, metrics, reference, exclude_samples, walk_file, assemblies_file, log)
    if retcode != 0:
        typer.Exit(code=retcode)
