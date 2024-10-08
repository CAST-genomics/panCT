#!/usr/bin/env python

from pathlib import Path
import typer
from typing_extensions import Annotated
from .complexity import AVAILABLE_METRICS

app = typer.Typer()


@app.command()
def complexity(
    output_file: Annotated[str, typer.Option("--out", help="Name of output file")],
    graph: Annotated[
        Path,
        typer.Argument(
            help="Path to the .gfa or .gbz file of a pangenome graph",
            exists=True,
            readable=True,
        ),
    ] = "",
    region: Annotated[
        str, typer.Option("--region", help="Region to compute complexity over")
    ] = "",
    region_file: Annotated[
        str,
        typer.Option(
            "--region-file", help="Bed file of regions to compute complexity over"
        ),
    ] = "",
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
        str, typer.Option("--reference", help="Which sequence to use as the reference")
    ] = "GRCh38",
    verbosity: str = "INFO",
):
    """
    Compute complexity scores
    """
    from .complexity import main as complexity_main
    from .logging import getLogger

    log = getLogger(name="complexity", level=verbosity)
    retcode = complexity_main(
        graph, output_file, region, region_file, metrics, reference, log
    )
    if retcode != 0:
        typer.Exit(code=retcode)


# Adding dummy command for now
# Removing this breaks Typer commands?
@app.command()
def map():
    print("Not implemented")


typer_click_object = typer.main.get_command(app)
if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panct
    app(prog_name="panct")
