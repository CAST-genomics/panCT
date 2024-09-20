#!/usr/bin/env python
import typer
from typing_extensions import Annotated
from .complexity import AVAILALBE_METRICS

app = typer.Typer()

@app.command()
def complexity(
    gbz_file: Annotated[
        str, typer.Option("--gbz", help="Path to .gbz file of the graph")
    ],
    output_file: Annotated[
        str, typer.Option("--out", help="Name of output file")
    ],
    region: Annotated[
        str, typer.Option("--region", help="Region to compute complexity over")
    ] = "",
    region_file: Annotated[
        str, typer.Option("--region-file", help="Bed file of regions to compute complexity over")
    ] = "",
    metrics: Annotated[
        str, typer.Option("--metrics", help="Comma-separated list of which complexity metrics to compute. "
                          "Options: " + ",".join(AVAILALBE_METRICS))
    ] = "sequniq",
    reference: Annotated[
        str, typer.Option("--reference", help="Which sequence to use as the reference")
    ] = "GRCh38"
):
    """
    Compute complexity scores
    """
    from .complexity import main as complexity_main
    complexity_main(gbz_file, output_file, region, region_file, 
        metrics, reference)

@app.command()
def map(gfa_file: str, output_file: str):
    """
    Generate a node-map file
    """
    typer.echo("Generating")


typer_click_object = typer.main.get_command(app)
if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panct
    app(prog_name="panct")
