#!/usr/bin/env python
import typer

app = typer.Typer()


@app.command()
def complexity(gfa_file: str, output_file: str):
    """
    Compute complexity scores
    """
    from .calculate_complexity_gfa import complexity_score as gfa_main

    print(gfa_main(gfa_file, output_file))


@app.command()
def map(gfa_file: str, output_file: str):
    """
    Generate a node-map file
    """
    typer.echo("Generating")


if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panct
    app(prog_name="panct")
