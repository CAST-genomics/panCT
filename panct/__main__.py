#!/usr/bin/env python
import typer
from rich import print

app = typer.Typer()


@app.command()
def main(gfa_file: str, output_file: str):
    # Import and run the function from gfa_to_complexity.py
    from .calculate_complexity_gfa import complexity_score as gfa_main
    print(gfa_main(gfa_file, output_file))


if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panct
    app(prog_name="panct")
