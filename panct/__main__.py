#!/usr/bin/env python
import typer

app = typer.Typer()


@app.command()
def main(gfa_file: str, output_file: str):
    # Import and run the function from gfa_to_complexity.py
    from .gfa_to_complexity import main as gfa_main
    import sys

    sys.argv = ["calculate_complexity_gfa.py", gfa_file, output_file]
    gfa_main()


if __name__ == "__main__":
    # Run the CLI when called from the command line or via python -m panCT
    app(prog_name="panCT")
