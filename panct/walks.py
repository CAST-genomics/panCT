"""
Extract walks (W lines) from a GFA file into an indexed tab-separated format
"""

import logging
import subprocess
from pathlib import Path

from .logging import getLogger


def extract_walks(
    graph: Path,
    output: Path = None,
    log: logging.Logger = None,
):
    """
    Creates a .walk file mapping nodes in the graph to sample IDs representing
    haplotypes

    Parameters
    ----------
    graph : Path
        The path to a pangenome graph in GFA file
    output : Path, optional
        The location to which to write output
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = getLogger(name="walks", level="ERROR")

    if output is None:
        output = graph.with_suffix(".walk")

    if output.suffix == ".gz":
        output = output.with_suffix("")

    # what is the path to the shell script build_node_sample_map.sh ?
    script_path = Path(__file__).parent / "build_node_sample_map.sh"

    result = subprocess.run(
        [script_path, graph, output],
        capture_output=True,
        text=True,
        check=True,
    )

    # TODO: bgzip and tabix index the resulting file

    return result
