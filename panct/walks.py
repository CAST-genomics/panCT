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
    Creates a VCF composed of haplotypes

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
        output = graph.with_suffix(".walks")

    result = subprocess.run(
        ["./build_node_sample_map.sh", graph, output], capture_output=True, text=True
    )

    return result
