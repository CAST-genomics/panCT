"""
Extract walks (W lines) from a GFA file into an indexed tab-separated format
"""

import logging
import subprocess
from pathlib import Path

from pysam import tabix_compress, tabix_index

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
        The location to which to write output. If not specified, we use the path to
        the graph, but with a .walk.gz file ending, instead.
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = getLogger(name="walks", level="ERROR")

    if output is None:
        if graph.suffix == ".gz":
            output = graph.with_suffix("").with_suffix(".walk.gz")
        else:
            output = graph.with_suffix(".walk.gz")

    also_index = False
    if output.suffix == ".gz":
        also_index = True
        output = output.with_suffix("")
    output_exists = output.exists()

    # what is the path to the shell script build_node_sample_map.sh ?
    script_path = Path(__file__).parent / "build_node_sample_map.sh"

    subprocess.run(
        [script_path, graph, output],
        capture_output=True,
        text=True,
        check=True,
    )

    # bgzip and tabix index the resulting file
    if also_index:
        gz_file = output.with_suffix(".walk.gz")
        tabix_compress(str(output), str(gz_file), force=True)
        if not output_exists:
            output.unlink()
        try:
            tabix_index(str(gz_file), preset="bed", force=True)
        except OSError as e:
            # check if the error message matches what we expect if the file is unsorted
            if str(e).startswith("building of index for "):
                log.error("Indexing failed. Are your walks properly formatted?")
            else:
                # otherwise, re-raise it
                raise
