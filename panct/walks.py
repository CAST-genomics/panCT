"""
Extract walks (W lines) from a GFA file into an indexed tab-separated format
"""

import logging
import subprocess
from shutil import which
from pathlib import Path
from tempfile import NamedTemporaryFile

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
    elif output == Path("/dev/stdout") or output == Path("-"):
        output = Path("")

    also_index = False
    tmp_walk = None
    gz_file = output
    if output.suffix == ".gz":
        also_index = True
        if which("bgzip") is None:
            # if bgzip isn't installed, we need to bgzip the output ourselves
            # so we first create an uncompressed file in $TMPDIR called tmp_walk
            gz_file = output.with_suffix("").with_suffix(".walk.gz")
            tmp_walk = NamedTemporaryFile(delete=False)
            output = Path(tmp_walk.name)

    # what is the path to the shell script build_node_sample_map.sh ?
    script_path = Path(__file__).parent / "build_node_sample_map.sh"

    args = [script_path, graph]
    # is the output stdout?
    if output != Path(""):
        args.append(output)

    log.info("Building a mapping of nodes to samples")
    subprocess.run(args, text=True, check=True)

    # bgzip and tabix index the resulting file
    if also_index:
        if tmp_walk is not None:
            log.info("Bgzipping the output file")
            tabix_compress(str(output), str(gz_file), force=True)
            # now, properly close and delete the temporary file
            tmp_walk.close()
            output.unlink()
        try:
            log.info("Indexing the output file")
            tabix_index(str(gz_file), seq_col=0, start_col=1, end_col=1, force=True)
        except OSError as e:
            # check if the error message matches what we expect if the file is unsorted
            if str(e).startswith("building of index for "):
                log.error("Indexing failed. Are your walks properly formatted?")
            else:
                # otherwise, re-raise it
                raise
