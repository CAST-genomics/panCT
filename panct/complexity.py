"""
Compute complexity scores for regions
of a pangenome graph
"""

import os
import sys
import time
import logging
from pathlib import Path
from typing import Optional

from . import utils as utils
from .logging import getLogger
from . import gbz_utils as gbz
from . import graph_utils as gutils

AVAILABLE_METRICS = ["sequniq-normwalk", "sequniq-normnode"]


def main(
    graph_file: Path,
    output_file: str,
    region: str,
    region_file: str,
    metrics: str,
    reference: str,
    log: logging.Logger = None,
):
    """
    Compute complexity scores for regions
    of a pangenome graph

    If a GFA file is given, compute complexity
    on the entire file.

    If a GBZ file is given, must specify a region
    (or file with list of regions)

    Parameters
    ----------
    graph_file : Path
        Path to GFA or GBZ file
    output_file : str
        Path to output file
    region : str
        chrom:start-end of region to process
    region_file : str
        Path to BED file of regions to process
    metrics : str
        Comma-separated list of metrics to compute
    reference : str
        Sample ID of reference
    log : logging.Logger
        logger object

    Returns
    -------
    retcode : int
        Return code of the program
    """
    if log is None:
        log = getLogger(name="complexity", level="ERROR")
    start_time = time.time()

    #### Check files and indices #####
    file_type = None
    if graph_file.suffix == ".gfa":
        file_type = "gfa"
    elif graph_file.suffix == ".gbz":
        file_type = "gbz"
        if not gbz.check_gbzbase_installed(log):
            return 1
        if not gbz.check_gbzfile(graph_file, log):
            return 1
    else:
        log.critical("Invalid graph type. Must be .gbz or .gfa")
        return 1

    #### Check requested metrics #####
    metrics_list = metrics.split(",")
    for m in metrics_list:
        if m not in AVAILABLE_METRICS:
            log.critical(f"Encountered invalid metric {m}")
            return 1

    ##### Set up output file #####
    outf = open(output_file, "w")
    header = []
    if file_type == "gbz":
        header = ["chrom", "start", "end"]
    header.extend(["numnodes", "total_length", "numwalks"] + metrics_list)
    outf.write("\t".join(header) + "\n")

    ##### If GFA, just process the whole graph #####
    if file_type == "gfa":
        if region != "" or region_file != "":
            log.warning("Regions are ignored when processing GFA")
        exclude = []
        if reference != "":
            exclude = [reference]
        node_table = gutils.NodeTable(str(graph_file), exclude)
        metric_results = []
        for m in metrics_list:
            metric_results.append(compute_complexity(node_table, m))
        items = [
            len(node_table.nodes.keys()),
            node_table.get_total_node_length(),
            node_table.numwalks,
        ] + metric_results
        outf.write("\t".join([str(item) for item in items]) + "\n")
        outf.flush()
        end_time = time.time()
        total_time = end_time - start_time
        sys.stderr.write(f"Total time: \t{total_time}\n")
        outf.close()
        return 0

    #### If GBZ: Set up list of regions to process #####
    regions = []
    if region != "":
        regions.append(utils.parse_region_string(region))
    if region_file != "":
        if not os.path.exists(region_file):
            log.critical(f"Could not find {region_file}")
            return 1
        regions.extend(utils.parse_regions_file(region_file))
    if len(regions) == 0:
        log.critical("Did not detect any regions")
        return 1

    ##### Process each region #####
    for region in regions:
        log.info(
            "Processing region {chrom}:{start}-{end}".format(
                chrom=region.chrom, start=region.start, end=region.end
            )
        )
        # Load node table for the region
        node_table = gbz.load_node_table_from_gbz(graph_file, region, reference)

        # Compute each requested complexity metric
        metric_results = []
        for m in metrics_list:
            metric_results.append(compute_complexity(node_table, m))

        # Output
        items = (
            [region.chrom, region.start, region.end]
            + [
                len(node_table.nodes.keys()),
                node_table.get_total_node_length(),
                node_table.numwalks,
            ]
            + metric_results
        )
        outf.write("\t".join([str(item) for item in items]) + "\n")
        outf.flush()

    ##### Cleanup #####
    end_time = time.time()
    time_per_region = (end_time - start_time) / len(regions)
    sys.stderr.write(f"Time per region\t{time_per_region}\n")
    outf.close()
    return 0


def compute_complexity(node_table: gutils.NodeTable, metric: str) -> Optional[float]:
    """
    Compute complexity for a node table. Options:

    sequniq-normwalk: sum_n  len(n)*p_n*(1-p_n)/L
       where L is the average walk length

    sequniq-normnode: sum_n  len(n)*p_n*(1-p_n)/L
       where L is the average node length

    Parameters
    ----------
    node_table : graph_utils.NodeTable
       Stores info on lengths/walks through each node
    metric : str
       Which metric to compute. See description above

    Returns
    -------
    complexity : float
       Complexity score

    Raises
    ------
    ValueError
       If invalid metric specified
    """
    if node_table.numwalks == 0:
        return None
    complexity = 0
    # Add up value for each node
    if metric in ("sequniq-normwalk", "sequniq-normnode"):
        for n in node_table.nodes.keys():
            length = node_table.nodes[n].length
            p = len(node_table.nodes[n].samples) / node_table.numwalks
            complexity += length * p * (1 - p)
        # Normalize
        if metric == "sequniq-normwalk":
            complexity = complexity / node_table.get_mean_walk_length()
        elif metric == "sequniq-normnode":
            complexity = complexity / node_table.get_mean_node_length()
        return complexity
    raise ValueError(f"Invalid metric {metric}")
