"""
Compute complexity scores for regions
of a pangenome graph
"""

import time
import logging
from pathlib import Path
from typing import Optional

import numpy as np

from .logging import getLogger
from . import gbz_utils as gbz
from .data import Region, Regions
from . import graph_utils as gutils

AVAILABLE_METRICS = ['sequniq-normwalk', 'sequniq-normnode','raw-percentage','sequniq-unnorm','sequniq-normdegree']


def main(
    graph_file: Path,
    output_file: Path = Path("/dev/stdout"),
    region_str: str | Path = None,
    metrics: str = "sequniq-normwalk",
    reference: str = "GRCh38",
    exclude_samples: str = "GRCh38,CHM13",
    walk_file : Path= None,
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
    output_file : str, optional
        Path to output file
    region_str : str|Path, optional
        chrom:start-end of region to process or a BED file of regions
    metrics : str, optional
        Comma-separated list of metrics to compute 
        sequniq-normwalk, sequniq-normnode, sequniq-normdegree: intended for computing complexity
        raw-percentage, sequniq-unnorm: intended for QC and debugging
    reference : str, optional
        Sample ID of reference
    walk_file : Path, optional
        Path to the associated walk file for assembly. 
        Needed for accurate walk computations.
    log : logging.Logger, optional
        Logger object

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
        # TODO: also handle .gfa.gz
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
    #to do: remove median metrics metrics
    header.extend(["numnodes", "total_length", "numwalks", "median_length", 
                   "median_length_drop_SNV","mean_length", "mean_length_drop_SNV", 
                   "ratio_SNV", "ratio_min50bp","mean_degree"] + metrics_list)
    outf.write("\t".join(header) + "\n")

    ##### If GFA, just process the whole graph #####
    if file_type == "gfa":
        if region_str is not None:
            log.warning("Regions are ignored when processing GFA")
        exclude_samples=str.split(exclude_samples,',')
        if reference != "":
            exclude_samples =set(exclude_samples+[reference])
            log.info(f'filtering out the following samples: {exclude_samples}')
        node_table = gutils.NodeTable(graph_file, exclude_samples,walk_file)
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
        log.debug(f"Total time: \t{total_time}\n")
        outf.close()
        return 0

    #### If GBZ: Set up list of regions to process #####
    regions = []
    if region_str is not None:
        if isinstance(region_str, Path):
            regions = Regions.read(region_str, log=log)
        else:
            region = Region.read(region_str)
            regions = Regions((region,), log=log)
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
        node_table = gbz.load_node_table_from_gbz(graph_file, region, reference, exclude_samples, walk_file)

        # Load link table for the region
        if (node_table.gfa_file!=None):
            link_table=gutils.LinkTable(node_table.gfa_file,reference)
        else:
            'Node table does not contain gfa file, using gbz file for link table.'
            link_table = gbz.load_link_table_from_gbz(graph_file, region, reference, exclude_samples, walk_file)
        # Compute each requested complexity metric
        metric_results = []
        for m in metrics_list:
            metric_results.append(compute_complexity(node_table,link_table, m))

        # Output
        items = (
            [region.chrom, region.start, region.end]
            + [
                len(node_table.nodes.keys()),
                node_table.get_total_node_length(),
                node_table.numwalks,
                node_table.get_median_node_length(),
                node_table.get_median_node_length_dropSNV(),
                (node_table.get_total_node_length()/len(node_table.nodes.keys())),
                node_table.get_mean_node_length_dropSNV(),
                (node_table.get_number_SNVs()/len(node_table.nodes.keys())),
                (node_table.get_number_min50bp()/len(node_table.nodes.keys())),
                2*len(list(link_table.links.keys()))/len(list(node_table.nodes.keys())),
            ]
            + metric_results
        )

        outf.write("\t".join([str(item) for item in items]) + "\n")
        outf.flush()

    ##### Cleanup #####
    end_time = time.time()
    time_per_region = (end_time - start_time) / len(regions)
    log.debug(f"Time per region\t{time_per_region}\n")
    outf.close()
    return 0

def compute_complexity(node_table: gutils.NodeTable,link_table: gutils.LinkTable, metric: str) -> Optional[float]:
    """
    Compute complexity for a node table. Options:

    sequniq-normwalk: sum_n  len(n)*p_n*(1-p_n)/L
       where L is the average walk length

    sequniq-normnode: sum_n  len(n)*p_n*(1-p_n)/L
       where L is the average node length

    raw-percentage: p_n(1-p_n)
        for troubleshooting
    sequniq-unnorm: len(n)*p_n*(1-p_n)
        for troubleshooting
        
    Parameters
    ----------
    node_table : graph_utils.NodeTable
       Stores info on lengths/walks through each node
    link_table : graph_utils.LinkTable
        Stores info on the links within the pangenome region
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
    if metric in ('sequniq-normwalk', 'sequniq-normnode','raw-percentage','sequniq-unnorm','sequniq-normdegree'):
        for n in node_table.nodes.keys():
            length = node_table.nodes[n].length
            p = len(node_table.nodes[n].samples) / node_table.numwalks
            if(metric=='raw-percentage'):
                complexity += p * (1 - p)
            elif(metric=='sequniq-normdegree'):
                degree=0
                for l in link_table.links.keys():
                    if ((n==link_table.links[l].node_1) | (n==link_table.links[l].node_2)):
                        degree+=1                
                complexity += degree * p * (1 - p)
            else:
                complexity += length * p * (1 - p)
        # Normalize
        if metric == "sequniq-normwalk":
            complexity = complexity / node_table.get_mean_walk_length()
        elif metric == "sequniq-normnode":
            complexity = complexity / node_table.get_mean_node_length()
        elif metric == 'raw-percentage':
            complexity = complexity
        elif metric == 'sequniq-unnorm':
            complexity = complexity
        elif metric == 'sequniq-normdegree':
            complexity=complexity/(2*len(list(link_table.links.keys()))/len(list(node_table.nodes.keys())))
            
        return complexity
    raise ValueError(f"Invalid metric {metric}")