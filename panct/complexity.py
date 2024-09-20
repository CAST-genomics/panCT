# Taken from grant proposal section C.2.1 in Research Strategy
# We will define sequence uniqueness as U = ∑ s in S (|s|*p_s*(1 − p_s))/L
# S is the set of nodes in a region,
# |s| is the length in bp of node s,
# p_s is the percent of sequences that go through node s, and
# L is the average length in bp of all paths traversing the subgraph of interest.

# This metric is meant to capture the relative amount of sequence in a region that is shared vs. polymorphic amongst haplotypes in a region.

import os
import sys
import time
import subprocess

from . import utils as utils
from . import gbz_utils as gbz
from . import graph_utils as gutils

AVAILALBE_METRICS = ["sequniq","sequniq2"]

def main(gbz_file: str, output_file: str, 
    region: str, region_file: str, metrics: str,
    reference: str):
    start_time = time.time()

    #### Check GBZ file and index #####
    if not gbz.CheckGBZBaseInstalled():
        return 1
    if not gbz.CheckGBZFile(gbz_file):
        return 1

    #### Check requested metrics #####
    metrics_list = metrics.split(",")
    for m in metrics_list:
        if m not in AVAILALBE_METRICS:
            utils.WARNING(f"Encountered invalid metric {m}")
            return 1

    #### Set up list of regions to process #####
    regions = []
    if region != "":
        region = utils.ParseRegionString(region)
        if region is None:
            return 1
        else: regions.append(region)
    if region_file != "":
        if not os.path.exists(region_file):
            utils.WARNING(f"Could not find {region_file}")
            return 1
        regions.extend(utils.ParseRegionsFile(region_file))

    ##### Set up output file #####
    outf = open(output_file, "w")
    outf.write("\t".join(["chrom","start","end"] + metrics_list)+"\n")

    ##### Process each region #####
    for region in regions:
        # Load node table
        #node_table = gutils.LoadNodeTableFromGBZ(gbz_file, region, reference)
        node_table = gutils.LoadNodeTableFromGFA("test.gfa") # TODO remove

        # Compute each requested complexity metric
        metric_results = []
        for m in metrics_list:
            metric_results.append(ComputeComplexity(node_table, m))

        # Output
        items = region + metric_results
        outf.write("\t".join([str(item) for item in items])+"\n")
    ##### Cleanup #####
    end_time = time.time()
    time_per_region = (end_time - start_time) / len(regions)
    sys.stderr.write(f"Time per region\t{time_per_region}\n")
    outf.close()

def ComputeComplexity(node_table, metric):
    """
    Compute complexity for a node table. Options:

    sequniq: sum_n |n|*p_n*(1-p_n)/|L| where |L| is the 
       average path length
    sequniq2: sum_n |n|*p_n*(1-p_n)/|L| where |L| is the 
       average node length
    """
    complexity = 0
    # Add up value for each node
    for n in node_table.nodes.keys():
        if metric == "sequniq" or metric == "sequniq2":
            length = node_table.nodes[n].length
            p = len(node_table.nodes[n].samples)/node_table.numwalks
            complexity += length*p*(1-p)
        else:
            utils.WARNING(f"Encountered invalid metric {metric}")
            return None
    # Normalize
    if metric == "sequniq":
        complexity = complexity/node_table.GetMeanPathLength()
    elif metric == "sequniq2":
        complexity = complexity/node_table.GetMeanNodeLength()   
    return complexity