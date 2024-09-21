"""
Utilities for panct package
"""

import re
import sys

def ParseRegionsFile(region_file):
    """
    Extract list of regions from BED file

    Parameters
    ----------
    region_file : str
       BED file of regions

    Returns
    -------
    regions_list : list of [[str, int, int]]
       with chrom, start, end of each region
    """
    regions_list = []
    with open(region_file, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            chrom = items[0]
            start = int(items[1])
            end = int(items[2])
            regions_list.append([chrom, start, end])
    return regions_list

def ParseRegionString(region, log):
    """
    Extract chrom, start, end from
    coordinate string

    Parameters
    ----------
    regions : str
       Coordinate string in the form
       chrom:start-end

    Returns
    -------
    region : list of [str, int, ind]
       Elements are chromosome, start, and end
       If we encounter an error parsing, then return None
    """
    if type(region) != str:
        return None
    if re.match(r"\w+:\d+-\d+", region) is None:
        return None
    chrom = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split(":")[1].split("-")[1])
    if start >= end:
        log.critical(r"Problem parsing coordinates {coords}. start>=end")
        return None
    return [chrom, start, end]