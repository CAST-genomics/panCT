"""
Utilities for panct package
"""

import re
from typing import List


class Region:
    """
    Store information about a genomic region

    Attributes
    ----------
    chrom : str
        Chromosome
    start : int
        Start coordinate
    end : int
        End coordinate
    """

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end


def parse_regions_file(region_file: str) -> List[Region]:
    """
    Extract list of regions from BED file

    Parameters
    ----------
    region_file : str
       BED file of regions

    Returns
    -------
    regions_list : list of Region

    Raises
    ------
    ValueError
       If a region line could not be parsed to
       chrom, start, end from the first 3 columns
    """
    regions_list = []
    with open(region_file, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            if len(items) < 3:
                raise ValueError(f"Improperly formatted region line: {line}")
            chrom = items[0]
            try:
                start = int(items[1])
            except ValueError:
                raise ValueError(f"Improper start coordinate on line: {line}")
            try:
                end = int(items[2])
            except ValueError:
                raise ValueError(f"Improper end coordinate on line: {line}")
            regions_list.append(Region(chrom, start, end))
    return regions_list


def parse_region_string(region_string) -> Region:
    """
    Extract chrom, start, end from
    coordinate string

    Parameters
    ----------
    region_string : str
       Coordinate string in the form
       chrom:start-end

    Returns
    -------
    region : Region
       Region object

    Raises
    ------
    ValueError
       If the region region string could not be parsed
    """
    if type(region_string) != str:
        raise ValueError(f"Problem parsing coordinates {region_string}. Invalid type")
    if re.match(r"\w+:\d+-\d+", region_string) is None:
        raise ValueError(f"Problem parsing coordinates {region_string}")
    chrom = region_string.split(":")[0]
    start = int(region_string.split(":")[1].split("-")[0])
    end = int(region_string.split(":")[1].split("-")[1])
    if start >= end:
        raise ValueError(f"Problem parsing coordinates {region_string}. start>=end")
    return Region(chrom, start, end)
