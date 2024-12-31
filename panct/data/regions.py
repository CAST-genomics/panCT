"""
Utilities for panct package
"""

from __future__ import annotations
import re
from pathlib import Path
from typing import Iterator, Type
from logging import getLogger, Logger

import numpy as np
import numpy.typing as npt

from .data import Data


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

    @classmethod
    def load(cls: Type[Region], region: str) -> Region:
        """
        Extract chrom, start, end from coordinate string

        Parameters
        ----------
        region : str
            Coordinate string in the form 'chrom:start-end'

        Returns
        -------
        region : Region
            Region object

        Raises
        ------
        ValueError
            If the region region string could not be parsed
        """
        if type(region) != str:
            raise ValueError(f"Problem parsing coordinates {region}. Invalid type")
        if re.match(r"\w+:\d+-\d+", region) is None:
            raise ValueError(f"Problem parsing coordinates {region}")
        chrom = region.split(":")[0]
        start = int(region.split(":")[1].split("-")[0])
        end = int(region.split(":")[1].split("-")[1])
        if start >= end:
            raise ValueError(f"Problem parsing coordinates {region}. start>=end")
        return cls(chrom, start, end)


class Regions(Data):
    """
    Store a bunch of Regions

    Attributes
    ----------
    data : tuple[Region]
        A bunch of Region objects
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(self, data: tuple[Region], log: Logger = None):
        super().__init__(log=log)
        self.data = data

    def __len__(self):
        return len(self.data)

    def __iter__(self) -> Iterator[Region]:
        return self.data.__iter__()

    def __getitem__(self, index):
        return self.data[index]

    @classmethod
    def load(cls: Type[Regions], fname: Path | str, log: Logger = None) -> Regions:
        """
        Extract list of regions from BED file

        Parameters
        ----------
        fname : Path | str
            BED file of regions

        Returns
        -------
        Regions
            A Regions object loaded with a bunch of regions

        Raises
        ------
        ValueError
            If a region line could not be parsed to chrom, start, end from the first
            3 columns
        """
        regions = []
        with open(fname, "r") as f:
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
                regions.append(Region(chrom, start, end))
        return cls(tuple(regions))
