"""
Utilities for processing .walk files
"""

from __future__ import annotations
from typing import Type
from pathlib import Path
from logging import Logger

from pysam import TabixFile

from .data import Data


class Node:
    """
    Stores metadata about a node in the graph

    Attributes
    ----------
    nodeid : int
        ID of the node
    samples : set of str
        IDs of samples (haplotypes) that pass through this node
    """

    def __init__(self, nodeid: int, samples: set[str] = set()):
        self.nodeid = nodeid
        self.samples = samples

    def add_sample(self, sampid: str):
        """
        Add a sample to the node

        Parameters
        ----------
        sampid : str
            ID of the sample (haplotype) to add
        """
        self.samples.add(sampid)


class Walks(Data):
    """
    Store walks from a .walk file

    Attributes
    ----------
    data : dict[str, Node]
        A bunch of Node objects, keyed by node ID
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(self, data: dict[str, Node], log: Logger = None):
        super().__init__(log=log)
        self.data = data

    def __len__(self):
        return len(self.data)

    @classmethod
    def read(
        cls: Type[Walks], fname: Path | str, region: str = None, log: Logger = None
    ) -> Walks:
        """
        Extract walks from a .walk file

        Parameters
        ----------
        fname: Path | str
            A .walk file of walks
        region: str, optional
            A region string denoting the start and end node IDs in the form
            of f'{start}-{end}'
        log: Logger, optional
            A Logger object to use for debugging statements

        Returns
        -------
        Walks
            A Walks object loaded with a bunch of Node objects
        """
        nodes = {}
        # Try to read the file with tabix
        if Path(fname).suffix == ".gz" and region is not None:
            # preprocess the region into a tabix region string
            region_str = ":" + region
            # iterate over the lines using tabix
            try:
                with TabixFile(filename=str(fname)) as f:
                    for line in f.fetch(region=region_str):
                        samples = line.strip().split("\t")
                        node = int(samples.pop(0))
                        nodes[node] = Node(node, set(samples))
                return cls(nodes, log)
            except ValueError:
                pass
        # If we couldn't parse with tabix, then fall back to slow loading
        # First, split the region into start and end coordinates
        start, end = -float("inf"), float("inf")
        if region is not None:
            start, end = tuple(
                (int(coord) if coord != "" else float("inf"))
                for coord in region.split("-")
            )
            if start == float("inf"):
                start = -start
        # Now iterate over the lines
        with open(fname, "r") as f:
            for line in f:
                samples = line.strip()
                node = int(samples.split("\t", maxsplit=1)[0])
                if node < start or node > end:
                    continue
                nodes[node] = Node(node, set(samples.split("\t")[1:]))
        return cls(nodes, log)
