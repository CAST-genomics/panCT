"""
Utilities for processing .walk files
"""

from __future__ import annotations
from typing import Type
from pathlib import Path
from logging import Logger
from collections import Counter

from pysam import TabixFile

from .data import Data


class Walks(Data):
    """
    Store walks from a .walk file

    Attributes
    ----------
    data : dict[int, Counter[tuple[str, int]]]
        A bunch of nodes, stored as a mapping of node IDs to tuples of (sample labels, haplotype ID)
    log: Logger
        A logging instance for recording debug statements.
    """

    def __init__(self, data: dict[int, Counter[tuple[str, int]]], log: Logger = None):
        super().__init__(log=log)
        self.data = data

    def __len__(self):
        return len(self.data)

    @classmethod
    def read(
        cls: Type[Walks],
        fname: Path | str,
        region: str = None,
        nodes: set[int] = None,
        exclude_samples: set[str] = set(),
        log: Logger = None,
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
        nodes: set[int], optional
            A subset of nodes to load. Defaults to all nodes.
        exclude_samples: set[str], optional
            If specifieed, we will not load these samples
        log: Logger, optional
            A Logger object to use for debugging statements

        Returns
        -------
        Walks
            A Walks object loaded with a bunch of Node objects
        """
        final_nodes = {}
        parse_samp = lambda samp: (samp[0], int(samp[1]))
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
                        if nodes is not None and node not in nodes:
                            continue
                        final_nodes[node] = Counter(
                            s
                            for samp in samples
                            if (s := parse_samp(samp.rsplit(":", 1)))[0]
                            not in exclude_samples
                        )
                if (
                    log is not None
                    and nodes is not None
                    and len(final_nodes) < len(nodes)
                ):
                    log.warning("Couldn't load all requested nodes")
                return cls(final_nodes, log)
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
        with cls.hook_compressed(fname, "r") as f:
            for line in f:
                samples = str(line.strip())
                node = int(samples.split("\t", maxsplit=1)[0])
                if (node < start or node > end) or (
                    nodes is not None and node not in nodes
                ):
                    continue
                final_nodes[node] = Counter(
                    s
                    for samp in samples.split("\t")[1:]
                    if (s := parse_samp(samp.rsplit(":", 1)))[0] not in exclude_samples
                )
        return cls(final_nodes, log)
