"""
Utilities for processing .walk files
"""

from __future__ import annotations
from typing import Type
from pathlib import Path
from logging import getLogger, Logger
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
        #super().__init__(log=log)
        self.log = log or getLogger(self.__class__.__name__)
        self.data = data

    def __len__(self):
        return len(self.data)

    @classmethod
    def read(
        cls: Type[Walks],
        fname: Path | str,
        region: str = None,
        nodes: set[str] = None,
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
            If specified, we will not load these samples
        log: Logger, optional
            A Logger object to use for debugging statements

        Returns
        -------
        Walks
            A Walks object loaded with a bunch of Node objects
        """
        if exclude_samples==None:
            exclude_samples=[]
        final_nodes = {}
        parse_samp = lambda samp: (samp[0], int(samp[1].replace('chr','')))
        log.info(f'processing nodes {region}')
        # Try to read the file with tabix
        if Path(fname).with_suffix(".gz.tbi").exists() and region is not None:
            # preprocess the region into a tabix region string
            region_str = ":" + region
            # iterate over the lines using tabix
            try:
                with TabixFile(filename=str(fname)) as f:
                    log.warning('Reading walks using tabix file.')
                    for line in f.fetch(region=region_str):
                        samples = line.strip().split("\t")
                        node = samples[0]
                        if ((nodes is not None) and (node not in nodes)):
                            continue

                        t=[]
                        for samp in samples:
                            s=samp.split('#')[0]
                            #get the sample ID (no haplotype)
                            if not s in list(exclude_samples)+[node]:
                                t.append(samp.split('#')[0]+'.'+samp.split('#')[1])
                        final_nodes[node]= (t)
                        #print(s)

                if ((log is not None)
                    and (nodes is not None)
                    and (len(final_nodes) < len(nodes))
                ):
                    log.warning(f"Couldn't load all requested nodes. "
                               f"{len(set(nodes).difference(final_nodes))/len(list(nodes))*100}% missing.")
            except ValueError:
                log.error('Failed to open tabix file.')
                pass
                
        # If we couldn't parse with tabix, then fall back to slow loading
        else:
            log.warning('Tabix file not found. Using walk file.'
                     ' This will be less efficient.')
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
                    samples = line.strip().split("\t")
                    node = samples[0]
                    if ((nodes is not None) and (node not in nodes)):
                        continue
    
                    t=[]
    
                    final_nodes[node]= (t)
                    if int(node) < start or int(node) > end:
                        continue
                    t=[]
                    for samp in samples:
                        s=samp.split('#')[0]
                        #print(s)
                        if not s in exclude_samples+[node]:
                            t.append(samp.split('#')[0]+'.'+samp.split('#')[1])
                    final_nodes[node]=t
            log.info('walk attribution calculation completed.')
        return cls(final_nodes, log)
