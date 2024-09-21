"""
Utilities for dealing with node tables
"""

import numpy as np
from typing import List
from . import gbz_utils as gbz


class Node:
    """
    Stores metadata about a node in the graph

    Attributes
    ----------
    length : int
        Length of the sequence at this node
    samples : set of str
        IDs of samples (haplotypes) that go
        through this node

    Methods
    -------
    add_sample(sampid)
        Add a sample to the node
    """

    def __init__(self, length=0):
        self.length = length
        self.samples = set()

    def add_sample(self, sampid):
        """
        Add a sample to the node

        Parameters
        ----------
        sampid : str
            ID of the sample (haplotype) to add
        """
        self.samples.add(sampid)


class NodeTable:
    """
    Table of nodes storing node metadata
    for a region

    Attributes
    ----------
    nodes : dict[str]->Node
        Dictionary of nodes, indexed by node ID
    numwalks : int
        Number of walks going through this region
    path_lengths : list of in
        List of lengths of paths through this region

    Methods
    -------
    get_path_length(nodelist=[])
        Get the total length of a walk
        through the given list of nodes
    """

    def __init__(self):
        self.nodes = {}  # node ID-> Node
        self.numwalks = 0
        self.path_lengths = []

    def get_path_length(self, nodelist: List[Node]) -> int:
        """
        Get the total length of a walk
        through the given list of nodes

        Parameters
        ----------
        nodelist : list of Node
            List of nodes of the walk

        Returns
        -------
        length : int
            Length (bp) of the walk

        Raises
        ------
        ValueError
            If we encounter a node ID not in the NodeTable
        """
        length = 0
        for n in nodelist:
            if n not in self.nodes.values():
                raise ValueError(f"Encountered unknown node {n}")
            length += self.nodes[n].length
        return length

    def GetMeanPathLength(self):
        return np.mean(self.path_lengths)

    def GetMeanNodeLength(self):
        return np.mean([n.length for n in self.nodes.values()])

    def GetTotalNodeLength(self):
        return np.sum([n.length for n in self.nodes.values()])


def CheckNodeSeq(seq):
    for char in seq:
        if char.upper() not in ["A", "C", "G", "T", "N"]:
            return False
    return True


def GetNodesFromWalk(walk_string):
    ws = walk_string.replace(">", ":").replace("<", ":").strip(":")
    return ws.split(":")


def LoadNodeTableFromGFA(gfa_file, log, exclude_samples=[]):
    nodetable = NodeTable()

    # First parse all the nodes
    with open(gfa_file, "r") as f:
        for line in f:
            linetype = line.split()[0]
            if linetype != "S":
                continue
            nodeid = line.split()[1]
            nodelen = 0
            nodeseq = line.strip().split()[2]
            if CheckNodeSeq(nodeseq):
                nodelen = len(nodeseq)
            else:
                for var in line.strip().split()[3:]:
                    if var.startswith("LN"):
                        length = int(var.split(":")[2])
            if nodelen == 0:
                raise ValueError(f"Could not determine node length for {nodeid}")
            nodetable.nodes[nodeid] = Node(length=nodelen)

    # Second pass to get the walks
    numwalks = 0
    path_lengths = []
    with open(gfa_file, "r") as f:
        for line in f:
            linetype = line.split()[0]
            if linetype != "W":
                continue
            sampid = line.split()[1]
            if sampid in exclude_samples:
                continue
            numwalks += 1
            hapid = line.split()[2]
            walk = line.split()[6]
            nodes = GetNodesFromWalk(walk)
            path_lengths.append(nodetable.get_path_length(nodes))
            for n in nodes:
                nodetable.nodes[n].add_sample(f"{sampid}:{hapid}")
    nodetable.numwalks = numwalks
    nodetable.path_lengths = path_lengths
    return nodetable


def LoadNodeTableFromGBZ(gbz_file, region, reference, log):
    gfa_file = gbz.ExtractRegionFromGBZ(gbz_file, region, reference)
    if gfa_file is None:
        return None
    return LoadNodeTableFromGFA(gfa_file.name, log, exclude_samples=[reference])
