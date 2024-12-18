"""
Utilities for dealing with node tables
"""

import numpy as np


class Node:
    """
    Stores metadata about a node in the graph

    Attributes
    ----------
    nodeid : str
        ID of the node
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

    def __init__(self, nodeid, length=0):
        self.nodeid = nodeid
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
    walk_lengths : list[int]
        List of lengths of walks through this region

    Methods
    -------
    load_from_gfa(gfafile)
        Generate NodeTable from GFA file
    add_node(node)
        Add node to the table
    add_walk(sampid, nodelist)
        Add a walk to the node table
    get_walk_length(nodelist=[])
        Get the total length of a walk
        through the given list of nodes
    get_mean_walk_length()
        Get mean length of all walks
    get_mean_node_length()
        Get mean length of all nodes
    get_total_node_length()
        Get total length of all nodes
    get_nodes_from_walk(walk_string)
        Get list of nodes from the walk
    """

    def __init__(self, gfa_file: str = None, exclude_samples: list[str] = []):
        self.nodes = {}  # node ID-> Node
        self.numwalks = 0
        self.walk_lengths = []
        if gfa_file is not None:
            self.load_from_gfa(gfa_file, exclude_samples)

    def add_node(self, node: Node):
        """
        Add a node to the node table

        Parameters
        ----------
        node : Node
            Node to add
        """
        self.nodes[node.nodeid] = node

    def add_walk(self, sampid: str, nodelist: list[str]):
        """
        Add a walk to the node table

        Parameters
        ----------
        sampid : str
            ID of the walk
        nodelist : list[str]
        """
        self.walk_lengths.append(self.get_walk_length(nodelist))
        for n in nodelist:
            self.nodes[n].add_sample(sampid)
        self.numwalks += 1

    def get_walk_length(self, nodelist: list[str]) -> int:
        """
        Get the total length of a walk
        through the given list of nodes

        Parameters
        ----------
        nodelist : list[str]
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
            if n not in self.nodes.keys():
                raise ValueError(f"Encountered unknown node {n}")
            length += self.nodes[n].length
        return length

    def get_mean_walk_length(self) -> float:
        """
        Get mean length of all walks

        Returns
        -------
        mean_walk_length : float
            Returns np.nan if there are no walks
        """
        if self.numwalks == 0:
            return np.nan
        return float(np.mean(self.walk_lengths))

    def get_mean_node_length(self) -> float:
        """
        Get mean length of all nodes

        Returns
        -------
        mean_node_length : float
            Returns np.nan if there are no walks
        """
        if len(self.nodes.keys()) == 0:
            return np.nan
        return float(np.mean([n.length for n in self.nodes.values()]))

    def get_total_node_length(self) -> int:
        """
        Get total length of all nodes

        Returns
        -------
        total_walk_length : int
        """
        return np.sum([n.length for n in self.nodes.values()])

    def get_nodes_from_walk(self, walk_string: str) -> list[str]:
        """
        Get list of nodes from a walk string

        Parameters
        ----------
        walk_string : str
            Walk string from a GFA file

        Returns
        -------
        list[str]
            List of node IDs
        """
        ws = walk_string.replace(">", ":").replace("<", ":").strip(":")
        return ws.split(":")

    def load_from_gfa(self, gfa_file: str, exclude_samples: list[str] = []):
        # First parse all the nodes
        with open(gfa_file, "r") as f:
            for line in f:
                linetype = line.split()[0]
                if linetype != "S":
                    continue
                nodeid = line.split()[1]
                nodelen = 0
                nodeseq = line.strip().split()[2]
                if nodeseq.strip() != "*":
                    nodelen = len(nodeseq)
                else:
                    for var in line.strip().split()[3:]:
                        if var.startswith("LN"):
                            nodelen = int(var.split(":")[2])
                if nodelen == 0:
                    raise ValueError(f"Could not determine node length for {nodeid}")
                self.add_node(Node(nodeid, length=nodelen))

        # Second pass to get the walks
        with open(gfa_file, "r") as f:
            for line in f:
                linetype = line.split()[0]
                if linetype != "W":
                    continue
                sampid = line.split()[1]
                if sampid in exclude_samples:
                    continue
                hapid = line.split()[2]
                walk = line.split()[6]
                nodes = self.get_nodes_from_walk(walk)
                self.add_walk(f"{sampid}:{hapid}", nodes)
