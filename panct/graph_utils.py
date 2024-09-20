"""
Utilities for extracting info from GFA
"""

from . import gbz_utils as gbz
from .utils import WARNING
import numpy as np

class Node:
    def __init__(self, length=0):
        self.length = length
        self.samples = set()

    def AddSample(self, sampid):
        self.samples.add(sampid)

class NodeTable:
    def __init__(self):
        self.nodes = {} # node ID-> Node
        self.numwalks = 0
        self.path_lengths = []

    def GetPathLength(self, nodelist):
        length = 0
        for n in nodelist:
            length += self.nodes[n].length
        return length

    def GetMeanPathLength(self):
        return np.mean(self.path_lengths)

    def GetMeanNodeLength(self):
        return np.mean([n.length for n in self.nodes.values()])

def CheckNodeSeq(seq):
    for char in seq:
        if char.upper() not in ["A","C","G","T","N"]:
            return False
    return True

def GetNodesFromWalk(walk_string):
    ws = walk_string.replace(">",":").replace("<",":").strip(":")
    return ws.split(":")

def LoadNodeTableFromGFA(gfa_file):
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
                WARNING(f"Could not determine node length for {nodeid}")
                return None
            nodetable.nodes[nodeid] = Node(length=nodelen)

    # Second pass to get the walks
    numwalks = 0
    path_lengths = []
    with open(gfa_file, "r") as f:
        for line in f:
            linetype = line.split()[0]
            if linetype != "W":                
                continue
            numwalks += 1
            sampid = line.split()[1]
            hapid = line.split()[2]
            walk = line.split()[6]
            nodes = GetNodesFromWalk(walk)
            path_lengths.append(nodetable.GetPathLength(nodes))
            for n in nodes:
                nodetable.nodes[n].AddSample(f"{sampid}:{hapid}")
    nodetable.numwalks = numwalks
    nodetable.path_lengths = path_lengths
    return nodetable

def LoadNodeTableFromGBZ(gbz_file, region, reference):
    gfa_file = gbz.ExtractRegionFromGBZ(gbz_file, region, reference)
    if gfa_file is None:
        return None
    return LoadNodeTableFromGFA(gfa_file.name)