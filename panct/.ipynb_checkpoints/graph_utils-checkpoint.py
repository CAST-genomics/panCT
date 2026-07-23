"""
Utilities for dealing with node tables
"""

from pathlib import Path
from logging import getLogger, Logger
import numpy as np

from .data import Walks
from collections import defaultdict

logger = getLogger(__name__)


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

    def __init__(self, nodeid, length=0, degree=0):
        self.nodeid = nodeid
        self.length = length
        self.samples = set()
        self.anc_count = dict()
        self.exp_het = dict()
        self.Fst = dict()
        self.degree= degree
    def add_sample(self, sampid):
        """
        Add a sample to the node

        Parameters
        ----------
        sampid : str
            ID of the sample (haplotype) to add
        """
        self.samples.add(sampid)
    def to_dict(self):
        return {
            "nodeid": self.nodeid,
            "length": self.length,
            "samples": self.samples,
            "anc_count":self.anc_count,
            "exp_het":self.exp_het,
            "Fst": self.Fst
        }



class NodeTable:
    """
    Table of nodes storing node metadata
    for a region

    Attributes
    ----------
    nodes : dict[str]->Node
        Dictionary of nodes, indexed by node ID.
    numwalks : int
        Number of walks going through this region.
    walk_lengths : list[int]
        List of lengths of walks through this region.
    walk_ids: list[str]
        List of the walk IDs through this region. 
        Maps by position in list to walk_lengths.
    gfa_filepath: str
        Path to the intermediate gfa file used to generate the node table. 
        Added to pass to other import table functions (i.e. Link_Table).
    sample_labels: bool
        Whether or not labelled walks will be used in the walk calculations.
        Labelled walks can be imported from walk files (passed in using walk_file variable).
        Walk files are calculated using ./walks.py, and should have a .walk.gz and a .walk.gz.tbi.
        If walk file not available, walks extracted from gbz, and samples will be unlabelled.
        If sample_labels==False using Human Pangenome v1.1, the walk stats will be inaccurate.
        If sample_labels==False for Human Pangenome v2.0, walks should be accurate.

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

    def __init__(self, gfa_file: Path = None, exclude_samples: list[str] = [], walk_file : Path = None):
        self.nodes = {}  # node ID-> Node
        self.numwalks = 0
        self.walk_lengths = []
        self.walk_ids = []
        if gfa_file is not None:
            self.load_from_gfa(gfa_file, exclude_samples,walk_file)
        self.gfa_file=gfa_file
        self.walk_file=walk_file
        self.sample_labels=False
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
        Add a walk to the node table if walk extracted from gbz file

        Parameters
        ----------
        sampid : str
            ID of the walk
        nodelist : list[str]
        """
        self.walk_lengths.append(self.get_walk_length(nodelist))
        self.walk_ids.append(sampid)
        for n in nodelist:
            self.nodes[n].add_sample(sampid)
        self.numwalks += 1
    def add_labelled_walks(self, walk_dict:dict, nodelist:list[str], exclude_samples:list[str]):
        """
        Add walks data to node table from walk file
        sampids: list[str]
            IDs of walks present in region
        nodelist : list[str]
            list of nodes in walk file
        exclude_samples: list[str]
            list of samples to exclude (typically reference)
        """
        for w in walk_dict.keys():
            if w in [exclude_samples]:  
                continue
            self.walk_lengths.append(self.get_walk_length(walk_dict[w]))
        self.walk_ids=list(walk_dict.keys())
        self.numwalks=len(list(walk_dict.keys()))

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
            else:
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
            Returns np.nan if there are no nodes
        mean_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan
        return float(np.mean([n.length for n in self.nodes.values()]))

    def get_mean_degree(self) -> float:
        """
        Get mean degree of all nodes

        Returns
        -------
            Returns np.nan if there are no nodes
        mean_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan
        return float(np.mean([n.degree for n in self.nodes.values()]))

    def get_total_node_length(self) -> int:
        """
        Get total length of all nodes
        Returns
        -------
        total_node_length : int
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
        
    def get_median_node_length(self) -> float:
        """
        Get median length of all nodes
        
        Returns
        -------
            Returns np.nan if there are no nodes
        median_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan
        return float(np.median([n.length for n in self.nodes.values()]))

    def get_median_node_length_dropSNV(self) -> float:
        """
        Get median length of all nodes longer than 1BP (aka remove SNVs)
        
        Returns
        -------
            Returns np.nan if there are no nodes
        median_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan        
        return float(np.median([n.length for n in self.nodes.values() if n.length!=1]))
                
    def get_mean_node_length_dropSNV(self) -> float:
        """
        Get mean length of all nodes longer than 1BP (aka remove SNVs)
        
        Returns
        -------
            Returns np.nan if there are no nodes
        mean_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan        
        return float(np.mean([n.length for n in self.nodes.values() if n.length!=1]))
        
    def get_number_SNVs (self) -> float:
        """
        Get the number of 1BP variants
        
        Returns
        -------
            Returns np.nan if there are no nodes
        mean_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan        
        return float(len([n.length for n in self.nodes.values() if n.length==1]))

    def get_number_min50bp (self) -> float:
        """
        Get the number of 1BP variants
        
        Returns
        -------
            Returns np.nan if there are no nodes
        mean_node_length : float
        """
        if len(self.nodes.keys()) == 0:
            return np.nan        
        return float(len([n.length for n in self.nodes.values() if n.length<=50]))
        
    def load_from_gfa(self, gfa_file: Path, exclude_samples: list[str] = [], walk_file: Path = None):
        #Does not automatically filter any samples out.
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

        if walk_file is None:
            if gfa_file.suffix == ".gz":
                walk_file = gfa_file.with_suffix("").with_suffix(".walk")
            else:
                walk_file = gfa_file.with_suffix(".walk")
            if not walk_file.exists():
                walk_file = walk_file.with_suffix(".walk.gz")
        walk_file=Path(walk_file)

        if walk_file.exists():
            logger.warning('Attempting to import walks from walk file.'
                          ' Recommended to use tabix file for efficient node-walk assignment.')
            node_set = set(self.nodes.keys())
            logger.info(f'{len(node_set)} nodes present in region')
        
            # find smallest and largest node for processing walks
            smallest_node = min(node_set, default="")
            largest_node = max(node_set, default="")
            # Get nodes from .walk file and add with self.add_walk()
            walks = Walks.read(
                walk_file,
                region=f"{smallest_node}-{largest_node}",
                nodes=node_set,
                exclude_samples=exclude_samples,
                log=logger)
            if all(not value for value in walks.data.values()):
                logger.warning('No walks imported from file. Consider using gbz file.')
            else:
                self.sample_labels=True
                for node in self.nodes:
                    if node in walks.data:
                        self.nodes[node].samples=list(set(walks.data[node]))
                        #assign walks to nodes
                walk_dict = defaultdict(list) #make a dictionary of walk:nodes
                for node, values in walks.data.items():
                    for value in values:
                        walk_dict[value].append(node)
                self.add_labelled_walks(walk_dict, node_set,exclude_samples)
        if (not walk_file.exists()):
            logger.warning('Walk file not found. Using walks present in gbz file. ' 
            'This may lead to inaccurate walk metrics.')
            #walk counts do not match for if you use walk file and gbz. 
            #confirmed walk file counts are correct using gfa file.
            self.sample_labels=False
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
                    if len(set(nodes).intersection(set(self.nodes.keys())))!=0:
                        self.add_walk(f"{sampid}:{hapid}", nodes)
        logger.info('node table construction complete.')

class Link:
    """
    Stores metadata about a link in the graph

    Attributes
    ----------
    link_index : int
        index value for the link. Arbitrary value, used for building link table, starting from 0 (by region)
    node_1 : str
        ID of node on one end of the link
    node_1_orient : str
        orientation of node 1
    node_2 : str
        ID of node on other end of the link
    node_2_orient : str
        orientation of node 2
    node_overlap : str
        CIGAR string indicating the overlap of the two nodes. Likely values== 0M, indicating 0 overlap
    
    nodeid : str
        ID of the node
    length : int
        Length of the sequence at this node
    samples : set of str
        IDs of samples (haplotypes) that go
        through this node

    Methods
    -------
    to_dict(Link)
        Convert the Link object into a dictionary object. For later conversion into dataframe.
    """

    def __init__(self, link_index, node_1, node_1_orient, node_2, node_2_orient, node_overlap='0M'):
        self.link_index=link_index
        self.node_1=node_1
        self.node_1_orient=node_1_orient
        self.node_2=node_2
        self.node_2_orient=node_2_orient
        
    def to_dict(self):
        return {
            "link_index": self.link_index,
            "node_1": self.node_1,
            "node_1_orient": self.node_1_orient,
            "node_2": self.node_2,
            "node_2_orient": self.node_2_orient,
            "node_overlap": getattr(self, "node_overlap", None),
        }

class LinkTable:
    """
    Table of links storing link metadata
    for a region

    Attributes
    ----------
    links : dict[str]->Link
        Dictionary of Links, indexed by arbitrary link_index (counts up from 0 for each region)
    node_1 : str
        ID of node on one end of the link
    node_1_orient : str
        orientation of node 1
    node_2 : str
        ID of node on other end of the link
    node_2_orient : str
        orientation of node 2
    node_overlap : str
        CIGAR string indicating the overlap of the two nodes. Likely values== 0M, indicating 0 overlap

    Methods
    -------
    load_links_from_gfa(gfafile)
        Generate LinkTable from GFA file
    add_link(link)
        Add link to the link_table

    """

    def __init__(self, gfa_file: Path = None, exclude_samples: list[str] = []):
        self.links =  {} 
        if gfa_file is not None:
            self.load_links_from_gfa(gfa_file, exclude_samples)

    def add_link(self, link: Link):
        """
        Add a node to the node table

        Parameters
        ----------
        node : Node
            Node to add
        """
        self.links[link.link_index] = link
        
    def load_links_from_gfa(self, gfa_file: Path, exclude_samples: list[str] = []):
        # First parse all the nodes
        link_index=0 # link index is a placeholder as a way to build the dictionary, relative and therefore meaningless
        with open(gfa_file, "r") as f:
            for line in f:
                linetype = line.split()[0]
                if linetype != "L":
                    continue
                node_1 = line.strip().split()[1]
                node_1_orient = line.strip().split()[2]
                node_2 = line.strip().split()[3]
                node_2_orient = line.strip().split()[4]
                node_overlap=line.strip().split()[5]
                self.add_link(Link(link_index,node_1,node_1_orient,node_2,node_2_orient,node_overlap))
                link_index+=1


#mostly for QC
class Walk:
    """
    Stores metadata about a walk in the graph

    Attributes
    ----------
    sampid : int
    SeqId : str
    SeqStart : str
    SeqEnd : str
    nodes : list
    Walk_str : 
    Methods
    -------
    to_dict(Walk)
        Convert the Walk object into a dictionary object. For later conversion into dataframe.
    """

    def __init__(self, sampid, SeqId, SeqStart, SeqEnd, Nodes, Walk_str):
        self.sampid=sampid
        self.SeqId=SeqId
        self.SeqStart=SeqStart
        self.SeqEnd=SeqEnd
        self.Nodes=Nodes
        self.Walk_str=Walk_str
        
    def to_dict(self):
        return {
            "SampID": self.sampid,
            "SeqId": self.SeqId,
            "SeqStart": self.SeqStart,
            "SeqEnd": self.SeqEnd,
            "Nodes": self.Nodes,
            "Walk_str": self.Walk_str,
            "N_Nodes":len(self.Nodes)
        }

class WalkTable:
    """
    Table of walks storing walk metadata
    for a region

    Attributes
    ----------
    sampid : int
    SeqId : str
    SeqStart : str
    SeqEnd : str
    Nodes : list
    Walk : str
    Methods
    -------
    load_walks_from_gfa(gfafile)
        Generate WalkTable from GFA file
    add_walk(walk)
        Add walk to the walk_table

    """

    def __init__(self, gfa_file: Path = None, exclude_samples: list[str] = []):
        self.walks =  {} 
        if gfa_file is not None:
            self.load_walks_from_gfa(gfa_file, exclude_samples)

    def add_walk(self, walk: Walk):
        """
        Add a node to the node table

        Parameters
        ----------
        node : Node
            Node to add
        """
        self.walks[walk.sampid] = walk
        
    def load_walks_from_gfa(self, gfa_file: Path, exclude_samples: list[str] = []):
        # First parse all the nodes
        with open(gfa_file, "r") as f:
            for line in f:
                linetype = line.split()[0]
                if linetype != "W":
                    continue
                SampID = line.strip().split()[1]
                if SampID in exclude_samples: #reason why number of walks doesn't match the total number of walks is that reference is excluded
                    continue
                SampID=SampID+':'+line.strip().split()[2]
                #hapid = line.strip().split()[2]
                SeqId = line.strip().split()[3]
                SeqStart = line.strip().split()[4]
                SeqEnd=line.strip().split()[5]
                Walk_str=line.strip().split()[6]
                Nodes=(line.strip().split()[6].replace('<','>').split('>')[1:])
                #as nodes is currently defined, the split will remove directionality.
                self.add_walk(Walk(SampID, SeqId, SeqStart, SeqEnd, Nodes, Walk_str))