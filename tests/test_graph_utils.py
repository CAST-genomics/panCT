import numpy as np
import os
from panct.graph_utils import *
from pathlib import Path
import pytest

DATADIR = Path(__file__).parent.joinpath("data")

def test_node():
	node = Node("n1")
	assert node.nodeid == "n1"
	assert node.length == 0

	node = Node("n2", 20)
	assert node.nodeid == "n2"
	assert node.length == 20
	node.add_sample("samp1")
	node.add_sample("samp2")
	assert len(node.samples) == 2


def test_node_table():
	nt = NodeTable()
	assert nt.numwalks == 0
	assert len(nt.walk_lengths) == 0
	assert nt.get_walk_length([]) == 0
	assert np.isnan(nt.get_mean_walk_length())
	assert np.isnan(nt.get_mean_node_length())
	assert nt.get_total_node_length() == 0
	walk_string = ">12438194>12438195>12438197"
	assert nt.get_nodes_from_walk(walk_string) == \
		["12438194", "12438195", "12438197"]

	walk_string = ">12438194<12438195>12438197"
	assert nt.get_nodes_from_walk(walk_string) == \
		["12438194", "12438195", "12438197"]

	nt.add_node(Node("n1", 200))
	assert nt.get_walk_length(["n1"]) == 200
	with pytest.raises(ValueError):
		nt.get_walk_length(["n1","n2"])
	nt.add_node(Node("n2", 50))
	assert nt.get_walk_length(["n1", "n2"]) == 250
	assert nt.get_mean_node_length() == 125

	nt.add_walk("samp:1", ["n1", "n2"])
	assert nt.numwalks == 1
	assert nt.get_mean_walk_length() == 250
	nt.add_walk("samp:1", ["n1"])
	assert nt.numwalks == 2
	assert nt.get_mean_walk_length() == 225

	# Load from GFA with no seqs
	nt = NodeTable(gfa_file=os.path.join(DATADIR, "test_noseq.gfa"))
	assert nt.get_mean_walk_length() == 95/4
	assert nt.numwalks == 4
	nt = NodeTable(gfa_file=os.path.join(DATADIR, "test_noseq.gfa"),
		exclude_samples=["GRCh38"])
	assert nt.get_mean_walk_length() == 70/3
	assert nt.numwalks == 3

	# Load from GFA with seqs
	nt = NodeTable(gfa_file=os.path.join(DATADIR, "test.gfa"))
	assert nt.get_mean_walk_length() == 38/4
	assert nt.get_mean_node_length() == 5

	# Load from GFA with no lengths
	with pytest.raises(ValueError):
		NodeTable(gfa_file=os.path.join(DATADIR, "test_nolen.gfa"))