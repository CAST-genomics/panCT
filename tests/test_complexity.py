import os
from pathlib import Path
from panct.complexity import complexity_score

DATADIR = Path(__file__).parent.joinpath("data")

def test_complexity():
    gfa_short_test_file = DATADIR / "short_test.gfa"
    tsv_short_test_file = DATADIR / "short_test_node_map.tsv"

    complexity = complexity_score(str(gfa_short_test_file), str(tsv_short_test_file))
    assert complexity == 0.875, f"Expected complexity score = 0.875, but got {complexity}"

def test_basic():
    assert True
