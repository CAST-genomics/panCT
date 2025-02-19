import os
from pathlib import Path
import pytest

DATADIR = Path(__file__).parent.joinpath("data")

from panct.complexity import *
from panct.graph_utils import *
from panct.logging import getLogger

# TODO add more tests of main once
# add gbz dependencies to test


def test_complexity_main(tmpdir):
    """
    panct complexity --out test.tab tests/data/basic.gfa
    """
    # Set up defaults
    graph_file = ""
    output_file = ""
    region = ""
    metrics = "sequniq-normwalk,sequniq-normnode"
    reference = ""
    log = getLogger(name="complexity", level="INFO")

    graph_file = Path("dummy")
    assert (
        main(
            graph_file,
            output_file,
            region,
            metrics,
            reference,
            log,
        )
        == 1
    )

    # Process GFA file
    graph_file = DATADIR / "basic.gfa"
    output_file = tmpdir / "test.tab"
    assert (
        main(
            graph_file,
            output_file,
            region,
            metrics,
            reference,
            log,
        )
        == 0
    )
    assert (
        main(
            graph_file,
            output_file,
            "chr1:1-2",
            metrics,
            reference,
            log,
        )
        == 0
    )

    # Invalid metric
    assert (
        main(graph_file, output_file, region, "xxx", reference, log) == 1
    )


def test_compute_complexity():
    # Set up node table
    nt = NodeTable()
    nt.add_node(Node("n1", 1))
    nt.add_node(Node("n2", 1))
    assert compute_complexity(nt, "sequniq-normwalk") is None

    # Add walks
    nt.add_walk("samp:1", ["n1", "n2"])

    # Compute different metrics
    assert compute_complexity(nt, "sequniq-normwalk") == 0
    assert compute_complexity(nt, "sequniq-normnode") == 0
    with pytest.raises(ValueError):
        compute_complexity(nt, "xxx")

    # Add another walk
    nt.add_walk("samp:2", ["n1"])
    assert compute_complexity(nt, "sequniq-normwalk") == 0.5 * (1 - 0.5) / 1.5
    assert compute_complexity(nt, "sequniq-normnode") == 0.5 * (1 - 0.5) / 1
