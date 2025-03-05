import os
from pathlib import Path
from logging import getLogger

import pytest
from typer.testing import CliRunner

from panct.__main__ import app
from panct.graph_utils import Node, NodeTable
from panct.complexity import main, compute_complexity

runner = CliRunner()

DATADIR = Path(__file__).parent.joinpath("data")

expected_basic_output = """numnodes\ttotal_length\tnumwalks\tsequniq-normwalk
2\t10\t3\t0.047619047619047616
"""


def prefix_expected_with_region(expected: str, region: list[str]):
    """
    Prefix the expected output with chrom, start, and end columns

    Parameters
    ----------
    expected: str
        The expected output without chrom, start, and end columns
    region: list[str]
        The region(s) to prefix the expected output with

    Returns
    -------
    str
        The expected output but with the chrom, start, and end columns added at the
        beginning
    """
    assert len(region) == len(expected.split("\n")) - 2
    region_bed_str = ["\t".join(map(str, i)) + "\t" for i in region]
    header = expected.split("\n", maxsplit=1)[0]
    header = "\t".join(("chrom", "start", "end")) + "\t" + header + "\n"
    non_header = filter(
        lambda e: not e.startswith("numnodes") and e != "",
        expected.split("\n"),
    )
    non_header = zip(region_bed_str, non_header)
    return header + "\n".join(s + e for s, e in non_header) + "\n"


def test_basic_stdout(capfd):
    """
    panct complexity tests/data/basic.gfa
    """
    in_file = DATADIR / "basic.gfa"
    expected = expected_basic_output

    cmd = f"complexity {in_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    # check that the output is the same as what we expect
    assert captured.out == expected
    assert result.exit_code == 0


def test_basic_stdout_region(capfd):
    """
    panct complexity --region chrTest:0-1 tests/data/basic.gbz
    """
    region = ("chrTest", 0, 1)
    region_str = f"{region[0]}:{region[1]}-{region[2]}"
    in_file = DATADIR / "basic.gbz"
    expected = expected_basic_output
    expected = prefix_expected_with_region(expected, (region,))

    cmd = f"complexity --region {region_str} {in_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    # check that the output is the same as what we expect
    assert captured.out == expected
    assert result.exit_code == 0


def test_basic_regions_bed(capfd):
    """
    panct complexity --out basic.tsv --region tests/data/basic.bed tests/data/basic.gbz
    """
    in_file = DATADIR / "basic.gbz"
    bed_file = DATADIR / "basic.bed"
    out_file = Path("basic.tsv")
    if out_file.exists():
        out_file.unlink()

    # create expected output
    expected = expected_basic_output.split("\n")
    expected.insert(2, expected[1])
    expected = "\n".join(expected)
    expected = prefix_expected_with_region(
        expected, [("chrTest", 0, 1), ("chrTest", 0, 2)]
    )

    cmd = f"complexity --out {out_file} --region {bed_file} {in_file}"
    result = runner.invoke(app, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    # check that the output is the same as what we expect
    with open(out_file, "r") as f:
        out_file_content = f.read()
    assert out_file_content == expected
    assert result.exit_code == 0

    out_file.unlink()


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
    log = getLogger()

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
    assert main(graph_file, output_file, region, "xxx", reference, log) == 1


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
